//! Separate somatic variant calls from library-preparation damage artifacts.
use std::process;

use std::path::{Path, PathBuf};

use anyhow::{Error, Result};
use chaff::classes::{validate_classes, DamageClass};
use chaff::copied_damage::CopiedDamage;
use chaff::evidence::PileupOptions;
use chaff::filter::{run_filter, FilterArgs, FilterKind, FilterOptions};
use chaff::model::{Distance, Model};
use chaff::read_end::{ATailing, EndRepairFillIn};
use clap::builder::styling::{AnsiColor, Effects, Style, Styles};
use clap::error::ErrorKind;
use clap::parser::ValueSource;
use clap::{ArgMatches, Command, CommandFactory, FromArgMatches, Parser};
use env_logger::Env;
use log::error;
use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

pub(crate) const HEADER: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
pub(crate) const USAGE: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
pub(crate) const LITERAL: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
pub(crate) const PLACEHOLDER: Style = AnsiColor::Cyan.on_default();
pub(crate) const ERROR: Style = AnsiColor::Red.on_default().effects(Effects::BOLD);
pub(crate) const VALID: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
pub(crate) const INVALID: Style = AnsiColor::Yellow.on_default().effects(Effects::BOLD);

/// The tool's one-line description at the top of the help: magenta, the same
/// color code as the [`FOOTER`], so the help opens and closes in one color.
pub(crate) const TITLE: Style = AnsiColor::Magenta.on_default();
/// Indented code examples and backtick-wrapped terms in help text: tertiary
/// color (yellow), a named ANSI slot that follows the user's terminal palette.
pub(crate) const CODE: Style = AnsiColor::Yellow.on_default();
/// The license/attribution footer at the bottom of the help: magenta, shared
/// with the [`TITLE`].
pub(crate) const FOOTER: Style = AnsiColor::Magenta.on_default();
/// The right-hand description column of an indented two-column table: a faded
/// gray (bright-black), so the left-hand code term (in [`CODE`]) stands out.
pub(crate) const FADED: Style = AnsiColor::BrightBlack.on_default();

/// Cargo's color style.
/// [source](https://github.com/crate-ci/clap-cargo/blob/master/src/style.rs)
pub(crate) const CARGO_STYLING: Styles = Styles::styled()
    .header(HEADER)
    .usage(USAGE)
    .literal(LITERAL)
    .placeholder(PLACEHOLDER)
    .error(ERROR)
    .valid(VALID)
    .invalid(INVALID);

/// Separate somatic variant calls from library-preparation damage artifacts.
///
/// Library preparation can turn DNA damage into a base change that both
/// strands of a duplex agree on. chaff reads a coordinate-sorted VCF/BCF of
/// somatic calls and the coordinate-sorted BAM of one of its samples,
/// merge-joins them without an index, and scores each SNV whose genotype is
/// heterozygous or missing by where its alternate molecules sit compared with
/// the reference molecules at the same site; homozygous, haploid, and indel
/// calls pass unscored. It writes the calls with INFO annotations and, past a
/// threshold, FILTERs.
///
/// MENTAL MODEL
///
///  1. Each template covering a call is one molecule; mates count once
///  2. A molecule's position is its distance from a template end
///  3. A true mutation's molecules sit where the reference molecules sit
///  4. An artifact's molecules crowd the end that made the artifact
///  5. Each call gets a likelihood ratio, artifact to mutation
///  6. EM over all calls learns each stratum's artifact fraction (the prior)
///     and each decay's scale
///  7. The posterior probability of a true mutation is written per call
///
/// MODELS
///
///   chaff   a learned prior, and decays with learned scales for end repair
///           fill-in and copied damage; A-tailing uses a window
///   fgbio   fgbio's (2 * maf)^2 prior and windows, as fgbio FilterSomaticVcf;
///           copied damage, which fgbio lacks, still decays
///
/// FILTERS
///
///   copied-damage        damage copied onto the other strand before strand
///                        tagging: alternates crowd the lesion strand's 5' end
///   end-repair-fill-in   errors in a filled-in recessed 3' end: alternates
///                        crowd the 3' end of the strand each template was
///                        copied from, or either end under fgbio (ERFAP)
///   a-tailing            adenines added to an over-digested 3' end: a T near
///                        the left end or an A near the right (fgbio ATAP)
///
/// EXAMPLES
///
///  1. Annotate every filter and write the per-sample metrics:
///
///   chaff -i calls.vcf.gz -b tumor.bam -r ref.fa -o out.vcf.gz \
///       --metrics tumor.chaff.tsv
///
///  2. Filter copied deamination in a Duplex Sequencing library:
///
///   chaff -i calls.vcf.gz -b tumor.bam -r ref.fa -s tumor \
///       -o calls.chaff.vcf.gz --metrics tumor.chaff.tsv \
///       --filters copied-damage --copied-damage-classes 'C>T' \
///       --copied-damage-threshold 0.05
///
///  3. Reproduce fgbio FilterSomaticVcf:
///
///   chaff -i calls.vcf -b tumor.bam -o out.vcf --model fgbio \
///       --filters end-repair-fill-in,a-tailing
#[derive(Debug, Parser)]
#[command(
    author,
    version,
    color = clap::ColorChoice::Always,
    rename_all = "kebab-case",
    verbatim_doc_comment,
    arg_required_else_help = true
)]
#[clap(styles = CARGO_STYLING)]
struct Cli {
    /// Input VCF/BCF of somatic calls, coordinate-sorted.
    ///
    /// Read twice (once to learn priors, once to write), so it must be a file.
    #[arg(short = 'i', long, value_name = "VCF", value_parser = file, verbatim_doc_comment)]
    input: PathBuf,

    /// Output VCF/BCF; the format follows the extension.
    ///
    ///   out.vcf      plain VCF
    ///   out.vcf.gz   BGZF-compressed VCF
    ///   out.bcf      BCF
    ///   -            plain VCF to standard output
    #[arg(short = 'o', long, value_name = "VCF", verbatim_doc_comment)]
    output: PathBuf,

    /// Coordinate-sorted BAM of the sample under test; no index needed.
    ///
    /// Distances from template ends count template bases through both reads'
    /// CIGARs, so a read whose mate maps to the same contig needs the mate
    /// CIGAR (`MC`) tag; the insert size (`TLEN`) is never read.
    #[arg(short = 'b', long, value_name = "BAM", verbatim_doc_comment)]
    bam: PathBuf,

    /// Indexed reference FASTA (`.fai` alongside) for CpG context.
    ///
    /// Required by the `copied-damage` filter and by `--spectrum`; giving it
    /// with neither is an error.
    #[arg(short = 'r', long = "ref", value_name = "FASTA", verbatim_doc_comment)]
    reference: Option<PathBuf>,

    /// The sample whose reads are in the BAM.
    ///
    /// Required when the VCF has more than one sample.
    #[arg(short = 's', long, value_name = "NAME", verbatim_doc_comment)]
    sample: Option<String>,

    /// Per-sample metrics TSV: one row per filter and stratum.
    #[arg(long, value_name = "TSV", value_parser = report_file, verbatim_doc_comment)]
    metrics: Option<PathBuf>,

    /// PDF of the sample's SNVs by trinucleotide context, before and after
    /// filtering, its panels on one scale.
    ///
    ///   SNVs                 every heterozygous SNV, whatever its FILTER
    ///   Expected Real SNVs   each weighed by the product of the posteriors
    ///                        of the filters --filters enables
    ///   Passing SNVs         those no filter flagged, when a filter has a
    ///                        threshold
    #[arg(long, value_name = "PDF", value_parser = report_file, verbatim_doc_comment)]
    spectrum: Option<PathBuf>,

    /// Minimum mapping quality of a read.
    #[arg(
        short = 'm',
        long,
        value_name = "MAPQ",
        default_value_t = 20,
        verbatim_doc_comment
    )]
    min_mapping_quality: u8,

    /// Minimum base quality at the call.
    #[arg(
        short = 'q',
        long,
        value_name = "QUAL",
        default_value_t = 20,
        verbatim_doc_comment
    )]
    min_base_quality: u8,

    /// Use only paired reads whose mate is also mapped.
    ///
    /// Duplicate, secondary, and supplementary reads are always left out.
    #[arg(short = 'p', long, verbatim_doc_comment)]
    paired_reads_only: bool,

    /// The filters to run, comma-separated.
    ///
    ///   --filters copied-damage                   only the copied damage filter
    ///   --filters end-repair-fill-in,a-tailing    only the fgbio filters
    #[arg(
        long,
        value_enum,
        value_delimiter = ',',
        value_name = "FILTER",
        default_value = "copied-damage,a-tailing,end-repair-fill-in",
        hide_possible_values = true,
        verbatim_doc_comment
    )]
    filters: Vec<FilterKind>,

    /// The model that scores calls: its prior and its distance models.
    ///
    ///   chaff   artifact fractions learned per sample and stratum by EM, and
    ///           decays for end repair fill-in and copied damage
    ///   fgbio   fgbio's per-call (2 * maf)^2 prior and windows, for parity
    #[arg(
        long,
        value_enum,
        value_name = "MODEL",
        default_value_t = Model::Chaff,
        hide_possible_values = true,
        verbatim_doc_comment
    )]
    model: Model,

    /// Damage classes as damaged base `>` read base, comma-separated.
    ///
    /// A class matches on either strand: `C>T` covers REF/ALT `C/T` (lesion on
    /// the forward strand) and `G/A` (lesion on the reverse strand).
    ///
    ///   C>T   cytosine or 5-methylcytosine deamination
    ///   G>T   guanine oxidation to 8-oxoguanine
    #[arg(
        long,
        value_delimiter = ',',
        value_name = "CLASS",
        default_value = "C>T,G>T",
        verbatim_doc_comment
    )]
    copied_damage_classes: Vec<DamageClass>,

    /// Decay scale in bases from the lesion strand's 5' end, or `learned`.
    ///
    /// The mean length over which a polymerase copies a lesion strand onto its
    /// partner, learned per library under either model unless a number fixes
    /// it. Molecules within it count as congruent in the metrics.
    #[arg(
        long,
        value_name = "BP",
        default_value = "learned",
        verbatim_doc_comment
    )]
    copied_damage_distance: Distance,

    /// Apply `CopiedDamageArtifact` at or below this posterior.
    #[arg(long, value_name = "P", value_parser = probability, verbatim_doc_comment)]
    copied_damage_threshold: Option<f64>,

    /// Distance in bases from a template end, or `learned`.
    ///
    ///   --model chaff   the decay scale from the 3' end of the strand each
    ///                   template was copied from, learned unless fixed
    ///   --model fgbio   fgbio's window from the nearest template end, 15
    ///                   unless fixed
    ///
    /// Molecules within it count as congruent in the metrics.
    #[arg(
        long,
        value_name = "BP",
        default_value = "learned",
        verbatim_doc_comment
    )]
    end_repair_fill_in_distance: Distance,

    /// Apply `EndRepairFillInArtifact` at or below this posterior.
    #[arg(
        long,
        alias = "end-repair-fill-in-p-value",
        value_name = "P",
        value_parser = probability,
        verbatim_doc_comment
    )]
    end_repair_fill_in_threshold: Option<f64>,

    /// Window in bases from the template end, under either model.
    #[arg(
        long,
        value_name = "BP",
        default_value_t = 2,
        value_parser = clap::value_parser!(u32).range(1..),
        verbatim_doc_comment
    )]
    a_tailing_distance: u32,

    /// Apply `ATailingArtifact` at or below this posterior.
    #[arg(
        long,
        alias = "a-tailing-p-value",
        value_name = "P",
        value_parser = probability,
        verbatim_doc_comment
    )]
    a_tailing_threshold: Option<f64>,
}

/// Parse a path to a file, refusing `-` for standard input.
fn file(text: &str) -> Result<PathBuf, String> {
    match text {
        "-" => Err("the input is read twice, so it must be a file, not standard input".into()),
        _ => Ok(PathBuf::from(text)),
    }
}

/// Parse a path to a report file, refusing `-`, since standard output carries
/// the VCF.
fn report_file(text: &str) -> Result<PathBuf, String> {
    match text {
        "-" => Err("standard output carries the VCF, so name a file".into()),
        _ => Ok(PathBuf::from(text)),
    }
}

/// The file a path names, its links and relative parts resolved, when the file
/// or its directory exists.
fn resolve(path: &Path) -> Option<PathBuf> {
    if let Ok(path) = path.canonicalize() {
        return Some(path);
    }
    let directory = match path.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent,
        _ => Path::new("."),
    };
    Some(directory.canonicalize().ok()?.join(path.file_name()?))
}

/// Parse a probability from zero to one.
fn probability(text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(value) if (0.0..=1.0).contains(&value) => Ok(value),
        _ => Err(format!("expected a probability from 0 to 1, found: {text}")),
    }
}

impl Cli {
    /// Reject, as usage errors of `cmd`, an option typed on the command line for
    /// a filter `--filters` leaves out, the copied damage filter or the
    /// spectrum without a reference, an input and outputs that name one file,
    /// and damage classes that repeat a change.
    fn validate(&self, matches: &ArgMatches, cmd: &mut Command) -> Result<(), clap::Error> {
        for kind in FilterKind::ALL {
            if self.filters.contains(&kind) {
                continue;
            }
            for id in kind.arguments() {
                if matches.value_source(id) != Some(ValueSource::CommandLine)
                    || (*id == "reference" && self.spectrum.is_some())
                {
                    continue;
                }
                let arg = cmd
                    .get_arguments()
                    .find(|arg| arg.get_id() == id)
                    .map(ToString::to_string)
                    .unwrap_or_else(|| id.to_string());
                let message = format!(
                    "the argument '{arg}' applies only to the {kind} filter, which '--filters' leaves out"
                );
                return Err(cmd.error(ErrorKind::ArgumentConflict, message));
            }
        }
        if self.model == Model::Fgbio
            && self.end_repair_fill_in_distance == Distance::Learned
            && matches.value_source("end_repair_fill_in_distance") == Some(ValueSource::CommandLine)
        {
            return Err(cmd.error(
                ErrorKind::ArgumentConflict,
                "'--end-repair-fill-in-distance learned' needs '--model chaff'; under '--model fgbio' it is a window of bases",
            ));
        }
        if self.filters.contains(&FilterKind::CopiedDamage) && self.reference.is_none() {
            return Err(cmd.error(
                ErrorKind::MissingRequiredArgument,
                "the copied-damage filter needs a reference FASTA: '--ref <FASTA>'",
            ));
        }
        if self.spectrum.is_some() && self.reference.is_none() {
            return Err(cmd.error(
                ErrorKind::MissingRequiredArgument,
                "'--spectrum' needs a reference FASTA: '--ref <FASTA>'",
            ));
        }
        let files = [
            ("--input", Some(&self.input)),
            (
                "--output",
                Some(&self.output).filter(|p| *p != Path::new("-")),
            ),
            ("--metrics", self.metrics.as_ref()),
            ("--spectrum", self.spectrum.as_ref()),
        ];
        for (i, (a, first)) in files.iter().enumerate() {
            for (b, second) in &files[i + 1..] {
                if let (Some(first), Some(second)) = (first, second) {
                    if resolve(first).is_some_and(|path| Some(path) == resolve(second)) {
                        let message = format!("'{a}' and '{b}' name the same file: {first:?}");
                        return Err(cmd.error(ErrorKind::ArgumentConflict, message));
                    }
                }
            }
        }
        let reads = [
            ("--bam", Some(&self.bam)),
            ("--ref", self.reference.as_ref()),
        ];
        for (a, written) in &files[1..] {
            for (b, read) in &reads {
                if let (Some(written), Some(read)) = (written, read) {
                    if resolve(written).is_some_and(|path| Some(path) == resolve(read)) {
                        let message = format!("'{a}' and '{b}' name the same file: {written:?}");
                        return Err(cmd.error(ErrorKind::ArgumentConflict, message));
                    }
                }
            }
        }
        validate_classes(&self.copied_damage_classes)
            .map_err(|error| cmd.error(ErrorKind::ValueValidation, error))
    }

    /// Gather the options into [`FilterArgs`].
    fn into_args(self) -> FilterArgs {
        let options = FilterOptions {
            sample: self.sample,
            filters: self.filters,
            model: self.model,
            copied_damage: CopiedDamage {
                classes: self.copied_damage_classes,
                distance: self.copied_damage_distance,
            },
            copied_damage_threshold: self.copied_damage_threshold,
            end_repair_fill_in: EndRepairFillIn {
                distance: self.end_repair_fill_in_distance,
            },
            end_repair_fill_in_threshold: self.end_repair_fill_in_threshold,
            a_tailing: ATailing {
                distance: self.a_tailing_distance,
            },
            a_tailing_threshold: self.a_tailing_threshold,
        };
        FilterArgs {
            input: self.input,
            output: self.output,
            bam: self.bam,
            reference: self.reference,
            metrics: self.metrics,
            spectrum: self.spectrum,
            pileup: PileupOptions {
                min_mapping_quality: self.min_mapping_quality,
                min_base_quality: self.min_base_quality,
                paired_reads_only: self.paired_reads_only,
            },
            options,
        }
    }
}

/// The ANSI escape that starts `style`, or an empty string when `color` is off
/// (honoring `NO_COLOR`). Its matching reset comes from [`esc_reset`].
fn esc(style: Style, color: bool) -> String {
    if color {
        style.render().to_string()
    } else {
        String::new()
    }
}

/// The reset escape for `style` (a full SGR reset), or empty when `color` is off.
fn esc_reset(style: Style, color: bool) -> String {
    if color {
        style.render_reset().to_string()
    } else {
        String::new()
    }
}

/// Paint the first line of `text` in the [`TITLE`] color (the tool's one-line
/// description), leaving the rest of the text untouched.
fn title_first_line(text: &str, color: bool) -> String {
    let title = esc(TITLE, color);
    let reset = esc_reset(TITLE, color);
    match text.split_once('\n') {
        Some((first, rest)) => format!("{title}{first}{reset}\n{rest}"),
        None => format!("{title}{text}{reset}"),
    }
}

/// Drop the backticks around terms (`` `like this` ``), painting each term with
/// the `code` escape. `after` is the escape that restores the surrounding text's
/// color once a term ends.
fn paint_backtick_terms(text: &str, code: &str, after: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(open) = rest.find('`') {
        out.push_str(&rest[..open]);
        let tail = &rest[open + 1..];
        match tail.find('`') {
            Some(close) => {
                out.push_str(code);
                out.push_str(&tail[..close]);
                out.push_str(after);
                rest = &tail[close + 1..];
            }
            None => {
                out.push('`');
                rest = tail;
                break;
            }
        }
    }
    out.push_str(rest);
    out
}

/// Split an indented two-column row into (left, gap, right) at the first run of
/// three or more spaces that follows the left-hand term.
fn split_two_column(content: &str) -> Option<(&str, &str, &str)> {
    let bytes = content.as_bytes();
    let mut seen_term = false;
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b' ' {
            let start = i;
            while i < bytes.len() && bytes[i] == b' ' {
                i += 1;
            }
            if seen_term && i - start >= 3 {
                return Some((&content[..start], &content[start..i], &content[i..]));
            }
        } else {
            seen_term = true;
            i += 1;
        }
    }
    None
}

/// Style one help string, dropping backticks throughout:
///
/// - indented two-column rows: the left code term in [`CODE`], the right
///   description column in [`FADED`];
/// - indented description continuations (deeper indent, no term): all [`FADED`];
/// - indented single-column examples (a bare command), and the lines a
///   trailing `\` continues: all [`CODE`];
/// - prose: default color, with backtick terms painted in [`CODE`].
fn style_help_text(text: &str, color: bool) -> String {
    let code = esc(CODE, color);
    let faded = esc(FADED, color);
    let reset = esc_reset(CODE, color);
    let mut continues = false;
    text.split_inclusive('\n')
        .map(|line| {
            let (content, newline) = match line.strip_suffix('\n') {
                Some(content) => (content, "\n"),
                None => (line, ""),
            };
            let continued = std::mem::replace(&mut continues, content.ends_with('\\'));
            if !content.starts_with("  ") || content.trim().is_empty() {
                return format!("{}{newline}", paint_backtick_terms(content, &code, &reset));
            }
            if continued {
                let painted = paint_backtick_terms(content, &code, &code);
                format!("{code}{painted}{reset}{newline}")
            } else if let Some((left, gap, right)) = split_two_column(content) {
                let left = paint_backtick_terms(left, &code, &code);
                let right = paint_backtick_terms(right, &code, &faded);
                format!("{code}{left}{reset}{gap}{faded}{right}{reset}{newline}")
            } else if content.len() - content.trim_start().len() > 2 {
                let painted = paint_backtick_terms(content, &code, &faded);
                format!("{faded}{painted}{reset}{newline}")
            } else {
                let painted = paint_backtick_terms(content, &code, &code);
                format!("{code}{painted}{reset}{newline}")
            }
        })
        .collect()
}

/// Style a command's help: paint the first line of the about text in the
/// [`TITLE`] color, apply [`style_help_text`] to the abouts and every option's
/// help, and add the license footer.
fn decorate_help(cmd: clap::Command, color: bool) -> clap::Command {
    let about = cmd.get_about().map(ToString::to_string);
    let long_about = cmd
        .get_long_about()
        .map(ToString::to_string)
        .or_else(|| about.clone());

    let footer = format!(
        "{}MIT License 2026 · Clint Valentine{}",
        esc(FOOTER, color),
        esc_reset(FOOTER, color),
    );

    let choice = if color {
        clap::ColorChoice::Always
    } else {
        clap::ColorChoice::Never
    };
    let mut cmd = cmd.color(choice);
    if let Some(about) = about {
        cmd = cmd.about(title_first_line(&style_help_text(&about, color), color));
    }
    if let Some(long_about) = long_about {
        cmd = cmd.long_about(title_first_line(
            &style_help_text(&long_about, color),
            color,
        ));
    }
    cmd = cmd
        .next_line_help(true)
        .after_help(footer.clone())
        .after_long_help(footer);
    cmd.mut_args(|mut arg| {
        if let Some(help) = arg.get_help().map(ToString::to_string) {
            arg = arg.help(style_help_text(&help, color));
        }
        if let Some(long_help) = arg.get_long_help().map(ToString::to_string) {
            arg = arg.long_help(style_help_text(&long_help, color));
        }
        arg
    })
}

/// Main binary entrypoint.
#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), Error> {
    let color = std::env::var_os("NO_COLOR").is_none();

    let env = Env::default().default_filter_or("info,usvg=error");
    let write_style = if color {
        env_logger::WriteStyle::Auto
    } else {
        env_logger::WriteStyle::Never
    };
    env_logger::Builder::from_env(env)
        .write_style(write_style)
        .init();

    // Help text is hand-wrapped, and clap would count the injected ANSI escapes as width.
    let mut cmd = decorate_help(Cli::command().term_width(usize::MAX), color);
    let matches = cmd.get_matches_mut();
    let cli = Cli::from_arg_matches(&matches).unwrap_or_else(|e| e.exit());
    cli.validate(&matches, &mut cmd)
        .unwrap_or_else(|e| e.exit());

    match run_filter(&cli.into_args()) {
        Ok(()) => process::exit(0),
        Err(e) => {
            error!("{e:#}");
            process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    fn args(extra: &[&str]) -> Result<FilterArgs, clap::Error> {
        let mut argv = vec!["chaff", "-b", "in.bam"];
        for (flag, default) in [("-i", "in.vcf"), ("-o", "out.vcf")] {
            if !extra.contains(&flag) {
                argv.extend([flag, default]);
            }
        }
        argv.extend(extra);
        let mut cmd = Cli::command().color(clap::ColorChoice::Never);
        let matches = cmd.try_get_matches_from_mut(argv)?;
        let cli = Cli::from_arg_matches(&matches)?;
        cli.validate(&matches, &mut cmd)?;
        Ok(cli.into_args())
    }

    #[rstest]
    #[case(&["--filters", "a-tailing", "--copied-damage-threshold", "0.05"], "the argument '--copied-damage-threshold <P>' applies only to the copied-damage filter")]
    #[case(&["--filters", "a-tailing", "--copied-damage-classes", "C>T"], "the argument '--copied-damage-classes <CLASS>' applies only to the copied-damage filter")]
    #[case(&["--filters", "copied-damage", "--a-tailing-distance", "2"], "the argument '--a-tailing-distance <BP>' applies only to the a-tailing filter")]
    #[case(&["--filters", "copied-damage", "--a-tailing-p-value", "0.01"], "the argument '--a-tailing-threshold <P>' applies only to the a-tailing filter")]
    #[case(&["--filters", "a-tailing", "--end-repair-fill-in-distance", "15"], "the argument '--end-repair-fill-in-distance <BP>' applies only to the end-repair-fill-in filter")]
    #[case(&["--filters", "a-tailing", "--copied-damage-distance", "20"], "the argument '--copied-damage-distance <BP>' applies only to the copied-damage filter")]
    #[case(&["--filters", "a-tailing,end-repair-fill-in", "--ref", "ref.fa"], "the argument '--ref <FASTA>' applies only to the copied-damage filter")]
    fn test_an_option_of_a_filter_left_out_is_a_usage_error(
        #[case] extra: &[&str],
        #[case] message: &str,
    ) {
        let error = args(extra).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::ArgumentConflict);
        assert_eq!(error.exit_code(), 2);
        assert!(error.to_string().contains(message), "{error}");
        assert!(error.to_string().contains("Usage: chaff"), "{error}");
    }

    #[test]
    fn test_the_copied_damage_filter_without_a_reference_is_a_usage_error() {
        let error = args(&["--filters", "copied-damage,a-tailing"]).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::MissingRequiredArgument);
        assert_eq!(error.exit_code(), 2);
        let message = "the copied-damage filter needs a reference FASTA: '--ref <FASTA>'";
        assert!(error.to_string().contains(message), "{error}");
    }

    #[test]
    fn test_the_spectrum_without_a_reference_is_a_usage_error() {
        let error = args(&["--filters", "a-tailing", "--spectrum", "out.pdf"]).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::MissingRequiredArgument);
        let message = "'--spectrum' needs a reference FASTA: '--ref <FASTA>'";
        assert!(error.to_string().contains(message), "{error}");
        let extra = [
            "--filters",
            "a-tailing",
            "--spectrum",
            "out.pdf",
            "--ref",
            "ref.fa",
        ];
        assert_eq!(
            args(&extra).unwrap().spectrum,
            Some(PathBuf::from("out.pdf"))
        );
    }

    #[test]
    fn test_standard_input_is_a_usage_error() {
        let error = args(&["-i", "-", "--ref", "ref.fa"]).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::ValueValidation);
        assert_eq!(error.exit_code(), 2);
        let message = "invalid value '-' for '--input <VCF>': the input is read twice";
        assert!(error.to_string().contains(message), "{error}");
    }

    #[rstest]
    #[case("--metrics", "<TSV>")]
    #[case("--spectrum", "<PDF>")]
    fn test_a_report_to_standard_output_is_a_usage_error(
        #[case] option: &str,
        #[case] value: &str,
    ) {
        let error = args(&["--ref", "ref.fa", option, "-"]).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::ValueValidation);
        let message =
            format!("invalid value '-' for '{option} {value}': standard output carries the VCF");
        assert!(error.to_string().contains(&message), "{error}");
    }

    #[test]
    fn test_outputs_that_name_the_input_or_each_other_are_a_usage_error() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("calls.vcf");
        std::fs::write(&input, "").unwrap();
        let link = dir.path().join("link.vcf");
        std::os::unix::fs::symlink(&input, &link).unwrap();
        let path = |p: &std::path::Path| p.to_str().unwrap().to_string();
        let (input, link) = (path(&input), path(&link));
        let metrics = path(&dir.path().join("out.tsv"));
        let filters = ["--filters", "a-tailing"];
        for (extra, message) in [
            (
                vec!["-o", &input],
                "'--input' and '--output' name the same file",
            ),
            (
                vec!["-o", &link],
                "'--input' and '--output' name the same file",
            ),
            (
                vec!["--metrics", &input],
                "'--input' and '--metrics' name the same file",
            ),
            (
                vec!["-o", &metrics, "--metrics", &metrics],
                "'--output' and '--metrics' name the same file",
            ),
            (
                vec![
                    "--metrics",
                    &metrics,
                    "--spectrum",
                    &metrics,
                    "--ref",
                    "ref.fa",
                ],
                "'--metrics' and '--spectrum' name the same file",
            ),
        ] {
            let error = args(&[&["-i", &input][..], &filters, &extra].concat()).unwrap_err();
            assert_eq!(error.kind(), ErrorKind::ArgumentConflict);
            assert_eq!(error.exit_code(), 2);
            assert!(error.to_string().contains(message), "{error}");
        }
        args(&[&["-i", &input, "-o", "-"][..], &filters].concat()).unwrap();
    }

    #[rstest]
    #[case("0")]
    #[case("-1")]
    #[case("wide")]
    fn test_a_distance_that_is_not_learned_or_positive_is_a_usage_error(#[case] value: &str) {
        let option = format!("--end-repair-fill-in-distance={value}");
        let error = args(&["--filters", "end-repair-fill-in", &option]).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::ValueValidation);
        assert_eq!(error.exit_code(), 2);
        let message = "expected 'learned' or a positive number of bases";
        assert!(error.to_string().contains(message), "{error}");
    }

    #[test]
    fn test_a_learned_end_repair_window_under_fgbio_is_a_usage_error() {
        let extra = ["--filters", "end-repair-fill-in", "--model", "fgbio"];
        let error = args(&[&extra[..], &["--end-repair-fill-in-distance", "learned"]].concat())
            .unwrap_err();
        assert_eq!(error.kind(), ErrorKind::ArgumentConflict);
        assert_eq!(error.exit_code(), 2);
        assert!(
            error.to_string().contains("needs '--model chaff'"),
            "{error}"
        );
        let window = args(&extra).unwrap().options.end_repair_fill_in.window();
        assert_eq!(window, 15.0);
    }

    #[test]
    fn test_a_damage_class_and_its_reverse_complement_are_a_usage_error() {
        let error = args(&["--ref", "ref.fa", "--copied-damage-classes", "C>T,G>A"]).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::ValueValidation);
        assert_eq!(error.exit_code(), 2);
        let message = "damage classes C>T and G>A describe the same change on opposite strands";
        assert!(error.to_string().contains(message), "{error}");
    }

    #[test]
    fn test_filter_arguments_are_every_option_named_for_a_filter() {
        let ids: Vec<String> = Cli::command()
            .get_arguments()
            .map(|arg| arg.get_id().to_string())
            .collect();
        for kind in FilterKind::ALL {
            let prefix = format!("{}_", kind.to_string().replace('-', "_"));
            for id in ids.iter().filter(|id| id.starts_with(&prefix)) {
                assert!(kind.arguments().contains(&id.as_str()), "{kind}: {id}");
            }
            for id in kind.arguments() {
                assert!(ids.iter().any(|known| known == id), "{kind}: {id}");
            }
        }
    }

    #[rstest]
    #[case(&["--filters", "a-tailing"])]
    #[case(&["--ref", "ref.fa"])]
    #[case(&["--filters", "copied-damage", "--ref", "ref.fa"])]
    #[case(&["--filters", "a-tailing", "--a-tailing-distance", "4", "--a-tailing-threshold", "0.001"])]
    #[case(&["--filters", "end-repair-fill-in", "--end-repair-fill-in-distance", "10"])]
    #[case(&["--filters", "end-repair-fill-in", "--model", "fgbio", "--end-repair-fill-in-distance", "10"])]
    #[case(&["--ref", "ref.fa", "--model", "fgbio", "--copied-damage-distance", "20"])]
    fn test_options_of_enabled_filters_and_defaults_are_accepted(#[case] extra: &[&str]) {
        args(extra).unwrap();
    }

    #[test]
    fn test_list_defaults_are_every_filter_and_both_classes_shown_comma_separated() {
        let options = args(&["--ref", "ref.fa"]).unwrap().options;
        assert_eq!(options.filters, FilterKind::ALL);
        let classes = [DamageClass::DEAMINATION, DamageClass::OXIDATION];
        assert_eq!(options.copied_damage.classes, classes);
        let help = Cli::command().render_long_help().to_string();
        assert!(help.contains("[default: copied-damage,a-tailing,end-repair-fill-in]"));
        assert!(help.contains("[default: C>T,G>T]"));
    }

    #[test]
    fn test_an_a_tailing_window_of_zero_bases_is_a_usage_error() {
        let error = args(&["--filters", "a-tailing", "--a-tailing-distance", "0"]).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::ValueValidation);
        assert_eq!(error.exit_code(), 2);
        let message = "invalid value '0' for '--a-tailing-distance <BP>'";
        assert!(error.to_string().contains(message), "{error}");
    }

    #[test]
    fn test_a_command_continued_by_a_backslash_is_code_on_every_line() {
        let text = "  chaff -i in.vcf \\\n      --metrics out.tsv\n      a description\n";
        let (code, faded, reset) = (esc(CODE, true), esc(FADED, true), esc_reset(CODE, true));
        let expected = format!(
            "{code}  chaff -i in.vcf \\{reset}\n{code}      --metrics out.tsv{reset}\n{faded}      a description{reset}\n"
        );
        assert_eq!(style_help_text(text, true), expected);
    }

    #[test]
    fn test_outputs_that_name_the_bam_or_the_reference_are_a_usage_error() {
        let filters = ["--filters", "a-tailing"];
        for (extra, message) in [
            (
                [&filters[..], &["-o", "in.bam"]].concat(),
                "'--output' and '--bam' name the same file",
            ),
            (
                [&filters[..], &["--metrics", "./in.bam"]].concat(),
                "'--metrics' and '--bam' name the same file",
            ),
            (
                vec!["--ref", "ref.fa", "--metrics", "ref.fa"],
                "'--metrics' and '--ref' name the same file",
            ),
        ] {
            let error = args(&extra).unwrap_err();
            assert_eq!(error.kind(), ErrorKind::ArgumentConflict);
            assert_eq!(error.exit_code(), 2);
            assert!(error.to_string().contains(message), "{error}");
        }
    }
}
