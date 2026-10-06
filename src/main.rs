//! Separate somatic variant calls from library-preparation damage artifacts.
use std::process;

use std::path::PathBuf;

use anyhow::{bail, Error, Result};
use chaff::classes::{validate_classes, DamageClass};
use chaff::copied_damage::CopiedDamage;
use chaff::filter::{run_filter, FilterArgs, FilterKind, FilterOptions};
use chaff::prior::PriorMode;
use chaff::read_end::{ATailing, EndRepairFillIn};
use chaff::template::ReadFilter;
use clap::builder::styling::{AnsiColor, Effects, Style, Styles};
use clap::parser::ValueSource;
use clap::{ArgMatches, CommandFactory, FromArgMatches, Parser};
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
/// merge-joins them without an index, and scores each call by where its
/// alternate molecules sit compared with the reference molecules at the same
/// site. It writes the calls with INFO annotations and, past a threshold,
/// FILTERs.
///
/// MENTAL MODEL
///
///  1. Each template covering a call is one molecule; mates count once
///  2. A molecule's position is its distance from a template end
///  3. A true mutation's molecules sit where the reference molecules sit
///  4. An artifact's molecules crowd the end that made the artifact
///  5. Each call gets a likelihood ratio, artifact to mutation
///  6. EM over all calls learns each stratum's artifact fraction (the prior)
///  7. The posterior probability of a true mutation is written per call
///
/// FILTERS
///
///   copied-damage        damage copied onto the other strand before strand
///                        tagging: alternates crowd the lesion strand's 5' end
///   end-repair-fill-in   damage copied into a filled-in recessed 3' end:
///                        alternates crowd either template end (fgbio ERFAP)
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
///  2. Filter an end-repaired, A-tailed duplex library with suspected deamination:
///
///   chaff -i calls.vcf.gz -b tumor.bam -r ref.fa -s tumor \
///       -o calls.chaff.vcf.gz --metrics tumor.chaff.tsv \
///       --filters copied-damage,end-repair-fill-in,a-tailing \
///       --copied-damage-classes C>T --copied-damage-threshold 0.05 \
///       --end-repair-fill-in-threshold 0.001 --a-tailing-threshold 0.001
///
///  3. Reproduce fgbio FilterSomaticVcf, including its prior:
///
///   chaff -i calls.vcf -b tumor.bam -o out.vcf --prior fgbio \
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
    #[arg(short = 'i', long, value_name = "VCF", verbatim_doc_comment)]
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
    /// A read whose mate maps to the same contig needs the mate CIGAR (`MC`)
    /// tag; the insert size (`TLEN`) is never read.
    #[arg(short = 'b', long, value_name = "BAM", verbatim_doc_comment)]
    bam: PathBuf,

    /// Indexed reference FASTA (`.fai` alongside) for CpG context.
    ///
    /// Required by the `copied-damage` filter.
    #[arg(short = 'r', long = "ref", value_name = "FASTA", verbatim_doc_comment)]
    reference: Option<PathBuf>,

    /// The sample whose reads are in the BAM.
    ///
    /// Required when the VCF has more than one sample.
    #[arg(short = 's', long, value_name = "NAME", verbatim_doc_comment)]
    sample: Option<String>,

    /// Per-sample metrics TSV: one row per filter and stratum.
    #[arg(long, value_name = "TSV", verbatim_doc_comment)]
    metrics: Option<PathBuf>,

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
        default_values_t = FilterKind::ALL,
        hide_possible_values = true,
        verbatim_doc_comment
    )]
    filters: Vec<FilterKind>,

    /// The prior that turns a likelihood ratio into a posterior.
    ///
    ///   learned   artifact fraction per sample and stratum, learned by EM
    ///   fgbio     fgbio's per-call mutation prior, (2 * maf)^2, for parity
    #[arg(
        long,
        value_enum,
        value_name = "PRIOR",
        default_value_t = PriorMode::Learned,
        hide_possible_values = true,
        verbatim_doc_comment
    )]
    prior: PriorMode,

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
        default_values = ["C>T", "G>T"],
        verbatim_doc_comment
    )]
    copied_damage_classes: Vec<DamageClass>,

    /// Mean length in bases of the resynthesis that copies a lesion.
    #[arg(long, value_name = "BP", default_value_t = 30.0, value_parser = positive, verbatim_doc_comment)]
    copied_damage_scale: f64,

    /// Apply `CopiedDamageArtifact` at or below this posterior.
    #[arg(long, value_name = "P", value_parser = probability, verbatim_doc_comment)]
    copied_damage_threshold: Option<f64>,

    /// Distance from a template end within which end repair fill-in acts.
    #[arg(long, value_name = "BP", default_value_t = 15, verbatim_doc_comment)]
    end_repair_fill_in_distance: u32,

    /// Score with a decay of this scale in bases from the nearest template end
    /// instead of the distance window.
    ///
    /// The distance window still sets the congruent molecules in the metrics.
    #[arg(long, value_name = "BP", value_parser = positive, verbatim_doc_comment)]
    end_repair_fill_in_scale: Option<f64>,

    /// Apply `EndRepairFillInArtifact` at or below this posterior.
    #[arg(
        long,
        alias = "end-repair-fill-in-p-value",
        value_name = "P",
        value_parser = probability,
        verbatim_doc_comment
    )]
    end_repair_fill_in_threshold: Option<f64>,

    /// Distance from a template end within which A-tailing acts.
    #[arg(long, value_name = "BP", default_value_t = 2, verbatim_doc_comment)]
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

/// Parse a finite length greater than zero.
fn positive(text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(value) if value.is_finite() && value > 0.0 => Ok(value),
        _ => Err(format!("expected a positive number, found: {text}")),
    }
}

/// Parse a probability from zero to one.
fn probability(text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(value) if (0.0..=1.0).contains(&value) => Ok(value),
        _ => Err(format!("expected a probability from 0 to 1, found: {text}")),
    }
}

/// The argument IDs of each filter's options.
fn filter_options(kind: FilterKind) -> &'static [&'static str] {
    match kind {
        FilterKind::CopiedDamage => &[
            "copied_damage_classes",
            "copied_damage_scale",
            "copied_damage_threshold",
        ],
        FilterKind::ATailing => &["a_tailing_distance", "a_tailing_threshold"],
        FilterKind::EndRepairFillIn => &[
            "end_repair_fill_in_distance",
            "end_repair_fill_in_scale",
            "end_repair_fill_in_threshold",
        ],
    }
}

impl Cli {
    /// Validate the options and gather them into [`FilterArgs`], rejecting an
    /// option given on the command line for a filter `--filters` leaves out.
    fn into_args(self, matches: &ArgMatches) -> Result<FilterArgs> {
        for kind in FilterKind::ALL {
            if self.filters.contains(&kind) {
                continue;
            }
            for id in filter_options(kind) {
                if matches.value_source(id) == Some(ValueSource::CommandLine) {
                    bail!(
                        "--{} applies only to the {kind} filter, which --filters leaves out",
                        id.replace('_', "-")
                    );
                }
            }
        }
        validate_classes(&self.copied_damage_classes)?;
        let options = FilterOptions {
            sample: self.sample,
            filters: self.filters,
            prior: self.prior,
            copied_damage: CopiedDamage {
                classes: self.copied_damage_classes,
                scale: self.copied_damage_scale,
            },
            copied_damage_threshold: self.copied_damage_threshold,
            end_repair_fill_in: EndRepairFillIn {
                distance: self.end_repair_fill_in_distance,
                scale: self.end_repair_fill_in_scale,
            },
            end_repair_fill_in_threshold: self.end_repair_fill_in_threshold,
            a_tailing: ATailing {
                distance: self.a_tailing_distance,
            },
            a_tailing_threshold: self.a_tailing_threshold,
        };
        Ok(FilterArgs {
            input: self.input,
            output: self.output,
            bam: self.bam,
            reference: self.reference,
            metrics: self.metrics,
            read_filter: ReadFilter {
                min_mapping_quality: self.min_mapping_quality,
                min_base_quality: self.min_base_quality,
                paired_reads_only: self.paired_reads_only,
            },
            options,
        })
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
/// - indented single-column examples (a bare command): all [`CODE`];
/// - prose: default color, with backtick terms painted in [`CODE`].
fn style_help_text(text: &str, color: bool) -> String {
    let code = esc(CODE, color);
    let faded = esc(FADED, color);
    let reset = esc_reset(CODE, color);
    text.split_inclusive('\n')
        .map(|line| {
            let (content, newline) = match line.strip_suffix('\n') {
                Some(content) => (content, "\n"),
                None => (line, ""),
            };
            if !content.starts_with("  ") || content.trim().is_empty() {
                return format!("{}{newline}", paint_backtick_terms(content, &code, &reset));
            }
            if let Some((left, gap, right)) = split_two_column(content) {
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

    let env = Env::default().default_filter_or("info");
    let write_style = if color {
        env_logger::WriteStyle::Auto
    } else {
        env_logger::WriteStyle::Never
    };
    env_logger::Builder::from_env(env)
        .write_style(write_style)
        .init();

    // Help text is hand-wrapped, and clap would count the injected ANSI escapes as width.
    let cmd = decorate_help(Cli::command().term_width(usize::MAX), color);
    let matches = cmd.get_matches();
    let cli = Cli::from_arg_matches(&matches).unwrap_or_else(|e| e.exit());

    match cli.into_args(&matches).and_then(|args| run_filter(&args)) {
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

    fn args(extra: &[&str]) -> Result<FilterArgs> {
        let base = ["chaff", "-i", "in.vcf", "-o", "out.vcf", "-b", "in.bam"];
        let matches = Cli::command().try_get_matches_from(base.iter().chain(extra))?;
        Cli::from_arg_matches(&matches)?.into_args(&matches)
    }

    #[rstest]
    #[case(&["--filters", "a-tailing", "--copied-damage-threshold", "0.05"], "--copied-damage-threshold applies only to the copied-damage filter")]
    #[case(&["--filters", "a-tailing", "--copied-damage-classes", "C>T"], "--copied-damage-classes applies only to the copied-damage filter")]
    #[case(&["--filters", "copied-damage", "--a-tailing-distance", "2"], "--a-tailing-distance applies only to the a-tailing filter")]
    #[case(&["--filters", "copied-damage", "--a-tailing-p-value", "0.01"], "--a-tailing-threshold applies only to the a-tailing filter")]
    #[case(&["--filters", "a-tailing", "--end-repair-fill-in-scale", "15"], "--end-repair-fill-in-scale applies only to the end-repair-fill-in filter")]
    fn test_an_option_of_a_filter_left_out_is_an_error(
        #[case] extra: &[&str],
        #[case] message: &str,
    ) {
        let error = args(extra).unwrap_err();
        assert!(error.to_string().contains(message), "{error}");
    }

    #[test]
    fn test_filter_options_are_every_option_named_for_a_filter() {
        let ids: Vec<String> = Cli::command()
            .get_arguments()
            .map(|arg| arg.get_id().to_string())
            .collect();
        for kind in FilterKind::ALL {
            let prefix = format!("{}_", kind.to_string().replace('-', "_"));
            let named: Vec<&str> = ids
                .iter()
                .map(String::as_str)
                .filter(|id| id.starts_with(&prefix))
                .collect();
            assert_eq!(named, filter_options(kind), "{kind}");
        }
    }

    #[rstest]
    #[case(&[])]
    #[case(&["--filters", "a-tailing"])]
    #[case(&["--filters", "a-tailing", "--a-tailing-distance", "4", "--a-tailing-threshold", "0.001"])]
    #[case(&["--filters", "end-repair-fill-in", "--end-repair-fill-in-distance", "10", "--end-repair-fill-in-scale", "15"])]
    fn test_options_of_enabled_filters_and_defaults_are_accepted(#[case] extra: &[&str]) {
        args(extra).unwrap();
    }
}
