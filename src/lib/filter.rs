//! Score each somatic call against library-preparation artifacts, learn the
//! artifact priors, and write the annotated calls.
//!
//! The VCF and the reads are merge-joined in coordinate order: each call asks
//! the [`Evidence`] for the molecules at its position, so neither file needs an
//! index. Learning a prior needs every call's likelihood ratio before any
//! posterior is known, so the VCF is read twice: once to score the calls and
//! once to write them.

use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context as _, Result};
use clap::ValueEnum;
use log::info;
use noodles::core::Position;
use noodles::vcf;
use noodles::vcf::header::record::value::map::info::{Number, Type};
use noodles::vcf::variant::record_buf::info::field::value::Array;
use noodles::vcf::variant::record_buf::info::field::Value;
use noodles::vcf::variant::record_buf::Filters;
use noodles::vcf::variant::RecordBuf;
use streampile::{RecordSource, StreamingPileupBuilder};

use crate::call::Genotype;
use crate::classes::{sbs6, Context};
use crate::copied_damage::{CopiedDamage, DamageSite};
use crate::evidence::{Evidence, Molecule, PileupEvidence, PileupOptions};
use crate::io::{add_filter, add_info, vcf_float, VariantReader, VariantWriter};
use crate::metrics::{write_metrics, StratumMetrics};
use crate::prior::{
    fgbio_artifact_prior, learn_artifact_fraction, posterior_mutation, PriorMode, PSEUDOCOUNT,
};
use crate::read_end::{is_filtered, ATailing, EndRepairFillIn, Score};
use crate::reference::Reference;

/// One of the artifact filters.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, ValueEnum)]
pub enum FilterKind {
    /// Damage copied onto the other strand before strand tagging.
    CopiedDamage,
    /// Adenines added to an over-digested recessed 3' end during A-tailing.
    ATailing,
    /// Damage copied into a filled-in recessed 3' end during end repair.
    EndRepairFillIn,
}

impl FilterKind {
    /// Every filter, in output order.
    pub const ALL: [FilterKind; 3] = [
        FilterKind::CopiedDamage,
        FilterKind::ATailing,
        FilterKind::EndRepairFillIn,
    ];

    /// The FILTER the filter applies.
    pub fn filter_id(self) -> &'static str {
        match self {
            FilterKind::CopiedDamage => CopiedDamage::FILTER,
            FilterKind::ATailing => ATailing::FILTER,
            FilterKind::EndRepairFillIn => EndRepairFillIn::FILTER,
        }
    }

    /// The INFO key of the posterior probability of a true mutation.
    pub fn posterior_id(self) -> &'static str {
        match self {
            FilterKind::CopiedDamage => CopiedDamage::INFO_POSTERIOR,
            FilterKind::ATailing => ATailing::INFO,
            FilterKind::EndRepairFillIn => EndRepairFillIn::INFO,
        }
    }

    /// Whether the filter scores a genotype.
    pub fn applies_to(self, gt: &Genotype) -> bool {
        match self {
            FilterKind::CopiedDamage => CopiedDamage::applies_to(gt),
            FilterKind::ATailing => ATailing::applies_to(gt),
            FilterKind::EndRepairFillIn => EndRepairFillIn::applies_to(gt),
        }
    }

    /// The posterior at or below which the filter's FILTER is applied.
    pub fn threshold(self, options: &FilterOptions) -> Option<f64> {
        match self {
            FilterKind::CopiedDamage => options.copied_damage_threshold,
            FilterKind::ATailing => options.a_tailing_threshold,
            FilterKind::EndRepairFillIn => options.end_repair_fill_in_threshold,
        }
    }

    /// The IDs of the command-line arguments that only this filter reads.
    pub fn arguments(self) -> &'static [&'static str] {
        match self {
            FilterKind::CopiedDamage => &[
                "reference",
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
}

impl fmt::Display for FilterKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FilterKind::CopiedDamage => write!(f, "copied-damage"),
            FilterKind::EndRepairFillIn => write!(f, "end-repair-fill-in"),
            FilterKind::ATailing => write!(f, "a-tailing"),
        }
    }
}

/// How calls are scored and filtered.
#[derive(Clone, Debug, PartialEq)]
pub struct FilterOptions {
    /// The sample whose reads are in the BAM; required with several samples.
    pub sample: Option<String>,
    /// The filters to run.
    pub filters: Vec<FilterKind>,
    /// The prior that turns likelihood ratios into posteriors.
    pub prior: PriorMode,
    /// The copied damage model.
    pub copied_damage: CopiedDamage,
    /// Apply `CopiedDamageArtifact` at or below this posterior.
    pub copied_damage_threshold: Option<f64>,
    /// The end repair fill-in model.
    pub end_repair_fill_in: EndRepairFillIn,
    /// Apply `EndRepairFillInArtifact` at or below this posterior.
    pub end_repair_fill_in_threshold: Option<f64>,
    /// The A-tailing model.
    pub a_tailing: ATailing,
    /// Apply `ATailingArtifact` at or below this posterior.
    pub a_tailing_threshold: Option<f64>,
}

impl Default for FilterOptions {
    fn default() -> Self {
        Self {
            sample: None,
            filters: FilterKind::ALL.to_vec(),
            prior: PriorMode::Learned,
            copied_damage: CopiedDamage::default(),
            copied_damage_threshold: None,
            end_repair_fill_in: EndRepairFillIn::default(),
            end_repair_fill_in_threshold: None,
            a_tailing: ATailing::default(),
            a_tailing_threshold: None,
        }
    }
}

impl FilterOptions {
    fn enabled(&self, kind: FilterKind) -> bool {
        self.filters.contains(&kind)
    }

    fn prior_text(&self) -> &'static str {
        match self.prior {
            PriorMode::Learned => "an artifact prior learned per sample and stratum",
            PriorMode::Fgbio => "fgbio's (2 * maf)^2 mutation prior",
        }
    }

    fn threshold_text(&self, kind: FilterKind) -> String {
        match kind.threshold(self) {
            Some(t) => format!("at or below a posterior of {t}"),
            None => "never applied without a threshold".to_string(),
        }
    }

    /// Add the INFO and FILTER lines of every enabled filter to a header.
    pub fn add_header_lines(&self, header: &mut vcf::Header) {
        let prior = self.prior_text();
        if self.enabled(FilterKind::CopiedDamage) {
            let classes: Vec<String> = self
                .copied_damage
                .classes
                .iter()
                .map(ToString::to_string)
                .collect();
            let model = format!(
                "damage classes {} and a {} bp copy scale",
                classes.join(","),
                self.copied_damage.scale
            );
            add_info(
                header,
                CopiedDamage::INFO_POSTERIOR,
                Number::Count(1),
                Type::Float,
                &format!("Posterior probability that the call is a true mutation rather than damage copied onto both strands, with {model} and {prior}."),
            );
            add_info(
                header,
                CopiedDamage::INFO_RATIO,
                Number::Count(1),
                Type::Float,
                "Log10 likelihood ratio of the copied damage artifact to a true mutation.",
            );
            add_info(
                header,
                CopiedDamage::INFO_ALT,
                Number::Count(2),
                Type::Integer,
                "Alternate molecules nearer the lesion strand's 5' end than its 3' end, and all alternate molecules measured.",
            );
            add_info(
                header,
                CopiedDamage::INFO_REF,
                Number::Count(2),
                Type::Integer,
                "Reference molecules nearer the lesion strand's 5' end than its 3' end, and all reference molecules measured.",
            );
            add_filter(
                header,
                CopiedDamage::FILTER,
                &format!(
                    "Call is likely damage copied onto both strands, with {model}, {}.",
                    self.threshold_text(FilterKind::CopiedDamage)
                ),
            );
        }
        if self.enabled(FilterKind::ATailing) {
            let distance = self.a_tailing.distance;
            add_info(
                header,
                ATailing::INFO,
                Number::Count(1),
                Type::Float,
                &format!("Posterior probability that the call is a true mutation rather than an A-tailing artifact, with a {distance} bp distance from the template end and {prior}."),
            );
            add_filter(
                header,
                ATailing::FILTER,
                &format!(
                    "Call is likely an A-tailing artifact, with a {distance} bp distance from the template end, {}.",
                    self.threshold_text(FilterKind::ATailing)
                ),
            );
        }
        if self.enabled(FilterKind::EndRepairFillIn) {
            let model = match self.end_repair_fill_in.scale {
                Some(scale) => format!("a {scale} bp decay from the nearest template end"),
                None => format!(
                    "a {} bp distance from the nearest template end",
                    self.end_repair_fill_in.distance
                ),
            };
            add_info(
                header,
                EndRepairFillIn::INFO,
                Number::Count(1),
                Type::Float,
                &format!("Posterior probability that the call is a true mutation rather than an end repair fill-in artifact, with {model} and {prior}."),
            );
            add_filter(
                header,
                EndRepairFillIn::FILTER,
                &format!(
                    "Call is likely an end repair fill-in artifact, with {model}, {}.",
                    self.threshold_text(FilterKind::EndRepairFillIn)
                ),
            );
        }
    }
}

/// The inputs and outputs of one run.
#[derive(Clone, Debug)]
pub struct FilterArgs {
    /// The somatic VCF or BCF, coordinate-sorted.
    pub input: PathBuf,
    /// The output VCF or BCF, or `-` for standard output.
    pub output: PathBuf,
    /// The coordinate-sorted BAM of the sample.
    pub bam: PathBuf,
    /// The indexed reference FASTA, needed by the copied damage filter.
    pub reference: Option<PathBuf>,
    /// The per-sample metrics TSV.
    pub metrics: Option<PathBuf>,
    /// The read and base floors.
    pub pileup: PileupOptions,
    /// The scoring and filtering options.
    pub options: FilterOptions,
}

/// One filter's score of one call.
#[derive(Clone, Debug, PartialEq)]
pub struct Annotation {
    /// The filter.
    pub kind: FilterKind,
    /// The stratum the artifact fraction is learned in.
    pub stratum: String,
    /// The likelihood ratio and molecule counts.
    pub score: Score,
    /// fgbio's per-call artifact prior.
    pub fgbio_prior: f64,
    /// The posterior probability of a true mutation, once known.
    pub posterior: Option<f64>,
}

/// The index of the sample under test, resolved as fgbio does: the named
/// sample, or the only sample when none is named.
pub fn resolve_sample(header: &vcf::Header, sample: Option<&str>) -> Result<(usize, String)> {
    let names = header.sample_names();
    match sample {
        Some(name) => match names.get_index_of(name) {
            Some(i) => Ok((i, name.to_string())),
            None => bail!("there is no genotype with the following sample in the input VCF/BCF: {name}"),
        },
        None => match names.len() {
            1 => Ok((0, names[0].clone())),
            0 => bail!("the input VCF/BCF has no samples"),
            n => bail!("the input VCF/BCF has {n} samples, so --sample must name the one whose reads are in the BAM"),
        },
    }
}

/// Rejects records that step backwards in coordinate order.
struct CoordinateOrder {
    ranks: HashMap<String, usize>,
    last: Option<(usize, usize, String)>,
}

impl CoordinateOrder {
    fn new(header: &vcf::Header) -> Self {
        let ranks = header
            .contigs()
            .keys()
            .enumerate()
            .map(|(i, name)| (name.clone(), i))
            .collect();
        Self { ranks, last: None }
    }

    fn check(&mut self, contig: &str, pos: usize) -> Result<()> {
        let next_rank = self.ranks.len();
        let rank = *self.ranks.entry(contig.to_string()).or_insert(next_rank);
        if let Some((last_rank, last_pos, last_contig)) = &self.last {
            if (rank, pos) < (*last_rank, *last_pos) {
                bail!("the input VCF/BCF is not coordinate sorted: {contig}:{pos} follows {last_contig}:{last_pos}");
            }
        }
        self.last = Some((rank, pos, contig.to_string()));
        Ok(())
    }
}

/// Score one call under every enabled filter that applies to it.
fn score_call(
    gt: &Genotype,
    contig: &str,
    pos: usize,
    options: &FilterOptions,
    evidence: &mut dyn Evidence,
    reference: &mut Option<&mut Reference>,
    cache: &mut Option<((String, usize), Vec<Molecule>)>,
) -> Result<Vec<Annotation>> {
    let kinds: Vec<FilterKind> = FilterKind::ALL
        .into_iter()
        .filter(|k| options.enabled(*k) && k.applies_to(gt))
        .collect();
    let (Some(ref_allele), Some(alt_allele)) = (gt.reference().bytes().next(), gt.first_alt())
    else {
        return Ok(Vec::new());
    };
    let Some(alt_allele) = alt_allele.bytes().next() else {
        return Ok(Vec::new());
    };
    if kinds.is_empty() {
        return Ok(Vec::new());
    }
    let (ref_base, alt_base) = (
        ref_allele.to_ascii_uppercase(),
        alt_allele.to_ascii_uppercase(),
    );
    let key = (contig.to_string(), pos);
    if cache.as_ref().map(|(k, _)| k) != Some(&key) {
        let position = Position::try_from(pos).context("a VCF position must be at least 1")?;
        let molecules = evidence
            .molecules(contig, position)
            .with_context(|| format!("failed to read molecules at {contig}:{pos}"))?;
        *cache = Some((key, molecules));
    }
    let molecules = &cache.as_ref().expect("filled above").1;
    let count = |base: u8| molecules.iter().filter(|m| m.base == base).count() as u32;
    let fgbio_prior =
        fgbio_artifact_prior(count(alt_base), count(ref_base), molecules.len() as u32);
    let substitution = sbs6(ref_base, alt_base).unwrap_or_else(|| "other".to_string());

    let mut annotations = Vec::with_capacity(kinds.len());
    for kind in kinds {
        let scored = match kind {
            FilterKind::EndRepairFillIn => Some((
                substitution.clone(),
                options
                    .end_repair_fill_in
                    .score(molecules, ref_base, alt_base),
            )),
            FilterKind::ATailing => Some((
                substitution.clone(),
                options.a_tailing.score(molecules, ref_base, alt_base),
            )),
            FilterKind::CopiedDamage => match (
                options.copied_damage.classify(ref_base, alt_base),
                reference.as_deref_mut(),
            ) {
                (None, _) | (_, None) => None,
                (Some((class, strand)), Some(reference)) => {
                    let (prev, base, next) = reference.context(contig, pos)?;
                    let damage = DamageSite {
                        class,
                        strand,
                        context: Context::of(prev, base, next),
                    };
                    let score = options
                        .copied_damage
                        .score(molecules, ref_base, alt_base, strand);
                    Some((damage.stratum(), score))
                }
            },
        };
        if let Some((stratum, score)) = scored {
            annotations.push(Annotation {
                kind,
                stratum,
                score,
                fgbio_prior,
                posterior: None,
            });
        }
    }
    Ok(annotations)
}

/// Learn the priors and fill in every annotation's posterior, returning the
/// learned artifact fraction of each filter and stratum.
fn assign_posteriors(
    calls: &mut [Vec<Annotation>],
    prior: PriorMode,
) -> BTreeMap<(FilterKind, String), f64> {
    let mut ratios: BTreeMap<(FilterKind, String), Vec<f64>> = BTreeMap::new();
    for annotation in calls.iter().flatten() {
        if let Some(llr) = annotation.score.log_likelihood_ratio {
            ratios
                .entry((annotation.kind, annotation.stratum.clone()))
                .or_default()
                .push(llr);
        }
    }
    let fractions: BTreeMap<(FilterKind, String), f64> = ratios
        .into_iter()
        .map(|(key, llrs)| (key, learn_artifact_fraction(&llrs, PSEUDOCOUNT)))
        .collect();
    for annotation in calls.iter_mut().flatten() {
        let Some(llr) = annotation.score.log_likelihood_ratio else {
            continue;
        };
        let artifact_prior = match prior {
            PriorMode::Learned => fractions[&(annotation.kind, annotation.stratum.clone())],
            PriorMode::Fgbio => annotation.fgbio_prior,
        };
        annotation.posterior = Some(posterior_mutation(llr, artifact_prior));
    }
    fractions
}

/// Write one call's annotations into its record.
fn annotate_record(record: &mut RecordBuf, annotations: &[Annotation], options: &FilterOptions) {
    let mut new_filters = Vec::new();
    for annotation in annotations {
        let (kind, score) = (annotation.kind, &annotation.score);
        let info = record.info_mut();
        if let (Some(posterior), Some(llr)) = (annotation.posterior, score.log_likelihood_ratio) {
            info.insert(
                kind.posterior_id().to_string(),
                Some(Value::Float(vcf_float(posterior))),
            );
            if kind == FilterKind::CopiedDamage {
                info.insert(
                    CopiedDamage::INFO_RATIO.to_string(),
                    Some(Value::Float(vcf_float(llr / std::f64::consts::LN_10))),
                );
            }
            if is_filtered(posterior, kind.threshold(options)) {
                new_filters.push(kind.filter_id().to_string());
            }
        }
        if kind == FilterKind::CopiedDamage {
            let pair = |a: u32, b: u32| {
                Some(Value::Array(Array::Integer(vec![
                    Some(a as i32),
                    Some(b as i32),
                ])))
            };
            info.insert(
                CopiedDamage::INFO_ALT.to_string(),
                pair(score.alt_congruent, score.alt_molecules),
            );
            info.insert(
                CopiedDamage::INFO_REF.to_string(),
                pair(score.ref_congruent, score.ref_molecules),
            );
        }
    }
    if new_filters.is_empty() {
        return;
    }
    let filters = record.filters_mut();
    if filters.is_pass() {
        *filters = Filters::default();
    }
    filters.as_mut().extend(new_filters);
}

/// Score, annotate, and filter the calls of `input` into `output` with
/// molecules from `evidence`, returning the per-stratum metrics.
pub fn filter_vcf(
    input: &Path,
    output: &Path,
    evidence: &mut dyn Evidence,
    mut reference: Option<&mut Reference>,
    options: &FilterOptions,
) -> Result<Vec<StratumMetrics>> {
    if options.enabled(FilterKind::CopiedDamage) && reference.is_none() {
        bail!("the copied damage filter needs a reference FASTA (--ref)");
    }
    if options.filters.is_empty() {
        log::warn!("every filter is disabled, so chaff will copy the input unchanged");
    }

    let mut reader = VariantReader::open(input)?;
    let header = reader
        .read_header()
        .context("failed to read the VCF/BCF header")?;
    let (sample_index, sample) = resolve_sample(&header, options.sample.as_deref())?;
    let mut order = CoordinateOrder::new(&header);
    let mut cache = None;
    let mut calls: Vec<Vec<Annotation>> = Vec::new();
    let mut record = RecordBuf::default();
    while reader.read_record(&header, &mut record)? != 0 {
        let contig = record.reference_sequence_name().to_string();
        let pos = record.variant_start().map(usize::from).unwrap_or(0);
        order.check(&contig, pos)?;
        let annotations = match Genotype::from_record(&record, sample_index) {
            Some(gt) if pos > 0 => score_call(
                &gt,
                &contig,
                pos,
                options,
                evidence,
                &mut reference,
                &mut cache,
            )?,
            _ => Vec::new(),
        };
        calls.push(annotations);
    }
    info!("scored {} calls of sample {sample}", calls.len());

    let fractions = assign_posteriors(&mut calls, options.prior);

    let mut out_header = header.clone();
    options.add_header_lines(&mut out_header);
    let mut writer = VariantWriter::create(output)?;
    writer.write_header(&out_header)?;
    let mut reader = VariantReader::open(input)?;
    let _ = reader.read_header()?;
    let mut index = 0;
    while reader.read_record(&header, &mut record)? != 0 {
        let annotations = calls
            .get(index)
            .context("the input VCF/BCF changed while it was read")?;
        annotate_record(&mut record, annotations, options);
        writer.write_record(&out_header, &record)?;
        index += 1;
    }
    writer.finish()?;

    let rows = metrics_rows(&sample, &calls, &fractions, options);
    for row in &rows {
        info!(
            "{} {}: {} calls, artifact fraction {}, {} filtered, {} of {} alternate molecules congruent",
            row.filter,
            row.stratum,
            row.calls,
            row.artifact_fraction
                .map_or_else(|| "per call".to_string(), |f| format!("{f:.4}")),
            row.filtered,
            row.alt_congruent,
            row.alt_molecules,
        );
    }
    Ok(rows)
}

/// Pool every annotation into one metrics row per filter and stratum.
fn metrics_rows(
    sample: &str,
    calls: &[Vec<Annotation>],
    fractions: &BTreeMap<(FilterKind, String), f64>,
    options: &FilterOptions,
) -> Vec<StratumMetrics> {
    let mut rows: BTreeMap<(FilterKind, String), StratumMetrics> = BTreeMap::new();
    for annotation in calls.iter().flatten() {
        let key = (annotation.kind, annotation.stratum.clone());
        let row = rows.entry(key.clone()).or_insert_with(|| StratumMetrics {
            sample: sample.to_string(),
            filter: annotation.kind.to_string(),
            stratum: annotation.stratum.clone(),
            calls: 0,
            filtered: 0,
            artifact_fraction: match options.prior {
                PriorMode::Learned => fractions.get(&key).copied(),
                PriorMode::Fgbio => None,
            },
            expected_artifacts: 0.0,
            alt_molecules: 0,
            alt_congruent: 0,
            alt_congruent_fraction: None,
            ref_molecules: 0,
            ref_congruent: 0,
            ref_congruent_fraction: None,
            asymmetry_p_value: None,
        });
        row.calls += 1;
        let score = &annotation.score;
        row.alt_molecules += u64::from(score.alt_molecules);
        row.alt_congruent += u64::from(score.alt_congruent);
        row.ref_molecules += u64::from(score.ref_molecules);
        row.ref_congruent += u64::from(score.ref_congruent);
        if let Some(posterior) = annotation.posterior {
            row.expected_artifacts += 1.0 - posterior;
            if is_filtered(posterior, annotation.kind.threshold(options)) {
                row.filtered += 1;
            }
        }
    }
    rows.into_values().map(StratumMetrics::finish).collect()
}

/// Filter the calls with molecules from `evidence`, writing the metrics when
/// asked.
pub fn run_filter_with(args: &FilterArgs, evidence: &mut dyn Evidence) -> Result<()> {
    let mut reference = args.reference.as_deref().map(Reference::open).transpose()?;
    let rows = filter_vcf(
        &args.input,
        &args.output,
        evidence,
        reference.as_mut(),
        &args.options,
    )?;
    if let Some(path) = &args.metrics {
        write_metrics(path, &rows)?;
    }
    Ok(())
}

/// Filter the calls with molecules piled up by `builder` from the BAM's
/// records, under the read and base floors of `args`.
pub fn run_filter_on<S: RecordSource>(
    args: &FilterArgs,
    builder: StreamingPileupBuilder<'_, S>,
) -> Result<()> {
    let mut evidence = PileupEvidence::new(builder, &args.pileup);
    run_filter_with(args, &mut evidence)
}

/// Filter the calls with the BAM named by `args`, streamed once through
/// streampile.
pub fn run_filter(args: &FilterArgs) -> Result<()> {
    let mut reader = noodles::bam::io::reader::Builder
        .build_from_path(&args.bam)
        .with_context(|| format!("failed to open BAM: {:?}", args.bam))?;
    let header = reader
        .read_header()
        .context("failed to read the BAM header")?;
    let builder = StreamingPileupBuilder::new(reader, &header).map_err(|error| match error {
        streampile::Error::NotCoordinateSorted { .. } => anyhow!(
            "the BAM must be coordinate sorted (@HD SO:coordinate): {:?}",
            args.bam
        ),
        error => error.into(),
    })?;
    run_filter_on(args, builder)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::MoleculeTable;
    use crate::testing::{gt, write_fasta, Variant, VcfBuilder};

    fn read_records(path: &Path) -> (vcf::Header, Vec<RecordBuf>) {
        let mut reader = VariantReader::open(path).unwrap();
        let header = reader.read_header().unwrap();
        let mut records = Vec::new();
        let mut record = RecordBuf::default();
        while reader.read_record(&header, &mut record).unwrap() != 0 {
            records.push(record.clone());
        }
        (header, records)
    }

    fn float(record: &RecordBuf, key: &str) -> Option<f32> {
        match record.info().get(key) {
            Some(Some(Value::Float(f))) => Some(*f),
            _ => None,
        }
    }

    #[test]
    fn test_resolve_sample() {
        let dir = tempfile::tempdir().unwrap();
        let path = VcfBuilder::new(&["tumor", "normal"]).write(&dir.path().join("a.vcf"));
        let mut reader = VariantReader::open(&path).unwrap();
        let header = reader.read_header().unwrap();
        assert_eq!(resolve_sample(&header, Some("normal")).unwrap().0, 1);
        assert!(resolve_sample(&header, None).is_err());
        assert!(resolve_sample(&header, Some("WhoDis")).is_err());
    }

    #[test]
    fn test_coordinate_order() {
        let mut header = vcf::Header::default();
        for name in ["chr1", "chr2"] {
            header.contigs_mut().insert(name.into(), Default::default());
        }
        let mut order = CoordinateOrder::new(&header);
        order.check("chr1", 5).unwrap();
        order.check("chr1", 5).unwrap();
        order.check("chr2", 1).unwrap();
        assert!(order.check("chr1", 9).is_err());
        let mut order = CoordinateOrder::new(&vcf::Header::default());
        order.check("chrB", 5).unwrap();
        order.check("chrA", 1).unwrap();
        assert!(order.check("chrB", 6).is_err());
    }

    #[test]
    fn test_copied_damage_end_to_end_on_a_molecule_table() {
        let dir = tempfile::tempdir().unwrap();
        let reference = write_fasta(dir.path(), "chr1", &"ACGTTCAA".repeat(250));
        let mut vcf = VcfBuilder::new(&["tumor"]);
        vcf.add(Variant::new(1002, &["C", "T"], vec![gt("tumor", "0/1")]));
        vcf.add(Variant::new(1006, &["C", "T"], vec![gt("tumor", "0/1")]));
        let input = vcf.write(&dir.path().join("in.vcf"));
        let output = dir.path().join("out.vcf");

        let mut table = MoleculeTable::new();
        let at = |base, d| Molecule::new(base, 90, d, 149 - d);
        for pos in [1002, 1006] {
            let mut molecules: Vec<Molecule> = (0..150).map(|d| at(b'C', d)).collect();
            if pos == 1002 {
                molecules.extend((0..6).map(|d| at(b'T', d)));
            } else {
                molecules.extend((0..150).step_by(25).map(|d| at(b'T', d)));
            }
            table.insert("chr1", pos, molecules);
        }
        let options = FilterOptions {
            filters: vec![FilterKind::CopiedDamage],
            copied_damage_threshold: Some(0.05),
            ..FilterOptions::default()
        };
        let mut reference = Reference::open(&reference).unwrap();
        let rows = filter_vcf(&input, &output, &mut table, Some(&mut reference), &options).unwrap();

        let (header, records) = read_records(&output);
        assert!(header.infos().contains_key(CopiedDamage::INFO_POSTERIOR));
        assert!(header.filters().contains_key(CopiedDamage::FILTER));
        assert!(!header.infos().contains_key(EndRepairFillIn::INFO));
        let artifact = &records[0];
        let mutation = &records[1];
        assert!(float(artifact, CopiedDamage::INFO_POSTERIOR).unwrap() < 0.05);
        assert!(artifact.filters().as_ref().contains(CopiedDamage::FILTER));
        assert!(float(mutation, CopiedDamage::INFO_POSTERIOR).unwrap() > 0.5);
        assert!(!mutation.filters().as_ref().contains(CopiedDamage::FILTER));
        assert_eq!(
            artifact.info().get(CopiedDamage::INFO_ALT),
            Some(Some(&Value::Array(Array::Integer(vec![Some(6), Some(6)]))))
        );

        assert_eq!(rows.len(), 2);
        let strata: Vec<&str> = rows.iter().map(|r| r.stratum.as_str()).collect();
        assert_eq!(strata, vec!["C>T:CpG", "C>T:non-CpG"]);
    }

    #[test]
    fn test_bcf_and_bgzf_outputs_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let mut vcf = VcfBuilder::new(&["tumor"]);
        vcf.add(Variant::new(10, &["G", "T"], vec![gt("tumor", "0/1")]));
        let input = vcf.write(&dir.path().join("in.vcf"));
        let mut table = MoleculeTable::new();
        let mut molecules: Vec<Molecule> =
            (0..20).map(|d| Molecule::new(b'G', 30, d, 90)).collect();
        molecules.push(Molecule::new(b'T', 30, 0, 90));
        table.insert("chr1", 10, molecules);
        let options = FilterOptions {
            filters: vec![FilterKind::ATailing, FilterKind::EndRepairFillIn],
            ..FilterOptions::default()
        };
        for name in ["out.bcf", "out.vcf.gz"] {
            let output = dir.path().join(name);
            filter_vcf(&input, &output, &mut table.clone(), None, &options).unwrap();
            let (header, records) = read_records(&output);
            assert!(header.infos().contains_key(EndRepairFillIn::INFO), "{name}");
            assert_eq!(records.len(), 1, "{name}");
            assert!(
                float(&records[0], EndRepairFillIn::INFO).is_some(),
                "{name}"
            );
            assert!(float(&records[0], ATailing::INFO).is_some(), "{name}");
        }
    }

    #[test]
    fn test_copied_damage_requires_a_reference() {
        let dir = tempfile::tempdir().unwrap();
        let input = VcfBuilder::new(&["tumor"]).write(&dir.path().join("in.vcf"));
        let error = filter_vcf(
            &input,
            &dir.path().join("out.vcf"),
            &mut MoleculeTable::new(),
            None,
            &FilterOptions::default(),
        )
        .unwrap_err();
        assert!(error.to_string().contains("reference"), "{error}");
    }
}
