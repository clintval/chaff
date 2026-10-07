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

use crate::call::{Genotype, Skip};
use crate::classes::{sbs6, Context};
use crate::copied_damage::{CopiedDamage, DamageSite};
use crate::evidence::{Evidence, Molecule, PileupEvidence, PileupOptions};
use crate::io::{add_filter, add_info, significant, vcf_float, VariantReader, VariantWriter};
use crate::metrics::{null_fraction, write_metrics, StratumMetrics};
use crate::model::{Distance, Model};
use crate::prior::{
    fgbio_artifact_prior, learn_artifact_fraction, learn_scale, posterior_mutation, BetaPrior,
    FILTER_PRIOR, STRATUM_PRIOR_STRENGTH,
};
use crate::read_end::{is_filtered, ATailing, Distances, EndRepairFillIn, ReferencePool, Score};
use crate::reference::Reference;
use crate::spectrum::{channel, After, Spectrum};

/// One of the artifact filters.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, ValueEnum)]
pub enum FilterKind {
    /// Damage copied onto the other strand before strand tagging.
    CopiedDamage,
    /// Adenines added to an over-digested recessed 3' end during A-tailing.
    ATailing,
    /// Errors in a recessed 3' end that end repair fills in.
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

    /// The IDs of every INFO field the filter writes.
    pub fn info_ids(self) -> &'static [&'static str] {
        match self {
            FilterKind::CopiedDamage => &[
                CopiedDamage::INFO_POSTERIOR,
                CopiedDamage::INFO_RATIO,
                CopiedDamage::INFO_ALT,
                CopiedDamage::INFO_REF,
            ],
            FilterKind::ATailing => &[ATailing::INFO],
            FilterKind::EndRepairFillIn => &[EndRepairFillIn::INFO],
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
                "copied_damage_distance",
                "copied_damage_threshold",
            ],
            FilterKind::ATailing => &["a_tailing_distance", "a_tailing_threshold"],
            FilterKind::EndRepairFillIn => &[
                "end_repair_fill_in_distance",
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
    /// The model: its prior and the shape of its distance models.
    pub model: Model,
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
            model: Model::Chaff,
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
        match self.model {
            Model::Chaff => "an artifact prior learned per sample and stratum",
            Model::Fgbio => "fgbio's (2 * maf)^2 mutation prior",
        }
    }

    fn threshold_text(&self, kind: FilterKind) -> String {
        match kind.threshold(self) {
            Some(t) => format!("at or below a posterior of {t}"),
            None => "never applied without a threshold".to_string(),
        }
    }

    /// Add the INFO and FILTER lines of every enabled filter to a header, with
    /// the distances the filters scored with.
    pub fn add_header_lines(&self, header: &mut vcf::Header, scales: &Scales) {
        let prior = self.prior_text();
        let decay = |distance: Distance, scale: f64, end: &str| {
            let learned = match distance {
                Distance::Learned => ", learned from the calls,",
                Distance::Bases(_) => "",
            };
            format!("a {} bp decay{learned} from {end}", significant(scale, 3))
        };
        if self.enabled(FilterKind::CopiedDamage) {
            let classes: Vec<String> = self
                .copied_damage
                .classes
                .iter()
                .map(ToString::to_string)
                .collect();
            let distance = significant(scales.copied_damage, 3);
            let model = format!(
                "damage classes {} and {}",
                classes.join(","),
                decay(
                    self.copied_damage.distance,
                    scales.copied_damage,
                    "the lesion strand's 5' end"
                ),
            );
            add_info(
                header,
                CopiedDamage::INFO_POSTERIOR,
                Number::Count(1),
                Type::Float,
                &format!("Posterior probability that the call is a real mutation rather than damage copied onto both strands, with {model} and {prior}."),
            );
            add_info(
                header,
                CopiedDamage::INFO_RATIO,
                Number::Count(1),
                Type::Float,
                "Log10 likelihood ratio of copied damage to a real mutation.",
            );
            add_info(
                header,
                CopiedDamage::INFO_ALT,
                Number::Count(2),
                Type::Integer,
                &format!("Alternate molecules within {distance} bp of the lesion strand's 5' end, and all alternate molecules measured."),
            );
            add_info(
                header,
                CopiedDamage::INFO_REF,
                Number::Count(2),
                Type::Integer,
                &format!("Reference molecules within {distance} bp of the lesion strand's 5' end, and all reference molecules measured."),
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
                &format!("Posterior probability that the call is a real mutation rather than an A-tailing artifact, with a {distance} bp window from the template end and {prior}."),
            );
            add_filter(
                header,
                ATailing::FILTER,
                &format!(
                    "Call is likely an A-tailing artifact, with a {distance} bp window from the template end, {}.",
                    self.threshold_text(FilterKind::ATailing)
                ),
            );
        }
        if self.enabled(FilterKind::EndRepairFillIn) {
            let model = match self.model {
                Model::Chaff => decay(
                    self.end_repair_fill_in.distance,
                    scales.end_repair_fill_in,
                    "the 3' end of the strand each template was copied from",
                ),
                Model::Fgbio => format!(
                    "a {} bp window from the nearest template end",
                    significant(scales.end_repair_fill_in, 3)
                ),
            };
            add_info(
                header,
                EndRepairFillIn::INFO,
                Number::Count(1),
                Type::Float,
                &format!("Posterior probability that the call is a real mutation rather than an end repair fill-in artifact, with {model} and {prior}."),
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
    /// The PDF of the sample's trinucleotide spectrum before and after
    /// filtering, which needs the reference FASTA.
    pub spectrum: Option<PathBuf>,
    /// The read and base floors.
    pub pileup: PileupOptions,
    /// The scoring and filtering options.
    pub options: FilterOptions,
}

/// A filter and one of its strata.
type Stratum = (FilterKind, String);

/// One filter's score of one call.
#[derive(Clone, Debug, PartialEq)]
pub struct Annotation {
    /// The filter.
    pub kind: FilterKind,
    /// The stratum the artifact fraction is learned in.
    pub stratum: String,
    /// The likelihood ratio and molecule counts.
    pub score: Score,
    /// The distances a decay scores the call by once its scale is known.
    pub distances: Option<Distances>,
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

/// Refuse a header that already declares an INFO or FILTER of an enabled
/// filter, from an earlier run of chaff or fgbio, whose values would otherwise
/// sit under this run's descriptions.
fn refuse_earlier_annotations(header: &vcf::Header, options: &FilterOptions) -> Result<()> {
    let mut found = Vec::new();
    for kind in options.filters.iter() {
        for id in kind.info_ids() {
            if header.infos().contains_key(*id) {
                found.push(format!("INFO/{id}"));
            }
        }
        if header.filters().contains_key(kind.filter_id()) {
            found.push(format!("FILTER/{}", kind.filter_id()));
        }
    }
    if !found.is_empty() {
        bail!(
            "the input VCF/BCF already has {}, from an earlier run; remove them first, e.g. with bcftools annotate -x {}",
            found.join(", ").replace('/', " "),
            found.join(",")
        );
    }
    Ok(())
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

/// The upper-cased REF and first alternate base of a heterozygous call whose
/// every called allele is one base: an SNV the filters can score.
fn snv(gt: &Genotype) -> Option<(u8, u8)> {
    if !(gt.is_het() && gt.calls_are_single_bases()) {
        return None;
    }
    let ref_base = gt.reference().bytes().next()?;
    let alt_base = gt.first_alt()?.bytes().next()?;
    Some((ref_base.to_ascii_uppercase(), alt_base.to_ascii_uppercase()))
}

/// The reference bases before, at, and after an SNV, refusing a REF that
/// disagrees with the FASTA.
fn snv_context(
    reference: &mut Reference,
    gt: &Genotype,
    contig: &str,
    pos: usize,
    ref_base: u8,
) -> Result<(Option<u8>, u8, Option<u8>)> {
    let context = reference.context(contig, pos)?;
    if context.1 != ref_base {
        bail!(
            "the call at {contig}:{pos} has REF {}, but the reference FASTA has {} there",
            gt.reference(),
            context.1 as char
        );
    }
    Ok(context)
}

/// Score one call, with its reference context when known, under every enabled
/// filter that applies to it.
fn score_call(
    gt: &Genotype,
    contig: &str,
    pos: usize,
    options: &FilterOptions,
    evidence: &mut dyn Evidence,
    context: Option<(Option<u8>, u8, Option<u8>)>,
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
            FilterKind::EndRepairFillIn => {
                let filter = &options.end_repair_fill_in;
                Some(match options.model {
                    Model::Chaff => (
                        substitution.clone(),
                        Score::default(),
                        Some(filter.distances(molecules, ref_base, alt_base)),
                    ),
                    Model::Fgbio => (
                        substitution.clone(),
                        filter.score(molecules, ref_base, alt_base),
                        None,
                    ),
                })
            }
            FilterKind::ATailing => {
                let mut score = options.a_tailing.score(molecules, ref_base, alt_base);
                if options.model == Model::Chaff && score.alt_molecules == 0 {
                    score.log_likelihood_ratio = None;
                }
                Some((substitution.clone(), score, None))
            }
            FilterKind::CopiedDamage => {
                match (options.copied_damage.classify(ref_base, alt_base), context) {
                    (Some((class, strand)), Some((prev, base, next))) => {
                        let damage = DamageSite {
                            class,
                            strand,
                            context: Context::of(prev, base, next),
                        };
                        let distances = options
                            .copied_damage
                            .distances(molecules, ref_base, alt_base, strand);
                        Some((damage.stratum(), Score::default(), Some(distances)))
                    }
                    _ => None,
                }
            }
        };
        if let Some((stratum, score, distances)) = scored {
            annotations.push(Annotation {
                kind,
                stratum,
                score,
                distances,
                fgbio_prior,
                posterior: None,
            });
        }
    }
    Ok(annotations)
}

/// The distance each filter scored with, in bases: a decay's scale, learned
/// or fixed, or a window.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Scales {
    /// The copied damage decay's scale.
    pub copied_damage: f64,
    /// The end repair fill-in decay's scale under the `chaff` model, or its
    /// window under `fgbio`.
    pub end_repair_fill_in: f64,
    /// The A-tailing window.
    pub a_tailing: f64,
}

impl Scales {
    /// The distance `kind` scored with.
    pub fn of(&self, kind: FilterKind) -> f64 {
        match kind {
            FilterKind::CopiedDamage => self.copied_damage,
            FilterKind::EndRepairFillIn => self.end_repair_fill_in,
            FilterKind::ATailing => self.a_tailing,
        }
    }
}

/// Learn or fix each decay's scale from all of its filter's calls, and score
/// every call the decay holds distances for at that scale, each call's
/// reference distances shrunk toward its stratum's pooled ones.
fn score_decays(calls: &mut [Vec<Annotation>], options: &FilterOptions) -> Scales {
    let mut scales = Scales {
        copied_damage: options
            .copied_damage
            .distance
            .bases(CopiedDamage::FALLBACK_SCALE),
        end_repair_fill_in: options.end_repair_fill_in.window(),
        a_tailing: f64::from(options.a_tailing.distance),
    };
    let decays = [
        (
            FilterKind::CopiedDamage,
            options.copied_damage.distance,
            CopiedDamage::FALLBACK_SCALE,
        ),
        (
            FilterKind::EndRepairFillIn,
            options.end_repair_fill_in.distance,
            EndRepairFillIn::FALLBACK_SCALE,
        ),
    ];
    for (kind, distance, fallback) in decays {
        let mut pools: BTreeMap<String, ReferencePool> = BTreeMap::new();
        let mut held: Vec<(&Distances, &str)> = Vec::new();
        for annotation in calls.iter().flatten().filter(|a| a.kind == kind) {
            if let Some(distances) = &annotation.distances {
                pools
                    .entry(annotation.stratum.clone())
                    .or_default()
                    .add(&distances.reference);
                held.push((distances, &annotation.stratum));
            }
        }
        if held.is_empty() {
            continue;
        }
        let scale = match distance {
            Distance::Bases(bases) => bases,
            Distance::Learned => learn_scale(
                |scale| {
                    held.iter()
                        .filter_map(|(d, stratum)| {
                            d.log_likelihood_ratio(scale, pools.get(*stratum))
                        })
                        .collect()
                },
                fallback,
            ),
        };
        match kind {
            FilterKind::CopiedDamage => scales.copied_damage = scale,
            _ => scales.end_repair_fill_in = scale,
        }
        if distance == Distance::Learned {
            info!(
                "{kind}: learned a decay scale of {scale:.2} bp from {}",
                plural(held.len(), "call")
            );
        }
        for annotation in calls.iter_mut().flatten().filter(|a| a.kind == kind) {
            if let Some(distances) = annotation.distances.take() {
                annotation.score = distances.score(scale, pools.get(&annotation.stratum));
            }
        }
    }
    scales
}

/// The learned artifact fractions of each filter and each of its strata.
#[derive(Clone, Debug, Default, PartialEq)]
struct Fractions {
    filters: BTreeMap<FilterKind, f64>,
    strata: BTreeMap<Stratum, f64>,
}

/// Learn the priors and fill in every annotation's posterior, returning the
/// learned artifact fractions.
fn assign_posteriors(calls: &mut [Vec<Annotation>], model: Model) -> Fractions {
    let mut filters: BTreeMap<FilterKind, Vec<f64>> = BTreeMap::new();
    let mut strata: BTreeMap<Stratum, Vec<f64>> = BTreeMap::new();
    for annotation in calls.iter().flatten() {
        if let Some(llr) = annotation.score.log_likelihood_ratio {
            filters.entry(annotation.kind).or_default().push(llr);
            strata
                .entry((annotation.kind, annotation.stratum.clone()))
                .or_default()
                .push(llr);
        }
    }
    let filters: BTreeMap<FilterKind, f64> = filters
        .into_iter()
        .map(|(kind, llrs)| (kind, learn_artifact_fraction(&llrs, FILTER_PRIOR)))
        .collect();
    let strata = strata
        .into_iter()
        .map(|(key, llrs)| {
            let prior = BetaPrior {
                mean: filters[&key.0],
                strength: STRATUM_PRIOR_STRENGTH,
            };
            (key, learn_artifact_fraction(&llrs, prior))
        })
        .collect();
    let fractions = Fractions { filters, strata };
    for annotation in calls.iter_mut().flatten() {
        let Some(llr) = annotation.score.log_likelihood_ratio else {
            continue;
        };
        let artifact_prior = match model {
            Model::Chaff => fractions.strata[&(annotation.kind, annotation.stratum.clone())],
            Model::Fgbio => annotation.fgbio_prior,
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
                    Some(Value::Float(CopiedDamage::log10_ratio(llr))),
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
    reference: Option<&mut Reference>,
    options: &FilterOptions,
) -> Result<Vec<StratumMetrics>> {
    let report = filter_vcf_report(input, output, evidence, reference, options, false)?;
    Ok(report.metrics)
}

/// What one run reports about the sample it filtered.
#[derive(Clone, Debug, PartialEq)]
pub struct Report {
    /// The sample under test.
    pub sample: String,
    /// One row per filter and stratum.
    pub metrics: Vec<StratumMetrics>,
    /// The trinucleotide spectrum of the scored SNVs, when asked for.
    pub spectrum: Option<Spectrum>,
}

/// As [`filter_vcf`], also tallying the sample's trinucleotide spectrum when
/// `spectrum` is set, which needs the reference.
pub fn filter_vcf_report(
    input: &Path,
    output: &Path,
    evidence: &mut dyn Evidence,
    mut reference: Option<&mut Reference>,
    options: &FilterOptions,
    spectrum: bool,
) -> Result<Report> {
    if options.enabled(FilterKind::CopiedDamage) && reference.is_none() {
        bail!("the copied damage filter needs a reference FASTA (--ref)");
    }
    if spectrum && reference.is_none() {
        bail!("the spectrum needs a reference FASTA (--ref)");
    }
    if options.filters.is_empty() {
        log::warn!("every filter is disabled, so chaff will copy the input unchanged");
    }

    let mut reader = VariantReader::open(input)?;
    let header = reader
        .read_header()
        .context("failed to read the VCF/BCF header")?;
    let (sample_index, sample) = resolve_sample(&header, options.sample.as_deref())?;
    refuse_earlier_annotations(&header, options)?;
    let mut order = CoordinateOrder::new(&header);
    let mut cache = None;
    let mut calls: Vec<Vec<Annotation>> = Vec::new();
    let mut channels: Vec<Option<usize>> = Vec::new();
    let mut skipped: BTreeMap<Skip, u64> = BTreeMap::new();
    let mut record = RecordBuf::default();
    while reader.read_record(&header, &mut record)? != 0 {
        let contig = record.reference_sequence_name().to_string();
        let pos = record.variant_start().map(usize::from).unwrap_or(0);
        order.check(&contig, pos)?;
        let gt = match Genotype::scored(&record, sample_index) {
            Ok(gt) if pos > 0 => gt,
            scored => {
                *skipped
                    .entry(scored.err().unwrap_or(Skip::NotSnv))
                    .or_default() += 1;
                calls.push(Vec::new());
                channels.push(None);
                continue;
            }
        };
        let bases = snv(&gt);
        let context = match (bases, reference.as_deref_mut()) {
            (Some((ref_base, _)), Some(reference)) => {
                Some(snv_context(reference, &gt, &contig, pos, ref_base)?)
            }
            _ => None,
        };
        channels.push(
            bases
                .zip(context)
                .and_then(|((_, alt_base), (prev, base, next))| {
                    channel(prev, base, alt_base, next)
                }),
        );
        let annotations = score_call(&gt, &contig, pos, options, evidence, context, &mut cache)?;
        if annotations.is_empty() {
            *skipped.entry(Skip::NoFilter).or_default() += 1;
        }
        calls.push(annotations);
    }
    log_scored(&sample, &calls, &skipped);

    let scales = score_decays(&mut calls, options);
    let fractions = assign_posteriors(&mut calls, options.model);

    let mut out_header = header.clone();
    options.add_header_lines(&mut out_header, &scales);
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

    let rows = metrics_rows(&sample, &calls, &fractions, &scales, options);
    for row in &rows {
        info!(
            "{} {}: {}, artifact fraction {}, {} filtered, {} of {} congruent",
            row.filter,
            row.stratum,
            plural(row.calls as usize, "call"),
            row.artifact_fraction
                .map_or_else(|| "per call".to_string(), |f| format!("{f:.4}")),
            row.filtered,
            row.alt_congruent,
            plural(row.alt_molecules as usize, "alternate molecule"),
        );
    }
    let spectrum = spectrum.then(|| tally_spectrum(&calls, &channels, options));
    Ok(Report {
        sample,
        metrics: rows,
        spectrum,
    })
}

/// Log how many of a sample's calls were scored, and why the others were
/// not, warning when none was.
fn log_scored(sample: &str, calls: &[Vec<Annotation>], skipped: &BTreeMap<Skip, u64>) {
    let scored = calls.iter().filter(|a| !a.is_empty()).count();
    info!(
        "scored {} of {} of sample {sample}",
        scored,
        plural(calls.len(), "call")
    );
    for (skip, count) in skipped {
        info!(
            "skipped {}: {}",
            plural(*count as usize, "call"),
            skip.reason()
        );
    }
    if scored == 0 {
        log::warn!(
            "no call of sample {sample} was scored: the filters score heterozygous SNVs, and SNVs without a genotype"
        );
    }
}

/// A count and a noun, the noun plural unless the count is one.
pub(crate) fn plural(count: usize, noun: &str) -> String {
    match count {
        1 => format!("1 {noun}"),
        n => format!("{n} {noun}s"),
    }
}

/// The trinucleotide spectrum of the calls with a channel: every one before
/// filtering, and after it the ones no filter flagged when any enabled filter
/// has a threshold, or else each weighed by the product of its posteriors.
fn tally_spectrum(
    calls: &[Vec<Annotation>],
    channels: &[Option<usize>],
    options: &FilterOptions,
) -> Spectrum {
    let thresholded = options
        .filters
        .iter()
        .any(|k| k.threshold(options).is_some());
    let mut spectrum = Spectrum::new(if thresholded {
        After::Passing
    } else {
        After::Weighted
    });
    for (annotations, channel) in calls.iter().zip(channels) {
        let Some(channel) = *channel else {
            continue;
        };
        let posteriors = annotations
            .iter()
            .filter_map(|a| a.posterior.map(|p| (a.kind, p)));
        let flagged = posteriors
            .clone()
            .any(|(kind, p)| is_filtered(p, kind.threshold(options)));
        let weight = posteriors.map(|(_, p)| p).product();
        spectrum.add(channel, flagged, weight);
    }
    spectrum
}

/// Pool every annotation into one metrics row per filter and stratum.
fn metrics_rows(
    sample: &str,
    calls: &[Vec<Annotation>],
    fractions: &Fractions,
    scales: &Scales,
    options: &FilterOptions,
) -> Vec<StratumMetrics> {
    let learned = options.model == Model::Chaff;
    let mut rows: BTreeMap<Stratum, (StratumMetrics, Vec<(u32, f64)>)> = BTreeMap::new();
    for annotation in calls.iter().flatten() {
        let key = (annotation.kind, annotation.stratum.clone());
        let (row, trials) = rows.entry(key.clone()).or_insert_with(|| {
            let row = StratumMetrics {
                sample: sample.to_string(),
                filter: annotation.kind.to_string(),
                stratum: annotation.stratum.clone(),
                artifact_fraction: learned
                    .then(|| fractions.strata.get(&key).copied())
                    .flatten(),
                filter_artifact_fraction: learned
                    .then(|| fractions.filters.get(&annotation.kind).copied())
                    .flatten(),
                distance: scales.of(annotation.kind),
                ..StratumMetrics::default()
            };
            (row, Vec::new())
        });
        row.calls += 1;
        let score = &annotation.score;
        row.alt_molecules += u64::from(score.alt_molecules);
        row.alt_congruent += u64::from(score.alt_congruent);
        row.ref_molecules += u64::from(score.ref_molecules);
        row.ref_congruent += u64::from(score.ref_congruent);
        if score.alt_molecules > 0 {
            let null = null_fraction(score.ref_congruent, score.ref_molecules);
            trials.push((score.alt_molecules, null));
        }
        match annotation.posterior {
            Some(posterior) => {
                row.expected_artifacts += 1.0 - posterior;
                row.expected_mutations += posterior;
                if is_filtered(posterior, annotation.kind.threshold(options)) {
                    row.filtered += 1;
                }
            }
            None => row.expected_mutations += 1.0,
        }
    }
    rows.into_values()
        .map(|(row, trials)| row.finish(&trials))
        .collect()
}

/// Filter the calls with molecules from `evidence`, writing the metrics and
/// the spectrum when asked.
pub fn run_filter_with(args: &FilterArgs, evidence: &mut dyn Evidence) -> Result<()> {
    let mut reference = args.reference.as_deref().map(Reference::open).transpose()?;
    let report = filter_vcf_report(
        &args.input,
        &args.output,
        evidence,
        reference.as_mut(),
        &args.options,
        args.spectrum.is_some(),
    )?;
    if let Some(path) = &args.metrics {
        write_metrics(path, &report.metrics)?;
    }
    if let (Some(path), Some(spectrum)) = (&args.spectrum, &report.spectrum) {
        spectrum.write_pdf(path, &report.sample)?;
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
            copied_damage: CopiedDamage {
                distance: Distance::Bases(30.0),
                ..CopiedDamage::default()
            },
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

    /// Without an alternate molecule, A-tailing has no evidence, so under the
    /// `chaff` model the call gets no posterior and is never filtered, while
    /// the `fgbio` model keeps fgbio's posterior from its prior alone.
    #[test]
    fn test_a_tailing_without_an_alternate_molecule_has_no_chaff_posterior() {
        let dir = tempfile::tempdir().unwrap();
        let mut vcf = VcfBuilder::new(&["tumor"]);
        vcf.add(Variant::new(10, &["G", "T"], vec![gt("tumor", "0/1")]));
        let input = vcf.write(&dir.path().join("in.vcf"));
        let mut table = MoleculeTable::new();
        let molecules: Vec<Molecule> = (0..20).map(|d| Molecule::new(b'G', 30, d, 90)).collect();
        table.insert("chr1", 10, molecules);
        for (model, scored) in [(Model::Chaff, false), (Model::Fgbio, true)] {
            let options = FilterOptions {
                filters: vec![FilterKind::ATailing],
                model,
                a_tailing_threshold: Some(1.0),
                ..FilterOptions::default()
            };
            let output = dir.path().join("out.vcf");
            filter_vcf(&input, &output, &mut table.clone(), None, &options).unwrap();
            let (_, records) = read_records(&output);
            assert_eq!(
                float(&records[0], ATailing::INFO).is_some(),
                scored,
                "{model}"
            );
            let filtered = records[0].filters().as_ref().contains(ATailing::FILTER);
            assert_eq!(filtered, scored, "{model}");
        }
    }

    /// An SNV whose genotype calls nothing is scored as heterozygous, while a
    /// homozygous alternate SNV is left as it was.
    #[test]
    fn test_an_snv_without_a_called_genotype_is_scored() {
        let dir = tempfile::tempdir().unwrap();
        let mut vcf = VcfBuilder::new(&["tumor"]);
        vcf.add(Variant::new(10, &["G", "T"], vec![gt("tumor", ".")]));
        vcf.add(Variant::new(20, &["G", "T"], vec![gt("tumor", "1/1")]));
        let input = vcf.write(&dir.path().join("in.vcf"));
        let mut table = MoleculeTable::new();
        for pos in [10, 20] {
            let mut molecules: Vec<Molecule> =
                (0..20).map(|d| Molecule::new(b'G', 30, d, 90)).collect();
            molecules.push(Molecule::new(b'T', 30, 0, 90));
            table.insert("chr1", pos, molecules);
        }
        let options = FilterOptions {
            filters: vec![FilterKind::EndRepairFillIn],
            ..FilterOptions::default()
        };
        for name in ["out.vcf", "out.bcf"] {
            let output = dir.path().join(name);
            filter_vcf(&input, &output, &mut table.clone(), None, &options).unwrap();
            let (_, records) = read_records(&output);
            assert!(
                float(&records[0], EndRepairFillIn::INFO).is_some(),
                "{name}"
            );
            assert!(
                float(&records[1], EndRepairFillIn::INFO).is_none(),
                "{name}"
            );
        }
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
    fn test_an_output_written_over_its_input_holds_every_call() {
        let dir = tempfile::tempdir().unwrap();
        let mut vcf = VcfBuilder::new(&["tumor"]);
        vcf.add(Variant::new(10, &["G", "T"], vec![gt("tumor", "0/1")]));
        vcf.add(Variant::new(20, &["G", "T"], vec![gt("tumor", "0/1")]));
        let path = vcf.write(&dir.path().join("calls.vcf"));
        let options = FilterOptions {
            filters: vec![FilterKind::EndRepairFillIn],
            ..FilterOptions::default()
        };
        filter_vcf(&path, &path, &mut MoleculeTable::new(), None, &options).unwrap();
        let (header, records) = read_records(&path);
        assert!(header.infos().contains_key(EndRepairFillIn::INFO));
        assert_eq!(records.len(), 2);
        let names: Vec<_> = std::fs::read_dir(dir.path()).unwrap().collect();
        assert_eq!(names.len(), 1, "{names:?}");
    }

    #[test]
    fn test_a_ref_that_differs_from_the_reference_fasta_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let reference = write_fasta(dir.path(), "chr1", &"ACGTTCAA".repeat(250));
        let mut vcf = VcfBuilder::new(&["tumor"]);
        vcf.add(Variant::new(1003, &["C", "T"], vec![gt("tumor", "0/1")]));
        let input = vcf.write(&dir.path().join("in.vcf"));
        let output = dir.path().join("out.vcf");
        let mut reference = Reference::open(&reference).unwrap();
        let options = FilterOptions::default();
        let error = filter_vcf(
            &input,
            &output,
            &mut MoleculeTable::new(),
            Some(&mut reference),
            &options,
        )
        .unwrap_err();
        let message = format!("{error:#}");
        assert!(message.contains("chr1:1003"), "{message}");
        assert!(message.contains("REF C"), "{message}");
        assert!(message.contains("has G"), "{message}");
    }

    /// A call whose molecules all sit far from the template ends takes no NaN
    /// into the prior, so the other call in its stratum keeps its posterior.
    #[test]
    fn test_a_call_far_from_every_end_leaves_its_stratum_finite() {
        let dir = tempfile::tempdir().unwrap();
        let mut vcf = VcfBuilder::new(&["tumor"]);
        vcf.add(Variant::new(10, &["G", "T"], vec![gt("tumor", "0/1")]));
        vcf.add(Variant::new(20, &["G", "T"], vec![gt("tumor", "0/1")]));
        let input = vcf.write(&dir.path().join("in.vcf"));
        let output = dir.path().join("out.vcf");
        let mut table = MoleculeTable::new();
        let mut near: Vec<Molecule> = (0..40).map(|d| Molecule::new(b'G', 30, d, 60)).collect();
        near.push(Molecule::new(b'T', 30, 3, 60));
        table.insert("chr1", 10, near);
        let mut far = vec![Molecule::new(b'G', 30, 900, 900); 40];
        far.push(Molecule::new(b'T', 30, 900, 900));
        table.insert("chr1", 20, far);
        let options = FilterOptions {
            filters: vec![FilterKind::EndRepairFillIn],
            end_repair_fill_in: EndRepairFillIn::new(1.0),
            ..FilterOptions::default()
        };
        let rows = filter_vcf(&input, &output, &mut table, None, &options).unwrap();
        let (_, records) = read_records(&output);
        for record in &records {
            let erfap = float(record, EndRepairFillIn::INFO).unwrap();
            assert!(erfap.is_finite(), "{erfap}");
        }
        assert!(rows[0].artifact_fraction.unwrap().is_finite());
    }

    /// Alternate molecules that copied damage cannot place, here without the
    /// template's far end, leave the call without a posterior or FILTER under
    /// either model; its counts still show that none was measured.
    #[test]
    fn test_copied_damage_without_a_measured_alternate_molecule_is_unscored() {
        let dir = tempfile::tempdir().unwrap();
        let reference = write_fasta(dir.path(), "chr1", &"ACATTCAA".repeat(250));
        let mut vcf = VcfBuilder::new(&["tumor"]);
        vcf.add(Variant::new(1002, &["C", "T"], vec![gt("tumor", "0/1")]));
        let input = vcf.write(&dir.path().join("in.vcf"));
        let output = dir.path().join("out.vcf");
        let mut table = MoleculeTable::new();
        let mut molecules: Vec<Molecule> = (0..100)
            .map(|d| Molecule::new(b'C', 40, d, 150 - d))
            .collect();
        molecules.extend((0..5).map(|d| Molecule {
            right: None,
            ..Molecule::new(b'T', 40, 40 + d, 0)
        }));
        table.insert("chr1", 1002, molecules);
        for model in [Model::Fgbio, Model::Chaff] {
            let options = FilterOptions {
                filters: vec![FilterKind::CopiedDamage],
                copied_damage_threshold: Some(0.05),
                model,
                ..FilterOptions::default()
            };
            let mut reference = Reference::open(&reference).unwrap();
            filter_vcf(&input, &output, &mut table, Some(&mut reference), &options).unwrap();
            let (_, records) = read_records(&output);
            let record = &records[0];
            assert_eq!(float(record, CopiedDamage::INFO_POSTERIOR), None, "{model}");
            assert_eq!(float(record, CopiedDamage::INFO_RATIO), None, "{model}");
            assert!(record.filters().as_ref().is_empty(), "{model}");
            assert_eq!(
                record.info().get(CopiedDamage::INFO_ALT),
                Some(Some(&Value::Array(Array::Integer(vec![Some(0), Some(0)])))),
            );
        }
    }

    /// Two calls whose alternate molecules each sit as their own reference
    /// molecules do are no asymmetry, however different the two calls are; one
    /// congruent alternate molecule with no congruent reference molecule is
    /// weak evidence, not certainty.
    #[test]
    fn test_the_asymmetry_test_compares_each_call_with_its_own_references() {
        let dir = tempfile::tempdir().unwrap();
        let reference = write_fasta(dir.path(), "chr1", &"ACATTCAA".repeat(250));
        let mut vcf = VcfBuilder::new(&["tumor"]);
        for pos in [1002, 1010, 1018] {
            vcf.add(Variant::new(pos, &["C", "T"], vec![gt("tumor", "0/1")]));
        }
        let input = vcf.write(&dir.path().join("in.vcf"));
        let output = dir.path().join("out.vcf");
        let near = |base| Molecule::new(base, 40, 10, 100);
        let far = |base| Molecule::new(base, 40, 100, 10);
        let molecules = |refs: (usize, usize), alts: (usize, usize)| -> Vec<Molecule> {
            let mut molecules = vec![near(b'C'); refs.0];
            molecules.extend(vec![far(b'C'); refs.1]);
            molecules.extend(vec![near(b'T'); alts.0]);
            molecules.extend(vec![far(b'T'); alts.1]);
            molecules
        };
        let options = FilterOptions {
            filters: vec![FilterKind::CopiedDamage],
            copied_damage: CopiedDamage {
                distance: Distance::Bases(30.0),
                ..CopiedDamage::default()
            },
            ..FilterOptions::default()
        };
        let mut table = MoleculeTable::new();
        table.insert("chr1", 1002, molecules((100, 900), (1, 9)));
        table.insert("chr1", 1010, molecules((90, 10), (90, 10)));
        let mut reference = Reference::open(&reference).unwrap();
        let rows = filter_vcf(&input, &output, &mut table, Some(&mut reference), &options).unwrap();
        let p = rows[0].asymmetry_p_value.unwrap();
        assert!(p > 0.1, "{p}");

        let mut table = MoleculeTable::new();
        table.insert("chr1", 1018, molecules((0, 50), (1, 0)));
        let rows = filter_vcf(&input, &output, &mut table, Some(&mut reference), &options).unwrap();
        let p = rows[0].asymmetry_p_value.unwrap();
        assert!((p - 1.0 / 52.0).abs() < 1e-12, "{p}");
    }

    /// Copied damage whose copies reach `d` with probability `exp(-d / 20)`
    /// teaches a scale near 20 bases, learned from its calls alongside true
    /// mutations whose alternate molecules sit anywhere.
    #[test]
    fn test_a_learned_scale_recovers_the_simulated_one() {
        let dir = tempfile::tempdir().unwrap();
        let reference = write_fasta(dir.path(), "chr1", &"ACGTTCAA".repeat(500));
        let mut state: u64 = 0x5eed;
        let mut uniform = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut vcf = VcfBuilder::new(&["tumor"]);
        let mut table = MoleculeTable::new();
        for call in 0..200 {
            let pos = 1002 + 8 * call;
            vcf.add(Variant::new(pos, &["C", "T"], vec![gt("tumor", "0/1")]));
            let at = |base, d: usize| Molecule::new(base, 40, d, 199 - d);
            let mut molecules: Vec<Molecule> = (0..200).map(|d| at(b'C', d)).collect();
            let mut alternates = 0;
            while alternates < 5 {
                let d = (uniform() * 200.0) as usize;
                if call % 2 == 1 || uniform() < (-(d as f64) / 20.0).exp() {
                    molecules.push(at(b'T', d));
                    alternates += 1;
                }
            }
            table.insert("chr1", pos, molecules);
        }
        let input = vcf.write(&dir.path().join("in.vcf"));
        let output = dir.path().join("out.vcf");
        let options = FilterOptions {
            filters: vec![FilterKind::CopiedDamage],
            ..FilterOptions::default()
        };
        let mut reference = Reference::open(&reference).unwrap();
        let rows = filter_vcf(&input, &output, &mut table, Some(&mut reference), &options).unwrap();
        let scale = rows[0].distance;
        assert!((scale - 20.0).abs() < 3.0, "{scale}");
        let fraction = rows[0].artifact_fraction.unwrap();
        assert!((fraction - 0.5).abs() < 0.1, "{fraction}");
    }

    /// Each stratum's expected mutations and expected artifacts are the sums of
    /// its calls' posteriors and their complements, a call without a
    /// posterior counting as a mutation, so together they count its calls,
    /// and the expected mutations match the calls' own `CDAP` values plus
    /// one per call without one.
    #[test]
    fn test_expected_mutations_sum_the_posteriors() {
        let dir = tempfile::tempdir().unwrap();
        let reference = write_fasta(dir.path(), "chr1", &"ACGTTCAA".repeat(250));
        let mut vcf = VcfBuilder::new(&["tumor"]);
        let mut table = MoleculeTable::new();
        let at = |base, d| Molecule::new(base, 40, d, 149 - d);
        let alternates = [[0, 1, 2, 3], [10, 50, 90, 130], [20, 60, 100, 140]];
        for (pos, distances) in [1002, 1010, 1018].into_iter().zip(alternates) {
            vcf.add(Variant::new(pos, &["C", "T"], vec![gt("tumor", "0/1")]));
            let mut molecules: Vec<Molecule> = (0..150).map(|d| at(b'C', d)).collect();
            molecules.extend(distances.map(|d| at(b'T', d)));
            table.insert("chr1", pos, molecules);
        }
        vcf.add(Variant::new(1026, &["C", "T"], vec![gt("tumor", "0/1")]));
        table.insert("chr1", 1026, (0..4).map(|d| at(b'T', d)).collect());
        let input = vcf.write(&dir.path().join("in.vcf"));
        let output = dir.path().join("out.vcf");
        let options = FilterOptions {
            filters: vec![FilterKind::CopiedDamage],
            ..FilterOptions::default()
        };
        let mut reference = Reference::open(&reference).unwrap();
        let rows = filter_vcf(&input, &output, &mut table, Some(&mut reference), &options).unwrap();
        let (_, records) = read_records(&output);
        let cdap: Vec<f64> = records
            .iter()
            .filter_map(|r| float(r, CopiedDamage::INFO_POSTERIOR).map(f64::from))
            .collect();
        assert_eq!(cdap.len(), 3);
        let row = &rows[0];
        assert_eq!(row.calls, 4);
        assert!(
            (row.expected_mutations + row.expected_artifacts - 4.0).abs() < 1e-9,
            "{row:?}"
        );
        let expected = cdap.iter().sum::<f64>() + 1.0;
        assert!(
            (row.expected_mutations - expected).abs() < 0.01,
            "{} vs {expected}",
            row.expected_mutations
        );
    }

    /// A second run on a first run's output would keep the first run's FILTERs
    /// under the second run's header lines, so chaff refuses a VCF that already
    /// declares an INFO or FILTER of a filter it is asked to run.
    #[test]
    fn test_a_vcf_annotated_by_a_filter_is_refused_by_it() {
        let dir = tempfile::tempdir().unwrap();
        let mut vcf = VcfBuilder::new(&["tumor"]);
        vcf.add(Variant::new(10, &["G", "T"], vec![gt("tumor", "0/1")]));
        let input = vcf.write(&dir.path().join("in.vcf"));
        let first = dir.path().join("first.vcf");
        let second = dir.path().join("second.vcf");
        let a_tailing = FilterOptions {
            filters: vec![FilterKind::ATailing],
            ..FilterOptions::default()
        };
        let mut table = MoleculeTable::new();
        filter_vcf(&input, &first, &mut table, None, &a_tailing).unwrap();
        let error = filter_vcf(&first, &second, &mut table, None, &a_tailing).unwrap_err();
        let message = format!("{error:#}");
        assert!(message.contains("INFO ATAP"), "{message}");
        assert!(message.contains("FILTER ATailingArtifact"), "{message}");
        assert!(!second.exists());
        let end_repair = FilterOptions {
            filters: vec![FilterKind::EndRepairFillIn],
            ..FilterOptions::default()
        };
        filter_vcf(&first, &second, &mut table, None, &end_repair).unwrap();
    }

    #[test]
    fn test_bcf_and_bgzf_inputs_score_as_a_plain_vcf_does() {
        let dir = tempfile::tempdir().unwrap();
        let mut vcf = VcfBuilder::new(&["tumor"]);
        vcf.add(Variant::new(10, &["G", "T"], vec![gt("tumor", "0/1")]));
        let plain = vcf.write(&dir.path().join("in.vcf"));
        let mut table = MoleculeTable::new();
        let mut molecules: Vec<Molecule> =
            (0..20).map(|d| Molecule::new(b'G', 30, d, 90)).collect();
        molecules.push(Molecule::new(b'T', 30, 0, 90));
        table.insert("chr1", 10, molecules);
        let copy = FilterOptions {
            filters: Vec::new(),
            ..FilterOptions::default()
        };
        let options = FilterOptions {
            filters: vec![FilterKind::ATailing, FilterKind::EndRepairFillIn],
            ..FilterOptions::default()
        };
        let scored = |input: &Path| {
            let output = dir.path().join("out.vcf");
            filter_vcf(input, &output, &mut table.clone(), None, &options).unwrap();
            let (_, records) = read_records(&output);
            let record = &records[0];
            (
                float(record, ATailing::INFO),
                float(record, EndRepairFillIn::INFO),
            )
        };
        let expected = scored(&plain);
        assert!(expected.0.is_some() && expected.1.is_some());
        for name in ["in.bcf", "in.vcf.gz"] {
            let input = dir.path().join(name);
            filter_vcf(&plain, &input, &mut table.clone(), None, &copy).unwrap();
            assert_eq!(scored(&input), expected, "{name}");
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

    /// The spectrum counts each heterozygous SNV in its channel, read from the
    /// pyrimidine, whatever its FILTER, and after filtering keeps the calls no
    /// threshold flags or, without a threshold, weighs each by its posterior.
    #[test]
    fn test_the_spectrum_counts_scored_snvs_before_and_after_filtering() {
        let dir = tempfile::tempdir().unwrap();
        let reference = write_fasta(dir.path(), "chr1", &"ACGTTCAA".repeat(250));
        let mut vcf = VcfBuilder::new(&["tumor"]);
        vcf.add(Variant::new(1002, &["C", "T"], vec![gt("tumor", "0/1")]));
        vcf.add(Variant {
            filters: vec!["LowQD".to_string()],
            ..Variant::new(1006, &["C", "T"], vec![gt("tumor", "0/1")])
        });
        vcf.add(Variant::new(1010, &["C", "T"], vec![gt("tumor", "1/1")]));
        vcf.add(Variant::new(1014, &["C", "CA"], vec![gt("tumor", "0/1")]));
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
        let mut reference = Reference::open(&reference).unwrap();
        let (acg, tca) = (16 * 2 + 2, 16 * 2 + 12);
        for threshold in [Some(0.05), None] {
            let options = FilterOptions {
                filters: vec![FilterKind::CopiedDamage],
                copied_damage: CopiedDamage {
                    distance: Distance::Bases(30.0),
                    ..CopiedDamage::default()
                },
                copied_damage_threshold: threshold,
                ..FilterOptions::default()
            };
            let report = filter_vcf_report(
                &input,
                &output,
                &mut table.clone(),
                Some(&mut reference),
                &options,
                true,
            )
            .unwrap();
            let spectrum = report.spectrum.unwrap();
            assert_eq!(spectrum.before.iter().sum::<f64>(), 2.0);
            assert_eq!((spectrum.before[acg], spectrum.before[tca]), (1.0, 1.0));
            let others: f64 =
                spectrum.after.iter().sum::<f64>() - spectrum.after[acg] - spectrum.after[tca];
            assert_eq!(others, 0.0);
            let (_, records) = read_records(&output);
            let cdap =
                |i: usize| f64::from(float(&records[i], CopiedDamage::INFO_POSTERIOR).unwrap());
            if threshold.is_some() {
                assert_eq!(spectrum.after_kind, After::Passing);
                assert_eq!((spectrum.after[acg], spectrum.after[tca]), (0.0, 1.0));
            } else {
                assert_eq!(spectrum.after_kind, After::Weighted);
                assert!((spectrum.after[acg] - cdap(0)).abs() < 1e-3, "{spectrum:?}");
                assert!((spectrum.after[tca] - cdap(1)).abs() < 1e-3, "{spectrum:?}");
            }
        }
    }

    #[test]
    fn test_the_spectrum_requires_a_reference() {
        let dir = tempfile::tempdir().unwrap();
        let input = VcfBuilder::new(&["tumor"]).write(&dir.path().join("in.vcf"));
        let options = FilterOptions {
            filters: vec![FilterKind::ATailing],
            ..FilterOptions::default()
        };
        let output = dir.path().join("out.vcf");
        let mut table = MoleculeTable::new();
        let error =
            filter_vcf_report(&input, &output, &mut table, None, &options, true).unwrap_err();
        assert!(
            error.to_string().contains("the spectrum needs a reference"),
            "{error}"
        );
    }
}
