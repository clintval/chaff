//! Read-end artifact filters ported from fgbio's `FilterSomaticVcf`: end
//! repair fill-in (`ERFAP`) and A-tailing (`ATAP`).
//!
//! Both ask whether a call's alternate molecules sit closer to a template end
//! than its reference molecules do. A molecule is congruent with the artifact
//! when the site lies within a distance of the relevant template end, counted
//! in the template's bases. Distances are 1-based, so the template's terminal
//! base is at distance 1, as in fgbio.
//!
//! The likelihoods are fgbio's. Let `f` be the congruent fraction of the
//! reference molecules and `e` an alternate molecule's base error probability:
//!
//! ```text
//!                       congruent alternate       incongruent alternate
//! P(. | artifact)       1 - e                     e
//! P(. | mutation)       (1 - e) f + e (1 - f)     (1 - e)(1 - f) + e f
//! ```
//!
//! These are the likelihoods of the `fgbio` model. Under the `chaff` model, end
//! repair fill-in instead uses the decay of [`tilt_log_likelihood_ratio`] on
//! the 0-based distance from the 3' end of the strand the template was copied
//! from, the end its polymerase extended, at a learned scale; A-tailing keeps
//! its window. Either way a molecule is congruent when the site lies within
//! the filter's distance of its end.

use std::collections::BTreeMap;

use crate::call::Genotype;
use crate::classes::Strand;
use crate::evidence::Molecule;
use crate::model::Distance;
use crate::prior::ln_add_exp;

/// What a filter extracts from one call's molecules.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Score {
    /// `ln P(molecules | artifact) - ln P(molecules | mutation)`, always
    /// finite, or `None` when alternate molecules are seen but no reference
    /// molecule can calibrate the null, or, for the decay models and for
    /// A-tailing under the `chaff` model, when no alternate molecule is
    /// measured.
    pub log_likelihood_ratio: Option<f64>,
    /// Alternate molecules measured.
    pub alt_molecules: u32,
    /// Alternate molecules congruent with the artifact.
    pub alt_congruent: u32,
    /// Reference molecules measured.
    pub ref_molecules: u32,
    /// Reference molecules congruent with the artifact.
    pub ref_congruent: u32,
}

/// The probability a base call is wrong, from its Phred quality, capped at the
/// error of a random base so a quality-zero base cannot zero a likelihood.
pub fn error_probability(quality: u8) -> f64 {
    10f64.powf(-f64::from(quality) / 10.0).min(0.75)
}

/// fgbio's windowed likelihood ratio over a call's molecules.
pub fn window_score(
    molecules: &[Molecule],
    ref_base: u8,
    alt_base: u8,
    congruent: impl Fn(&Molecule) -> bool,
) -> Score {
    let mut score = Score::default();
    for m in molecules.iter().filter(|m| m.base == ref_base) {
        score.ref_molecules += 1;
        if congruent(m) {
            score.ref_congruent += 1;
        }
    }
    let f = f64::from(score.ref_congruent) / f64::from(score.ref_molecules);
    let mut ll_artifact = 0.0;
    let mut ll_mutation = 0.0;
    for m in molecules.iter().filter(|m| m.base == alt_base) {
        score.alt_molecules += 1;
        let e = error_probability(m.quality);
        if congruent(m) {
            score.alt_congruent += 1;
            ll_artifact += (-e).ln_1p();
            ll_mutation += ((1.0 - e) * f + e * (1.0 - f)).ln();
        } else {
            ll_artifact += e.ln();
            ll_mutation += ((1.0 - e) * (1.0 - f) + e * f).ln();
        }
    }
    score.log_likelihood_ratio = if score.alt_molecules > 0 && score.ref_molecules == 0 {
        None
    } else {
        Some(ll_artifact - ll_mutation)
    };
    score
}

/// The reference molecules a decay's calls hold, pooled per stratum by
/// distance and span, toward which each call's own reference molecules are
/// shrunk.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReferencePool {
    distances: BTreeMap<usize, u64>,
    spans: BTreeMap<usize, u64>,
    molecules: u64,
}

impl ReferencePool {
    /// Add one call's reference distances and spans.
    pub fn add(&mut self, distances: &[usize], spans: &[usize]) {
        for &d in distances {
            *self.distances.entry(d).or_default() += 1;
        }
        for &s in spans {
            *self.spans.entry(s).or_default() += 1;
        }
        self.molecules += distances.len() as u64;
    }

    /// `ln` of the mean decay weight `exp(-d / scale)` over the pooled
    /// molecules, `ln W_p`, shrunk toward their [`ln_spread_weight`] by
    /// [`REFERENCE_PSEUDO_MOLECULES`], or `None` without any.
    pub fn ln_mean_weight(&self, scale: f64) -> Option<f64> {
        let terms = self
            .distances
            .iter()
            .map(|(&d, &n)| (n as f64).ln() - d as f64 / scale);
        let n = self.molecules as f64;
        let spans = self.spans.iter().map(|(&s, &n)| (s, n as f64));
        Some(ln_shrunk(
            ln_sum_exp(terms)? - n.ln(),
            n,
            ln_spread_weight(spans, scale),
        ))
    }
}

/// The pseudo-molecules each mean reference weight is shrunk toward the one
/// above it by: ten, as a stratum's fraction is shrunk toward its filter's by
/// ten pseudo-calls.
pub const REFERENCE_PSEUDO_MOLECULES: f64 = 10.0;

/// `ln U`, the mean decay weight `exp(-d / scale)` of a site spread evenly
/// along its molecules' templates, from their spans and counts: a span of `S`
/// distances gives `(1 - exp(-S / scale)) / (S (1 - exp(-1 / scale)))`, about
/// `scale / S`. `None` without a span.
pub fn ln_spread_weight(spans: impl Iterator<Item = (usize, f64)>, scale: f64) -> Option<f64> {
    let spans: Vec<(f64, f64)> = spans.map(|(s, n)| (s as f64, n)).collect();
    let count: f64 = spans.iter().map(|(_, n)| n).sum();
    let terms = spans
        .iter()
        .map(|&(s, n)| n.ln() + (-(-s / scale).exp_m1()).ln() - s.ln());
    Some(ln_sum_exp(terms)? - count.ln() - (-(-1.0 / scale).exp_m1()).ln())
}

/// `ln((n W + k U) / (n + k))`: a mean weight `W` over `n` molecules shrunk
/// toward `U` by [`REFERENCE_PSEUDO_MOLECULES`] `k`, in log space, or `W`
/// itself without a `U`.
fn ln_shrunk(ln_weight: f64, n: f64, ln_toward: Option<f64>) -> f64 {
    match ln_toward {
        Some(ln_toward) => {
            ln_add_exp(
                n.ln() + ln_weight,
                REFERENCE_PSEUDO_MOLECULES.ln() + ln_toward,
            ) - (n + REFERENCE_PSEUDO_MOLECULES).ln()
        }
        None => ln_weight,
    }
}

/// `ln(sum(exp(x)))` without overflow, or `None` for no terms.
fn ln_sum_exp(terms: impl Iterator<Item = f64>) -> Option<f64> {
    let terms: Vec<f64> = terms.collect();
    let max = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return None;
    }
    Some(max + terms.iter().map(|t| (t - max).exp()).sum::<f64>().ln())
}

/// A continuous artifact model: the artifact's alternate molecules follow the
/// reference molecules' distance distribution tilted by `w(d) = exp(-d /
/// scale)`, the chance a copy reaches distance `d`. With `W` the mean of `w`
/// over the reference molecules and `e` a base error probability, each
/// alternate molecule contributes
///
/// ```text
/// ln((1 - e) w(d) / W + e)
/// ```
///
/// since a base error lands anywhere a reference molecule could. A site with
/// few reference molecules measures `W` poorly, so `W` is shrunk by
/// [`REFERENCE_PSEUDO_MOLECULES`] `k` toward `ln_pool_weight`, the
/// [`ReferencePool::ln_mean_weight`] at `scale` of its stratum's reference
/// molecules, or without a pool toward the [`ln_spread_weight`] of its own
/// reference molecules' `ref_spans`: `W = (n W_site + k W_p) / (n + k)` for
/// `n` reference molecules at the site. The pool is itself shrunk toward its
/// molecules' spread weight, which bounds `W` from below: reference molecules
/// that all sit far from the end, at a site or across its pool, can't make
/// alternates merely nearer the end than they are look like a copy. The
/// terms are summed in log space, so the ratio is finite however far a
/// molecule sits from the end. Returns `None` without both a reference and an
/// alternate distance, as such a call carries no evidence either way.
pub fn tilt_log_likelihood_ratio(
    ref_distances: &[usize],
    ref_spans: &[usize],
    alt: &[(usize, u8)],
    scale: f64,
    ln_pool_weight: Option<f64>,
) -> Option<f64> {
    if alt.is_empty() || ref_distances.is_empty() {
        return None;
    }
    let ln_w = |d: usize| -(d as f64) / scale;
    let n = ref_distances.len() as f64;
    let ln_site = ln_sum_exp(ref_distances.iter().map(|&d| ln_w(d)))? - n.ln();
    let ln_toward =
        ln_pool_weight.or_else(|| ln_spread_weight(ref_spans.iter().map(|&s| (s, 1.0)), scale));
    let ln_mean = ln_shrunk(ln_site, n, ln_toward);
    Some(
        alt.iter()
            .map(|&(d, q)| {
                let e = error_probability(q);
                ln_add_exp((-e).ln_1p() + ln_w(d) - ln_mean, e.ln())
            })
            .sum(),
    )
}

/// The 0-based distances a decay scores one call by: its measured reference
/// molecules', and its measured alternate molecules' with their base
/// qualities.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Distances {
    /// The reference molecules' distances.
    pub reference: Vec<usize>,
    /// The spans of the reference molecules whose templates' two ends are
    /// known: how many distances a site anywhere on the template could take.
    pub spans: Vec<usize>,
    /// The alternate molecules' distances and base qualities.
    pub alternate: Vec<(usize, u8)>,
    /// Whether most of the call's molecules are a duplex consensus, which
    /// holds both strands.
    pub duplex: bool,
}

impl Distances {
    /// The distances `distance` gives a call's reference and alternate
    /// molecules, leaving out the molecules it gives none, and the spans
    /// `span` gives its measured reference molecules.
    pub fn of(
        molecules: &[Molecule],
        ref_base: u8,
        alt_base: u8,
        distance: impl Fn(&Molecule) -> Option<usize>,
        span: impl Fn(&Molecule) -> Option<usize>,
    ) -> Self {
        let mut distances = Self::default();
        for m in molecules {
            match distance(m) {
                Some(d) if m.base == ref_base => {
                    distances.reference.push(d);
                    distances.spans.extend(span(m));
                }
                Some(d) if m.base == alt_base => distances.alternate.push((d, m.quality)),
                _ => {}
            }
        }
        let duplex = molecules.iter().filter(|m| m.origin.is_none()).count();
        distances.duplex = 2 * duplex > molecules.len();
        distances
    }

    /// The [`tilt_log_likelihood_ratio`] at `scale`, shrunk toward a pool's
    /// mean weight.
    pub fn log_likelihood_ratio(&self, scale: f64, ln_pool_weight: Option<f64>) -> Option<f64> {
        tilt_log_likelihood_ratio(
            &self.reference,
            &self.spans,
            &self.alternate,
            scale,
            ln_pool_weight,
        )
    }

    /// The score at `scale`, shrunk toward a pool's mean weight: its ratio,
    /// and the molecules less than `scale` bases from the end as congruent.
    pub fn score(&self, scale: f64, ln_pool_weight: Option<f64>) -> Score {
        let near = |d: usize| (d as f64) < scale;
        let count = |n: usize| u32::try_from(n).unwrap_or(u32::MAX);
        Score {
            log_likelihood_ratio: self.log_likelihood_ratio(scale, ln_pool_weight),
            alt_molecules: count(self.alternate.len()),
            alt_congruent: count(self.alternate.iter().filter(|(d, _)| near(*d)).count()),
            ref_molecules: count(self.reference.len()),
            ref_congruent: count(self.reference.iter().filter(|&&d| near(d)).count()),
        }
    }
}

/// The 1-based distance of the site from the nearest template end the
/// molecule knows.
fn nearest_end(m: &Molecule) -> Option<usize> {
    m.left.into_iter().chain(m.right).min().map(|d| d + 1)
}

/// The end repair fill-in artifact filter.
///
/// End repair blunts a fragment: polymerase extends a recessed 3' end across
/// the opposite strand's 5' overhang, and an exonuclease trims a 3' overhang.
/// A base the polymerase misincorporates, or copies from a damaged overhang,
/// sits on the extended strand near its 3' end.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EndRepairFillIn {
    /// The distance in bases: the decay scale under the `chaff` model, learned
    /// by default, and the window from the nearest template end under `fgbio`,
    /// fgbio's [`EndRepairFillIn::FGBIO_WINDOW`] by default.
    pub distance: Distance,
}

impl EndRepairFillIn {
    /// The INFO key, fgbio's.
    pub const INFO: &'static str = "ERFAP";
    /// The FILTER name, fgbio's.
    pub const FILTER: &'static str = "EndRepairFillInArtifact";

    /// fgbio's default window, in bases.
    pub const FGBIO_WINDOW: f64 = 15.0;

    /// The decay scale without a call to learn it from, in bases.
    pub const FALLBACK_SCALE: f64 = 15.0;

    /// The decay scale of a duplex consensus without a call to learn it from,
    /// in bases: the changes both of its strands agree on crowd the last few
    /// bases of a template.
    pub const DUPLEX_FALLBACK_SCALE: f64 = 5.0;

    /// The decay scale to learn from, and to shrink a learned scale toward,
    /// for calls that are mostly a duplex consensus or not.
    pub fn fallback_scale(duplex: bool) -> f64 {
        if duplex {
            Self::DUPLEX_FALLBACK_SCALE
        } else {
            Self::FALLBACK_SCALE
        }
    }

    /// A filter at a fixed `distance` in bases.
    pub fn new(distance: f64) -> Self {
        Self {
            distance: Distance::Bases(distance),
        }
    }

    /// The window of the `fgbio` model.
    pub fn window(&self) -> f64 {
        self.distance.bases(Self::FGBIO_WINDOW)
    }

    /// The 0-based distance of the site from the 3' end of the strand the
    /// template was copied from, which end repair extended: the rightmost base
    /// for the forward strand, the leftmost for the reverse, and the nearer for
    /// a duplex consensus, where either strand's end can hold the error.
    pub fn three_prime_distance(m: &Molecule) -> Option<usize> {
        match m.origin {
            Some(Strand::Forward) => m.right,
            Some(Strand::Reverse) => m.left,
            None => m.left.into_iter().chain(m.right).min(),
        }
    }

    /// How many distances [`Self::three_prime_distance`] could give a site
    /// anywhere on the molecule's template, when both its ends are known: the
    /// template's length, or half of it, rounded up, from the nearer end of a
    /// duplex consensus.
    pub fn three_prime_span(m: &Molecule) -> Option<usize> {
        let length = m.length()?;
        Some(match m.origin {
            Some(_) => length,
            None => length.div_ceil(2),
        })
    }

    /// The distances the `chaff` model's decay scores a call by.
    pub fn distances(&self, molecules: &[Molecule], ref_base: u8, alt_base: u8) -> Distances {
        Distances::of(
            molecules,
            ref_base,
            alt_base,
            Self::three_prime_distance,
            Self::three_prime_span,
        )
    }

    /// Heterozygous calls whose every called allele is one base.
    pub fn applies_to(gt: &Genotype) -> bool {
        gt.is_het() && gt.calls_are_single_bases()
    }

    /// Whether the site is within the window of the nearest known template
    /// end.
    pub fn is_congruent(&self, m: &Molecule) -> bool {
        nearest_end(m).is_some_and(|d| d as f64 <= self.window())
    }

    /// The call's score under the `fgbio` model, from fgbio's window.
    pub fn score(&self, molecules: &[Molecule], ref_base: u8, alt_base: u8) -> Score {
        window_score(molecules, ref_base, alt_base, |m| self.is_congruent(m))
    }
}

/// The A-tailing artifact filter.
///
/// End repair can over-digest a 3' end and leave it recessed; A-tailing then
/// fills it with adenines. The template gains `A`s at a 3' end, which read as
/// `T`s at the 5' end of the complementary strand: on the forward strand, a
/// `T` near the leftmost template end or an `A` near the rightmost one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ATailing {
    /// The distance from a template end within which a molecule is congruent.
    pub distance: u32,
}

impl Default for ATailing {
    fn default() -> Self {
        Self { distance: 2 }
    }
}

impl ATailing {
    /// The INFO key, fgbio's.
    pub const INFO: &'static str = "ATAP";
    /// The FILTER name, fgbio's.
    pub const FILTER: &'static str = "ATailingArtifact";

    /// Heterozygous, not heterozygous non-reference, one-base calls with an
    /// `A` or `T` alternate allele.
    pub fn applies_to(gt: &Genotype) -> bool {
        gt.is_het()
            && !gt.is_het_non_ref()
            && gt.calls_are_single_bases()
            && gt
                .alts()
                .any(|a| a.eq_ignore_ascii_case("A") || a.eq_ignore_ascii_case("T"))
    }

    /// Whether a molecule, of either allele, is congruent: the site is within
    /// the distance of its nearest template end, and that end is where the
    /// alternate allele would appear by A addition (the leftmost end for a
    /// forward-strand `T`, the rightmost for an `A`). At a tie the nearest end
    /// is the kept mate's own, as fgbio has it.
    pub fn is_congruent(&self, alt_base: u8, m: &Molecule) -> bool {
        let left = m.left.map(|d| d + 1);
        let right = m.right.map(|d| d + 1);
        let within = |d: Option<usize>| d.is_some_and(|d| d <= self.distance as usize);
        let (near_left, near_right) = match (left, right) {
            (Some(l), Some(r)) if l == r => (!m.reverse, m.reverse),
            (Some(l), Some(r)) => (l < r, r < l),
            (Some(_), None) => (true, false),
            (None, Some(_)) => (false, true),
            (None, None) => (false, false),
        };
        match alt_base.to_ascii_uppercase() {
            b'T' => near_left && within(left),
            b'A' => near_right && within(right),
            _ => false,
        }
    }

    /// The call's score from its molecules.
    pub fn score(&self, molecules: &[Molecule], ref_base: u8, alt_base: u8) -> Score {
        window_score(molecules, ref_base, alt_base, |m| {
            self.is_congruent(alt_base, m)
        })
    }
}

/// Whether a threshold flags a posterior: at or below it, as fgbio does.
pub fn is_filtered(posterior_mutation: f64, threshold: Option<f64>) -> bool {
    threshold.is_some_and(|t| posterior_mutation <= t)
}

#[cfg(test)]
mod tests {
    use noodles::core::Position;
    use noodles::sam::alignment::RecordBuf;
    use streampile::testing::{Frag, Pair, SamBuilder, Strand};

    use super::*;
    use crate::classes::Strand as Origin;
    use crate::evidence::{Evidence, PileupEvidence, PileupOptions};
    use crate::io::vcf_float;
    use crate::prior::{
        fgbio_artifact_prior, learn_artifact_fraction, posterior_mutation, BetaPrior, FILTER_PRIOR,
        STRATUM_PRIOR_STRENGTH,
    };

    const A: u8 = b'A';
    const C: u8 = b'C';
    const G: u8 = b'G';
    const T: u8 = b'T';

    fn single_genotype(alleles: &[&str]) -> Genotype {
        Genotype::new(alleles, alleles)
    }

    /// The molecule one read alone shows at `pos`, like one of fgbio's
    /// `BaseEntry` values.
    fn entry(record: &RecordBuf, pos: usize) -> Molecule {
        let mut reads = SamBuilder::new();
        reads.extend([record.clone()]);
        let molecules = molecules_at(&reads, pos);
        assert_eq!(molecules.len(), 1, "{record:?}");
        molecules[0]
    }

    /// The posterior fgbio reports for one call, with its `(2 * maf)^2` prior.
    /// The exact values asserted against it were produced by fgbio 4.1.1
    /// `FilterSomaticVcf` on the same reads.
    fn fgbio_posterior(molecules: &[Molecule], score: &Score) -> f64 {
        let depth = molecules.len() as u32;
        let prior = fgbio_artifact_prior(score.alt_molecules, score.ref_molecules, depth);
        posterior_mutation(score.log_likelihood_ratio.unwrap(), prior)
    }

    /// The posterior chaff reports for a VCF holding only this call.
    fn learned_posterior(score: &Score) -> f64 {
        let llr = score.log_likelihood_ratio.unwrap();
        let filter = BetaPrior {
            mean: learn_artifact_fraction(&[llr], FILTER_PRIOR),
            strength: STRATUM_PRIOR_STRENGTH,
        };
        let prior = learn_artifact_fraction(&[llr], filter);
        posterior_mutation(llr, prior)
    }

    fn molecules_at(reads: &SamBuilder, pos: usize) -> Vec<Molecule> {
        let mut evidence =
            PileupEvidence::new(reads.to_pileup_builder(), &PileupOptions::default());
        evidence
            .molecules("chr1", Position::try_from(pos).unwrap())
            .unwrap()
    }

    /// fgbio: "filters should only set the filter if a threshold was provide
    /// and the pvalue is <= threshold".
    #[test]
    fn test_filters_only_set_with_a_threshold_and_a_value_at_or_below_it() {
        for p in [1.0, 1e-3, 1e-4, 1e-5, 0.0] {
            assert!(!is_filtered(p, None));
        }
        let threshold = Some(0.001);
        assert!(!is_filtered(1.0, threshold));
        assert!(!is_filtered(0.01, threshold));
        assert!(is_filtered(0.001, threshold));
        assert!(is_filtered(0.00099, threshold));
        assert!(is_filtered(1e-20, threshold));
        assert!(is_filtered(0.0, threshold));
    }

    /// fgbio: "EndRepairFillInArtifactLikelihoodFilter.appliesTo should return
    /// false for any event that is not a SNP".
    #[test]
    fn test_end_repair_fill_in_does_not_apply_to_non_snvs() {
        assert!(!EndRepairFillIn::applies_to(&single_genotype(&["G", "GT"])));
        assert!(!EndRepairFillIn::applies_to(&single_genotype(&["GT", "G"])));
        assert!(!EndRepairFillIn::applies_to(&single_genotype(&[
            "CT", "GC"
        ])));
    }

    /// fgbio: "... should return false for homozygous genotypes".
    #[test]
    fn test_end_repair_fill_in_does_not_apply_to_homozygous_genotypes() {
        assert!(!EndRepairFillIn::applies_to(&Genotype::new(
            &["A"],
            &["A", "A"]
        )));
    }

    /// fgbio: "... should return false for events with multiple alt alleles if
    /// any alts are not SNVs".
    #[test]
    fn test_end_repair_fill_in_does_not_apply_when_any_alt_is_not_an_snv() {
        let gt = Genotype::new(&["A", "G"], &["A", "G", "TT"]);
        assert!(!EndRepairFillIn::applies_to(&gt));
    }

    /// fgbio: "... should return true for any event that is a SNP".
    #[test]
    fn test_end_repair_fill_in_applies_to_every_snv() {
        for r in ["A", "C", "G", "T"] {
            for a in ["A", "C", "G", "T"] {
                if r != a {
                    assert!(EndRepairFillIn::applies_to(&single_genotype(&[r, a])));
                }
            }
        }
    }

    /// fgbio: "... should return true for events with multiple alt alleles if
    /// they are all SNVs".
    #[test]
    fn test_end_repair_fill_in_applies_when_every_alt_is_an_snv() {
        let gt = Genotype::new(&["A", "G"], &["A", "G", "T"]);
        assert!(EndRepairFillIn::applies_to(&gt));
    }

    /// fgbio: "EndRepairFillInArtifactLikelihoodFilter.isArtifactCongruent
    /// should return false for any base that is not within the defined read end".
    #[test]
    fn test_end_repair_fill_in_incongruent_away_from_the_ends() {
        let filter = EndRepairFillIn::new(15.0);
        let mut builder = SamBuilder::new().read_length(50);
        let recs = builder.add_pair(
            Pair::at(101, 101)
                .bases1("ACACA".repeat(10))
                .bases2("ACACA".repeat(10)),
        );
        for r in &recs {
            for pos in 116..=135 {
                assert!(!filter.is_congruent(&entry(r, pos)));
            }
        }
    }

    /// fgbio: "... should return true for any base within the defined read end".
    #[test]
    fn test_end_repair_fill_in_congruent_near_the_ends() {
        let filter = EndRepairFillIn::new(15.0);
        let mut builder = SamBuilder::new().read_length(50);
        let recs = builder.add_pair(
            Pair::at(101, 101)
                .bases1("ACACA".repeat(10))
                .bases2("ACACA".repeat(10)),
        );
        for r in &recs {
            for pos in (101..=115).chain(136..=150) {
                assert!(filter.is_congruent(&entry(r, pos)));
            }
        }
    }

    /// A deletion between the site and the template's leftmost base brings the
    /// site nearer that end, in template bases, than its reference positions.
    #[test]
    fn test_end_repair_fill_in_counts_template_bases_across_a_deletion() {
        let filter = EndRepairFillIn::new(15.0);
        let mut builder = SamBuilder::new().read_length(50);
        let plain = builder.add_frag(Frag::at(101));
        let deleted = builder.add_frag(Frag::at(101).cigar("10M5D40M"));
        let (plain, deleted) = (entry(&plain[0], 120), entry(&deleted[0], 120));
        assert_eq!((plain.left, deleted.left), (Some(19), Some(14)));
        assert!(!filter.is_congruent(&plain));
        assert!(filter.is_congruent(&deleted));
    }

    /// The reads of fgbio's "distributed throughout the reads" annotation tests:
    /// ten reference fragments per start and strand over starts 1 to 25, and
    /// one alternate fragment per start and strand.
    fn distributed_reads() -> SamBuilder {
        let mut builder = SamBuilder::new().read_length(50).base_quality(40);
        for start in 1..=25 {
            for _ in 1..=10 {
                for strand in [Strand::Plus, Strand::Minus] {
                    builder.add_frag(Frag::at(start).strand(strand).bases("G".repeat(50)));
                }
            }
        }
        for start in 1..=25 {
            for strand in [Strand::Plus, Strand::Minus] {
                builder.add_frag(Frag::at(start).strand(strand).bases("T".repeat(50)));
            }
        }
        builder
    }

    /// Five reference fragments per start and strand over starts 1 to 25, and
    /// one forward alternate fragment per start in `alt_starts`.
    fn biased_reads(alt_starts: std::ops::RangeInclusive<usize>) -> SamBuilder {
        let mut builder = SamBuilder::new().read_length(50).base_quality(40);
        for start in 1..=25 {
            for _ in 1..=5 {
                for strand in [Strand::Plus, Strand::Minus] {
                    builder.add_frag(Frag::at(start).strand(strand).bases("G".repeat(50)));
                }
            }
        }
        for start in alt_starts {
            builder.add_frag(Frag::at(start).bases("T".repeat(50)));
        }
        builder
    }

    /// fgbio: "EndRepairFillInArtifactLikelihoodFilter.annotations should
    /// compute a non-significant p-value when data is distributed throughout
    /// the reads".
    #[test]
    fn test_end_repair_fill_in_not_significant_when_distributed() {
        let filter = EndRepairFillIn::new(15.0);
        let molecules = molecules_at(&distributed_reads(), 25);
        let score = filter.score(&molecules, G, T);
        assert!(score.log_likelihood_ratio.is_some());
        assert!(fgbio_posterior(&molecules, &score) > 0.5);
        assert_eq!(vcf_float(fgbio_posterior(&molecules, &score)), 1.0);
        assert!(learned_posterior(&score) > 0.5);
    }

    /// fgbio: "... should compute a significant p-value when data is heavily
    /// biased".
    #[test]
    fn test_end_repair_fill_in_significant_when_biased() {
        let filter = EndRepairFillIn::new(15.0);
        let molecules = molecules_at(&biased_reads(11..=25), 25);
        let score = filter.score(&molecules, G, T);
        assert!(score.log_likelihood_ratio.is_some());
        assert!(fgbio_posterior(&molecules, &score) < 1e-6);
        assert_eq!(vcf_float(fgbio_posterior(&molecules, &score)), 1.869e-10);
        assert!(learned_posterior(&score) < 1e-6);
    }

    /// fgbio: "ATailingArtifactLikelihoodFilter.appliesTo should return false
    /// for any event that is not a SNP".
    #[test]
    fn test_a_tailing_does_not_apply_to_non_snvs() {
        assert!(!ATailing::applies_to(&single_genotype(&["A", "AT"])));
        assert!(!ATailing::applies_to(&single_genotype(&["TA", "A"])));
        assert!(!ATailing::applies_to(&single_genotype(&["AT", "GC"])));
    }

    /// fgbio: "... should return true for *>A and *>T and false for all other
    /// SNVs".
    #[test]
    fn test_a_tailing_applies_only_to_a_and_t_alternates() {
        let applies = |r, a| ATailing::applies_to(&single_genotype(&[r, a]));
        assert!(!applies("A", "C"));
        assert!(!applies("A", "G"));
        assert!(applies("A", "T"));
        assert!(applies("C", "A"));
        assert!(!applies("C", "G"));
        assert!(applies("C", "T"));
        assert!(applies("G", "A"));
        assert!(!applies("G", "C"));
        assert!(applies("G", "T"));
        assert!(applies("T", "A"));
        assert!(!applies("T", "C"));
        assert!(!applies("T", "G"));
    }

    /// fgbio: "ATailingArtifactLikelihoodFilter.isArtifactCongruent should
    /// return false for any base in a read that is not near the ends".
    #[test]
    fn test_a_tailing_incongruent_away_from_the_ends() {
        let filter = ATailing { distance: 5 };
        let mut builder = SamBuilder::new().read_length(50);
        let recs = builder.add_pair(
            Pair::at(101, 101)
                .bases1("ACACA".repeat(10))
                .bases2("ACACA".repeat(10)),
        );
        for r in &recs {
            for pos in 106..=145 {
                assert!(!filter.is_congruent(A, &entry(r, pos)));
            }
        }
    }

    /// fgbio: "... should return false for any base at the ends of the insert
    /// that is not allele-matched".
    #[test]
    fn test_a_tailing_incongruent_at_the_ends_when_not_allele_matched() {
        let filter = ATailing { distance: 5 };
        let mut builder = SamBuilder::new().read_length(50);
        let bases = format!("CACAC{}GTGTG", "N".repeat(40));
        let recs: Vec<_> = builder
            .add_pair(Pair::at(101, 101))
            .into_iter()
            .map(|r| SamBuilder::with_bases(r, &bases))
            .collect();
        for r in &recs {
            for pos in 101..=105 {
                let m = entry(r, pos);
                assert!(m.base == C || m.base == A);
                assert!(!filter.is_congruent(A, &m));
            }
            for pos in 146..=150 {
                let m = entry(r, pos);
                assert!(m.base == G || m.base == T);
                assert!(!filter.is_congruent(T, &m));
            }
        }
    }

    /// fgbio: "... should return true for bases at the ends of the insert that
    /// are allele-matched".
    #[test]
    fn test_a_tailing_congruent_at_the_ends_when_allele_matched() {
        let filter = ATailing { distance: 5 };
        let mut builder = SamBuilder::new().read_length(50);
        let bases = format!("GTGTG{}CACAC", "N".repeat(40));
        let recs: Vec<_> = builder
            .add_pair(Pair::at(101, 101))
            .into_iter()
            .map(|r| SamBuilder::with_bases(r, &bases))
            .collect();
        for r in &recs {
            for pos in 101..=105 {
                assert!(filter.is_congruent(T, &entry(r, pos)));
            }
            for pos in 146..=150 {
                assert!(filter.is_congruent(A, &entry(r, pos)));
            }
        }
    }

    /// Soft-clipped bases are template bases and hard-clipped bases are not.
    #[test]
    fn test_a_tailing_counts_soft_clips_and_not_hard_clips() {
        let filter = ATailing::default();
        let mut builder = SamBuilder::new().read_length(50);
        let soft = builder.add_frag(Frag::at(101).bases("T".repeat(50)).cigar("2S48M"));
        let hard = builder.add_frag(Frag::at(101).bases("T".repeat(48)).cigar("2H48M"));
        let (soft, hard) = (entry(&soft[0], 101), entry(&hard[0], 101));
        assert_eq!((soft.left, hard.left), (Some(2), Some(0)));
        assert!(!filter.is_congruent(T, &soft));
        assert!(filter.is_congruent(T, &hard));
    }

    /// fgbio: "ATailingArtifactLikelihoodFilter.annotations should compute a
    /// non-significant p-value when data is distributed throughout the reads".
    #[test]
    fn test_a_tailing_not_significant_when_distributed() {
        let filter = ATailing { distance: 5 };
        let molecules = molecules_at(&distributed_reads(), 25);
        let score = filter.score(&molecules, G, T);
        assert!(score.log_likelihood_ratio.is_some());
        assert!(fgbio_posterior(&molecules, &score) > 0.5);
        assert_eq!(vcf_float(fgbio_posterior(&molecules, &score)), 1.0);
        assert!(learned_posterior(&score) > 0.5);
    }

    /// fgbio: "... should compute a significant p-value when data is heavily
    /// biased". Intended difference: with a prior learned from this lone call
    /// (about 0.7) rather than `(2 * maf)^2` (about 0.0015), five congruent
    /// molecules at a congruent reference fraction of 0.1 give about 4e-6, not
    /// fgbio's 1.5e-8.
    #[test]
    fn test_a_tailing_significant_when_biased() {
        let filter = ATailing { distance: 5 };
        let molecules = molecules_at(&biased_reads(21..=25), 25);
        let score = filter.score(&molecules, G, T);
        assert!(score.log_likelihood_ratio.is_some());
        assert!(fgbio_posterior(&molecules, &score) < 1e-6);
        assert_eq!(vcf_float(fgbio_posterior(&molecules, &score)), 1.547e-8);
        let learned = learned_posterior(&score);
        assert!(learned > 1e-6 && learned < 1e-5, "{learned}");
    }

    /// fgbio: "... should compute an intermediate p-value when data is heavily
    /// biased for both ref and alt". Intended difference: with the learned
    /// prior the posterior stays intermediate, about 0.005, where fgbio's prior
    /// drives it below 1e-4.
    #[test]
    fn test_a_tailing_intermediate_when_both_alleles_are_biased() {
        let filter = ATailing { distance: 5 };
        let mut builder = SamBuilder::new().read_length(50).base_quality(40);
        for start in 20..=25 {
            for _ in 1..=10 {
                for strand in [Strand::Plus, Strand::Minus] {
                    builder.add_frag(Frag::at(start).strand(strand).bases("G".repeat(50)));
                }
            }
        }
        for start in 21..=25 {
            builder.add_frag(Frag::at(start).bases("T".repeat(50)));
        }
        let molecules = molecules_at(&builder, 25);
        let score = filter.score(&molecules, G, T);
        assert!(score.log_likelihood_ratio.is_some());
        assert!(fgbio_posterior(&molecules, &score) < 1e-4);
        assert_eq!(vcf_float(fgbio_posterior(&molecules, &score)), 8.094e-5);
        let learned = learned_posterior(&score);
        assert!(learned > 1e-3 && learned < 1e-2, "{learned}");
    }

    /// fgbio: "... should not throw an exception if there is no alt allele
    /// coverage or the pileup is empty".
    #[test]
    fn test_a_tailing_annotates_without_alt_coverage_or_reads() {
        let filter = ATailing { distance: 10 };
        let mut builder = SamBuilder::new().read_length(50).base_quality(40);

        let molecules = molecules_at(&builder, 25);
        let score = filter.score(&molecules, G, T);
        assert_eq!(score.log_likelihood_ratio, Some(0.0));
        assert!((fgbio_posterior(&molecules, &score) - 0.9999).abs() < 1e-12);

        for start in 1..=20 {
            for _ in 1..=10 {
                for strand in [Strand::Plus, Strand::Minus] {
                    builder.add_frag(Frag::at(start).strand(strand).bases("G".repeat(50)));
                }
            }
        }
        let molecules = molecules_at(&builder, 25);
        let score = filter.score(&molecules, G, T);
        assert_eq!(score.log_likelihood_ratio, Some(0.0));
        assert!((fgbio_posterior(&molecules, &score) - 0.000025).abs() < 1e-12);
    }

    #[test]
    fn test_window_score_without_reference_molecules_has_no_ratio() {
        let alt = Molecule::new(T, 30, 1, 50);
        let score = window_score(&[alt], G, T, |_| true);
        assert_eq!(score.log_likelihood_ratio, None);
        assert_eq!(score.alt_molecules, 1);
    }

    #[test]
    fn test_window_score_matches_hand_computation() {
        let refs = (0..4).map(|i| Molecule::new(G, 30, i * 10, 200));
        let alt = Molecule::new(T, 30, 100, 200);
        let molecules: Vec<_> = refs.chain([alt]).collect();
        let score = window_score(&molecules, G, T, |m| {
            m.left == Some(0) || m.left == Some(100)
        });
        let e: f64 = 1e-3;
        let f = 0.25;
        let expected = (1.0 - e).ln() - ((1.0 - e) * f + e * (1.0 - f)).ln();
        assert!((score.log_likelihood_ratio.unwrap() - expected).abs() < 1e-12);
        assert_eq!((score.ref_congruent, score.ref_molecules), (1, 4));
        assert_eq!((score.alt_congruent, score.alt_molecules), (1, 1));
    }

    #[test]
    fn test_tilt_ratio_rewards_alternates_nearer_the_end_than_the_references() {
        let refs: Vec<usize> = (0..100).collect();
        let spans = vec![100; 100];
        let ratio = |alt: &[(usize, u8)]| tilt_log_likelihood_ratio(&refs, &spans, alt, 15.0, None);
        let near = ratio(&[(1, 90), (3, 90)]).unwrap();
        let far = ratio(&[(80, 90), (95, 90)]).unwrap();
        assert!(near > 0.0, "{near}");
        assert!(far < 0.0, "{far}");
        assert_eq!(ratio(&[]), None);
        assert_eq!(
            tilt_log_likelihood_ratio(&[], &[], &[(1, 30)], 15.0, None),
            None
        );
    }

    /// Weights of distances hundreds of scales from the end underflow in
    /// linear space; the ratio stays finite and keeps its sign. Without a span
    /// to bound `W`, the tilt is relative, so shifting every distance alike
    /// leaves the ratio as it was, and an alternate molecule at the end of
    /// templates whose references all sit 800 bases in counts as hundreds of
    /// nats; their spans bound it to about `ln(S / scale)`.
    #[test]
    fn test_tilt_ratio_is_finite_far_from_the_end() {
        let refs = [800, 900];
        let ratio = |refs: &[usize], spans: &[usize], d: usize| {
            tilt_log_likelihood_ratio(refs, spans, &[(d, 30)], 1.0, None).unwrap()
        };
        let flat = ratio(&refs, &[], 850);
        let near = ratio(&refs, &[], 0);
        let far = ratio(&refs, &[], 5000);
        assert!(flat.is_finite() && near.is_finite() && far.is_finite());
        assert!(near > 700.0, "{near}");
        assert!(far < 0.0, "{far}");
        let shifted = ratio(&[0, 100], &[], 50);
        assert!((flat - shifted).abs() < 1e-9, "{flat} vs {shifted}");
        let bounded = ratio(&refs, &[2000, 2000], 0);
        assert!(bounded > 6.0 && bounded < 8.0, "{bounded}");
        assert!(ratio(&refs, &[2000, 2000], 5000) < 0.0);
    }

    #[test]
    fn test_tilt_ratio_is_flat_when_alternates_match_the_references() {
        let refs = vec![5, 5, 5];
        let llr = tilt_log_likelihood_ratio(&refs, &[], &[(5, 255)], 30.0, None).unwrap();
        assert!(llr.abs() < 1e-12, "{llr}");
    }

    /// A site spread evenly along a template takes each distance of its span
    /// equally often: each distance from one end once, and each distance of
    /// the first half twice from the nearer end of a duplex consensus.
    #[test]
    fn test_a_span_holds_the_distances_a_site_anywhere_on_its_template_takes() {
        let scale = 7.0;
        for (length, origin) in [
            (40, Some(Origin::Forward)),
            (40, None),
            (41, Some(Origin::Reverse)),
        ] {
            let sites: Vec<Molecule> = (0..length)
                .map(|left| {
                    let m = Molecule::new(G, 40, left, length - 1 - left);
                    origin.map_or(m, |strand| m.from_strand(strand))
                })
                .collect();
            let weights: f64 = sites
                .iter()
                .map(|m| {
                    (-(EndRepairFillIn::three_prime_distance(m).unwrap() as f64) / scale).exp()
                })
                .sum();
            let span = EndRepairFillIn::three_prime_span(&sites[0]).unwrap();
            let spread = ln_spread_weight([(span, 3.0)].into_iter(), scale).unwrap();
            assert!(
                (spread - (weights / length as f64).ln()).abs() < 1e-12,
                "{length} {origin:?}"
            );
        }
        assert_eq!(
            EndRepairFillIn::three_prime_span(&Molecule::new(G, 40, 20, 20)),
            Some(21)
        );
        let half = Molecule {
            left: None,
            ..Molecule::new(G, 40, 20, 20)
        };
        assert_eq!(EndRepairFillIn::three_prime_span(&half), None);
        let short = ln_spread_weight([(1000, 1.0)].into_iter(), 5.0).unwrap();
        assert!((short - (5.0f64 / 1000.0).ln()).abs() < 0.1, "{short}");
        assert_eq!(ln_spread_weight(std::iter::empty(), 5.0), None);
    }

    /// Reference molecules that all sit farther from the end than the
    /// alternate molecules, at the site and across its pool alike, set `W` by
    /// the nearest of them, so without spans two alternates 31 bases from the
    /// end read as a copy at a 5-base scale, though a copy reaches 31 bases
    /// with a weight of `exp(-31 / 5)`. Shrunk toward the weight molecules
    /// spread evenly along their 200-base templates would give, the same
    /// alternates favor a mutation, with or without a pool, while two at the
    /// end still favor the artifact.
    #[test]
    fn test_references_all_farther_from_the_end_than_distant_alternates_bound_the_weight() {
        let refs = [35, 40, 45, 50, 55, 60];
        let spans = [200; 6];
        let distant = [(31, 40), (31, 40)];
        let mut pool = ReferencePool::default();
        pool.add(&refs, &spans);
        let weight = pool.ln_mean_weight(5.0);
        let ratio = |spans: &[usize], alt: &[(usize, u8)], pool: Option<f64>| {
            tilt_log_likelihood_ratio(&refs, spans, alt, 5.0, pool).unwrap()
        };

        let mut unbounded = ReferencePool::default();
        unbounded.add(&refs, &[]);
        let unbounded = unbounded.ln_mean_weight(5.0);
        assert!(ratio(&[], &distant, None) > 4.0);
        assert!(ratio(&[], &distant, unbounded) > 4.0);

        let alone = ratio(&spans, &distant, None);
        let pooled = ratio(&spans, &distant, weight);
        assert!(alone < -4.0, "{alone}");
        assert!(pooled < -3.0, "{pooled}");
        let near = ratio(&spans, &[(0, 40), (1, 40)], weight);
        assert!(near > 7.0, "{near}");
    }

    /// A pool that holds reference molecules at nearly every distance, here
    /// none within three bases of the end, keeps nearly its own mean weight:
    /// the even spread it is shrunk toward weighs as ten of its 3,000
    /// molecules.
    #[test]
    fn test_a_pool_at_every_distance_barely_moves_toward_the_even_spread() {
        let refs: Vec<usize> = (3..150).cycle().take(3000).collect();
        let mut spread = ReferencePool::default();
        spread.add(&refs, &vec![150; refs.len()]);
        let mut unbounded = ReferencePool::default();
        unbounded.add(&refs, &[]);
        let alt = [(0, 40), (1, 40), (60, 40)];
        let site: Vec<usize> = refs.iter().copied().step_by(15).collect();
        for scale in [5.0, 30.0] {
            let ratio = |pool: &ReferencePool| {
                tilt_log_likelihood_ratio(&site, &[], &alt, scale, pool.ln_mean_weight(scale))
                    .unwrap()
            };
            let (spread, unbounded) = (ratio(&spread), ratio(&unbounded));
            assert!(
                (spread - unbounded).abs() < 0.01,
                "{scale}: {spread} {unbounded}"
            );
        }
    }

    /// Reference molecules spread along 150-base templates copied from
    /// `origin`, and four alternate molecules `d` bases from the leftmost end.
    fn copied_from(origin: Option<Origin>, alt_left: usize) -> Vec<Molecule> {
        let at = |base, d: usize| {
            let m = Molecule::new(base, 40, d, 149 - d);
            match origin {
                Some(strand) => m.from_strand(strand),
                None => m,
            }
        };
        let mut molecules: Vec<Molecule> = (0..150).map(|d| at(G, d)).collect();
        molecules.extend((0..4).map(|i| at(T, alt_left + i)));
        molecules
    }

    /// One reference molecule far from the end can't make three alternates
    /// near it look like a copy: the site's reference weight is shrunk toward
    /// its stratum's pool, or without one toward the even spread its own
    /// template gives, while a site with hundreds of reference molecules
    /// keeps nearly its own.
    #[test]
    fn test_a_site_with_few_reference_molecules_is_shrunk_toward_its_pool() {
        let mut pool = ReferencePool::default();
        pool.add(
            &(0..150).cycle().take(3000).collect::<Vec<_>>(),
            &[150; 3000],
        );
        let weight = pool.ln_mean_weight(30.0);
        let alt = [(0, 40), (1, 40), (2, 40)];
        let ratio = |refs: &[usize], spans: &[usize], pool| {
            tilt_log_likelihood_ratio(refs, spans, &alt, 30.0, pool).unwrap()
        };
        let lone = ratio(&[100], &[], None);
        let shrunk = ratio(&[100], &[150], weight);
        let alone = ratio(&[100], &[150], None);
        assert!(lone > 9.5, "{lone}");
        assert!(shrunk < 5.5 && shrunk > 0.0, "{shrunk}");
        assert!((alone - shrunk).abs() < 1e-9, "{alone} {shrunk}");
        let many: Vec<usize> = (0..150).cycle().take(300).collect();
        let own = ratio(&many, &[150; 300], None);
        let pooled = ratio(&many, &[150; 300], weight);
        assert!((own - pooled).abs() < 0.01, "{own} {pooled}");
        assert_eq!(
            tilt_log_likelihood_ratio(&[], &[], &alt, 30.0, weight),
            None
        );
        assert!(ReferencePool::default().ln_mean_weight(30.0).is_none());
    }

    /// Fill-in errors sit near the 3' end of the strand a template was copied
    /// from: the rightmost end for the forward strand, the leftmost for the
    /// reverse, and either for a duplex consensus. fgbio's window counts
    /// either end for every template.
    #[test]
    fn test_end_repair_fill_in_decays_from_the_copied_strand_s_three_prime_end() {
        let filter = EndRepairFillIn::new(15.0);
        let decay = |origin, alt_left| {
            let molecules = copied_from(origin, alt_left);
            filter.distances(&molecules, G, T).score(15.0, None)
        };
        let window = |origin, alt_left| {
            let molecules = copied_from(origin, alt_left);
            filter.score(&molecules, G, T).log_likelihood_ratio.unwrap()
        };
        let (left, right) = (1, 145);
        let score = decay(Some(Origin::Reverse), left);
        assert!(score.log_likelihood_ratio.unwrap() > 5.0, "{score:?}");
        assert_eq!((score.alt_congruent, score.ref_congruent), (4, 15));
        assert!(
            decay(Some(Origin::Forward), right)
                .log_likelihood_ratio
                .unwrap()
                > 5.0
        );
        assert!(decay(None, left).log_likelihood_ratio.unwrap() > 5.0);
        assert!(decay(None, right).log_likelihood_ratio.unwrap() > 5.0);

        let wrong = decay(Some(Origin::Forward), left);
        assert!(wrong.log_likelihood_ratio.unwrap() < -5.0, "{wrong:?}");
        assert_eq!(wrong.alt_congruent, 0);
        assert!(
            decay(Some(Origin::Reverse), right)
                .log_likelihood_ratio
                .unwrap()
                < -5.0
        );
        assert!(window(Some(Origin::Forward), left) > 5.0);
    }

    /// The window counts an alternate molecule 15 bases from the end as the
    /// artifact and one 16 bases away as a mutation; the decay weighs the two
    /// almost alike.
    #[test]
    fn test_end_repair_fill_in_decay_has_no_cliff_at_its_distance() {
        let filter = EndRepairFillIn::new(15.0);
        let refs: Vec<Molecule> = (0..150).map(|d| Molecule::new(G, 40, d, 149 - d)).collect();
        let molecules = |d: usize| {
            let mut molecules = refs.clone();
            molecules.push(Molecule::new(T, 40, d - 1, 150 - d));
            molecules
        };
        let decay = |d| {
            let distances = filter.distances(&molecules(d), G, T);
            distances.log_likelihood_ratio(15.0, None).unwrap()
        };
        let window = |d| {
            let score = filter.score(&molecules(d), G, T);
            score.log_likelihood_ratio.unwrap()
        };
        let (near, far) = (decay(15), decay(16));
        assert!((near - far).abs() < 0.1, "{near} {far}");
        let (near, far) = (window(15), window(16));
        assert!((near - far).abs() > 5.0, "{near} {far}");
    }

    /// fgbio takes the other end only when it is strictly nearer than the kept
    /// mate's own 5' end.
    #[test]
    fn test_a_tailing_tie_goes_to_the_kept_read_s_own_end() {
        let filter = ATailing::default();
        let forward = Molecule::new(A, 30, 1, 1);
        assert!(filter.is_congruent(T, &forward));
        assert!(!filter.is_congruent(A, &forward));
        let reverse = Molecule {
            reverse: true,
            ..forward
        };
        assert!(filter.is_congruent(A, &reverse));
        assert!(!filter.is_congruent(T, &reverse));
        assert!(!filter.is_congruent(C, &reverse));
    }
}
