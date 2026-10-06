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
//! With `--end-repair-fill-in-scale`, end repair fill-in instead uses the
//! continuous model of [`tilt_log_likelihood_ratio`] on the 0-based distance
//! from the nearest template end.

use crate::call::Genotype;
use crate::evidence::Molecule;

/// What a filter extracts from one call's molecules.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Score {
    /// `ln P(molecules | artifact) - ln P(molecules | mutation)`, always
    /// finite, or `None` when alternate molecules are seen but no reference
    /// molecule can calibrate the null, or, for the decay models, when no
    /// alternate molecule is measured.
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
/// since a base error lands anywhere a reference molecule could. The weights
/// are taken relative to the nearest reference molecule's and the terms summed
/// in log space, so the ratio is finite however far a molecule sits from the
/// end. Returns `None` without both a reference and an alternate distance, as
/// such a call carries no evidence either way.
pub fn tilt_log_likelihood_ratio(
    ref_distances: &[usize],
    alt: &[(usize, u8)],
    scale: f64,
) -> Option<f64> {
    if alt.is_empty() {
        return None;
    }
    let nearest = *ref_distances.iter().min()? as f64;
    let ln_w = |d: usize| -(d as f64 - nearest) / scale;
    let sum: f64 = ref_distances.iter().map(|&d| ln_w(d).exp()).sum();
    let ln_mean = (sum / ref_distances.len() as f64).ln();
    Some(
        alt.iter()
            .map(|&(d, q)| {
                let e = error_probability(q);
                ln_add_exp((-e).ln_1p() + ln_w(d) - ln_mean, e.ln())
            })
            .sum(),
    )
}

/// `ln(exp(a) + exp(b))`, without overflow.
fn ln_add_exp(a: f64, b: f64) -> f64 {
    a.max(b) + (-(a - b).abs()).exp().ln_1p()
}

/// The [`tilt_log_likelihood_ratio`] of a call's reference and alternate
/// molecules at the 0-based distances `distance` gives them, leaving out the
/// molecules it gives none.
pub fn molecule_tilt_ratio(
    molecules: &[Molecule],
    ref_base: u8,
    alt_base: u8,
    scale: f64,
    distance: impl Fn(&Molecule) -> Option<usize>,
) -> Option<f64> {
    let mut ref_distances = Vec::new();
    let mut alt = Vec::new();
    for m in molecules {
        match distance(m) {
            Some(d) if m.base == ref_base => ref_distances.push(d),
            Some(d) if m.base == alt_base => alt.push((d, m.quality)),
            _ => {}
        }
    }
    tilt_log_likelihood_ratio(&ref_distances, &alt, scale)
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
/// A damaged base in the single-stranded overhang, such as 8-oxoguanine, is
/// copied into the extended strand, so after amplification both strands carry
/// the change near the template end.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EndRepairFillIn {
    /// The distance from a template end within which a molecule is congruent.
    pub distance: u32,
    /// The scale of the continuous model, when it replaces the window.
    pub scale: Option<f64>,
}

impl Default for EndRepairFillIn {
    fn default() -> Self {
        Self {
            distance: 15,
            scale: None,
        }
    }
}

impl EndRepairFillIn {
    /// The INFO key, fgbio's.
    pub const INFO: &'static str = "ERFAP";
    /// The FILTER name, fgbio's.
    pub const FILTER: &'static str = "EndRepairFillInArtifact";

    /// A window filter at `distance`.
    pub fn new(distance: u32) -> Self {
        Self {
            distance,
            scale: None,
        }
    }

    /// Heterozygous calls whose every called allele is one base.
    pub fn applies_to(gt: &Genotype) -> bool {
        gt.is_het() && gt.calls_are_single_bases()
    }

    /// Whether the site is within the distance of the nearest known template end.
    pub fn is_congruent(&self, m: &Molecule) -> bool {
        nearest_end(m).is_some_and(|d| d <= self.distance as usize)
    }

    /// The call's score from its molecules.
    pub fn score(&self, molecules: &[Molecule], ref_base: u8, alt_base: u8) -> Score {
        let mut score = window_score(molecules, ref_base, alt_base, |m| self.is_congruent(m));
        if let Some(scale) = self.scale {
            score.log_likelihood_ratio =
                molecule_tilt_ratio(molecules, ref_base, alt_base, scale, |m| {
                    nearest_end(m).map(|d| d - 1)
                });
        }
        score
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
    /// forward-strand `T`, the rightmost for an `A`). At a tie both ends are
    /// candidates.
    pub fn is_congruent(&self, alt_base: u8, m: &Molecule) -> bool {
        let left = m.left.map(|d| d + 1);
        let right = m.right.map(|d| d + 1);
        let within = |d: Option<usize>| d.is_some_and(|d| d <= self.distance as usize);
        let (near_left, near_right) = match (left, right) {
            (Some(l), Some(r)) => (l <= r, r <= l),
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
        let filter = EndRepairFillIn::new(15);
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
        let filter = EndRepairFillIn::new(15);
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
        let filter = EndRepairFillIn::new(15);
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
        let filter = EndRepairFillIn::new(15);
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
        let filter = EndRepairFillIn::new(15);
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
        let near = tilt_log_likelihood_ratio(&refs, &[(1, 90), (3, 90)], 15.0).unwrap();
        let far = tilt_log_likelihood_ratio(&refs, &[(80, 90), (95, 90)], 15.0).unwrap();
        assert!(near > 0.0, "{near}");
        assert!(far < 0.0, "{far}");
        assert_eq!(tilt_log_likelihood_ratio(&refs, &[], 15.0), None);
        assert_eq!(tilt_log_likelihood_ratio(&[], &[(1, 30)], 15.0), None);
    }

    /// Weights of distances hundreds of scales from the end underflow in
    /// linear space; the ratio stays finite and keeps its sign.
    #[test]
    fn test_tilt_ratio_is_finite_far_from_the_end() {
        let refs = [800, 900];
        let flat = tilt_log_likelihood_ratio(&refs, &[(850, 30)], 1.0).unwrap();
        let near = tilt_log_likelihood_ratio(&refs, &[(0, 30)], 1.0).unwrap();
        let far = tilt_log_likelihood_ratio(&refs, &[(5000, 30)], 1.0).unwrap();
        assert!(flat.is_finite() && near.is_finite() && far.is_finite());
        assert!(near > 700.0, "{near}");
        assert!(far < 0.0, "{far}");
        let shifted = tilt_log_likelihood_ratio(&[0, 100], &[(50, 30)], 1.0).unwrap();
        assert!((flat - shifted).abs() < 1e-9, "{flat} vs {shifted}");
    }

    #[test]
    fn test_tilt_ratio_is_flat_when_alternates_match_the_references() {
        let refs = vec![5, 5, 5];
        let llr = tilt_log_likelihood_ratio(&refs, &[(5, 255)], 30.0).unwrap();
        assert!(llr.abs() < 1e-12, "{llr}");
    }

    #[test]
    fn test_continuous_end_repair_fill_in_uses_the_tilt_model() {
        let filter = EndRepairFillIn {
            distance: 15,
            scale: Some(15.0),
        };
        let molecules = molecules_at(&biased_reads(11..=25), 25);
        let score = filter.score(&molecules, G, T);
        assert!(score.log_likelihood_ratio.unwrap() > 0.0);
        assert_eq!(score.alt_congruent, 15);
    }

    #[test]
    fn test_the_distance_window_sets_the_congruent_counts_under_a_scale() {
        let molecules = molecules_at(&biased_reads(11..=25), 25);
        let score = |distance| {
            EndRepairFillIn {
                distance,
                scale: Some(15.0),
            }
            .score(&molecules, G, T)
        };
        let (wide, narrow) = (score(15), score(5));
        assert_eq!(wide.log_likelihood_ratio, narrow.log_likelihood_ratio);
        assert_eq!((wide.alt_congruent, narrow.alt_congruent), (15, 5));
    }

    #[test]
    fn test_a_tailing_tie_is_congruent_for_either_end() {
        let filter = ATailing::default();
        let m = Molecule::new(A, 30, 1, 1);
        assert!(filter.is_congruent(A, &m));
        assert!(filter.is_congruent(T, &m));
        assert!(!filter.is_congruent(C, &m));
    }
}
