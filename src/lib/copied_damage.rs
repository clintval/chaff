//! The copied damage filter: damage that polymerase copied onto the other strand
//! before the strands were tagged.
//!
//! A lesion on one strand, such as a deaminated cytosine or an 8-oxoguanine,
//! is templated into the complementary strand when polymerase resynthesizes it
//! during end repair fill-in, nick translation, or gap filling. Both strands
//! then carry the change, so duplex consensus agrees on it. Resynthesis runs
//! 5' to 3' along the new strand, so it copies the lesion strand from a nick
//! or recessed end toward the lesion strand's 5' end: copied lesions sit near
//! the lesion strand's 5' end and are depleted near its 3' end, over tens to
//! more than a hundred bases.
//!
//! The lesion strand comes from the substitution class: a `C>T` class puts the
//! lesion on the strand that carries the reference `C`, a `G>T` class on the
//! strand that carries the reference `G`. For each molecule with both template
//! ends known, `d` is the 0-based distance of the site from the lesion strand's
//! 5' end in the template's bases: from the leftmost base for a forward-strand
//! lesion, from the rightmost for a reverse-strand one.
//!
//! Under a true mutation, the alternate molecules' distances follow the
//! reference molecules' distances at the same site, which absorbs capture and
//! fragment-length skew. Under the artifact, a lesion at distance `d` survives
//! into a copy with probability `w(d) = exp(-d / distance)`, an exponentially
//! distributed resynthesis length with mean `distance`; the alternate distances
//! follow the reference distances tilted by `w`. See
//! [`crate::read_end::tilt_log_likelihood_ratio`] for the per-molecule terms.
//! The decay holds under both models, since fgbio has no such filter; its
//! scale is learned per library unless fixed, and a molecule is congruent when
//! `d` is less than it.
//!
//! The model sees only copies that start at a 5' end. A lesion copied from an
//! internal nick, by nick translation or strand displacement, or across the
//! gap an abasic site leaves, can sit anywhere in the template, so its
//! alternate molecules look like a mutation's: per-call scores cannot separate
//! them, and they show only as an excess of the damage class across a library.

use crate::call::Genotype;
use crate::classes::{Context, DamageClass, Strand};
use crate::evidence::Molecule;
use crate::model::{CopiedDamagePrior, Distance};
use crate::read_end::{Distances, Score};

/// The copied damage filter.
#[derive(Clone, Debug, PartialEq)]
pub struct CopiedDamage {
    /// The damage classes to test.
    pub classes: Vec<DamageClass>,
    /// The decay scale in bases from the lesion strand's 5' end, the mean
    /// resynthesis length: learned by default.
    pub distance: Distance,
    /// Where each call's artifact prior comes from under the `chaff` model:
    /// the library's chance model by default.
    pub prior: CopiedDamagePrior,
}

impl Default for CopiedDamage {
    fn default() -> Self {
        Self {
            classes: vec![DamageClass::DEAMINATION, DamageClass::OXIDATION],
            distance: Distance::Learned,
            prior: CopiedDamagePrior::Chance,
        }
    }
}

/// The damage class, lesion strand, and CpG context one call is tested under.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DamageSite {
    /// The class that explains the call.
    pub class: DamageClass,
    /// The strand carrying the lesion.
    pub strand: Strand,
    /// Whether the lesion base sits in a CpG.
    pub context: Context,
}

impl DamageSite {
    /// The stratum label, e.g. `C>T:CpG`.
    pub fn stratum(&self) -> String {
        format!("{}:{}", self.class, self.context)
    }
}

impl CopiedDamage {
    /// The INFO key for the posterior probability of a true mutation.
    pub const INFO_POSTERIOR: &'static str = "CDAP";
    /// The INFO key for the log10 likelihood ratio, artifact to mutation.
    pub const INFO_RATIO: &'static str = "CDLR";
    /// The INFO key for the alternate molecules within the distance of the
    /// lesion strand's 5' end, and all alternate molecules measured.
    pub const INFO_ALT: &'static str = "CDAC";
    /// The INFO key for the reference molecules within the distance of the
    /// lesion strand's 5' end, and all reference molecules measured.
    pub const INFO_REF: &'static str = "CDRC";
    /// The FILTER name.
    pub const FILTER: &'static str = "CopiedDamageArtifact";

    /// The decay scale without a call to learn it from, in bases.
    pub const FALLBACK_SCALE: f64 = 30.0;

    /// The `CDLR` value of a natural-log likelihood ratio: its log10 at four
    /// significant digits, whatever its sign.
    pub fn log10_ratio(log_likelihood_ratio: f64) -> f32 {
        crate::io::significant(log_likelihood_ratio / std::f64::consts::LN_10, 4)
    }

    /// Heterozygous calls whose every called allele is one base.
    pub fn applies_to(gt: &Genotype) -> bool {
        gt.is_het() && gt.calls_are_single_bases()
    }

    /// The first configured class that explains `ref_base>alt_base`, with the
    /// strand it puts the lesion on.
    pub fn classify(&self, ref_base: u8, alt_base: u8) -> Option<(DamageClass, Strand)> {
        self.classes
            .iter()
            .find_map(|class| class.lesion_strand(ref_base, alt_base).map(|s| (*class, s)))
    }

    /// The 0-based distance of the site from the lesion strand's 5' end, in the
    /// template's bases, for a molecule with both ends known.
    pub fn five_prime_distance(m: &Molecule, strand: Strand) -> Option<usize> {
        m.length()?;
        match strand {
            Strand::Forward => m.left,
            Strand::Reverse => m.right,
        }
    }

    /// The distances from the lesion strand's 5' end that the decay scores a
    /// call by.
    pub fn distances(
        &self,
        molecules: &[Molecule],
        ref_base: u8,
        alt_base: u8,
        strand: Strand,
    ) -> Distances {
        Distances::of(
            molecules,
            ref_base,
            alt_base,
            |m| Self::five_prime_distance(m, strand),
            Molecule::length,
        )
    }

    /// The call's score at a decay scale of `scale` bases.
    pub fn score(
        &self,
        molecules: &[Molecule],
        ref_base: u8,
        alt_base: u8,
        strand: Strand,
        scale: f64,
    ) -> Score {
        self.distances(molecules, ref_base, alt_base, strand)
            .score(scale, None)
    }
}

#[cfg(test)]
mod tests {
    use noodles::core::Position;
    use streampile::testing::{Pair, SamBuilder};

    use super::*;
    use crate::evidence::{Evidence, PileupEvidence, PileupOptions};

    const C: u8 = b'C';
    const G: u8 = b'G';
    const T: u8 = b'T';
    const A: u8 = b'A';

    /// A molecule 150 bp long whose leftmost base is `d` bases before the site.
    fn at(base: u8, d: usize) -> Molecule {
        Molecule::new(base, 90, d, 149 - d)
    }

    /// Reference molecules spread evenly around a site, 150 bp long.
    fn spread_references(base: u8) -> Vec<Molecule> {
        (0..150).map(|d| at(base, d)).collect()
    }

    #[test]
    fn test_classify_picks_the_class_and_strand() {
        let filter = CopiedDamage::default();
        assert_eq!(
            filter.classify(C, T),
            Some((DamageClass::DEAMINATION, Strand::Forward))
        );
        assert_eq!(
            filter.classify(G, A),
            Some((DamageClass::DEAMINATION, Strand::Reverse))
        );
        assert_eq!(
            filter.classify(G, T),
            Some((DamageClass::OXIDATION, Strand::Forward))
        );
        assert_eq!(
            filter.classify(C, A),
            Some((DamageClass::OXIDATION, Strand::Reverse))
        );
        assert_eq!(filter.classify(T, C), None);
    }

    #[test]
    fn test_distance_from_the_lesion_strand_five_prime_end() {
        let first = Molecule::new(C, 30, 0, 99);
        assert_eq!(
            CopiedDamage::five_prime_distance(&first, Strand::Forward),
            Some(0)
        );
        assert_eq!(
            CopiedDamage::five_prime_distance(&first, Strand::Reverse),
            Some(99)
        );
        let last = Molecule::new(C, 30, 99, 0);
        assert_eq!(
            CopiedDamage::five_prime_distance(&last, Strand::Reverse),
            Some(0)
        );
        let half = Molecule {
            right: None,
            ..first
        };
        assert_eq!(
            CopiedDamage::five_prime_distance(&half, Strand::Forward),
            None
        );
    }

    /// A deletion between the site and the lesion strand's 5' end brings the
    /// site nearer that end, in template bases, and within the distance.
    #[test]
    fn test_distance_counts_template_bases_across_a_deletion() {
        let mut builder = SamBuilder::new().read_length(50);
        for cigar1 in ["50M", "20M4D30M"] {
            builder.add_pair(Pair::at(101, 141).cigar1(cigar1));
        }
        let mut evidence =
            PileupEvidence::new(builder.to_pileup_builder(), &PileupOptions::default());
        let position = Position::try_from(132).unwrap();
        let molecules = evidence.molecules("chr1", position).unwrap();
        let distances: Vec<_> = molecules
            .iter()
            .map(|m| CopiedDamage::five_prime_distance(m, Strand::Forward).unwrap())
            .collect();
        assert_eq!(distances, vec![31, 27]);
    }

    #[test]
    fn test_congruent_molecules_are_within_the_scale() {
        let filter = CopiedDamage::default();
        let molecules: Vec<Molecule> = [0, 2, 3, 29, 30].map(|d| at(T, d)).into();
        let mut with_refs = spread_references(C);
        with_refs.extend(molecules);
        let count = |scale| {
            filter
                .score(&with_refs, C, T, Strand::Forward, scale)
                .alt_congruent
        };
        assert_eq!((count(30.0), count(2.5)), (4, 2));
    }

    #[test]
    fn test_alternates_near_the_lesion_five_prime_end_favor_the_artifact() {
        let filter = CopiedDamage::default();
        let mut molecules = spread_references(C);
        for d in [2, 5, 9, 14, 20] {
            molecules.push(at(T, d));
        }
        let score = filter.score(&molecules, C, T, Strand::Forward, 30.0);
        assert!(score.log_likelihood_ratio.unwrap() > 3.0, "{score:?}");
        assert_eq!((score.alt_congruent, score.alt_molecules), (5, 5));
        assert_eq!((score.ref_congruent, score.ref_molecules), (30, 150));
    }

    #[test]
    fn test_the_same_alternates_favor_a_mutation_on_the_other_strand() {
        let filter = CopiedDamage::default();
        let mut molecules = spread_references(G);
        for d in [2, 5, 9, 14, 20] {
            molecules.push(at(A, d));
        }
        let score = filter.score(&molecules, G, A, Strand::Reverse, 30.0);
        assert!(score.log_likelihood_ratio.unwrap() < -3.0, "{score:?}");
        assert_eq!(score.alt_congruent, 0);
    }

    #[test]
    fn test_alternates_spread_like_the_references_favor_a_mutation() {
        let filter = CopiedDamage::default();
        let mut molecules = spread_references(C);
        for d in (0..150).step_by(10) {
            molecules.push(at(T, d));
        }
        let score = filter.score(&molecules, C, T, Strand::Forward, 30.0);
        assert!(score.log_likelihood_ratio.unwrap() < -5.0, "{score:?}");
        assert_eq!((score.alt_congruent, score.alt_molecules), (3, 15));
    }

    /// A site deep enough to outweigh the even spread its `W` is shrunk toward
    /// keeps its own skew.
    #[test]
    fn test_capture_skew_shared_by_both_alleles_is_absorbed() {
        let filter = CopiedDamage::default();
        let skewed = |base| -> Vec<Molecule> { (0..20).map(|d| at(base, d)).collect() };
        let mut molecules: Vec<Molecule> = skewed(C).into_iter().cycle().take(400).collect();
        molecules.extend(skewed(T).into_iter().step_by(4));
        let score = filter.score(&molecules, C, T, Strand::Forward, 30.0);
        assert!(score.log_likelihood_ratio.unwrap().abs() < 0.5, "{score:?}");
    }

    #[test]
    fn test_molecules_without_both_ends_are_not_measured() {
        let filter = CopiedDamage::default();
        let molecules = vec![
            Molecule {
                right: None,
                ..Molecule::new(T, 90, 10, 0)
            },
            Molecule::new(C, 90, 10, 100),
        ];
        let score = filter.score(&molecules, C, T, Strand::Forward, 30.0);
        assert_eq!(score.alt_molecules, 0);
        assert_eq!(score.ref_molecules, 1);
        assert_eq!(score.log_likelihood_ratio, None);
        let distances = filter.distances(&molecules, C, T, Strand::Forward);
        assert_eq!(distances.spans, vec![111]);
    }

    #[test]
    fn test_log10_ratio_keeps_four_significant_digits_of_either_sign() {
        let ln = std::f64::consts::LN_10;
        assert_eq!(CopiedDamage::log10_ratio(1.234_56 * ln), 1.235);
        assert_eq!(CopiedDamage::log10_ratio(-1.234_56 * ln), -1.235);
        assert_eq!(CopiedDamage::log10_ratio(0.012_345_6 * ln), 0.012_35);
        assert_eq!(CopiedDamage::log10_ratio(123.456 * ln), 123.5);
        assert_eq!(CopiedDamage::log10_ratio(0.0), 0.0);
    }

    #[test]
    fn test_damage_site_stratum() {
        let site = DamageSite {
            class: DamageClass::DEAMINATION,
            strand: Strand::Reverse,
            context: Context::CpG,
        };
        assert_eq!(site.stratum(), "C>T:CpG");
    }
}
