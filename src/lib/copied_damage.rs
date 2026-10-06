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
//! into a copy with probability `w(d) = exp(-d / scale)`, an exponentially
//! distributed resynthesis length with mean `scale`; the alternate distances
//! follow the reference distances tilted by `w`. See
//! [`crate::read_end::tilt_log_likelihood_ratio`] for the per-molecule terms.

use crate::call::Genotype;
use crate::classes::{Context, DamageClass, Strand};
use crate::evidence::Molecule;
use crate::read_end::{tilt_log_likelihood_ratio, Score};

/// The copied damage filter.
#[derive(Clone, Debug, PartialEq)]
pub struct CopiedDamage {
    /// The damage classes to test.
    pub classes: Vec<DamageClass>,
    /// The mean resynthesis length in bases.
    pub scale: f64,
}

impl Default for CopiedDamage {
    fn default() -> Self {
        Self {
            classes: vec![DamageClass::DEAMINATION, DamageClass::OXIDATION],
            scale: 30.0,
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
    /// The INFO key for the alternate molecules nearer the lesion strand's 5'
    /// end, and all alternate molecules measured.
    pub const INFO_ALT: &'static str = "CDAC";
    /// The INFO key for the reference molecules nearer the lesion strand's 5'
    /// end, and all reference molecules measured.
    pub const INFO_REF: &'static str = "CDRC";
    /// The FILTER name.
    pub const FILTER: &'static str = "CopiedDamageArtifact";

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

    /// The 0-based distance of the site from the lesion strand's 5' end and the
    /// template length, in the template's bases, for a molecule with both ends
    /// known.
    pub fn distance(m: &Molecule, strand: Strand) -> Option<(usize, usize)> {
        let d = match strand {
            Strand::Forward => m.left?,
            Strand::Reverse => m.right?,
        };
        Some((d, m.length()?))
    }

    /// Whether a distance is nearer the lesion strand's 5' end than its 3' end.
    pub fn is_five_prime_proximal(d: usize, length: usize) -> bool {
        2 * d + 1 < length
    }

    /// The call's score: the tilt likelihood ratio and the 5'-proximal counts.
    pub fn score(
        &self,
        molecules: &[Molecule],
        ref_base: u8,
        alt_base: u8,
        strand: Strand,
    ) -> Score {
        let mut score = Score {
            depth: molecules.len() as u32,
            ..Score::default()
        };
        let mut ref_distances = Vec::new();
        let mut alt = Vec::new();
        for m in molecules {
            let Some((d, length)) = Self::distance(m, strand) else {
                continue;
            };
            let proximal = Self::is_five_prime_proximal(d, length);
            if m.base == ref_base {
                score.ref_molecules += 1;
                score.ref_congruent += u32::from(proximal);
                ref_distances.push(d);
            } else if m.base == alt_base {
                score.alt_molecules += 1;
                score.alt_congruent += u32::from(proximal);
                alt.push((d, m.quality));
            }
        }
        score.log_likelihood_ratio = tilt_log_likelihood_ratio(&ref_distances, &alt, self.scale);
        score
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
            CopiedDamage::distance(&first, Strand::Forward),
            Some((0, 100))
        );
        assert_eq!(
            CopiedDamage::distance(&first, Strand::Reverse),
            Some((99, 100))
        );
        let last = Molecule::new(C, 30, 99, 0);
        assert_eq!(
            CopiedDamage::distance(&last, Strand::Reverse),
            Some((0, 100))
        );
        let half = Molecule {
            right: None,
            ..first
        };
        assert_eq!(CopiedDamage::distance(&half, Strand::Forward), None);
    }

    /// A deletion between the site and the lesion strand's 5' end brings the
    /// site nearer that end, in template bases, past the template's midpoint.
    #[test]
    fn test_distance_counts_template_bases_across_a_deletion() {
        let mut builder = SamBuilder::new().read_length(50);
        for cigar1 in ["50M", "20M4D30M"] {
            builder.add_pair(Pair {
                cigar1: Some(cigar1.into()),
                ..Pair::at(101, 141)
            });
        }
        let mut evidence =
            PileupEvidence::new(builder.to_pileup_builder(), &PileupOptions::default());
        let position = Position::try_from(146).unwrap();
        let molecules = evidence.molecules("chr1", position).unwrap();
        let distances: Vec<_> = molecules
            .iter()
            .map(|m| CopiedDamage::distance(m, Strand::Forward).unwrap())
            .collect();
        assert_eq!(distances, vec![(45, 90), (41, 86)]);
        assert!(!CopiedDamage::is_five_prime_proximal(45, 90));
        assert!(CopiedDamage::is_five_prime_proximal(41, 86));
    }

    #[test]
    fn test_five_prime_proximal_splits_the_molecule_in_half() {
        assert!(CopiedDamage::is_five_prime_proximal(0, 100));
        assert!(CopiedDamage::is_five_prime_proximal(49, 100));
        assert!(!CopiedDamage::is_five_prime_proximal(50, 100));
        assert!(!CopiedDamage::is_five_prime_proximal(2, 5));
        assert!(CopiedDamage::is_five_prime_proximal(1, 5));
    }

    #[test]
    fn test_alternates_near_the_lesion_five_prime_end_favor_the_artifact() {
        let filter = CopiedDamage::default();
        let mut molecules = spread_references(C);
        for d in [2, 5, 9, 14, 20] {
            molecules.push(at(T, d));
        }
        let score = filter.score(&molecules, C, T, Strand::Forward);
        assert!(score.log_likelihood_ratio.unwrap() > 3.0, "{score:?}");
        assert_eq!((score.alt_congruent, score.alt_molecules), (5, 5));
        assert_eq!((score.ref_congruent, score.ref_molecules), (75, 150));
    }

    #[test]
    fn test_the_same_alternates_favor_a_mutation_on_the_other_strand() {
        let filter = CopiedDamage::default();
        let mut molecules = spread_references(G);
        for d in [2, 5, 9, 14, 20] {
            molecules.push(at(A, d));
        }
        let score = filter.score(&molecules, G, A, Strand::Reverse);
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
        let score = filter.score(&molecules, C, T, Strand::Forward);
        assert!(score.log_likelihood_ratio.unwrap() < -5.0, "{score:?}");
        assert_eq!((score.alt_congruent, score.alt_molecules), (8, 15));
    }

    #[test]
    fn test_capture_skew_shared_by_both_alleles_is_absorbed() {
        let filter = CopiedDamage::default();
        let skewed = |base| -> Vec<Molecule> { (0..20).map(|d| at(base, d)).collect() };
        let mut molecules = skewed(C);
        molecules.extend(skewed(T).into_iter().step_by(4));
        let score = filter.score(&molecules, C, T, Strand::Forward);
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
        let score = filter.score(&molecules, C, T, Strand::Forward);
        assert_eq!(score.alt_molecules, 0);
        assert_eq!(score.ref_molecules, 1);
        assert_eq!(score.log_likelihood_ratio, Some(0.0));
        assert_eq!(score.depth, 2);
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
