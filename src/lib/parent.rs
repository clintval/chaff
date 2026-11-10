//! The parent allele of a call: the base its alternate allele arose from,
//! where the sample does not carry the reference base.
//!
//! Each filter compares a call's alternate molecules with the molecules of
//! the base the change happened on, and classes the change by that base. A
//! VCF names the reference base as REF, so at a site where the sample is
//! homozygous for another germline allele a call is classed from a base the
//! sample lacks and has no reference molecules to compare with: a `G>T`
//! against the reference is a `C>T` on a sample homozygous for `C` there, at
//! a CpG when the next reference base is `G`. Nothing need be phased. Where
//! nearly every molecule that does not hold the alternate base holds one
//! other base, that base is what the alternate arose from, so the filters
//! count its molecules as the reference molecules and class the call from it.
//! A heterozygous site, where REF holds about half of those molecules, keeps
//! REF, since the alternate's parent there is unknown.

use crate::evidence::Molecule;

/// The INFO key of the parent allele, written on a call scored from a base
/// other than REF.
pub const INFO: &str = "PA";

/// The fewest molecules of the parent base, so a handful of molecules at a
/// thinly covered site cannot name one.
pub const MIN_PARENT_MOLECULES: usize = 10;

/// The most reference molecules, as a share of the molecules that are not
/// the alternate base: well above a consensus's error rate at a site the
/// reference allele is absent from, and far below the half a heterozygous
/// site holds.
pub const MAX_REFERENCE_SHARE: f64 = 0.02;

/// The least share of the molecules that are not the alternate base the
/// parent must hold, which leaves room for a few mismeasured bases.
pub const MIN_PARENT_SHARE: f64 = 0.95;

/// The base the alternate allele arose from when it is not `ref_base`: the
/// one base, other than `alt_base`, that at least [`MIN_PARENT_SHARE`] of the
/// molecules not holding `alt_base` hold, [`MIN_PARENT_MOLECULES`] or more of
/// them, with `ref_base` on at most [`MAX_REFERENCE_SHARE`] of them. `None`
/// keeps `ref_base`. Only `A`, `C`, `G`, and `T` molecules count.
pub fn parent_allele(molecules: &[Molecule], ref_base: u8, alt_base: u8) -> Option<u8> {
    const BASES: [u8; 4] = *b"ACGT";
    let count = |base: u8| molecules.iter().filter(|m| m.base == base).count();
    let counts = BASES.map(count);
    let others = BASES
        .iter()
        .zip(counts)
        .filter(|(&base, _)| base != alt_base)
        .map(|(_, n)| n)
        .sum::<usize>() as f64;
    let (&parent, n) = BASES
        .iter()
        .zip(counts)
        .filter(|(&base, _)| base != alt_base && base != ref_base)
        .max_by_key(|(_, n)| *n)?;
    let enough = n >= MIN_PARENT_MOLECULES
        && n as f64 >= MIN_PARENT_SHARE * others
        && count(ref_base) as f64 <= MAX_REFERENCE_SHARE * others;
    enough.then_some(parent)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::evidence::DELETION;

    const A: u8 = b'A';
    const C: u8 = b'C';
    const G: u8 = b'G';
    const T: u8 = b'T';

    /// Molecules holding each base, `n` of them, 150 bp long.
    fn molecules(bases: &[(u8, usize)]) -> Vec<Molecule> {
        bases
            .iter()
            .flat_map(|&(base, n)| {
                (0..n).map(move |d| Molecule::new(base, 90, d % 150, 149 - d % 150))
            })
            .collect()
    }

    /// A sample homozygous for `C` where the reference has `G` names `C` as
    /// the parent of a `T`, with or without a stray reference molecule.
    #[rstest]
    #[case(&[(C, 183), (T, 2)], Some(C))]
    #[case(&[(C, 183), (G, 1), (T, 2)], Some(C))]
    #[case(&[(C, 183), (G, 3), (T, 2)], Some(C))]
    #[case(&[(C, 183), (G, 4), (T, 2)], None)]
    fn test_a_sample_without_ref_names_the_base_it_carries(
        #[case] bases: &[(u8, usize)],
        #[case] expected: Option<u8>,
    ) {
        assert_eq!(parent_allele(&molecules(bases), G, T), expected);
    }

    /// A heterozygous site keeps REF: the alternate may have arisen from
    /// either allele.
    #[rstest]
    #[case(&[(C, 90), (G, 90), (T, 2)])]
    #[case(&[(C, 120), (G, 60), (T, 2)])]
    #[case(&[(C, 5), (G, 5), (T, 2)])]
    fn test_a_heterozygous_site_keeps_ref(#[case] bases: &[(u8, usize)]) {
        assert_eq!(parent_allele(&molecules(bases), G, T), None);
    }

    /// A site homozygous for REF keeps it, whatever else is there.
    #[test]
    fn test_a_site_homozygous_for_ref_keeps_it() {
        assert_eq!(parent_allele(&molecules(&[(G, 183), (T, 2)]), G, T), None);
        assert_eq!(
            parent_allele(&molecules(&[(G, 183), (C, 1), (T, 2)]), G, T),
            None
        );
    }

    #[rstest]
    #[case(9, None)]
    #[case(10, Some(C))]
    fn test_the_parent_needs_enough_molecules(#[case] n: usize, #[case] expected: Option<u8>) {
        assert_eq!(parent_allele(&molecules(&[(C, n), (T, 2)]), G, T), expected);
    }

    /// A few molecules of a third base are mismeasured bases; many are
    /// another allele.
    #[rstest]
    #[case(&[(C, 95), (A, 5), (T, 2)], Some(C))]
    #[case(&[(C, 94), (A, 6), (T, 2)], None)]
    fn test_a_third_base_within_the_share_is_tolerated(
        #[case] bases: &[(u8, usize)],
        #[case] expected: Option<u8>,
    ) {
        assert_eq!(parent_allele(&molecules(bases), G, T), expected);
    }

    /// Disagreeing mates and deletions are neither allele, so they count
    /// toward no share.
    #[test]
    fn test_n_and_deletions_do_not_count() {
        let bases = [(C, 20), (b'N', 50), (DELETION, 50), (T, 2)];
        assert_eq!(parent_allele(&molecules(&bases), G, T), Some(C));
    }

    /// With every molecule alternate, or none at all, no base can be the
    /// parent.
    #[test]
    fn test_no_other_molecules_name_no_parent() {
        assert_eq!(parent_allele(&molecules(&[(T, 100)]), G, T), None);
        assert_eq!(parent_allele(&[], G, T), None);
    }

    /// The parent is read from the molecules, whatever REF and ALT are, so a
    /// reference `A` with a sample homozygous for `G` gives a `G>A` to a `T`
    /// alternate its `G` parent.
    #[test]
    fn test_the_parent_can_be_any_base_but_ref_and_alt() {
        assert_eq!(
            parent_allele(&molecules(&[(G, 300), (A, 2), (T, 2)]), A, T),
            Some(G)
        );
    }
}
