//! The parent alleles of a call: the bases its alternate allele may have
//! arisen from, read from the molecules at the site.
//!
//! Each filter compares a call's alternate molecules with the molecules of
//! the base the change happened on, and classes the change by that base. A
//! VCF gives the reference base as REF, so at a site where the sample is
//! homozygous for another germline allele a call is classed from a base the
//! sample lacks and has no reference molecules to compare with: a `G>T`
//! against the reference is a `C>T` on a sample homozygous for `C` there, at
//! a CpG when the next reference base is `G`. Nothing need be phased. Where
//! the molecules that do not hold the alternate base hold almost none of REF
//! and ten or more of one other base, that base is what the alternate arose
//! from, so the filters count its molecules as the reference molecules and
//! class the call from it.
//!
//! Where the sample carries REF and another base both, each on more than a
//! sliver of those molecules, as at a heterozygous germline site or under
//! allelic imbalance, the alternate arose from either, and which one is
//! unknown.
//! Each filter then scores the call under each parent, with that parent's
//! molecules as the reference molecules and the class the parent makes, and
//! [`combine`] folds the two into one posterior: an artifact arises on a
//! parent's molecules at its class's damage rate, which the artifact share
//! of the class's stratum gives as odds, and a mutation on any molecule at
//! one rate, so a parent's artifact hypothesis weighs its molecules times its
//! stratum's artifact odds. A parent whose change no filter scores can only
//! have given a mutation.

use crate::evidence::Molecule;
use crate::prior::ln_add_exp;

/// The INFO key of the base other than REF the alternate may have arisen
/// from, written in place of REF where the sample does not carry REF and
/// beside it where the sample carries both.
pub const INFO: &str = "PA";

/// The INFO key of the weight the copied-damage posterior gave the artifact
/// hypothesis of the base [`INFO`] holds, where the sample carries REF too.
pub const INFO_WEIGHT: &str = "CDPW";

/// The fewest molecules of a parent base other than REF, so a handful of
/// molecules at a thinly covered site cannot make one.
pub const MIN_PARENT_MOLECULES: usize = 10;

/// The share of the molecules that are not the alternate base at or below
/// which REF or another base drops out as a parent, leaving the other alone;
/// above it both are parents. Two percent sits well above a consensus's
/// error rate, so a stray molecule of either base makes no parent.
pub const MIN_ALLELE_SHARE: f64 = 0.02;

/// A base the alternate allele may have arisen from, with its molecules at
/// the site.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Parent {
    /// The base.
    pub base: u8,
    /// Its molecules at the site.
    pub molecules: u32,
}

/// The bases the alternate allele may have arisen from at a site.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Parents {
    /// REF alone, as the VCF has it.
    Reference(Parent),
    /// Another base in place of REF, which the sample does not carry.
    Other(Parent),
    /// REF and another base, both of which the sample carries, REF first.
    Both(Parent, Parent),
}

impl Parents {
    /// Each parent, REF first where the sample carries it.
    pub fn each(&self) -> Vec<Parent> {
        match self {
            Parents::Reference(parent) | Parents::Other(parent) => vec![*parent],
            Parents::Both(reference, other) => vec![*reference, *other],
        }
    }

    /// The base other than REF, if any: the value of [`INFO`].
    pub fn other(&self) -> Option<u8> {
        match self {
            Parents::Reference(_) => None,
            Parents::Other(other) | Parents::Both(_, other) => Some(other.base),
        }
    }

    /// The base standing in place of REF, where the sample does not carry
    /// REF.
    pub fn in_place_of_reference(&self) -> Option<u8> {
        match self {
            Parents::Other(other) => Some(other.base),
            _ => None,
        }
    }
}

/// The bases the alternate allele may have arisen from. Only `A`, `C`, `G`,
/// and `T` molecules count, and the base other than `alt_base` and `ref_base`
/// that the most of them hold is the candidate, with [`MIN_PARENT_MOLECULES`]
/// or more of them. Of the molecules not holding `alt_base`, `ref_base` on at
/// most [`MIN_ALLELE_SHARE`] leaves the candidate the one parent, in place of
/// `ref_base`; the candidate on at most that share leaves `ref_base` the one
/// parent; and both above it are the two parents. A third base takes no part.
pub fn parents(molecules: &[Molecule], ref_base: u8, alt_base: u8) -> Parents {
    const BASES: [u8; 4] = *b"ACGT";
    let count = |base: u8| molecules.iter().filter(|m| m.base == base).count();
    let counts = BASES.map(count);
    let others = BASES
        .iter()
        .zip(counts)
        .filter(|(&base, _)| base != alt_base)
        .map(|(_, n)| n)
        .sum::<usize>() as f64;
    let reference = Parent {
        base: ref_base,
        molecules: count(ref_base) as u32,
    };
    let candidate = BASES
        .iter()
        .zip(counts)
        .filter(|(&base, _)| base != alt_base && base != ref_base)
        .max_by_key(|(_, n)| *n)
        .filter(|(_, n)| *n >= MIN_PARENT_MOLECULES);
    let Some((&base, n)) = candidate else {
        return Parents::Reference(reference);
    };
    let other = Parent {
        base,
        molecules: n as u32,
    };
    let (n, r) = (n as f64, f64::from(reference.molecules));
    if r <= MIN_ALLELE_SHARE * others {
        Parents::Other(other)
    } else if n <= MIN_ALLELE_SHARE * others {
        Parents::Reference(reference)
    } else {
        Parents::Both(reference, other)
    }
}

/// The base standing in place of `ref_base`, where the sample does not carry
/// it: see [`parents`].
pub fn parent_allele(molecules: &[Molecule], ref_base: u8, alt_base: u8) -> Option<u8> {
    parents(molecules, ref_base, alt_base).in_place_of_reference()
}

/// One parent's artifact hypothesis as [`combine`] weighs it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Weighed {
    /// The parent's share of the parents' molecules.
    pub share: f64,
    /// The artifact prior of the parent's stratum, the share of its calls
    /// that are artifacts, or 0 where the filter scores no such change or
    /// chance explains none of the stratum's positions.
    pub prior: f64,
    /// The hypothesis's natural-log likelihood ratio, artifact to mutation,
    /// where it was scored.
    pub ratio: Option<f64>,
}

/// What [`combine`] gives one call.
#[derive(Clone, Debug, PartialEq)]
pub struct Combined {
    /// The natural-log likelihood ratio, artifact to mutation, of the
    /// parents' artifact hypotheses averaged by their weights.
    pub ratio: f64,
    /// The artifact prior: the parents' artifact odds summed, each its share
    /// times its stratum's odds, as a probability.
    pub prior: f64,
    /// Each parent's weight: its share times its stratum's artifact odds over
    /// their sum, 0 for a parent without a scored artifact hypothesis.
    pub weights: Vec<f64>,
}

/// Fold the parents' hypotheses into one ratio and prior, which
/// [`crate::prior::posterior_mutation`] turns into the call's posterior.
///
/// Calls arise on parent `p`'s molecules as artifacts at its class's damage
/// rate `ρ_p` and as mutations at a rate `μ`, taken to be the same on every
/// base, so with `π_p` the artifact share of the calls in `p`'s stratum,
/// `ρ_p / μ = π_p / (1 - π_p)`. At a site whose parents hold shares `q_p` of
/// the molecules, an artifact from `p` arises at `q_p ρ_p` against mutations
/// at `μ` over all of them, so the prior odds of an artifact are
/// `A = Σ q_p π_p / (1 - π_p)`, and with `l_p` each hypothesis's likelihood
/// ratio the posterior odds are `Σ q_p (π_p / (1 - π_p)) exp(l_p)`, which a
/// prior `π = A / (1 + A)` and a ratio `l = ln Σ ω_p exp(l_p)`,
/// `ω_p = q_p π_p / ((1 - π_p) A)`, reproduce. A parent whose change the
/// filter does not score has `π_p = 0`: it adds mutations and no artifact. A
/// lone parent gives its own ratio and prior back, a prior of 1 makes the
/// call an artifact outright, a prior of 0 on every scored parent makes it a
/// mutation outright, and `None` is returned when no parent holds a scored
/// hypothesis.
pub fn combine(parents: &[Weighed]) -> Option<Combined> {
    let scored: Vec<(usize, f64, f64)> = parents
        .iter()
        .enumerate()
        .filter(|(_, p)| p.share > 0.0)
        .filter_map(|(i, p)| {
            p.ratio
                .map(|ratio| (i, p.share * p.prior / (1.0 - p.prior), ratio))
        })
        .collect();
    if scored.is_empty() {
        return None;
    }
    let mut weights = vec![0.0; parents.len()];
    if let Some((i, _, ratio)) = scored.iter().find(|(_, odds, _)| !odds.is_finite()) {
        weights[*i] = 1.0;
        return Some(Combined {
            ratio: *ratio,
            prior: 1.0,
            weights,
        });
    }
    if let [(i, _, ratio)] = scored.as_slice() {
        if parents[*i].share == 1.0 {
            weights[*i] = 1.0;
            return Some(Combined {
                ratio: *ratio,
                prior: parents[*i].prior,
                weights,
            });
        }
    }
    let odds: f64 = scored.iter().map(|(_, odds, _)| odds).sum();
    if odds <= 0.0 {
        let (i, _, ratio) = scored[0];
        weights[i] = 1.0;
        return Some(Combined {
            ratio,
            prior: 0.0,
            weights,
        });
    }
    let ln_sum = scored
        .iter()
        .map(|(_, odds, ratio)| odds.ln() + ratio)
        .fold(f64::NEG_INFINITY, ln_add_exp);
    for (i, each, _) in &scored {
        weights[*i] = each / odds;
    }
    Some(Combined {
        ratio: ln_sum - odds.ln(),
        prior: odds / (1.0 + odds),
        weights,
    })
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::evidence::DELETION;
    use crate::prior::{logit, posterior_mutation};

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

    fn parent(base: u8, molecules: u32) -> Parent {
        Parent { base, molecules }
    }

    /// A sample homozygous for `C` where the reference has `G` gives `C` as
    /// the parent of a `T`, with or without a stray reference molecule, until
    /// REF passes 2% of the molecules that are not alternate.
    #[rstest]
    #[case(&[(C, 183), (T, 2)], Some(C))]
    #[case(&[(C, 183), (G, 1), (T, 2)], Some(C))]
    #[case(&[(C, 183), (G, 3), (T, 2)], Some(C))]
    #[case(&[(C, 183), (G, 4), (T, 2)], None)]
    fn test_a_sample_without_ref_gives_the_base_it_carries(
        #[case] bases: &[(u8, usize)],
        #[case] expected: Option<u8>,
    ) {
        assert_eq!(parent_allele(&molecules(bases), G, T), expected);
    }

    /// A site carrying REF and another base, each above 2% of the molecules
    /// that are not alternate, gives both, REF first, so nothing stands in
    /// place of REF; a third base takes no part.
    #[rstest]
    #[case(&[(C, 90), (G, 90), (T, 2)], parent(G, 90), parent(C, 90))]
    #[case(&[(C, 120), (G, 60), (T, 2)], parent(G, 60), parent(C, 120))]
    #[case(&[(C, 20), (G, 80), (T, 2)], parent(G, 80), parent(C, 20))]
    #[case(&[(C, 80), (G, 20), (T, 2)], parent(G, 20), parent(C, 80))]
    #[case(&[(C, 11), (G, 489), (T, 2)], parent(G, 489), parent(C, 11))]
    #[case(&[(C, 489), (G, 11), (T, 2)], parent(G, 11), parent(C, 489))]
    #[case(&[(C, 183), (G, 4), (T, 2)], parent(G, 4), parent(C, 183))]
    #[case(&[(C, 40), (G, 40), (A, 20), (T, 2)], parent(G, 40), parent(C, 40))]
    fn test_a_site_with_both_bases_gives_both(
        #[case] bases: &[(u8, usize)],
        #[case] reference: Parent,
        #[case] other: Parent,
    ) {
        let found = parents(&molecules(bases), G, T);
        assert_eq!(found, Parents::Both(reference, other));
        assert_eq!(found.each(), [reference, other]);
        assert_eq!(found.other(), Some(other.base));
        assert_eq!(found.in_place_of_reference(), None);
        assert_eq!(parent_allele(&molecules(bases), G, T), None);
    }

    /// At or under 2% of the molecules that are not alternate a base drops
    /// out: another base leaves REF alone, and REF leaves the other base
    /// alone in its place, whatever a third base holds.
    #[rstest]
    #[case(&[(C, 10), (G, 490), (T, 2)], Parents::Reference(parent(G, 490)))]
    #[case(&[(C, 490), (G, 10), (T, 2)], Parents::Other(parent(C, 490)))]
    #[case(&[(C, 60), (A, 40), (T, 2)], Parents::Other(parent(C, 60)))]
    fn test_a_base_at_the_floor_drops_out(
        #[case] bases: &[(u8, usize)],
        #[case] expected: Parents,
    ) {
        assert_eq!(parents(&molecules(bases), G, T), expected);
    }

    /// A site homozygous for REF keeps it, whatever else is there.
    #[test]
    fn test_a_site_homozygous_for_ref_keeps_it() {
        assert_eq!(
            parents(&molecules(&[(G, 183), (T, 2)]), G, T),
            Parents::Reference(parent(G, 183))
        );
        assert_eq!(
            parents(&molecules(&[(G, 183), (C, 1), (T, 2)]), G, T),
            Parents::Reference(parent(G, 183))
        );
    }

    #[rstest]
    #[case(&[(C, 9), (T, 2)], None)]
    #[case(&[(C, 10), (T, 2)], Some(C))]
    #[case(&[(C, 9), (G, 11), (T, 2)], None)]
    fn test_a_parent_needs_enough_molecules(
        #[case] bases: &[(u8, usize)],
        #[case] expected: Option<u8>,
    ) {
        assert_eq!(parents(&molecules(bases), G, T).other(), expected);
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
    fn test_no_other_molecules_give_no_parent() {
        assert_eq!(
            parents(&molecules(&[(T, 100)]), G, T),
            Parents::Reference(parent(G, 0))
        );
        assert_eq!(parents(&[], G, T), Parents::Reference(parent(G, 0)));
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

    fn weighed(share: f64, prior: f64, ratio: Option<f64>) -> Weighed {
        Weighed {
            share,
            prior,
            ratio,
        }
    }

    /// A lone parent gives its own ratio and prior back, so the posterior is
    /// the one its filter would write alone.
    #[test]
    fn test_a_lone_parent_combines_to_itself() {
        let combined = combine(&[weighed(1.0, 0.3, Some(2.5))]).unwrap();
        assert_eq!((combined.ratio, combined.prior), (2.5, 0.3));
        assert_eq!(combined.weights, [1.0]);
        assert_eq!(
            posterior_mutation(combined.ratio, combined.prior),
            posterior_mutation(2.5, 0.3)
        );
    }

    /// A parent with no molecules drops out: all REF gives REF's result,
    /// and no REF gives the other parent's.
    #[rstest]
    #[case(1.0, 0.0, 2.0, 0.3)]
    #[case(0.0, 1.0, -1.0, 0.6)]
    fn test_a_parent_without_molecules_drops_out(
        #[case] share_ref: f64,
        #[case] share_other: f64,
        #[case] ratio: f64,
        #[case] prior: f64,
    ) {
        let parents = [
            weighed(share_ref, 0.3, Some(2.0)),
            weighed(share_other, 0.6, Some(-1.0)),
        ];
        let combined = combine(&parents).unwrap();
        assert!((combined.ratio - ratio).abs() < 1e-12, "{combined:?}");
        assert!((combined.prior - prior).abs() < 1e-12, "{combined:?}");
        let expected = if share_ref > 0.0 {
            [1.0, 0.0]
        } else {
            [0.0, 1.0]
        };
        assert_eq!(combined.weights, expected);
    }

    /// Two parents with equal molecules and evidence sum their artifact
    /// odds, and the one whose class is far more often damage takes nearly
    /// all of the artifact hypothesis.
    #[test]
    fn test_the_high_rate_parent_dominates() {
        let parents = [weighed(0.5, 0.005, Some(3.0)), weighed(0.5, 0.5, Some(3.0))];
        let combined = combine(&parents).unwrap();
        let odds = 0.5 * 0.005 / 0.995 + 0.5;
        assert!(
            (combined.prior - odds / (1.0 + odds)).abs() < 1e-12,
            "{combined:?}"
        );
        assert!((combined.ratio - 3.0).abs() < 1e-12, "{combined:?}");
        assert!(combined.weights[1] > 0.99, "{combined:?}");
        assert!((combined.weights.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    /// The posterior odds of an artifact are the sum of each parent's
    /// artifact odds, its share times its stratum's odds times its likelihood
    /// ratio, against mutation odds of one.
    #[test]
    fn test_the_posterior_odds_sum_the_parents_artifact_odds() {
        let parents = [weighed(0.6, 0.2, Some(1.0)), weighed(0.4, 0.5, Some(-0.5))];
        let combined = combine(&parents).unwrap();
        let artifact = 0.6 * (0.2 / 0.8) * 1f64.exp() + 0.4 * (0.5 / 0.5) * (-0.5f64).exp();
        let posterior = posterior_mutation(combined.ratio, combined.prior);
        assert!((posterior - 1.0 / (artifact + 1.0)).abs() < 1e-12);
    }

    /// A parent whose change the filter does not score can only have given a
    /// mutation: it takes no weight, and at equal shares it halves the odds
    /// of the artifact, since mutations arise on its molecules as on the
    /// other parent's.
    #[test]
    fn test_an_unscored_parent_only_dilutes_the_prior() {
        let parents = [weighed(0.5, 0.0, None), weighed(0.5, 0.4, Some(2.0))];
        let combined = combine(&parents).unwrap();
        assert!((combined.prior - 0.25).abs() < 1e-12, "{combined:?}");
        assert!((logit(combined.prior) - (0.5 * 0.4f64 / 0.6).ln()).abs() < 1e-12);
        assert!((combined.ratio - 2.0).abs() < 1e-12);
        assert_eq!(combined.weights, [0.0, 1.0]);
        assert_eq!(combine(&[weighed(1.0, 0.0, None)]), None);
        assert_eq!(combine(&[weighed(1.0, 0.4, None)]), None);
    }

    /// A scored parent whose stratum holds no artifact at all makes the call
    /// a mutation outright, with its posterior of one written, as a lone
    /// parent at a prior of 0 always was.
    #[test]
    fn test_a_prior_of_zero_makes_the_call_a_mutation() {
        let alone = combine(&[weighed(1.0, 0.0, Some(5.0))]).unwrap();
        assert_eq!((alone.ratio, alone.prior), (5.0, 0.0));
        assert_eq!(alone.weights, [1.0]);
        assert_eq!(posterior_mutation(alone.ratio, alone.prior), 1.0);
        let both = combine(&[weighed(0.5, 0.0, Some(5.0)), weighed(0.5, 0.0, Some(2.0))]).unwrap();
        assert_eq!(both.prior, 0.0);
        assert_eq!(both.weights, [1.0, 0.0]);
        assert_eq!(posterior_mutation(both.ratio, both.prior), 1.0);
    }

    /// A stratum whose every call is an artifact makes the call one, as it
    /// does for a lone parent, whatever the other parent says.
    #[test]
    fn test_a_prior_of_one_makes_the_call_an_artifact() {
        let parents = [weighed(0.5, 0.3, Some(1.0)), weighed(0.5, 1.0, Some(0.0))];
        let combined = combine(&parents).unwrap();
        assert_eq!(combined.prior, 1.0);
        assert_eq!(combined.weights, [0.0, 1.0]);
        assert_eq!(posterior_mutation(combined.ratio, combined.prior), 0.0);
    }
}
