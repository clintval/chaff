//! Priors on whether a call is an artifact, and the posterior they give.
//!
//! The default prior is learned. Each call `i` in a stratum carries a
//! likelihood ratio `exp(l_i) = P(molecules | artifact) / P(molecules |
//! mutation)`. The calls form a two-component mixture with an unknown artifact
//! fraction `pi`, which expectation-maximization estimates per sample and per
//! stratum, the way GATK's LearnReadOrientationModel learns its artifact
//! priors:
//!
//! ```text
//! E-step: r_i = 1 / (1 + exp(-(l_i + logit(pi))))
//! M-step: pi  = (sum_i r_i + c) / (n + 2c)
//! ```
//!
//! The pseudocount `c` is a Beta(c + 1, c + 1) prior on `pi`, so `pi` is the
//! maximum a posteriori estimate and stays strictly between 0 and 1 when a
//! stratum holds few calls. The objective is concave in `pi`, so the fixed
//! point is unique.
//!
//! fgbio's prior is kept for parity: a mutation prior of `min((2 * maf)^2,
//! 0.9999)`, where `maf` is the call's alternate molecule fraction, or one over
//! the depth when no alternate molecule is seen. At the low allele fractions of
//! duplex sequencing it is near zero, which makes any call whose alternate
//! molecules all sit inside the window an artifact.

use std::fmt;

use clap::ValueEnum;

/// The pseudocount on each side of the learned artifact fraction.
pub const PSEUDOCOUNT: f64 = 1.0;

/// Which prior turns a likelihood ratio into a posterior.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum PriorMode {
    /// Learn the artifact fraction per sample and stratum by EM.
    #[default]
    Learned,
    /// fgbio's per-call `(2 * maf)^2` mutation prior, for parity.
    Fgbio,
}

impl fmt::Display for PriorMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PriorMode::Learned => write!(f, "learned"),
            PriorMode::Fgbio => write!(f, "fgbio"),
        }
    }
}

/// The logistic function, evaluated without overflow.
pub fn sigmoid(x: f64) -> f64 {
    if x >= 0.0 {
        1.0 / (1.0 + (-x).exp())
    } else {
        let e = x.exp();
        e / (1.0 + e)
    }
}

/// The log odds of a probability.
pub fn logit(p: f64) -> f64 {
    (p / (1.0 - p)).ln()
}

/// The posterior probability that a call is a true mutation, given its
/// artifact-to-mutation log likelihood ratio and the artifact prior.
pub fn posterior_mutation(log_likelihood_ratio: f64, artifact_prior: f64) -> f64 {
    sigmoid(-(log_likelihood_ratio + logit(artifact_prior)))
}

/// fgbio's artifact prior for one call: one minus `min((2 * m)^2, 0.9999)`,
/// where `m` is the alternate molecule fraction among reference and alternate
/// molecules, or `1 / depth` when it is zero.
pub fn fgbio_artifact_prior(alt_molecules: u32, ref_molecules: u32, depth: u32) -> f64 {
    let total = alt_molecules + ref_molecules;
    let maf = if total > 0 {
        f64::from(alt_molecules) / f64::from(total)
    } else {
        0.0
    };
    let m = if maf != 0.0 {
        maf
    } else {
        1.0 / f64::from(depth)
    };
    let prior_mutation = (m * 2.0).powi(2).min(0.9999);
    1.0 - prior_mutation
}

/// The maximum a posteriori artifact fraction of a stratum, by EM over the
/// calls' log likelihood ratios. An empty stratum gets the prior mean, one half.
pub fn learn_artifact_fraction(log_likelihood_ratios: &[f64], pseudocount: f64) -> f64 {
    let n = log_likelihood_ratios.len() as f64;
    let mut pi = 0.5;
    for _ in 0..10_000 {
        let odds = logit(pi);
        let responsibility: f64 = log_likelihood_ratios
            .iter()
            .map(|l| sigmoid(l + odds))
            .sum();
        let next = (responsibility + pseudocount) / (n + 2.0 * pseudocount);
        let converged = (next - pi).abs() < 1e-12;
        pi = next;
        if converged {
            break;
        }
    }
    pi
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tolerance: f64) -> bool {
        (a - b).abs() <= tolerance
    }

    #[test]
    fn test_sigmoid_is_stable_at_extremes() {
        assert_eq!(sigmoid(0.0), 0.5);
        assert!(close(sigmoid(1000.0), 1.0, 0.0));
        assert!(close(sigmoid(-1000.0), 0.0, 0.0));
        assert!(sigmoid(-1000.0).is_finite());
    }

    #[test]
    fn test_posterior_mutation_is_the_prior_without_evidence() {
        assert!(close(posterior_mutation(0.0, 0.25), 0.75, 1e-12));
        assert!(posterior_mutation(10.0, 0.5) < 1e-4);
        assert!(posterior_mutation(-10.0, 0.5) > 1.0 - 1e-4);
    }

    #[test]
    fn test_fgbio_artifact_prior() {
        assert!(close(fgbio_artifact_prior(1, 99, 100), 1.0 - 0.0004, 1e-12));
        assert!(close(
            fgbio_artifact_prior(50, 50, 100),
            1.0 - 0.9999,
            1e-12
        ));
        assert!(close(
            fgbio_artifact_prior(0, 400, 400),
            1.0 - 0.000025,
            1e-12
        ));
        assert!(close(fgbio_artifact_prior(0, 0, 0), 1.0 - 0.9999, 1e-12));
    }

    #[test]
    fn test_learned_fraction_of_an_empty_stratum_is_one_half() {
        assert_eq!(learn_artifact_fraction(&[], PSEUDOCOUNT), 0.5);
    }

    #[test]
    fn test_learned_fraction_ignores_uninformative_calls() {
        let informative = [20.0, 20.0, -20.0, -20.0, -20.0, -20.0];
        let mut with_flat = informative.to_vec();
        with_flat.extend([0.0; 50]);
        let a = learn_artifact_fraction(&informative, PSEUDOCOUNT);
        let b = learn_artifact_fraction(&with_flat, PSEUDOCOUNT);
        assert!(close(a, 3.0 / 8.0, 1e-6), "{a}");
        assert!(close(a, b, 1e-6), "{a} vs {b}");
    }

    #[test]
    fn test_learned_fraction_recovers_a_simulated_mixture() {
        let mut llrs = Vec::new();
        for i in 0..1000 {
            let artifact = i % 5 == 0;
            let strength = 2.0 + f64::from(i % 7) * 0.5;
            llrs.push(if artifact { strength } else { -strength });
        }
        let pi = learn_artifact_fraction(&llrs, PSEUDOCOUNT);
        assert!(close(pi, 0.2, 0.05), "{pi}");
    }

    #[test]
    fn test_learned_fraction_of_one_strong_artifact_call() {
        let pi = learn_artifact_fraction(&[50.0], PSEUDOCOUNT);
        assert!(close(pi, 2.0 / 3.0, 1e-9), "{pi}");
    }
}
