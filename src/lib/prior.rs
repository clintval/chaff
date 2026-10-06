//! Priors on whether a call is an artifact, and the posterior they give.
//!
//! The default prior is learned. Each call `i` carries a likelihood ratio
//! `exp(l_i) = P(molecules | artifact) / P(molecules | mutation)`. A set of
//! calls forms a two-component mixture with an unknown artifact fraction
//! `pi`, estimated per sample the way GATK's LearnReadOrientationModel learns
//! its artifact priors: as the fixed point of expectation-maximization,
//!
//! ```text
//! E-step: r_i = 1 / (1 + exp(-(l_i + logit(pi))))
//! M-step: pi  = (sum_i r_i + k m) / (n + k)
//! ```
//!
//! which is the maximum a posteriori `pi` under a [`BetaPrior`] of `k`
//! pseudo-calls at mean `m`, `Beta(k m + 1, k (1 - m) + 1)`. The objective is
//! concave in `pi`, so the fixed point is unique, and chaff solves for it
//! directly rather than iterating.
//!
//! The prior is learned twice. Each filter first learns one fraction from all
//! of its calls under the weak [`FILTER_PRIOR`], `Beta(2, 2)`. Each stratum
//! then learns its own under [`STRATUM_PRIOR_STRENGTH`] pseudo-calls at its
//! filter's fraction, so a stratum of one or two calls mostly inherits the
//! filter's fraction rather than moving its prior toward its own calls, and a
//! stratum of hundreds keeps nearly its own.
//!
//! fgbio's prior is kept for parity: a mutation prior of `min((2 * maf)^2,
//! 0.9999)`, where `maf` is the call's alternate molecule fraction, or one over
//! the depth when no alternate molecule is seen. At the low allele fractions of
//! duplex sequencing it is near zero, which makes any call whose alternate
//! molecules all sit inside the window an artifact.

use std::fmt;

use clap::ValueEnum;

/// A Beta prior on an artifact fraction: `strength` pseudo-calls at `mean`
/// over a flat `Beta(1, 1)`, `Beta(strength mean + 1, strength (1 - mean) + 1)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BetaPrior {
    /// The fraction the pseudo-calls hold.
    pub mean: f64,
    /// The number of pseudo-calls.
    pub strength: f64,
}

/// The prior on a filter's fraction: two pseudo-calls at one half, `Beta(2,
/// 2)`, which keeps the fraction strictly between 0 and 1.
pub const FILTER_PRIOR: BetaPrior = BetaPrior {
    mean: 0.5,
    strength: 2.0,
};

/// The pseudo-calls a stratum's prior holds at its filter's fraction: ten, so
/// a stratum's own calls outweigh its filter's once it holds more than ten, and
/// a lone call moves its stratum at most one eleventh of the way to itself.
pub const STRATUM_PRIOR_STRENGTH: f64 = 10.0;

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

/// The maximum a posteriori artifact fraction of a set of calls from their log
/// likelihood ratios under `prior`. The fixed point EM converges to solves
/// `sum_i r_i + k m = pi (n + k)`, whose left side less its right falls as
/// `pi` rises, so bisection finds it to machine precision however slowly EM
/// would creep there. Without calls it is the prior's mean.
pub fn learn_artifact_fraction(log_likelihood_ratios: &[f64], prior: BetaPrior) -> f64 {
    let n = log_likelihood_ratios.len() as f64;
    let excess = |pi: f64| {
        let odds = logit(pi);
        let responsibility: f64 = log_likelihood_ratios
            .iter()
            .map(|l| sigmoid(l + odds))
            .sum();
        responsibility + prior.strength * prior.mean - pi * (n + prior.strength)
    };
    let (mut low, mut high) = (0.0f64, 1.0f64);
    loop {
        let pi = 0.5 * (low + high);
        if pi <= low || pi >= high {
            return pi;
        }
        if excess(pi) > 0.0 {
            low = pi;
        } else {
            high = pi;
        }
    }
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
        assert_eq!(learn_artifact_fraction(&[], FILTER_PRIOR), 0.5);
    }

    #[test]
    fn test_learned_fraction_ignores_uninformative_calls() {
        let informative = [20.0, 20.0, -20.0, -20.0, -20.0, -20.0];
        let mut with_flat = informative.to_vec();
        with_flat.extend([0.0; 50]);
        let a = learn_artifact_fraction(&informative, FILTER_PRIOR);
        let b = learn_artifact_fraction(&with_flat, FILTER_PRIOR);
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
        let pi = learn_artifact_fraction(&llrs, FILTER_PRIOR);
        assert!(close(pi, 0.2, 0.05), "{pi}");
    }

    /// Uninformative calls slow EM to a crawl: 100,000 calls at a ratio of 0
    /// and 10 at 20 have their fixed point at `(10 + 1) / (10 + 2)`.
    #[test]
    fn test_learned_fraction_reaches_its_fixed_point_among_many_flat_calls() {
        let mut llrs = vec![0.0; 100_000];
        llrs.extend([20.0; 10]);
        let pi = learn_artifact_fraction(&llrs, FILTER_PRIOR);
        assert!(close(pi, 11.0 / 12.0, 1e-7), "{pi}");
    }

    /// Ten pseudo-calls at the filter's fraction hold a lone call's stratum
    /// near it, and a stratum of a thousand calls near its own fraction.
    #[test]
    fn test_a_stratum_prior_shrinks_toward_its_filter_by_ten_calls() {
        let prior = BetaPrior {
            mean: 0.1,
            strength: STRATUM_PRIOR_STRENGTH,
        };
        assert!(close(learn_artifact_fraction(&[], prior), 0.1, 1e-12));
        let lone = learn_artifact_fraction(&[50.0], prior);
        assert!(close(lone, 2.0 / 11.0, 1e-9), "{lone}");
        let mut llrs = vec![50.0; 500];
        llrs.extend([-50.0; 500]);
        let many = learn_artifact_fraction(&llrs, prior);
        assert!(close(many, 501.0 / 1010.0, 1e-9), "{many}");
    }

    #[test]
    fn test_learned_fraction_of_one_strong_artifact_call() {
        let pi = learn_artifact_fraction(&[50.0], FILTER_PRIOR);
        assert!(close(pi, 2.0 / 3.0, 1e-9), "{pi}");
    }
}
