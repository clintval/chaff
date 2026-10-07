//! Priors on whether a call is an artifact, and the posterior they give.
//!
//! The `chaff` model learns its prior. Each call `i` carries a likelihood ratio
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
//! concave in `pi`, so the fixed point is unique, and `chaff` solves for it
//! directly rather than iterating.
//!
//! The prior is learned twice. Each filter first learns one fraction from all
//! of its calls under the weak [`FILTER_PRIOR`], `Beta(2, 2)`. Each stratum
//! then learns its own under [`STRATUM_PRIOR_STRENGTH`] pseudo-calls at its
//! filter's fraction, so a stratum of one or two calls mostly inherits the
//! filter's fraction rather than moving its prior toward its own calls, and a
//! stratum of hundreds keeps nearly its own.
//!
//! The `chaff` model learns each decay filter's scale `s` with its fraction.
//! The maximum likelihood scale maximizes the filter's marginal likelihood
//! over `(pi, s)`, a profile over `s` with `pi` solved exactly at each, and is
//! then shrunk toward the filter's default the way a stratum's fraction is
//! shrunk toward its filter's: the two are averaged in log space, weighted by
//! the calls' expected artifacts `sum_i r_i` and [`SCALE_PRIOR_STRENGTH`]
//! pseudo-calls at the default,
//!
//! ```text
//! ln s = (sum_i r_i ln s_mle + k ln s_default) / (sum_i r_i + k)
//! ```
//!
//! The weight counts artifacts rather than calls because only an artifact's
//! molecules carry its scale, so a library of a few artifacts keeps nearly the
//! default and one of hundreds keeps nearly its own. The scale is one per
//! filter, shared by its strata, since how far a polymerase copies is a
//! property of the library's enzymes and not of the substitution, and pooling
//! the strata gives the fit the most artifact calls.
//!
//! The `fgbio` model keeps fgbio's prior for parity: a mutation prior of `min((2 * maf)^2,
//! 0.9999)`, where `maf` is the call's alternate molecule fraction, or one over
//! the depth when no alternate molecule is seen. At the low allele fractions of
//! Duplex Sequencing it is near zero, which makes any call whose alternate
//! molecules all sit inside the window an artifact.

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

/// The decay scales, in bases, a scale is learned between.
pub const SCALE_RANGE: (f64, f64) = (1.0, 1000.0);

/// The pseudo-calls at its default that a learned scale is shrunk toward: ten,
/// as many as a stratum's prior holds at its filter's fraction.
pub const SCALE_PRIOR_STRENGTH: f64 = STRATUM_PRIOR_STRENGTH;

/// `ln(exp(a) + exp(b))`, without overflow.
pub(crate) fn ln_add_exp(a: f64, b: f64) -> f64 {
    a.max(b) + (-(a - b).abs()).exp().ln_1p()
}

/// The log marginal likelihood of a set of calls at artifact fraction `pi`,
/// relative to every call being a mutation, plus the log density of `prior`:
/// `sum_i ln(1 - pi + pi exp(l_i)) + k m ln(pi) + k (1 - m) ln(1 - pi)`.
pub fn log_marginal_likelihood(log_likelihood_ratios: &[f64], pi: f64, prior: BetaPrior) -> f64 {
    let (ln_pi, ln_rest) = (pi.ln(), (-pi).ln_1p());
    let calls: f64 = log_likelihood_ratios
        .iter()
        .map(|l| ln_add_exp(ln_rest, ln_pi + l))
        .sum();
    calls + prior.strength * (prior.mean * ln_pi + (1.0 - prior.mean) * ln_rest)
}

/// A set of calls' decay scale: the maximum likelihood scale of
/// [`max_likelihood_scale`] and `default` averaged in log space, weighted by
/// the calls' expected artifacts at that scale and [`SCALE_PRIOR_STRENGTH`]
/// pseudo-calls at `default`. `log_likelihood_ratios` gives the calls' ratios
/// at a scale; without a call to learn from, the scale is `default`.
pub fn learn_scale(
    log_likelihood_ratios: impl Fn(f64) -> Vec<f64>,
    prior: BetaPrior,
    default: f64,
) -> f64 {
    if log_likelihood_ratios(default).is_empty() {
        return default;
    }
    let mle = max_likelihood_scale(&log_likelihood_ratios, prior);
    let llrs = log_likelihood_ratios(mle);
    let odds = logit(learn_artifact_fraction(&llrs, prior));
    let artifacts: f64 = llrs.iter().map(|l| sigmoid(l + odds)).sum();
    let k = SCALE_PRIOR_STRENGTH;
    ((artifacts * mle.ln() + k * default.ln()) / (artifacts + k)).exp()
}

/// The decay scale that maximizes a set of calls' marginal likelihood, with
/// their artifact fraction under `prior` solved exactly at each scale. The
/// profile is searched over a log grid of [`SCALE_RANGE`] and refined by
/// golden-section search around the grid's best.
pub fn max_likelihood_scale(
    log_likelihood_ratios: impl Fn(f64) -> Vec<f64>,
    prior: BetaPrior,
) -> f64 {
    const STEPS: usize = 48;
    let profile = |ln_scale: f64| {
        let llrs = log_likelihood_ratios(ln_scale.exp());
        let pi = learn_artifact_fraction(&llrs, prior);
        log_marginal_likelihood(&llrs, pi, prior)
    };
    let (low, high) = (SCALE_RANGE.0.ln(), SCALE_RANGE.1.ln());
    let grid: Vec<f64> = (0..=STEPS)
        .map(|i| low + (high - low) * i as f64 / STEPS as f64)
        .collect();
    let values: Vec<f64> = grid.iter().map(|&x| profile(x)).collect();
    let best = (0..=STEPS)
        .max_by(|&i, &j| values[i].total_cmp(&values[j]))
        .expect("the grid has points");
    let ratio = (5f64.sqrt() - 1.0) / 2.0;
    let (mut a, mut b) = (grid[best.saturating_sub(1)], grid[(best + 1).min(STEPS)]);
    let (mut c, mut d) = (b - ratio * (b - a), a + ratio * (b - a));
    let (mut fc, mut fd) = (profile(c), profile(d));
    while b - a > 1e-6 {
        if fc >= fd {
            (b, d, fd) = (d, c, fc);
            c = b - ratio * (b - a);
            fc = profile(c);
        } else {
            (a, c, fc) = (c, d, fd);
            d = a + ratio * (b - a);
            fd = profile(d);
        }
    }
    let refined = 0.5 * (a + b);
    if profile(refined) >= values[best] {
        refined.exp()
    } else {
        grid[best].exp()
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
    fn test_log_marginal_likelihood_weighs_each_call_by_the_fraction() {
        let prior = BetaPrior {
            mean: 0.5,
            strength: 0.0,
        };
        let pi: f64 = 0.25;
        let expected = (0.75 + 0.25 * 3f64.exp()).ln() + (0.75 + 0.25 * (-2f64).exp()).ln();
        let got = log_marginal_likelihood(&[3.0, -2.0], pi, prior);
        assert!(close(got, expected, 1e-12), "{got}");
        let weak = log_marginal_likelihood(&[], pi, FILTER_PRIOR);
        assert!(close(weak, pi.ln() + (1.0 - pi).ln(), 1e-12), "{weak}");
    }

    /// Calls whose ratio peaks at a scale of 20 bases have that maximum
    /// likelihood scale.
    #[test]
    fn test_max_likelihood_scale_finds_the_peak_of_the_profile() {
        let llrs = |scale: f64| vec![4.0 - (scale.ln() - 20f64.ln()).powi(2); 3];
        let scale = max_likelihood_scale(llrs, FILTER_PRIOR);
        assert!(close(scale, 20.0, 1e-3), "{scale}");
    }

    /// Five hundred artifact calls peaking at 20 bases keep nearly that scale,
    /// two barely move it from the default of 30, and as many calls that are
    /// mutations, whose ratios also peak at 20, do not move it at all.
    #[test]
    fn test_learn_scale_weighs_the_default_as_ten_artifact_calls() {
        let peaked =
            |n, height: f64| move |scale: f64| vec![height - (scale.ln() - 20f64.ln()).powi(2); n];
        let many = learn_scale(peaked(500, 8.0), FILTER_PRIOR, 30.0);
        let expected = ((500.0 * 20f64.ln() + 10.0 * 30f64.ln()) / 510.0).exp();
        assert!(close(many, expected, 0.01), "{many} vs {expected}");
        let two = learn_scale(peaked(2, 8.0), FILTER_PRIOR, 30.0);
        assert!(two > 27.5 && two < 30.0, "{two}");
        let mutations = learn_scale(peaked(500, -8.0), FILTER_PRIOR, 30.0);
        assert!(close(mutations, 30.0, 0.1), "{mutations}");
        assert_eq!(learn_scale(|_| Vec::new(), FILTER_PRIOR, 30.0), 30.0);
    }

    #[test]
    fn test_learned_fraction_of_one_strong_artifact_call() {
        let pi = learn_artifact_fraction(&[50.0], FILTER_PRIOR);
        assert!(close(pi, 2.0 / 3.0, 1e-9), "{pi}");
    }
}
