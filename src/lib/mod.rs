//! `chaff`: separate somatic variant calls from library-preparation damage
//! artifacts in Duplex Sequencing and other UMI sequencing.
//!
//! # Models
//!
//! Each template covering a call is one molecule (see [`evidence`]), with the
//! base it holds at the site, that base's quality, and the site's 0-based
//! distances from the template's leftmost and rightmost bases, counted in
//! template bases. A filter compares where a call's alternate molecules sit
//! with where its reference molecules sit, and writes the posterior
//! probability that the call is a real mutation. A base's error probability
//! is `e = min(10^(-q / 10), 0.75)` for base quality `q`.
//!
//! ## Likelihood ratios
//!
//! Each call `i` gets a natural-log likelihood ratio of artifact to mutation,
//! `l_i = ln P(molecules | artifact) - ln P(molecules | mutation)`.
//!
//! - **Decays**, for copied damage under either model and end repair fill-in
//!   under the `chaff` model (see [`copied_damage`] and [`read_end`]). A copy
//!   reaches distance `d` from the artifact's end with probability
//!   `w(d) = exp(-d / s)`, so the artifact's alternate molecules follow the
//!   reference molecules' distances tilted by `w`, while a mutation's follow
//!   them untilted. With `W` the mean of `w(d)` over the call's reference
//!   molecules, `l_i = sum ln((1 - e) w(d) / W + e)` over its alternate
//!   molecules, since a base error lands wherever a reference molecule could.
//!   Copied damage measures `d` from the lesion strand's 5' end: the leftmost
//!   base for a lesion on the forward strand (a `C>T` or `G>T`), the rightmost
//!   for one on the reverse strand (a `G>A` or `C>A`). End repair fill-in
//!   measures it from the 3' end of the strand each template was copied from:
//!   read 1's strand, the rightmost base of an F1R2 template and the leftmost
//!   of an F2R1 one, or the nearer end of a duplex consensus whose reads carry
//!   fgbio's `aD` and `bD` depths of both strands.
//! - **Windows**, for A-tailing under either model and end repair fill-in
//!   under the `fgbio` model, are fgbio's. A molecule is congruent when the
//!   site lies within the window of the relevant end, at a 1-based distance of
//!   at most the window: the nearest end for end repair fill-in, and for
//!   A-tailing the leftmost end for a `T` or the rightmost for an `A`. With `f`
//!   the congruent fraction of the call's reference molecules, a congruent
//!   alternate molecule contributes `ln(1 - e) - ln((1 - e) f + e (1 - f))`
//!   and an incongruent one `ln(e) - ln((1 - e)(1 - f) + e f)`.
//!
//! A call with alternate but no reference molecules, or for a decay without a
//! measured alternate molecule, gets no ratio and no posterior.
//!
//! ## Priors and posteriors
//!
//! The posterior probability of a real mutation is `σ(-(l_i + logit π))` for
//! an artifact prior `π`, and a filter applies its FILTER where the posterior
//! is at or below its threshold.
//!
//! Under the `chaff` model the prior is learned per sample (see [`prior`]).
//! Each filter's calls form a two-component mixture whose artifact fraction
//! is the fixed point of expectation-maximization, solved exactly by
//! bisection: `π_f = (sum r_i + 1) / (n + 2)` with
//! `r_i = σ(l_i + logit π_f)`, the maximum a posteriori fraction under
//! `Beta(2, 2)`. Each stratum, the damage class and CpG context for copied
//! damage and the substitution for the others, then learns
//! `π = (sum r_i + k π_f) / (n + k)` from its own calls and `k = 10`
//! pseudo-calls at its filter's fraction, so a stratum of one or two calls
//! mostly inherits `π_f` and one of hundreds keeps nearly its own.
//!
//! On a duplex BAM whose reads carry each strand's single-strand consensus,
//! copied damage under the `chaff` model takes its prior from the library
//! instead (see [`simplex`]): for a call with `k` alternate molecules, the
//! share `E(k) / S(k)` of the library's positions with `k` duplex changes of
//! its stratum that chance explains, as `π = (min(E(k), S(k)) + k' π_s) /
//! (S(k) + k')` with `k' = 10` pseudo-positions at the stratum's learned
//! fraction `π_s`.
//!
//! Under the `fgbio` model the prior is fgbio's per call: an artifact prior of
//! `1 - min((2 m)^2, 0.9999)`, where `m` is the call's alternate molecule
//! fraction, or one over its depth when no alternate molecule is seen.
//!
//! ## Decay scales
//!
//! Each decay filter's scale is one per filter, shared by its strata, and
//! learned unless fixed. The maximum likelihood scale `s_mle` maximizes the
//! filter's marginal likelihood `sum ln(1 - π + π exp(l_i(s)))`, with `π`
//! solved exactly at each `s`, over a log grid of 1 to 1,000 bases refined by
//! golden-section search. It is then shrunk toward the filter's default `s_0`,
//! 30 bases for copied damage and 15 for end repair fill-in, as
//! `ln s = (sum r_i ln s_mle + k ln s_0) / (sum r_i + k)` with `k = 10`: only
//! artifact calls carry a scale, so their expected count weighs the data
//! against the pseudo-calls.
//!
//! ## Metrics
//!
//! The metrics have one row per filter and stratum (see [`metrics`]). A
//! molecule is congruent when the site lies within the filter's distance of
//! the end its artifact favors: a 0-based distance less than a decay's scale,
//! or a 1-based distance of at most a window. The asymmetry test asks whether
//! more alternate molecules are congruent than each call's own reference
//! molecules predict. Call `i` contributes its `a_i` alternate molecules, each
//! congruent under the null with probability
//! `f_i = (ref_congruent_i + 1) / (ref_molecules_i + 2)`, and the p-value is
//! `P(X >= alt_congruent)` for `X` the sum of the calls' `Binomial(a_i, f_i)`,
//! the exact Poisson-binomial upper tail. Because each call is measured
//! against its own reference molecules, skew that a site's alleles share,
//! such as capture or fragment length, stays out of the test.
//!
//! # Differences from fgbio
//!
//! The `chaff` crate ports the filters, likelihoods, and tests of fgbio's
//! `FilterSomaticVcf`, and with `--model fgbio` it writes fgbio 4.1.1's values
//! and FILTERs wherever overlapping mates agree in base and quality, no read
//! has an indel or soft clip between the call and its mate's 5' end, and every
//! base is Q2 or better. It matches fgbio where fgbio's choices are arbitrary:
//! a deletion at the site counts in the depth of its prior, a spanning
//! deletion `*` is no called allele, and an A-tailing site equally far from
//! both template ends is nearer the kept read's own end.
//!
//! It differs on purpose here:
//!
//! - **The `chaff` model, the default.** fgbio's mutation prior is near zero
//!   at duplex allele fractions, so alternate molecules inside the window make
//!   a call an artifact however many reference molecules sit there too; the
//!   `chaff` model learns the prior per library instead (see [`prior`]).
//!   fgbio's end repair window counts either template end and stops at a fixed
//!   distance; the `chaff` model measures from the 3' end of the strand each
//!   template was copied from and decays with distance at a learned scale (see
//!   [`model`]).
//! - **End repair.** A polymerase extends a recessed 3' end across a 5'
//!   overhang; fgbio's docs describe filling in a 3' overhang.
//! - **Template ends.** Distances from both template ends count template
//!   bases: the `chaff` crate walks both reads' CIGARs, the mate's from its
//!   `MC` tag, as
//!   [fgbio #1172](https://github.com/fulcrumgenomics/fgbio/pull/1172) does
//!   for clipping, so an indel counts by its length. Soft clips count and hard
//!   clips do not. fgbio measures the far end by insert size, while the
//!   `chaff` crate never reads `TLEN`.
//! - **Overlapping mates.** They are called into one base: mates that agree
//!   keep the higher quality, and mates that disagree count as neither allele.
//!   fgbio keeps the first read of each name, so values differ where
//!   overlapping mates differ in base or quality.
//! - **Base errors.** A base's error probability is capped at 0.75, a random
//!   base's, so a Q0 or Q1 base cannot zero a likelihood.
//! - **Missing evidence.** A call with alternate but no reference molecules
//!   gets no INFO value; fgbio writes `NaN`.
//! - **Streaming.** The BAM is always streamed, never queried by index.
//! - **Number format.** Values keep htsjdk's rounding but are written in
//!   decimal: `0.00003218` for fgbio's `3.218e-05`.
#![warn(missing_docs)]

pub mod call;
pub mod classes;
pub mod copied_damage;
pub mod evidence;
pub mod filter;
pub mod io;
pub mod metrics;
pub mod model;
pub mod prior;
pub mod read_end;
pub mod reference;
pub mod simplex;
pub mod spectrum;
pub mod testing;
