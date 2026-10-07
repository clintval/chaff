//! What a duplex BAM's single-strand consensus says about a library's damage.
//!
//! A duplex consensus from fgbio's or fgumi's `CallDuplexConsensusReads`
//! carries, per base, the single-strand consensus of each of the molecule's two
//! strands: `ac` and `bc` hold their bases and `ad` and `bd` the raw reads
//! behind them. Once aligned, as by `ZipperBams`, the bases are reverse
//! complemented with the consensus, so they sit in the reference's
//! orientation. A consensus counts only with all four tags, read past any
//! hard clip they still hold. Where mates overlap, only the one that sorts
//! first counts, so each molecule counts once.
//!
//! A library profile reads every consensus once and counts, at each reference base
//! a damage class can change, two kinds of molecule:
//!
//! - **Single-strand changes:** one strand holds the reference base and the
//!   other the class's damaged base, each read by at least
//!   [`MIN_STRAND_READS`] raw reads. These are lesions polymerase never
//!   copied.
//! - **Duplex changes:** the consensus base, at the base quality floor, and
//!   both strands' bases are the damaged base. These are copied lesions and
//!   real mutations alike.
//!
//! Their rates per molecule, and the conversion ratio of the duplex rate to the
//! single-strand rate, describe the library's damage; the rates leave out
//! positions whose changes are at least 3 and 1% of their molecules, as
//! germline or clonal.
//!
//! Duplex changes also give the chance of a call. Copied lesions fall on a
//! position independently, molecule by molecule, so most positions with one
//! hold one, while a real mutation shared by a clone puts several on one
//! position. A lesion's rate varies from position to position, as with
//! methylation at CpG, which makes chance put several on one position more
//! often than a uniform rate would. Single-strand changes measure that
//! variation free of real mutations, which only duplex changes carry: their
//! counts at each position follow a negative binomial, a Poisson whose rate
//! varies as a gamma of shape `a`, fitted so that at their rate per molecule
//! the positions expected with no change match those observed. With `n_j`
//! molecules at position `j` and `S(k)` positions holding `k` duplex
//! changes, the duplex rate `r` solves `sum_j NB(0; n_j r, a) = S(0)`, so
//! that chance explains every position without a change and, but for the
//! few real mutations on one molecule, every position with one, and chance
//! puts `k` changes on
//!
//! ```text
//! E(k) = sum_j NB(k; n_j r, a)
//! ```
//!
//! positions. `E(k) / S(k)` is the share of positions with `k` changes that
//! chance explains: few of a clean library's positions with two or more, and
//! most of a damaged library's. Chance puts `k` changes on a deep position
//! far more often than on a shallow one, so `E(k)` and `S(k)` are counted
//! within depth bins, one per doubling of `n_j`. Chance cannot explain more
//! positions than there are, so when `a` has `E(k)` exceed `S(k)` in some
//! bin, for some `k` of two or more, by more than [`FIT_TOLERANCE`] standard
//! deviations of a Poisson count, the model raises `a` toward a Poisson's
//! until it no longer does, and reports what still exceeds `S(k)` as its
//! misfit. Deep counts are sparse, so from [`MAX_CHANGES`] on a count takes
//! the positions with it or more. Positions whose changes are at least 2 and
//! 20% of their molecules are germline and count toward neither model, nor
//! do positions so shallow that 2 changes would make them germline, where a
//! germline variant shows as one change, as damage would.

use std::collections::{BTreeMap, HashMap, VecDeque};
use std::fs::File;
use std::num::NonZero;
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{bail, Context as _, Result};
use log::{info, warn};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::value::Array;
use noodles::sam::alignment::record::data::field::Value;
use noodles_bgzf as bgzf;

use crate::classes::{complement, Context, DamageClass};
use crate::evidence::PileupOptions;
use crate::reference::Reference;

/// The fewest raw reads behind each strand's base for a molecule to count
/// toward single-strand changes, so a single read's error does not.
pub const MIN_STRAND_READS: i64 = 2;

/// The changes from which the chance model takes positions with a count or
/// more rather than that count alone, as deep counts are sparse: a call
/// with this many or more takes the share of positions with as many or
/// more, and the fit checks positions with this many or more together.
pub const MAX_CHANGES: u32 = 8;

/// The fewest changes at a position the chance model gives a call's prior
/// for: a lone change is as likely a copied lesion as a mutation on one
/// molecule, so its count says nothing, and the rate takes nearly all of
/// them as chance.
pub const MIN_CHANCE_CHANGES: u32 = 2;

/// The records a library profile reads before it decides the BAM has no
/// single-strand consensus.
const TAG_PROBE: u64 = 1_000;

/// The most threads that decompress the BAM a library profile reads, as many
/// as the machine offers up to this.
const MAX_DECOMPRESSION_THREADS: NonZero<usize> = NonZero::new(4).unwrap();

/// The reference bases a library profile holds in memory at once.
const REFERENCE_CHUNK: usize = 1 << 20;

/// The longest reference span of a consensus the profile counts, beyond which it
/// is skipped rather than held.
const MAX_SPAN: usize = 100_000;

/// The standard deviations of a Poisson count by which the positions chance
/// expects with a count of changes may exceed those observed before the
/// chance model takes its gamma shape to vary too much.
pub const FIT_TOLERANCE: f64 = 3.0;

/// The deepest positions the chance model keeps at their own depth.
const EXACT_DEPTH: u32 = 64;

/// The steps per doubling within which deeper positions pool.
const DEPTH_STEPS: f64 = 64.0;

/// The least and greatest finite gamma shapes the models search.
const DISPERSIONS: (f64, f64) = (0.05, 1e4);

/// The fewest changes that, with [`CLONAL_FRACTION`], make a position
/// germline or clonal, left out of the rates.
const CLONAL_CHANGES: u32 = 3;

/// The least share of a position's molecules that, with [`CLONAL_CHANGES`],
/// makes a position germline or clonal.
const CLONAL_FRACTION: f64 = 0.01;

/// The least share of a position's molecules whose changes make it germline,
/// left out of the chance model too.
const GERMLINE_FRACTION: f64 = 0.2;

/// The lowest agreement of the strands' bases with the consensus base, over
/// consensus on the reverse strand, for the tags to be taken as aligned with
/// their consensus.
const MIN_ORIENTATION_AGREEMENT: f64 = 0.9;

/// Positions by their molecules and their changes among them.
pub type Positions = BTreeMap<(u32, u32), u64>;

/// One damage stratum's counts over a library.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct StratumProfile {
    /// Molecules whose consensus base passes the quality floor at a base of
    /// the class.
    pub molecules: u64,
    /// Those whose consensus base is the damaged base.
    pub changes: u64,
    /// Molecules with both strands' bases called, each by at least
    /// [`MIN_STRAND_READS`] raw reads.
    pub strand_molecules: u64,
    /// Those with a single-strand change.
    pub single_strand_changes: u64,
    /// The positions by their molecules and duplex changes.
    pub positions: Positions,
    /// The positions by their molecules with both strands called and their
    /// single-strand changes.
    pub strand_positions: Positions,
}

impl StratumProfile {
    /// Duplex changes per molecule.
    pub fn change_rate(&self) -> Option<f64> {
        ratio(self.changes, self.molecules)
    }

    /// Single-strand changes per molecule with both strands called.
    pub fn single_strand_rate(&self) -> Option<f64> {
        ratio(self.single_strand_changes, self.strand_molecules)
    }

    /// The duplex change rate over the single-strand change rate.
    pub fn conversion_ratio(&self) -> Option<f64> {
        let single = self.single_strand_rate().filter(|&rate| rate > 0.0)?;
        Some(self.change_rate()? / single)
    }

    /// The stratum's chance model: the dispersion fitted to its single-strand
    /// changes, raised toward a Poisson's until chance expects no more
    /// positions with any count of two or more than were observed, the rate
    /// of duplex changes it matches to the positions with none, and the
    /// positions expected and observed with each count.
    pub fn chance(&self) -> Chance {
        let fitted = self.single_strand_dispersion();
        let fit = |dispersion: f64| Chance {
            fitted,
            ..self.chance_at(dispersion)
        };
        let chance = fit(fitted);
        if chance.fits() || fitted.is_infinite() {
            return chance;
        }
        if !fit(DISPERSIONS.1).fits() {
            return fit(f64::INFINITY);
        }
        let (mut low, mut high) = (fitted.ln(), DISPERSIONS.1.ln());
        for _ in 0..20 {
            let mid = (low + high) / 2.0;
            if fit(mid.exp()).fits() {
                high = mid;
            } else {
                low = mid;
            }
        }
        fit(high.exp())
    }

    /// The gamma shape of the single-strand change rate's variation across
    /// positions, infinite when it varies no more than a Poisson allows.
    pub fn single_strand_dispersion(&self) -> f64 {
        fitted_dispersion(&self.strand_positions)
    }

    /// The chance model at a given dispersion.
    pub fn chance_at(&self, dispersion: f64) -> Chance {
        let depths = depths(&self.positions);
        let zeros = observed(&self.positions, 0) as f64;
        let rate = matched_rate(&depths, zeros, dispersion);
        let mut bins: BTreeMap<u32, Tally> = BTreeMap::new();
        for (&(n, k), &count) in &self.positions {
            bins.entry(depth_bin(n)).or_default().observe(k, count);
        }
        for (n, count) in depths {
            let tally = bins.entry(depth_bin(n)).or_default();
            let mean = f64::from(n) * rate;
            for (k, p) in (0..).zip(probabilities(mean, dispersion)) {
                if germline(k, n) || (k >= MAX_CHANGES && f64::from(k) > mean && p < 1e-16) {
                    break;
                }
                tally.expect(k, count * p);
            }
        }
        let mut pooled = Tally::default();
        let (mut beyond, mut seen) = (0.0, 0);
        for tally in bins.values() {
            for (k, &expected) in (0..).zip(&tally.expected) {
                pooled.expect(k, expected);
            }
            for (k, &observed) in (0..).zip(&tally.observed) {
                pooled.observe(k, observed);
            }
            let (excess, observed) = tally.excess();
            beyond += excess;
            seen += observed;
        }
        Chance {
            rate,
            dispersion,
            fitted: dispersion,
            pooled,
            bins,
            excess: beyond / seen.max(1) as f64,
        }
    }
}

/// The depth bin of a position or call with `molecules` molecules, one per
/// doubling, within which the chance model compares positions.
pub fn depth_bin(molecules: u32) -> u32 {
    molecules.max(1).ilog2()
}

/// Whether `changes` of a position's or call's `molecules` make it
/// germline: at least 2 and 20% of them.
pub fn germline(changes: u32, molecules: u32) -> bool {
    changes >= 2 && f64::from(changes) / f64::from(molecules.max(1)) >= GERMLINE_FRACTION
}

/// Whether a position of `molecules` molecules counts toward the chance
/// model: deep enough that [`MIN_CHANCE_CHANGES`] changes would not make it
/// germline, so a call there could take its prior from the model. A
/// shallower position shows a germline variant as one change, or none, as
/// damage would, and off-target consensus leaves many of them.
fn priced(molecules: u32) -> bool {
    !germline(MIN_CHANCE_CHANGES, molecules)
}

/// The probabilities of 0, 1, 2, and more changes at a position whose
/// molecules expect `mean` of them, as [`negative_binomial`] gives them, by
/// the ratio of each to the last.
fn probabilities(mean: f64, dispersion: f64) -> impl Iterator<Item = f64> {
    let first = negative_binomial(0, mean, dispersion);
    std::iter::successors(Some((0.0, first)), move |&(k, p): &(f64, f64)| {
        let ratio = if dispersion.is_infinite() {
            mean / (k + 1.0)
        } else {
            (dispersion + k) / (k + 1.0) * mean / (dispersion + mean)
        };
        Some((k + 1.0, p * ratio))
    })
    .map(|(_, p)| p)
}

/// Each depth's positions, whatever their changes. Depths above
/// [`EXACT_DEPTH`] pool, within steps of 1/[`DEPTH_STEPS`] of a doubling, at
/// their mean, which moves a position's chance of a change by far less than
/// its noise and spares the fits most of their depths.
fn depths(positions: &Positions) -> BTreeMap<u32, f64> {
    let mut steps: BTreeMap<(bool, u32), (f64, f64)> = BTreeMap::new();
    for (&(n, _), &count) in positions {
        let step = if n <= EXACT_DEPTH {
            (false, n)
        } else {
            (true, (DEPTH_STEPS * f64::from(n).log2()) as u32)
        };
        let (molecules, positions) = steps.entry(step).or_default();
        *molecules += f64::from(n) * count as f64;
        *positions += count as f64;
    }
    let mut depths = BTreeMap::new();
    for (molecules, positions) in steps.into_values() {
        *depths
            .entry((molecules / positions).round() as u32)
            .or_default() += positions;
    }
    depths
}

/// The positions chance expects, `E(k)`, and those observed, `S(k)`, with
/// each count of changes `k`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Tally {
    /// `E(k)`, from `k` of 0.
    pub expected: Vec<f64>,
    /// `S(k)`, from `k` of 0.
    pub observed: Vec<u64>,
}

impl Tally {
    fn expect(&mut self, changes: u32, positions: f64) {
        let k = changes as usize;
        if self.expected.len() <= k {
            self.expected.resize(k + 1, 0.0);
        }
        self.expected[k] += positions;
    }

    fn observe(&mut self, changes: u32, positions: u64) {
        let k = changes as usize;
        if self.observed.len() <= k {
            self.observed.resize(k + 1, 0);
        }
        self.observed[k] += positions;
    }

    /// The expected and observed positions with `changes` changes, or with
    /// at least as many from [`MAX_CHANGES`] on.
    pub fn at(&self, changes: u32) -> (f64, u64) {
        let k = changes as usize;
        if changes < MAX_CHANGES {
            let expected = self.expected.get(k).copied().unwrap_or(0.0);
            return (expected, self.observed.get(k).copied().unwrap_or(0));
        }
        let expected = self.expected.iter().skip(k).sum();
        (expected, self.observed.iter().skip(k).sum())
    }

    /// The counts of changes the fit checks: each from 2 below
    /// [`MAX_CHANGES`], and those from it on together.
    fn checked(&self) -> impl Iterator<Item = (f64, u64)> + '_ {
        (MIN_CHANCE_CHANGES..=MAX_CHANGES).map(|k| self.at(k))
    }

    /// The positions with two or more changes that chance expects beyond
    /// those observed and [`FIT_TOLERANCE`] standard deviations of a Poisson
    /// count, and those observed.
    fn excess(&self) -> (f64, u64) {
        self.checked()
            .fold((0.0, 0), |(beyond, seen), (expected, observed)| {
                let allowed = observed as f64 + FIT_TOLERANCE * expected.sqrt();
                (beyond + (expected - allowed).max(0.0), seen + observed)
            })
    }

    /// Whether chance expects no more positions with any count of two or
    /// more than were observed, within [`FIT_TOLERANCE`] standard deviations
    /// of a Poisson count.
    fn fits(&self) -> bool {
        self.checked().all(|(expected, observed)| {
            expected <= observed as f64 + FIT_TOLERANCE * expected.sqrt()
        })
    }
}

/// A stratum's chance model, solved once.
#[derive(Clone, Debug, PartialEq)]
pub struct Chance {
    /// Duplex changes per molecule that fall by chance.
    pub rate: f64,
    /// The gamma shape of the rate's variation across positions, infinite
    /// when it does not vary.
    pub dispersion: f64,
    /// The gamma shape fitted to the single-strand changes, which
    /// `dispersion` raises when it makes chance explain more positions than
    /// there are.
    pub fitted: f64,
    /// The positions over every depth.
    pub pooled: Tally,
    /// The positions of each [`depth_bin`].
    pub bins: BTreeMap<u32, Tally>,
    /// The positions with two or more changes that chance expects beyond
    /// those observed at their depths and their noise, [`FIT_TOLERANCE`]
    /// standard deviations of a Poisson count, as a share of those observed,
    /// or of one when none are: zero when the model fits.
    pub excess: f64,
}

impl Chance {
    /// The expected and observed positions with `changes` changes over every
    /// depth, or with at least as many from [`MAX_CHANGES`] on.
    pub fn at(&self, changes: u32) -> (f64, u64) {
        self.pooled.at(changes)
    }

    /// The expected and observed positions with `changes` changes in the
    /// depth bin of `molecules`, or with at least as many from
    /// [`MAX_CHANGES`] on.
    pub fn at_depth(&self, molecules: u32, changes: u32) -> (f64, u64) {
        self.bins
            .get(&depth_bin(molecules))
            .map_or((0.0, 0), |tally| tally.at(changes))
    }

    /// Whether, at every depth, chance expects no more positions with any
    /// count of two or more than were observed, within [`FIT_TOLERANCE`]
    /// standard deviations of a Poisson count.
    pub fn fits(&self) -> bool {
        self.bins.values().all(Tally::fits)
    }
}

/// `numerator / denominator`, or `None` without a denominator.
fn ratio(numerator: u64, denominator: u64) -> Option<f64> {
    (denominator > 0).then(|| numerator as f64 / denominator as f64)
}

/// The positions with exactly `changes` changes, or with at least
/// [`MAX_CHANGES`] at that.
fn observed(positions: &Positions, changes: u32) -> u64 {
    let k = changes.min(MAX_CHANGES);
    positions
        .iter()
        .filter(|((_, c), _)| (*c).min(MAX_CHANGES) == k)
        .map(|(_, count)| count)
        .sum()
}

/// The probability of `k` changes at a position whose molecules expect
/// `mean` of them, when the rate varies across positions as a gamma of
/// shape `dispersion` and mean one: a negative binomial, or a Poisson when
/// the dispersion is infinite.
pub fn negative_binomial(k: u32, mean: f64, dispersion: f64) -> f64 {
    if mean <= 0.0 {
        return f64::from(u8::from(k == 0));
    }
    let k_f = f64::from(k);
    if dispersion.is_infinite() {
        let ln_factorial: f64 = (1..=k).map(|i| f64::from(i).ln()).sum();
        return (k_f * mean.ln() - mean - ln_factorial).exp();
    }
    let coefficient: f64 = (0..k)
        .map(|i| ((dispersion + f64::from(i)) / f64::from(i + 1)).ln())
        .sum();
    let zero = -dispersion * (mean / dispersion).ln_1p();
    let rest = k_f * (mean / (dispersion + mean)).ln();
    (coefficient + zero + rest).exp()
}

/// The positions of `depths` expected with `k` changes at `rate` per
/// molecule.
fn expected_at(depths: &BTreeMap<u32, f64>, k: u32, rate: f64, dispersion: f64) -> f64 {
    depths
        .iter()
        .map(|(&n, &count)| count * negative_binomial(k, f64::from(n) * rate, dispersion))
        .sum()
}

/// The rate per molecule, up to one, at which the positions of `depths`
/// expected with no change match `zeros`; zero when every position is
/// without one. A position's chance of no change falls as the rate rises,
/// so the rate is unique, where the positions expected with one change, which
/// rise with the rate until each position expects one and fall past that,
/// would meet their count at a low rate and again at a high one.
fn matched_rate(depths: &BTreeMap<u32, f64>, zeros: f64, dispersion: f64) -> f64 {
    let positions: f64 = depths.values().sum();
    let molecules: f64 = depths.iter().map(|(&n, &count)| f64::from(n) * count).sum();
    if positions - zeros <= 0.0 || molecules <= 0.0 {
        return 0.0;
    }
    let none = |rate: f64| expected_at(depths, 0, rate, dispersion);
    if none(1.0) >= zeros {
        return 1.0;
    }
    let (mut low, mut high) = (((positions - zeros) / molecules).ln(), 0.0);
    for _ in 0..60 {
        let mid = (low + high) / 2.0;
        if none(mid.exp()) > zeros {
            low = mid;
        } else {
            high = mid;
        }
    }
    ((low + high) / 2.0).exp()
}

/// The gamma shape of a rate's variation across positions, fitted so the
/// positions expected with no change match those observed at the rate the
/// changes give per molecule: infinite when they vary no more than a Poisson
/// allows. Single-strand changes give it free of real mutations, which only
/// duplex changes carry. At a given mean, a position's chance of no change
/// rises as the shape falls from a Poisson's, so the shape is where the
/// expected zeros reach those observed.
fn fitted_dispersion(positions: &Positions) -> f64 {
    let depths = depths(positions);
    let zeros = observed(positions, 0) as f64;
    let (molecules, changes) =
        positions
            .iter()
            .fold((0.0, 0.0), |(molecules, changes), (&(n, k), &count)| {
                let count = count as f64;
                (
                    molecules + f64::from(n) * count,
                    changes + f64::from(k) * count,
                )
            });
    if changes == 0.0 {
        return f64::INFINITY;
    }
    let rate = changes / molecules;
    let excess = |ln_dispersion: f64| expected_at(&depths, 0, rate, ln_dispersion.exp()) - zeros;
    if expected_at(&depths, 0, rate, f64::INFINITY) >= zeros {
        return f64::INFINITY;
    }
    let (mut low, mut high) = (DISPERSIONS.0.ln(), DISPERSIONS.1.ln());
    if excess(low) < 0.0 {
        return DISPERSIONS.0;
    }
    for _ in 0..40 {
        let mid = (low + high) / 2.0;
        if excess(mid) < 0.0 {
            high = mid;
        } else {
            low = mid;
        }
    }
    ((low + high) / 2.0).exp()
}

/// The single-strand and duplex changes of a library, per damage stratum,
/// keyed by stratum label such as `C>T:CpG`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LibraryProfile {
    /// Each stratum's counts.
    pub strata: BTreeMap<String, StratumProfile>,
}

impl LibraryProfile {
    /// The counts of a stratum.
    pub fn stratum(&self, stratum: &str) -> Option<&StratumProfile> {
        self.strata.get(stratum)
    }
}

/// One reference position's molecules while consensus still covers it.
#[derive(Clone, Copy, Debug, Default)]
struct Site {
    reference: u8,
    context: Option<Context>,
    molecules: u32,
    changes: [u32; 4],
    strand_molecules: u32,
    single_strand: [u32; 4],
}

/// The index of an upper-case base in `ACGT`.
fn index(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// One consensus's single-strand bases and raw-read depths, from its first
/// base on.
struct Strands<'a> {
    a_bases: &'a [u8],
    b_bases: &'a [u8],
    a_reads: &'a [i64],
    b_reads: &'a [i64],
}

/// Where a consensus's first base sits in its per-base tags of `length`: at
/// 0 when they match its bases, past its leading hard clip when they still
/// hold the clipped bases, or `None` when they fit neither.
fn tag_offset(record: &bam::Record, length: usize) -> Result<Option<usize>> {
    let bases = record.sequence().len();
    if length == bases {
        return Ok(Some(0));
    }
    let kinds = record
        .cigar()
        .iter()
        .map(|op| op.map(|op| (op.kind(), op.len())))
        .collect::<std::io::Result<Vec<_>>>()?;
    let hard = |op: Option<&(Kind, usize)>| match op {
        Some(&(Kind::HardClip, len)) => len,
        _ => 0,
    };
    let leading = hard(kinds.first());
    let trailing = if kinds.len() > 1 {
        hard(kinds.last())
    } else {
        0
    };
    Ok((length == bases + leading + trailing).then_some(leading))
}

/// A contig's reference bases, read a chunk at a time as the consensus advances.
struct Window {
    contig: String,
    start: usize,
    bases: Vec<u8>,
}

impl Window {
    fn base(&self, pos: usize) -> Option<u8> {
        pos.checked_sub(self.start)
            .and_then(|i| self.bases.get(i))
            .copied()
    }
}

/// Counts the molecules of a coordinate-sorted BAM into a [`LibraryProfile`].
struct Scanner<'a> {
    classes: &'a [DamageClass],
    options: PileupOptions,
    reference: Reference,
    window: Option<Window>,
    contig: Option<(usize, String)>,
    head: usize,
    sites: VecDeque<Site>,
    profile: LibraryProfile,
    tagged: u64,
    agreement: [(u64, u64); 2],
    mates: Mates,
}

/// The spans counted by consensus whose mates sort after them, by name, so
/// each mate leaves its mate's span to it.
#[derive(Default)]
struct Mates {
    spans: HashMap<Vec<u8>, Span>,
    prune_at: usize,
}

/// A counted consensus's 0-based, half-open reference span and where its
/// mate starts.
#[derive(Clone, Copy)]
struct Span {
    start: usize,
    end: usize,
    mate_start: usize,
}

impl Mates {
    /// The span a consensus leaves to its mate, after recording its own for
    /// a mate yet to come: its mate's span when its mate was counted before
    /// it, or `None`.
    fn overlap(
        &mut self,
        record: &bam::Record,
        contig_id: usize,
        start: usize,
        end: usize,
    ) -> Option<(usize, usize)> {
        let flags = record.flags();
        if !flags.is_segmented() || flags.is_mate_unmapped() {
            return None;
        }
        let (Some(Ok(mate_contig)), Some(Ok(mate_start))) = (
            record.mate_reference_sequence_id(),
            record.mate_alignment_start(),
        ) else {
            return None;
        };
        let (name, mate_start) = (record.name()?, usize::from(mate_start) - 1);
        if mate_contig != contig_id || mate_start >= end {
            return None;
        }
        if let Some(span) = self.spans.remove(name.as_ref() as &[u8]) {
            return Some((span.start, span.end));
        }
        if mate_start >= start {
            if self.spans.len() >= self.prune_at {
                self.spans.retain(|_, span| span.mate_start >= start);
                self.prune_at = (2 * self.spans.len()).max(1 << 10);
            }
            let span = Span {
                start,
                end,
                mate_start,
            };
            self.spans.insert(name.to_vec(), span);
        }
        None
    }
}

impl<'a> Scanner<'a> {
    fn new(classes: &'a [DamageClass], options: PileupOptions, reference: Reference) -> Self {
        Self {
            classes,
            options,
            reference,
            window: None,
            contig: None,
            head: 0,
            sites: VecDeque::new(),
            profile: LibraryProfile::default(),
            tagged: 0,
            agreement: [(0, 0); 2],
            mates: Mates::default(),
        }
    }

    /// Whether a damage class changes the reference base `base`.
    fn tracked(&self, base: u8) -> bool {
        self.classes
            .iter()
            .any(|c| c.lesion == base || complement(c.lesion) == base)
    }

    /// The reference base at a 0-based position of a contig, reading the
    /// chunk holding it when needed.
    fn reference_base(&mut self, contig: &str, pos: usize) -> Result<Option<u8>> {
        if self.reference.length(contig).is_none() {
            return Ok(None);
        }
        let current = self
            .window
            .as_ref()
            .is_some_and(|w| w.contig == contig && w.base(pos).is_some());
        if !current {
            let start = pos.saturating_sub(1);
            let bases = self
                .reference
                .bases(contig, start, start + REFERENCE_CHUNK)?;
            self.window = Some(Window {
                contig: contig.to_string(),
                start,
                bases,
            });
        }
        Ok(self.window.as_ref().and_then(|w| w.base(pos)))
    }

    /// Fold every position before `pos` into the profile.
    fn flush_before(&mut self, pos: usize) {
        while self.head < pos {
            let Some(site) = self.sites.pop_front() else {
                self.head = pos;
                break;
            };
            self.head += 1;
            self.fold(site);
        }
    }

    fn fold(&mut self, site: Site) {
        let Some(context) = site.context else {
            return;
        };
        for class in self.classes {
            let alt = if site.reference == class.lesion {
                class.reads_as
            } else if site.reference == complement(class.lesion) {
                complement(class.reads_as)
            } else {
                continue;
            };
            let Some(a) = index(alt) else {
                continue;
            };
            let stratum = self
                .profile
                .strata
                .entry(format!("{class}:{context}"))
                .or_default();
            let changes = site.changes[a];
            let share = f64::from(changes) / f64::from(site.molecules.max(1));
            if changes < CLONAL_CHANGES || share < CLONAL_FRACTION {
                stratum.molecules += u64::from(site.molecules);
                stratum.changes += u64::from(changes);
                stratum.strand_molecules += u64::from(site.strand_molecules);
                stratum.single_strand_changes += u64::from(site.single_strand[a]);
            }
            if priced(site.molecules) && !germline(changes, site.molecules) {
                *stratum
                    .positions
                    .entry((site.molecules, changes))
                    .or_default() += 1;
            }
            let single = site.single_strand[a];
            if priced(site.strand_molecules) && !germline(single, site.strand_molecules) {
                *stratum
                    .strand_positions
                    .entry((site.strand_molecules, single))
                    .or_default() += 1;
            }
        }
    }

    fn site(&mut self, pos: usize) -> &mut Site {
        let offset = pos - self.head;
        if offset >= self.sites.len() {
            self.sites.resize(offset + 1, Site::default());
        }
        &mut self.sites[offset]
    }

    /// Count one consensus.
    fn add(&mut self, record: &bam::Record, header: &noodles::sam::Header) -> Result<()> {
        let data = record.data();
        let strand_bases = |tag: &[u8; 2]| match data.get(tag) {
            Some(Ok(Value::String(s))) => Some(s.to_vec()),
            _ => None,
        };
        let strand_reads = |tag: &[u8; 2]| data.get(tag).and_then(Result::ok).and_then(integers);
        let (Some(a_bases), Some(b_bases), Some(a_reads), Some(b_reads)) = (
            strand_bases(b"ac"),
            strand_bases(b"bc"),
            strand_reads(b"ad"),
            strand_reads(b"bd"),
        ) else {
            return Ok(());
        };
        self.tagged += 1;

        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || (self.options.paired_reads_only
                && (!flags.is_segmented() || flags.is_mate_unmapped()))
        {
            return Ok(());
        }
        if record
            .mapping_quality()
            .is_some_and(|q| q.get() < self.options.min_mapping_quality)
        {
            return Ok(());
        }
        let length = a_bases.len();
        if [b_bases.len(), a_reads.len(), b_reads.len()] != [length; 3] {
            return Ok(());
        }
        let Some(offset) = tag_offset(record, length)? else {
            return Ok(());
        };

        let (Some(Ok(contig_id)), Some(Ok(start))) =
            (record.reference_sequence_id(), record.alignment_start())
        else {
            return Ok(());
        };
        let start = usize::from(start) - 1;
        if self.contig.as_ref().map(|(id, _)| *id) != Some(contig_id) {
            self.flush_before(usize::MAX);
            self.sites.clear();
            self.mates.spans.clear();
            let name = header
                .reference_sequences()
                .get_index(contig_id)
                .map(|(name, _)| name.to_string())
                .context("a record's contig is not in the BAM header")?;
            self.contig = Some((contig_id, name));
            self.head = start;
        }
        if start < self.head {
            bail!("the BAM must be coordinate sorted (@HD SO:coordinate) for its single-strand consensus to be profiled");
        }
        self.flush_before(start);
        let span: usize = record
            .cigar()
            .iter()
            .filter_map(Result::ok)
            .filter(|op| op.kind().consumes_reference())
            .map(|op| op.len())
            .sum();
        if span > MAX_SPAN {
            return Ok(());
        }
        let contig = self
            .contig
            .as_ref()
            .map(|(_, name)| name.clone())
            .unwrap_or_default();
        let mate_span = self.mates.overlap(record, contig_id, start, start + span);
        let sequence = record.sequence();
        let qualities = record.quality_scores();
        let qualities = qualities.as_ref();
        let reverse = flags.is_reverse_complemented();
        let strands = Strands {
            a_bases: &a_bases[offset..],
            b_bases: &b_bases[offset..],
            a_reads: &a_reads[offset..],
            b_reads: &b_reads[offset..],
        };

        let (mut query, mut pos) = (0usize, start);
        for op in record.cigar().iter() {
            let op = op?;
            let len = op.len();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for i in 0..len {
                        let (q, p) = (query + i, pos + i);
                        if mate_span.is_some_and(|(s, e)| s <= p && p < e) {
                            continue;
                        }
                        self.add_base(&contig, p, q, &sequence, qualities, &strands, reverse)?;
                    }
                    query += len;
                    pos += len;
                }
                Kind::Insertion | Kind::SoftClip => query += len,
                Kind::Deletion | Kind::Skip => pos += len,
                Kind::HardClip | Kind::Pad => {}
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn add_base(
        &mut self,
        contig: &str,
        pos: usize,
        query: usize,
        sequence: &bam::record::Sequence<'_>,
        qualities: &[u8],
        strands: &Strands<'_>,
        reverse: bool,
    ) -> Result<()> {
        let Some(reference) = self.reference_base(contig, pos)? else {
            return Ok(());
        };
        let seq = sequence.get(query).map(|b| b.to_ascii_uppercase());
        let a = strands.a_bases.get(query).map(u8::to_ascii_uppercase);
        let b = strands.b_bases.get(query).map(u8::to_ascii_uppercase);
        if let (Some(seq), Some(a)) = (seq, a) {
            if index(seq).is_some() && index(a).is_some() {
                let tally = &mut self.agreement[usize::from(reverse)];
                tally.0 += u64::from(seq == a);
                tally.1 += 1;
            }
        }
        if index(reference).is_none() || !self.tracked(reference) {
            return Ok(());
        }
        let prev = if pos > 0 {
            self.reference_base(contig, pos - 1)?
        } else {
            None
        };
        let next = self.reference_base(contig, pos + 1)?;
        let min_quality = self.options.min_base_quality;
        let site = self.site(pos);
        site.reference = reference;
        site.context = Some(Context::of(prev, reference, next));
        if let (Some(seq), Some(&quality)) = (seq, qualities.get(query)) {
            if let (Some(i), true) = (index(seq), quality >= min_quality) {
                site.molecules += 1;
                if seq != reference && a == Some(seq) && b == Some(seq) {
                    site.changes[i] += 1;
                }
            }
        }
        let (Some(a), Some(b)) = (a, b) else {
            return Ok(());
        };
        let (Some(ia), Some(ib)) = (index(a), index(b)) else {
            return Ok(());
        };
        let reads = |r: &[i64]| r.get(query).is_some_and(|&n| n >= MIN_STRAND_READS);
        if !reads(strands.a_reads) || !reads(strands.b_reads) {
            return Ok(());
        }
        site.strand_molecules += 1;
        if b == reference && a != reference {
            site.single_strand[ia] += 1;
        } else if a == reference && b != reference {
            site.single_strand[ib] += 1;
        }
        Ok(())
    }

    fn finish(mut self) -> Option<LibraryProfile> {
        self.flush_before(usize::MAX);
        if self.tagged == 0 {
            return None;
        }
        let (agree, total) = self.agreement[1];
        if total > 0 && (agree as f64) < MIN_ORIENTATION_AGREEMENT * total as f64 {
            warn!(
                "the single-strand consensus bases of reverse-strand consensus agree with its base at only {agree} of {total} bases, \
                 so they were not reverse complemented when the consensus was aligned and chaff leaves them out"
            );
            return None;
        }
        Some(self.profile)
    }
}

/// The values of an integer array field.
fn integers(value: Value<'_>) -> Option<Vec<i64>> {
    let Value::Array(array) = value else {
        return None;
    };
    let values: std::io::Result<Vec<i64>> = match array {
        Array::Int8(v) => v.iter().map(|x| x.map(i64::from)).collect(),
        Array::UInt8(v) => v.iter().map(|x| x.map(i64::from)).collect(),
        Array::Int16(v) => v.iter().map(|x| x.map(i64::from)).collect(),
        Array::UInt16(v) => v.iter().map(|x| x.map(i64::from)).collect(),
        Array::Int32(v) => v.iter().map(|x| x.map(i64::from)).collect(),
        Array::UInt32(v) => v.iter().map(|x| x.map(i64::from)).collect(),
        Array::Float(_) => return None,
    };
    values.ok()
}

/// Profile the single-strand and duplex changes of a coordinate-sorted BAM,
/// decompressed on up to four threads, or `None`
/// when its records carry no single-strand consensus, carry it unaligned, or
/// `stop` is set before the last is read.
pub fn profile_library(
    bam: &Path,
    reference: &Path,
    classes: &[DamageClass],
    options: &PileupOptions,
    stop: &AtomicBool,
) -> Result<Option<LibraryProfile>> {
    let file = File::open(bam).with_context(|| format!("failed to open BAM: {bam:?}"))?;
    let workers = std::thread::available_parallelism()
        .map_or(NonZero::<usize>::MIN, |n| n.min(MAX_DECOMPRESSION_THREADS));
    let decoder = bgzf::io::MultithreadedReader::with_worker_count(workers, file);
    let mut reader = bam::io::Reader::from(decoder);
    let header = reader
        .read_header()
        .context("failed to read the BAM header")?;
    let mut scanner = Scanner::new(classes, *options, Reference::open(reference)?);
    let mut record = bam::Record::default();
    let mut records = 0u64;
    while reader.read_record(&mut record)? != 0 {
        if stop.load(Ordering::Relaxed) {
            return Ok(None);
        }
        records += 1;
        scanner.add(&record, &header)?;
        if records == TAG_PROBE && scanner.tagged == 0 {
            break;
        }
    }
    let profile = scanner.finish();
    match &profile {
        Some(profile) => info!(
            "profiled the single-strand consensus of {records} records in {} strata",
            profile.strata.len()
        ),
        None => info!("the BAM's records carry no single-strand consensus (ac, bc, ad and bd), so chaff learns its priors from the calls alone"),
    }
    Ok(profile)
}

#[cfg(test)]
mod tests {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::value::Array as BufArray;
    use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
    use streampile::testing::{Frag, Pair, SamBuilder, Strand};

    use super::*;
    use crate::testing::{poisson_stratum, write_fasta};

    const UNIT: &str = "AACGTTCA";

    fn reference() -> String {
        UNIT.repeat(50)
    }

    /// A consensus of the reference's first 40 bases with both strands' bases,
    /// and any changes `(0-based offset, consensus base, quality, a, b, a raw
    /// reads)`.
    fn frag(changes: &[(usize, u8, u8, u8, u8, i16)], strand: Strand) -> Frag {
        let mut seq = reference().as_bytes()[..40].to_vec();
        let mut quals = vec![30u8; 40];
        let (mut a, mut b) = (seq.clone(), seq.clone());
        let mut a_reads = vec![3i16; 40];
        for &(i, base, quality, sa, sb, reads) in changes {
            seq[i] = base;
            quals[i] = quality;
            a[i] = sa;
            b[i] = sb;
            a_reads[i] = reads;
        }
        let text = |v: &[u8]| BufValue::from(String::from_utf8(v.to_vec()).unwrap());
        Frag::at(1)
            .bases(String::from_utf8(seq).unwrap())
            .quals(quals)
            .strand(strand)
            .attr(Tag::new(b'a', b'c'), text(&a))
            .attr(Tag::new(b'b', b'c'), text(&b))
            .attr(
                Tag::new(b'a', b'd'),
                BufValue::Array(BufArray::Int16(a_reads)),
            )
            .attr(Tag::new(b'b', b'd'), depths(40))
    }

    /// Raw-read depths of 3 at each of `length` bases.
    fn depths(length: usize) -> BufValue {
        BufValue::Array(BufArray::Int16(vec![3; length]))
    }

    /// The four tags of a consensus whose strands both read `bases`, but for
    /// strand A's `(offset, base)` change.
    fn tags(bases: &str, change: Option<(usize, u8)>) -> [(Tag, BufValue); 4] {
        let mut a = bases.as_bytes().to_vec();
        if let Some((i, base)) = change {
            a[i] = base;
        }
        [
            (
                Tag::new(b'a', b'c'),
                BufValue::from(String::from_utf8(a).unwrap()),
            ),
            (Tag::new(b'b', b'c'), BufValue::from(bases.to_string())),
            (Tag::new(b'a', b'd'), depths(bases.len())),
            (Tag::new(b'b', b'd'), depths(bases.len())),
        ]
    }

    /// A record with the four tags of strands that both read `bases`.
    fn tagged_record(
        mut record: noodles::sam::alignment::RecordBuf,
        bases: &str,
    ) -> noodles::sam::alignment::RecordBuf {
        for (tag, value) in tags(bases, None) {
            record.data_mut().insert(tag, value);
        }
        record
    }

    fn tagged(mut frag: Frag, tags: [(Tag, BufValue); 4]) -> Frag {
        for (tag, value) in tags {
            frag = frag.attr(tag, value);
        }
        frag
    }

    fn profile(reads: &SamBuilder) -> Option<LibraryProfile> {
        let dir = tempfile::tempdir().unwrap();
        let fasta = write_fasta(dir.path(), "chr1", &reference());
        let bam = dir.path().join("reads.bam");
        reads.write_bam(&bam).unwrap();
        let classes = [DamageClass::DEAMINATION, DamageClass::OXIDATION];
        let stop = AtomicBool::new(false);
        profile_library(&bam, &fasta, &classes, &PileupOptions::default(), &stop).unwrap()
    }

    /// Of eleven consensus over positions 1 to 40, one holds a duplex C>T at
    /// the C of a CpG, one a C>T on one strand at a C outside CpG, and one
    /// the same change read by a single raw read, which does not count. The
    /// chance model holds the positions, deep enough for two changes not to
    /// be germline, which ten consensus would not be.
    #[test]
    fn test_a_profile_counts_duplex_and_single_strand_changes_per_stratum() {
        let mut reads = SamBuilder::new().read_length(40);
        reads.add_frag(frag(&[(2, b'T', 30, b'T', b'T', 3)], Strand::Plus));
        reads.add_frag(frag(&[(6, b'N', 2, b'T', b'C', 3)], Strand::Plus));
        reads.add_frag(frag(&[(6, b'N', 2, b'T', b'C', 1)], Strand::Plus));
        for _ in 0..8 {
            reads.add_frag(frag(&[], Strand::Plus));
        }
        let profile = profile(&reads).unwrap();
        let cpg = profile.stratum("C>T:CpG").unwrap();
        let other = profile.stratum("C>T:non-CpG").unwrap();
        // Positions 1 to 40 hold five units: a CpG C and G and one other C each.
        assert_eq!(cpg.molecules, 110);
        assert_eq!(cpg.changes, 1);
        assert_eq!(cpg.strand_molecules, 110);
        assert_eq!(cpg.single_strand_changes, 0);
        assert_eq!(cpg.positions.get(&(11, 1)), Some(&1));
        assert_eq!(cpg.positions.get(&(11, 0)), Some(&9));
        assert_eq!(other.molecules, 53);
        assert_eq!(other.changes, 0);
        assert_eq!(other.strand_molecules, 54);
        assert_eq!(other.single_strand_changes, 1);
        assert_eq!(other.single_strand_rate(), Some(1.0 / 54.0));
        assert_eq!(cpg.conversion_ratio(), None);
        let oxidation = profile.stratum("G>T:CpG").unwrap();
        assert_eq!((oxidation.molecules, oxidation.changes), (110, 0));

        let mut reads = SamBuilder::new().read_length(40);
        for _ in 0..10 {
            reads.add_frag(frag(&[], Strand::Plus));
        }
        let shallow = super::tests::profile(&reads).unwrap();
        let cpg = shallow.stratum("C>T:CpG").unwrap();
        assert_eq!((cpg.molecules, cpg.positions.len()), (100, 0));
        assert!(cpg.strand_positions.is_empty());
    }

    /// A consensus base that only one strand carries is no duplex change,
    /// even at the quality floor, but a single-strand one.
    #[test]
    fn test_a_duplex_change_needs_both_strands() {
        let mut reads = SamBuilder::new().read_length(40);
        reads.add_frag(frag(&[(2, b'T', 30, b'T', b'C', 3)], Strand::Plus));
        reads.add_frag(frag(&[(10, b'T', 30, b'T', b'T', 3)], Strand::Plus));
        let profile = profile(&reads).unwrap();
        let cpg = profile.stratum("C>T:CpG").unwrap();
        assert_eq!(
            (cpg.molecules, cpg.changes, cpg.single_strand_changes),
            (20, 1, 1)
        );
    }

    #[test]
    fn test_a_bam_without_single_strand_consensus_has_no_profile() {
        let mut reads = SamBuilder::new().read_length(40);
        for _ in 0..3 {
            reads.add_frag(Frag::at(1).bases(&reference()[..40]));
        }
        assert_eq!(profile(&reads), None);
    }

    /// A profile told to stop reads no further and gives none.
    #[test]
    fn test_a_stopped_profile_has_none() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = write_fasta(dir.path(), "chr1", &reference());
        let bam = dir.path().join("reads.bam");
        let mut reads = SamBuilder::new().read_length(40);
        reads.add_frag(frag(&[], Strand::Plus));
        reads.write_bam(&bam).unwrap();
        let classes = [DamageClass::DEAMINATION];
        let options = PileupOptions::default();
        let stop = AtomicBool::new(true);
        let profile = profile_library(&bam, &fasta, &classes, &options, &stop).unwrap();
        assert_eq!(profile, None);
    }

    /// Strand bases left in the sequencing orientation disagree with a reverse
    /// strand consensus's base, so the profile leaves them out.
    #[test]
    fn test_unaligned_strand_bases_of_reverse_consensus_have_no_profile() {
        let mut reads = SamBuilder::new().read_length(40);
        reads.add_frag(frag(&[], Strand::Minus));
        assert!(profile(&reads).is_some());
        let mut reads = SamBuilder::new().read_length(40);
        let complemented: String = reference()[..40]
            .bytes()
            .map(|b| complement(b) as char)
            .collect();
        reads.add_frag(tagged(
            Frag::at(1).bases(&reference()[..40]).strand(Strand::Minus),
            tags(&complemented, None),
        ));
        assert_eq!(profile(&reads), None);
    }

    /// Mates over the same positions are one molecule there: the mate that
    /// sorts second leaves its mate's span to it, with or without the mates'
    /// CIGARs, and counts alone where its mate is left out. Of positions 1 to
    /// 50, mates over 1 to 40 and 11 to 50 cover all 12 CpG positions.
    #[test]
    fn test_overlapping_mates_count_once() {
        let molecules = |mate_cigars: bool, first_mapq: u8| {
            let (mut scratch, mut reads) = (
                SamBuilder::new().read_length(40),
                SamBuilder::new().read_length(40),
            );
            let pair = Pair::at(1, 11)
                .bases1(&reference()[..40])
                .bases2(&reference()[10..50])
                .mapq1(first_mapq);
            for record in scratch.add_pair(pair) {
                let bases = String::from_utf8(record.sequence().as_ref().to_vec()).unwrap();
                let mut record = tagged_record(record, &bases);
                if !mate_cigars {
                    record = SamBuilder::without_mate_cigar(record);
                }
                reads.extend([record]);
            }
            let profile = profile(&reads).unwrap();
            profile.stratum("C>T:CpG").unwrap().molecules
        };
        assert_eq!(molecules(true, 60), 12);
        assert_eq!(molecules(false, 60), 12);
        assert_eq!(molecules(true, 5), 10);
    }

    /// A strand change on a reverse consensus, its tags reverse complemented
    /// with it, counts at the reference base it sits on: a G>A at the G of a
    /// CpG is a C>T there.
    #[test]
    fn test_a_reverse_consensus_counts_its_strand_changes_where_they_align() {
        let mut reads = SamBuilder::new().read_length(40);
        reads.add_frag(frag(&[(3, b'N', 2, b'A', b'G', 3)], Strand::Minus));
        let profile = profile(&reads).unwrap();
        assert_eq!(profile.stratum("C>T:CpG").unwrap().single_strand_changes, 1);
        assert_eq!(profile.stratum("G>T:CpG").unwrap().single_strand_changes, 0);
    }

    /// A strand change after an insertion and a deletion counts at its own
    /// reference base, the C at position 31, outside CpG.
    #[test]
    fn test_strand_changes_follow_insertions_and_deletions() {
        let reference = reference();
        let bases = format!(
            "{}GG{}{}",
            &reference[..10],
            &reference[10..20],
            &reference[23..41]
        );
        let mut reads = SamBuilder::new();
        reads.add_frag(tagged(
            Frag::at(1).bases(&bases).cigar("10M2I10M3D18M"),
            tags(&bases, Some((29, b'T'))),
        ));
        let profile = profile(&reads).unwrap();
        let other = profile.stratum("C>T:non-CpG").unwrap();
        assert_eq!(other.single_strand_changes, 1);
        assert_eq!(profile.stratum("C>T:CpG").unwrap().single_strand_changes, 0);
    }

    /// A hard-clipped consensus whose tags still hold the clipped bases reads
    /// them past the clip, as does one whose tags were clipped with it, while
    /// tags that fit neither leave the consensus out, and a BAM whose tags
    /// lack the raw-read depths has no profile.
    #[test]
    fn test_hard_clipped_consensus_reads_its_tags_past_the_clip() {
        let reference = reference();
        let full = &reference[..40];
        let clipped = |tags: [(Tag, BufValue); 4]| {
            let mut reads = SamBuilder::new();
            reads.add_frag(tagged(
                Frag::at(6).bases(&full[5..38]).cigar("5H33M2H"),
                tags,
            ));
            profile(&reads).unwrap()
        };
        for tags in [
            tags(full, Some((6, b'T'))),
            tags(&full[5..38], Some((1, b'T'))),
        ] {
            let other = clipped(tags);
            let other = other.stratum("C>T:non-CpG").unwrap();
            assert_eq!(
                (other.single_strand_changes, other.strand_molecules),
                (1, 4)
            );
        }
        let misfit = clipped(tags(&full[..39], Some((6, b'T'))));
        assert!(misfit.strata.is_empty(), "{misfit:?}");

        let mut reads = SamBuilder::new().read_length(40);
        let [ac, bc, _, _] = tags(full, None);
        reads.add_frag(Frag::at(1).bases(full).attr(ac.0, ac.1).attr(bc.0, bc.1));
        assert_eq!(profile(&reads), None);
    }

    /// Secondary, supplementary, duplicate, and poorly mapped consensus are
    /// left out, and a BAM that opens with more poorly mapped consensus than
    /// the probe reads is still profiled.
    #[test]
    fn test_only_primary_well_mapped_consensus_counts() {
        let full = &reference()[..40];
        let mut reads = SamBuilder::new().read_length(40);
        for _ in 0..TAG_PROBE {
            reads.add_frag(tagged(Frag::at(1).bases(full).mapq(5), tags(full, None)));
        }
        let mut scratch = SamBuilder::new().read_length(40);
        for bits in [0x100, 0x400, 0x800, 0] {
            let record = scratch
                .add_frag(tagged(Frag::at(2).bases(full), tags(full, None)))
                .remove(0);
            reads.extend([SamBuilder::with_flags(record, bits)]);
        }
        let profile = profile(&reads).unwrap();
        assert_eq!(profile.stratum("C>T:CpG").unwrap().molecules, 10);
    }

    /// Single-strand changes that vary no more than a Poisson leave the
    /// chance model a Poisson, whose rate reproduces the positions with no
    /// change.
    #[test]
    fn test_the_chance_model_matches_the_positions_with_no_change() {
        let mut stratum = StratumProfile::default();
        for (k, count) in [(0, 900), (1, 90), (2, 8), (9, 2)] {
            stratum.positions.insert((100, k), count);
        }
        stratum.strand_positions.insert((100, 0), 990);
        stratum.strand_positions.insert((100, 1), 10);
        assert!(stratum.single_strand_dispersion().is_infinite());
        assert!(stratum.chance().dispersion.is_infinite());
        let chance = stratum.chance_at(f64::INFINITY);
        let mean = 100.0 * chance.rate;
        assert!((1000.0 * (-mean).exp() - 900.0).abs() < 1e-6, "{chance:?}");
        assert!((chance.at(0).0 - 900.0).abs() < 1e-6);
        assert!((chance.at(1).0 - 1000.0 * mean * (-mean).exp()).abs() < 1e-6);
        assert!((chance.at(2).0 - 1000.0 * mean * mean / 2.0 * (-mean).exp()).abs() < 1e-6);
        assert_eq!(chance.at(2).1, 8);
        assert_eq!(chance.at(MAX_CHANGES).1, 2);
        assert_eq!(chance.at(9).1, 2);
        assert_eq!(chance.at(10).1, 0);
        assert!(chance.at(MAX_CHANGES).0 < 1e-6);
        assert_eq!(StratumProfile::default().chance().rate, 0.0);
    }

    /// Single-strand changes drawn from a gamma-varying rate give back its
    /// shape, and duplex changes sharing it fall together more often than a
    /// Poisson allows, as the chance model takes them to.
    #[test]
    fn test_the_dispersion_of_single_strand_changes_widens_the_chance_model() {
        let mut stratum = StratumProfile::default();
        for u in 0..12u32 {
            let count = (100_000.0 * negative_binomial(u, 0.3, 1.5)).round() as u64;
            stratum.strand_positions.insert((1000, u), count);
        }
        for k in 0..=MAX_CHANGES {
            let count = (100_000.0 * negative_binomial(k, 0.3, 1.5)).round() as u64;
            stratum
                .positions
                .insert((1000, k), count + if k == 2 { 500 } else { 0 });
        }
        let dispersion = stratum.single_strand_dispersion();
        assert!((dispersion - 1.5).abs() < 0.05, "{dispersion}");
        let chance = stratum.chance_at(dispersion);
        let poisson = stratum.chance_at(f64::INFINITY);
        assert!(
            chance.at(2).0 > 1.3 * poisson.at(2).0,
            "{chance:?} {poisson:?}"
        );
        assert!((chance.at(0).0 - chance.at(0).1 as f64).abs() < 1e-3);
        assert!(chance.fits(), "{chance:?}");
        assert_eq!(stratum.chance(), chance);
    }

    /// Poisson damage at 1e-4 per molecule over 100,000 positions of 1,000
    /// molecules, with 200 real mutations of 2 molecules each, leaves chance
    /// explaining 72% of the positions with 2 changes: 69% are its, and the
    /// real positions, without a change of chance's, draw its rate up by 2%.
    #[test]
    fn test_the_chance_model_leaves_real_mutations_to_the_positions_beyond_chance() {
        let stratum = poisson_stratum(&[1000], 100_000.0, 1e-4, 200);
        let chance = stratum.chance();
        assert!((chance.rate / 1e-4 - 1.02).abs() < 0.005, "{chance:?}");
        let (expected, observed) = chance.at(2);
        let share = expected / observed as f64;
        assert!((share - 0.722).abs() < 0.01, "{share} {chance:?}");
        assert!(chance.fits() && chance.excess < 0.01, "{chance:?}");
    }

    /// Chance puts 2 changes on a position of 2,000 molecules a hundred times
    /// as often as on one of 200, so with Poisson damage at 1e-4 over 50,000
    /// positions of each depth and 100 real mutations of 2 molecules at each,
    /// chance explains 9% of the shallow positions with 2 changes and 93% of
    /// the deep ones, where over both depths it would explain 84% of either.
    #[test]
    fn test_the_chance_model_compares_positions_of_like_depth() {
        let stratum = poisson_stratum(&[200, 2000], 50_000.0, 1e-4, 100);
        let chance = stratum.chance();
        let share = |(expected, observed): (f64, u64)| expected / observed as f64;
        let shallow = share(chance.at_depth(200, 2));
        assert!((shallow - 0.093).abs() < 0.005, "{shallow} {chance:?}");
        let deep = share(chance.at_depth(2000, 2));
        assert!((deep - 0.926).abs() < 0.005, "{deep} {chance:?}");
        let pooled = share(chance.at(2));
        assert!((pooled - 0.837).abs() < 0.005, "{pooled} {chance:?}");
        assert_eq!(chance.at_depth(150, 2), chance.at_depth(200, 2));
        assert_eq!(chance.at_depth(5000, 2), (0.0, 0));
    }

    /// Chance puts changes on a position only where they would not make it
    /// germline, as the observed positions are: two of ten molecules would,
    /// as would four of twenty, but not two or three of twenty.
    #[test]
    fn test_chance_expects_no_position_it_would_make_germline() {
        let mut stratum = StratumProfile::default();
        for n in [10, 20] {
            stratum.positions.insert((n, 0), 900);
            stratum.positions.insert((n, 1), 90);
        }
        let chance = stratum.chance_at(f64::INFINITY);
        assert_eq!(chance.at_depth(10, 2).0, 0.0);
        assert!(chance.at_depth(20, 2).0 > 1.0, "{chance:?}");
        assert!(chance.at_depth(20, 3).0 > 0.01, "{chance:?}");
        assert_eq!(chance.at_depth(20, 4).0, 0.0);
    }

    /// From [`MAX_CHANGES`] on, a count takes the positions with it or more,
    /// so a call with 40 changes, far beyond chance, finds only the clonal
    /// positions with as many.
    #[test]
    fn test_deep_counts_take_the_positions_with_as_many_or_more() {
        let tally = Tally {
            expected: vec![
                90.0, 9.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.25, 0.125,
            ],
            observed: vec![90, 9, 2, 1, 0, 0, 0, 0, 3, 1, 0, 0, 2],
        };
        assert_eq!(tally.at(2), (1.0, 2));
        assert_eq!(tally.at(7), (0.03125, 0));
        assert_eq!(tally.at(8), (0.375, 6));
        assert_eq!(tally.at(9), (0.125, 3));
        assert_eq!(tally.at(12), (0.0, 2));
        assert_eq!(tally.at(40), (0.0, 0));

        let mut stratum = poisson_stratum(&[1000], 100_000.0, 1e-4, 0);
        stratum.positions.insert((1000, 40), 3);
        let (expected, observed) = stratum.chance().at_depth(1000, 40);
        assert!(expected < 1e-12 && observed == 3, "{expected} {observed}");
    }

    /// Single-strand changes that vary more than the duplex changes allow
    /// would have chance explain more positions than there are, so the model
    /// raises their shape until it does not, and reports the misfit of any
    /// shape that still would.
    #[test]
    fn test_the_chance_model_varies_no_more_than_the_duplex_changes_allow() {
        let mut stratum = StratumProfile::default();
        for k in 0..=MAX_CHANGES {
            let poisson = (100_000.0 * negative_binomial(k, 0.1, f64::INFINITY)).round();
            stratum.positions.insert((1000, k), poisson as u64);
            let varied = (100_000.0 * negative_binomial(k, 0.1, 0.2)).round();
            stratum.strand_positions.insert((1000, k), varied as u64);
        }
        let fitted = stratum.single_strand_dispersion();
        assert!((fitted - 0.2).abs() < 0.02, "{fitted}");
        let wide = stratum.chance_at(fitted);
        assert!(!wide.fits() && wide.excess > 1.0, "{wide:?}");
        let chance = stratum.chance();
        assert!(chance.fits(), "{chance:?}");
        assert_eq!(chance.fitted, fitted);
        assert!(chance.dispersion > 10.0 * fitted, "{chance:?}");
        let (expected, observed) = chance.at(2);
        assert!(expected <= observed as f64 + FIT_TOLERANCE * expected.sqrt());
    }

    /// The rate is where the positions expected with no change meet those
    /// observed, which a count of ones, met at a low rate and again at a
    /// high one, could not fix, and it is at most one change per molecule.
    #[test]
    fn test_the_matched_rate_reproduces_the_positions_with_no_change() {
        let depths = BTreeMap::from([(1, 100_000.0), (10_000, 1_000.0)]);
        let at = |k, rate| expected_at(&depths, k, rate, f64::INFINITY);
        for truth in [3e-5, 1e-3, 0.5] {
            let rate = matched_rate(&depths, at(0, truth), f64::INFINITY);
            assert!((rate / truth - 1.0).abs() < 1e-6, "{rate} {truth}");
        }
        assert_eq!(matched_rate(&depths, 101_000.0, f64::INFINITY), 0.0);
        assert_eq!(matched_rate(&depths, 0.0, f64::INFINITY), 1.0);
        assert_eq!(matched_rate(&BTreeMap::new(), 0.0, f64::INFINITY), 0.0);
    }

    /// Positions that expect two changes each meet their count of ones at
    /// a rate below the true one too, where a Poisson expects more zeros than
    /// there are: the chance rate and the single-strand shape both hold,
    /// at a Poisson and at a shape near which the ones barely move.
    #[test]
    fn test_the_fits_hold_where_positions_expect_more_than_one_change() {
        for (dispersion, tolerance) in [(f64::INFINITY, 0.0), (1.5, 0.02), (0.3, 0.02)] {
            let mut stratum = StratumProfile::default();
            for n in [500, 1000, 3000] {
                let mean = f64::from(n) * 2e-3;
                for k in 0..200 {
                    let count = (100_000.0 * negative_binomial(k, mean, dispersion)).round();
                    if count > 0.0 && !germline(k, n) {
                        stratum.positions.insert((n, k), count as u64);
                        stratum.strand_positions.insert((n, k), count as u64);
                    }
                }
            }
            let fitted = stratum.single_strand_dispersion();
            assert!(
                fitted == dispersion || (fitted - dispersion).abs() <= tolerance * dispersion,
                "{fitted} {dispersion}"
            );
            let chance = stratum.chance();
            assert!((chance.rate / 2e-3 - 1.0).abs() < 0.01, "{chance:?}");
            assert!(chance.fits() && chance.excess == 0.0, "{chance:?}");
        }
    }

    /// Positions deeper than 64 molecules pool at their mean depth within
    /// steps of a 64th of a doubling, keeping every position.
    #[test]
    fn test_deep_positions_pool_at_their_mean_depth() {
        let positions = Positions::from([
            ((10, 0), 3),
            ((10, 1), 1),
            ((1000, 0), 1),
            ((1001, 2), 3),
            ((2000, 0), 2),
        ]);
        let expected = BTreeMap::from([(10, 4.0), (1001, 4.0), (2000, 2.0)]);
        assert_eq!(super::depths(&positions), expected);
    }

    #[test]
    fn test_the_negative_binomial_sums_to_one_and_tends_to_a_poisson() {
        for (mean, dispersion) in [(0.3, 1.5), (2.0, 0.5), (0.01, 10.0)] {
            let total: f64 = (0..200)
                .map(|k| negative_binomial(k, mean, dispersion))
                .sum();
            assert!((total - 1.0).abs() < 1e-9, "{mean} {dispersion}");
        }
        let poisson = negative_binomial(2, 0.5, f64::INFINITY);
        assert!((poisson - 0.125 * (-0.5f64).exp()).abs() < 1e-12);
        assert!((negative_binomial(2, 0.5, 1e9) - poisson).abs() < 1e-8);
        assert_eq!(negative_binomial(0, 0.0, 1.0), 1.0);
        assert_eq!(negative_binomial(1, 0.0, 1.0), 0.0);
    }
}
