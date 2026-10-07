//! What a duplex BAM's single-strand consensus says about a library's damage.
//!
//! A duplex consensus from fgbio's or fgumi's `CallDuplexConsensusReads`
//! carries, per base, the single-strand consensus of each of the molecule's two
//! strands: `ac` and `bc` hold their bases and `ad` and `bd` the raw reads
//! behind them. Once aligned, as by `ZipperBams`, the bases are reverse
//! complemented with the consensus, so they sit in the reference's
//! orientation. A consensus counts only with all four tags, read past any
//! hard clip they still hold. A consensus whose mate overlaps it counts only
//! outside its mate's span, from the mate's `MC` tag, so each molecule counts
//! once.
//!
//! A library profile reads every consensus once and counts, at each reference base
//! a damage class can change, two kinds of molecule:
//!
//! - **Single-strand changes:** one strand holds the reference base and the
//!   other the class's damaged base, read by at least [`MIN_STRAND_READS`] raw
//!   reads. These are lesions polymerase never copied.
//! - **Duplex changes:** the consensus base, at the base quality floor, is the
//!   damaged base. These are copied lesions and real mutations alike.
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
//! varies as a gamma of shape `a`, fitted so the positions expected with no
//! change and with one both match those observed. Copying adds variation of
//! its own that single-strand changes cannot show, since a lesion is copied
//! only near where a fragment ends or a nick opens, so the chance model
//! takes `a` at most [`MAX_DISPERSION`], a rate at least as varied as an
//! exponential's. With `n_j` molecules at
//! position `j` and `S(k)` positions holding `k` duplex changes, the duplex
//! rate `r` solves `sum_j NB(1; n_j r, a) = S(1)`, and chance puts `k`
//! changes on
//!
//! ```text
//! E(k) = sum_j NB(k; n_j r, a)
//! ```
//!
//! positions. `E(k) / S(k)` is the share of positions with `k` changes that
//! chance explains: few of a clean library's positions with two or more, and
//! most of a damaged library's. Positions whose changes are at least 2 and
//! 20% of their molecules are germline and count toward neither model.

use std::collections::{BTreeMap, VecDeque};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{bail, Context as _, Result};
use log::{info, warn};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::value::Array;
use noodles::sam::alignment::record::data::field::Value;

use crate::classes::{complement, Context, DamageClass};
use crate::evidence::PileupOptions;
use crate::reference::Reference;

/// The fewest raw reads behind a strand's base for its change to count as a
/// single-strand change, so a single read's error does not.
pub const MIN_STRAND_READS: i64 = 2;

/// The most changes at a position the chance model tells apart; positions
/// with more pool with it.
pub const MAX_CHANGES: u32 = 8;

/// The fewest changes at a position the chance model gives a call's prior
/// for: its rate matches the positions with one, so it explains them all.
pub const MIN_CHANCE_CHANGES: u32 = 2;

/// The records a library profile reads before it decides the BAM has no
/// single-strand consensus.
const TAG_PROBE: u64 = 1_000;

/// The reference bases a library profile holds in memory at once.
const REFERENCE_CHUNK: usize = 1 << 20;

/// The longest reference span of a consensus the profile counts, beyond which it
/// is skipped rather than held.
const MAX_SPAN: usize = 100_000;

/// The largest gamma shape the chance model takes: copying a lesion varies
/// across positions at least as much as an exponential does, more than the
/// single-strand changes alone show, since it also depends on where
/// fragments end.
pub const MAX_DISPERSION: f64 = 1.0;

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
    /// Molecules with both strands' bases called.
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
        let single = self.single_strand_rate()?;
        (single > 0.0).then(|| self.change_rate().unwrap_or(0.0) / single)
    }

    /// The stratum's chance model: the dispersion fitted to its single-strand
    /// changes, at most [`MAX_DISPERSION`], the rate of duplex changes it
    /// matches to the positions with one, and the positions expected and
    /// observed with each count.
    pub fn chance(&self) -> Chance {
        self.chance_at(self.single_strand_dispersion().min(MAX_DISPERSION))
    }

    /// The gamma shape of the single-strand change rate's variation across
    /// positions, infinite when it varies no more than a Poisson allows.
    pub fn single_strand_dispersion(&self) -> f64 {
        fitted_dispersion(&self.strand_positions)
    }

    /// The chance model at a given dispersion.
    pub fn chance_at(&self, dispersion: f64) -> Chance {
        let rate = matched_rate(&self.positions, dispersion);
        let mut expected = vec![0.0; MAX_CHANGES as usize + 1];
        for (&(n, _), &count) in &self.positions {
            let mean = f64::from(n) * rate;
            let mut below = 0.0;
            for (k, e) in expected.iter_mut().enumerate().take(MAX_CHANGES as usize) {
                let p = negative_binomial(k as u32, mean, dispersion);
                below += p;
                *e += count as f64 * p;
            }
            expected[MAX_CHANGES as usize] += count as f64 * (1.0 - below).max(0.0);
        }
        let observed = (0..=MAX_CHANGES)
            .map(|k| observed(&self.positions, k))
            .collect();
        Chance {
            rate,
            dispersion,
            expected,
            observed,
        }
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
    /// `E(k)` for `k` from 0 to [`MAX_CHANGES`], the last pooling deeper ones.
    pub expected: Vec<f64>,
    /// `S(k)` likewise.
    pub observed: Vec<u64>,
}

impl Chance {
    /// The expected and observed positions with `changes` changes.
    pub fn at(&self, changes: u32) -> (f64, u64) {
        let k = changes.min(MAX_CHANGES) as usize;
        (self.expected[k], self.observed[k])
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

/// The positions expected with `k` changes at `rate` per molecule.
fn expected_at(positions: &Positions, k: u32, rate: f64, dispersion: f64) -> f64 {
    positions
        .iter()
        .map(|(&(n, _), &count)| {
            count as f64 * negative_binomial(k, f64::from(n) * rate, dispersion)
        })
        .sum()
}

/// The rate per molecule at which the positions expected with one change
/// match those observed, on the rising side of that count, or where it
/// peaks when it never rises that high; zero without such positions.
fn matched_rate(positions: &Positions, dispersion: f64) -> f64 {
    let ones = observed(positions, 1) as f64;
    let deepest = positions.keys().map(|&(n, _)| n).max().unwrap_or(0);
    if ones == 0.0 || deepest == 0 {
        return 0.0;
    }
    let singles = |ln_rate: f64| expected_at(positions, 1, ln_rate.exp(), dispersion);
    let (mut low, mut high) = ((1e-12f64).ln(), (100.0 / f64::from(deepest)).ln());
    let mut peak = (low, high);
    for _ in 0..100 {
        let third = (peak.1 - peak.0) / 3.0;
        if singles(peak.0 + third) < singles(peak.1 - third) {
            peak.0 += third;
        } else {
            peak.1 -= third;
        }
    }
    high = (peak.0 + peak.1) / 2.0;
    if singles(high) <= ones {
        return high.exp();
    }
    for _ in 0..100 {
        let mid = (low + high) / 2.0;
        if singles(mid) < ones {
            low = mid;
        } else {
            high = mid;
        }
    }
    ((low + high) / 2.0).exp()
}

/// The gamma shape of a rate's variation across positions, fitted so the
/// positions expected with no change match those observed once the rate
/// matches the positions with one: infinite when they vary no more than a
/// Poisson allows. Single-strand changes give it free of real mutations,
/// which only duplex changes carry. Expected zeros fall as the shape falls
/// from a Poisson's until the rate can no longer match the ones, so the
/// shape is the first crossing on the way down.
fn fitted_dispersion(positions: &Positions) -> f64 {
    const RANGE: (f64, f64) = (0.05, 1e4);
    const STEPS: usize = 48;
    let zeros = observed(positions, 0) as f64;
    if observed(positions, 1) == 0 {
        return f64::INFINITY;
    }
    let excess = |ln_dispersion: f64| {
        let dispersion = ln_dispersion.exp();
        expected_at(
            positions,
            0,
            matched_rate(positions, dispersion),
            dispersion,
        ) - zeros
    };
    let (top, bottom) = (RANGE.1.ln(), RANGE.0.ln());
    if excess(top) <= 0.0 {
        return f64::INFINITY;
    }
    let step = (top - bottom) / STEPS as f64;
    let mut high = top;
    for i in 1..=STEPS {
        let low = top - step * i as f64;
        if excess(low) < 0.0 {
            let (mut low, mut high) = (low, high);
            for _ in 0..50 {
                let mid = (low + high) / 2.0;
                if excess(mid) < 0.0 {
                    low = mid;
                } else {
                    high = mid;
                }
            }
            return ((low + high) / 2.0).exp();
        }
        high = low;
    }
    RANGE.0
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
            if site.molecules > 0 && (changes < 2 || share < GERMLINE_FRACTION) {
                *stratum
                    .positions
                    .entry((site.molecules, changes))
                    .or_default() += 1;
            }
            let single = site.single_strand[a];
            let single_share = f64::from(single) / f64::from(site.strand_molecules.max(1));
            if site.strand_molecules > 0 && (single < 2 || single_share < GERMLINE_FRACTION) {
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
        let mate_span = mate_span(record, contig_id)?;
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
                if seq != reference {
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
        site.strand_molecules += 1;
        let reads = |r: &[i64]| r.get(query).is_some_and(|&n| n >= MIN_STRAND_READS);
        if b == reference && a != reference && reads(strands.a_reads) {
            site.single_strand[ia] += 1;
        } else if a == reference && b != reference && reads(strands.b_reads) {
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

/// The 0-based, half-open reference span of a second of pair's mate, from
/// its position and `MC` tag, which the consensus leaves to its mate so a
/// molecule counts once where its mates overlap.
fn mate_span(record: &bam::Record, contig_id: usize) -> Result<Option<(usize, usize)>> {
    let flags = record.flags();
    if !flags.is_segmented() || flags.is_mate_unmapped() || !flags.is_last_segment() {
        return Ok(None);
    }
    let (Some(Ok(mate_contig)), Some(Ok(mate_start))) = (
        record.mate_reference_sequence_id(),
        record.mate_alignment_start(),
    ) else {
        return Ok(None);
    };
    if mate_contig != contig_id {
        return Ok(None);
    }
    let Some(Ok(Value::String(cigar))) = record.data().get(b"MC") else {
        return Ok(None);
    };
    let start = usize::from(mate_start) - 1;
    Ok(Some((start, start + reference_length(cigar))))
}

/// The reference bases a CIGAR string spans.
fn reference_length(cigar: &[u8]) -> usize {
    let mut length = 0;
    let mut number = 0;
    for &c in cigar {
        if c.is_ascii_digit() {
            number = number * 10 + usize::from(c - b'0');
        } else {
            if matches!(c, b'M' | b'D' | b'N' | b'=' | b'X') {
                length += number;
            }
            number = 0;
        }
    }
    length
}

/// Profile the single-strand and duplex changes of a coordinate-sorted BAM,
/// or `None` when its records carry no single-strand consensus, carry it
/// unaligned, or `stop` is set before the last is read.
pub fn profile_library(
    bam: &Path,
    reference: &Path,
    classes: &[DamageClass],
    options: &PileupOptions,
    stop: &AtomicBool,
) -> Result<Option<LibraryProfile>> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(bam)
        .with_context(|| format!("failed to open BAM: {bam:?}"))?;
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
    use crate::testing::write_fasta;

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

    /// Of ten consensus over positions 1 to 40, one holds a duplex C>T at the C of
    /// a CpG, one a C>T on one strand at a C outside CpG, and one the same
    /// change read by a single raw read, which does not count.
    #[test]
    fn test_a_profile_counts_duplex_and_single_strand_changes_per_stratum() {
        let mut reads = SamBuilder::new().read_length(40);
        reads.add_frag(frag(&[(2, b'T', 30, b'T', b'T', 3)], Strand::Plus));
        reads.add_frag(frag(&[(6, b'N', 2, b'T', b'C', 3)], Strand::Plus));
        reads.add_frag(frag(&[(6, b'N', 2, b'T', b'C', 1)], Strand::Plus));
        for _ in 0..7 {
            reads.add_frag(frag(&[], Strand::Plus));
        }
        let profile = profile(&reads).unwrap();
        let cpg = profile.stratum("C>T:CpG").unwrap();
        let other = profile.stratum("C>T:non-CpG").unwrap();
        // Positions 1 to 40 hold five units: a CpG C and G and one other C each.
        assert_eq!(cpg.molecules, 100);
        assert_eq!(cpg.changes, 1);
        assert_eq!(cpg.strand_molecules, 100);
        assert_eq!(cpg.single_strand_changes, 0);
        assert_eq!(cpg.positions.get(&(10, 1)), Some(&1));
        assert_eq!(cpg.positions.get(&(10, 0)), Some(&9));
        assert_eq!(other.molecules, 48);
        assert_eq!(other.changes, 0);
        assert_eq!(other.strand_molecules, 50);
        assert_eq!(other.single_strand_changes, 1);
        assert_eq!(other.single_strand_rate(), Some(0.02));
        assert_eq!(cpg.conversion_ratio(), None);
        let oxidation = profile.stratum("G>T:CpG").unwrap();
        assert_eq!((oxidation.molecules, oxidation.changes), (100, 0));
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

    /// Mates over the same positions are one molecule there: the second of pair
    /// leaves its mate's span to it.
    #[test]
    fn test_overlapping_mates_count_once() {
        let mut reads = SamBuilder::new().read_length(40);
        let mut pair = Pair::at(1, 1)
            .bases1(&reference()[..40])
            .bases2(&reference()[..40]);
        for (tag, value) in tags(&reference()[..40], None) {
            pair = pair.attr(tag, value);
        }
        reads.add_pair(pair);
        let profile = profile(&reads).unwrap();
        assert_eq!(profile.stratum("C>T:CpG").unwrap().molecules, 10);
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
    /// chance model at the cap, and a Poisson chance model's rate reproduces
    /// the positions with one change.
    #[test]
    fn test_the_chance_model_matches_the_positions_with_one_change() {
        let mut stratum = StratumProfile::default();
        for (k, count) in [(0, 900), (1, 90), (2, 8), (9, 2)] {
            stratum.positions.insert((100, k), count);
        }
        stratum.strand_positions.insert((100, 0), 990);
        stratum.strand_positions.insert((100, 1), 10);
        assert!(stratum.single_strand_dispersion().is_infinite());
        assert_eq!(stratum.chance().dispersion, MAX_DISPERSION);
        let chance = stratum.chance_at(f64::INFINITY);
        let mean = 100.0 * chance.rate;
        assert!(
            (1000.0 * mean * (-mean).exp() - 90.0).abs() < 1e-6,
            "{chance:?}"
        );
        assert!((chance.at(1).0 - 90.0).abs() < 1e-6);
        assert!((chance.at(2).0 - 1000.0 * mean * mean / 2.0 * (-mean).exp()).abs() < 1e-6);
        assert_eq!(chance.at(2).1, 8);
        assert_eq!(chance.at(12), chance.at(MAX_CHANGES));
        assert_eq!(chance.at(MAX_CHANGES).1, 2);
        assert!(chance.at(MAX_CHANGES).0 < 1e-6);
        assert_eq!(StratumProfile::default().chance().rate, 0.0);
    }

    /// Single-strand changes drawn from a gamma-varying rate give back its
    /// shape, duplex changes sharing it fall together more often than a
    /// Poisson allows, and the chance model varies at least as much as the
    /// cap.
    #[test]
    fn test_the_dispersion_of_single_strand_changes_widens_the_chance_model() {
        let mut stratum = StratumProfile::default();
        for u in 0..12u32 {
            let count = (100_000.0 * negative_binomial(u, 0.3, 1.5)).round() as u64;
            stratum.strand_positions.insert((1000, u), count);
        }
        for (k, count) in [(0, 70_000), (1, 21_000), (2, 6_000), (3, 3_000)] {
            stratum.positions.insert((1000, k), count);
        }
        let dispersion = stratum.single_strand_dispersion();
        assert!((dispersion - 1.5).abs() < 0.05, "{dispersion}");
        let chance = stratum.chance_at(dispersion);
        let poisson = stratum.chance_at(f64::INFINITY);
        assert!(
            chance.at(2).0 > 1.3 * poisson.at(2).0,
            "{chance:?} {poisson:?}"
        );
        assert!((chance.at(1).0 - 21_000.0).abs() < 1e-3);
        let capped = stratum.chance();
        assert_eq!(capped.dispersion, MAX_DISPERSION);
        assert!(capped.at(2).0 > chance.at(2).0);
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
