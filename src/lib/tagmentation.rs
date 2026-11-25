//! The tagmentation filter: errors a transposase library carries at fixed
//! positions from a template end.
//!
//! A transposase cuts a fragment and joins its adapters in one step, leaving a
//! short gap at each end that a polymerase fills in. The bases at fixed
//! positions from a template's end, the first above all, can then read as
//! another base more often than the bases inside it, on both strands alike, so
//! a duplex consensus agrees on the change. Which positions and which
//! substitutions differs between chemistries, so the filter fixes none: it
//! learns them from the library.
//!
//! # The end profile
//!
//! An [`EndProfile`] reads every record once and counts, at each of the first
//! [`PROFILE_LENGTH`] bases from the record's own 5' end, the reference base
//! against the base read, in the read's orientation, with a `C` or a `G` kept
//! apart by CpG context. A record counts when its first [`PROFILE_LENGTH`]
//! bases align as one block without a clip, so that its first base is the
//! template's. A substitution counts once per site, strand, position, and base
//! read, so an allele that many molecules share weighs what one molecule's
//! error does.
//!
//! # The transposase's preference
//!
//! A transposase prefers some bases at each position near its cut, and it
//! reads the sample's base, not the reference's. A real allele `r>a` is
//! therefore found at position `p` more often than the reference base it
//! replaces, by
//!
//! ```text
//! fold(p) = preference(p, a) / preference(p, r)
//! preference(p, b) = share of b among the reference bases at p
//!                  / share of b among those of the interior
//! ```
//!
//! with no error involved, where the interior is the positions past
//! [`WINDOW`]. A substitution's rate at a position is thus expected to be
//! `interior rate * max(1, fold(p))`: every interior substitution is taken
//! for a real allele the preference can enrich, and a substitution the
//! preference depletes is never an excess, so the excess is biased low.
//!
//! # Learned positions
//!
//! A position within [`WINDOW`] is learned for a substitution when it has at
//! least [`MIN_EVENTS`] substitutions, at [`MIN_FOLD`] times the expected rate
//! or more, and [`MIN_Z`] standard deviations of a Poisson count above the
//! expected count. Its weight is its excess rate over the expected one. A
//! library without such a position learns nothing, and none of its calls is
//! scored.
//!
//! # Scoring a call
//!
//! A molecule whose site lies `left` bases from its template's leftmost base
//! and `right` from its rightmost reads the call `r>a` at the left end and its
//! reverse complement at the right, so it takes the excess `x` and the
//! preference fold `s` of both:
//!
//! ```text
//! x = excess(r>a, left) + excess(reverse complement, right)
//! s = fold(r>a, left) * fold(reverse complement, right)
//! ```
//!
//! each zero or one outside [`WINDOW`]. Under the artifact a call's alternate
//! molecules follow its reference molecules tilted by `x`, and under a
//! mutation tilted by `s`, so with `X` and `S` the means of `x` and `s` over
//! the reference molecules and `e` a base error probability, each alternate
//! molecule contributes
//!
//! ```text
//! ln((1 - e) (x / X) / (s / S) + e)
//! ```
//!
//! `X` and `S` are shrunk by [`REFERENCE_PSEUDO_MOLECULES`] toward the pooled
//! reference molecules of the calls of the same substitution, and those toward
//! a site spread evenly along its templates, so a call without a reference
//! molecule is scored from the pool: where a template ends does not depend on
//! the site.
//!
//! A call none of whose alternate molecules sits at a learned position has no
//! molecule the artifact could have made, so its posterior is one whatever its
//! prior.

use std::collections::{BTreeMap, HashSet};
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::num::NonZero;
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{bail, Context as _, Result};
use log::info;
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles_bgzf as bgzf;

use crate::call::Genotype;
use crate::classes::{complement, sbs6, Context};
use crate::evidence::{Molecule, PileupOptions};
use crate::io::StagedFile;
use crate::read_end::{error_probability, Score, REFERENCE_PSEUDO_MOLECULES};
use crate::reference::Reference;

/// The bases from a record's 5' end that a profile counts.
pub const PROFILE_LENGTH: usize = 20;

/// The positions from a template end that can be learned; those past it are
/// the interior the others are compared with.
pub const WINDOW: usize = 12;

/// The fewest substitutions at a position for it to be learned.
pub const MIN_EVENTS: u64 = 50;

/// The least ratio of a position's rate to its expected rate for it to be
/// learned.
pub const MIN_FOLD: f64 = 1.5;

/// The fewest standard deviations of a Poisson count by which a position's
/// substitutions exceed those expected for it to be learned.
pub const MIN_Z: f64 = 8.0;

/// The least preference fold a molecule is scored with, so a base never seen
/// at a position cannot divide a ratio by zero.
const MIN_PREFERENCE: f64 = 1e-3;

/// The template length a site is spread along when no molecule gives one.
const NOMINAL_LENGTH: f64 = 200.0;

/// The threads that decompress the BAM a profile reads.
const DECOMPRESSION_THREADS: NonZero<usize> = NonZero::new(2).unwrap();

/// The reference bases a profile reads at once.
const REFERENCE_CHUNK: usize = 1 << 20;

const BASES: [u8; 4] = *b"ACGT";

fn index(base: u8) -> Option<usize> {
    BASES.iter().position(|&b| b == base)
}

/// A substitution in a read's orientation, with the CpG context of a `C>T` or
/// a `G>A`, whose real alleles are far more common at a CpG.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct EndClass {
    /// The reference base, as the read's strand reads it.
    pub reference: u8,
    /// The base read.
    pub read: u8,
    /// The context of a `C>T` or a `G>A`, and `None` for any other.
    pub context: Option<Context>,
}

impl EndClass {
    /// The class of `reference>read` at a base in `context`.
    pub fn of(reference: u8, read: u8, context: Context) -> Self {
        let split = matches!((reference, read), (b'C', b'T') | (b'G', b'A'));
        Self {
            reference,
            read,
            context: split.then_some(context),
        }
    }

    /// The same change as the other strand reads it.
    pub fn reverse_complement(self) -> Self {
        Self {
            reference: complement(self.reference),
            read: complement(self.read),
            context: self.context,
        }
    }

    /// The stratum its calls learn their artifact fraction in: the
    /// substitution from its pyrimidine, with the context of a `C>T`.
    pub fn stratum(self) -> String {
        let substitution = sbs6(self.reference, self.read).unwrap_or_else(|| "other".to_string());
        match self.context {
            Some(context) => format!("{substitution}:{context}"),
            None => substitution,
        }
    }

    /// Every class, in order.
    pub fn all() -> Vec<EndClass> {
        let mut classes = Vec::new();
        for reference in BASES {
            for read in BASES.into_iter().filter(|&read| read != reference) {
                for context in Context::ALL {
                    let class = Self::of(reference, read, context);
                    if !classes.contains(&class) {
                        classes.push(class);
                    }
                }
            }
        }
        classes
    }
}

impl fmt::Display for EndClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}>{}", self.reference as char, self.read as char)?;
        match self.context {
            Some(context) => write!(f, ":{context}"),
            None => Ok(()),
        }
    }
}

/// Every stratum the filter can score a call in.
pub fn strata() -> Vec<String> {
    let mut strata: Vec<String> = EndClass::all().into_iter().map(EndClass::stratum).collect();
    strata.sort();
    strata.dedup();
    strata
}

/// A library's bases by position from a record's 5' end: the reference base
/// against the base read, in the read's orientation, by CpG context.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EndProfile {
    counts: Vec<u64>,
}

impl Default for EndProfile {
    fn default() -> Self {
        Self {
            counts: vec![0; PROFILE_LENGTH * 4 * 4 * 2],
        }
    }
}

impl EndProfile {
    fn slot(position: usize, reference: u8, read: u8, context: Context) -> Option<usize> {
        if position >= PROFILE_LENGTH {
            return None;
        }
        let context = usize::from(context == Context::CpG);
        Some(((position * 4 + index(reference)?) * 4 + index(read)?) * 2 + context)
    }

    /// Count `count` bases at a 0-based `position` from a 5' end where
    /// `reference` read as `read`; bases other than `ACGT` and positions past
    /// [`PROFILE_LENGTH`] are left out.
    pub fn add(&mut self, position: usize, reference: u8, read: u8, context: Context, count: u64) {
        if let Some(slot) = Self::slot(position, reference, read, context) {
            self.counts[slot] += count;
        }
    }

    /// The bases counted at a position where `reference` read as `read`.
    pub fn count(&self, position: usize, reference: u8, read: u8, context: Context) -> u64 {
        Self::slot(position, reference, read, context).map_or(0, |slot| self.counts[slot])
    }

    /// The bases counted in all.
    pub fn bases(&self) -> u64 {
        self.counts.iter().sum()
    }

    fn class_count(&self, position: usize, reference: u8, read: u8, class: EndClass) -> u64 {
        match class.context {
            Some(context) => self.count(position, reference, read, context),
            None => Context::ALL
                .iter()
                .map(|&context| self.count(position, reference, read, context))
                .sum(),
        }
    }

    /// Learn each class's interior rate, preference folds, and excess rates.
    pub fn learn(&self) -> Weights {
        let depth = |position: usize, base: u8| -> u64 {
            BASES
                .iter()
                .flat_map(|&read| Context::ALL.map(|context| (read, context)))
                .map(|(read, context)| self.count(position, base, read, context))
                .sum()
        };
        let total = |position: usize| -> u64 { BASES.iter().map(|&b| depth(position, b)).sum() };
        let interior_total: u64 = (WINDOW..PROFILE_LENGTH).map(total).sum();
        let preference = |position: usize, base: u8| -> f64 {
            let inside: u64 = (WINDOW..PROFILE_LENGTH).map(|p| depth(p, base)).sum();
            let here = total(position);
            if inside == 0 || here == 0 || interior_total == 0 {
                return 0.0;
            }
            (depth(position, base) as f64 / here as f64) / (inside as f64 / interior_total as f64)
        };
        let mut classes = BTreeMap::new();
        for class in EndClass::all() {
            let events = |p: usize| self.class_count(p, class.reference, class.read, class);
            let bases = |p: usize| -> u64 {
                BASES
                    .iter()
                    .map(|&read| self.class_count(p, class.reference, read, class))
                    .sum()
            };
            let interior_bases: u64 = (WINDOW..PROFILE_LENGTH).map(bases).sum();
            if interior_bases == 0 {
                continue;
            }
            let interior_events: u64 = (WINDOW..PROFILE_LENGTH).map(events).sum();
            let background = interior_events as f64 / interior_bases as f64;
            let mut weights = ClassWeights {
                background,
                ..ClassWeights::default()
            };
            for position in 0..WINDOW {
                let (n, k) = (bases(position), events(position));
                let from = preference(position, class.reference);
                if n == 0 || from == 0.0 {
                    continue;
                }
                let fold = preference(position, class.read) / from;
                let rate = k as f64 / n as f64;
                let expected = background * fold.max(1.0);
                let count = expected * n as f64;
                let exceeds = if count > 0.0 {
                    (k as f64 - count) / count.sqrt() >= MIN_Z
                } else {
                    true
                };
                weights.preference[position] = fold;
                weights.rate[position] = rate;
                weights.expected[position] = expected;
                weights.events[position] = k;
                if k >= MIN_EVENTS && rate >= MIN_FOLD * expected && exceeds {
                    weights.excess[position] = rate - expected;
                }
            }
            classes.insert(class, weights);
        }
        Weights { classes }
    }

    /// Write the counts as a tab-separated file that [`EndProfile::read`]
    /// reads back.
    pub fn write(&self, path: &Path) -> Result<()> {
        let staged = StagedFile::create(path)?;
        let mut writer = BufWriter::new(staged.writer()?);
        writeln!(writer, "{}", PROFILE_COLUMNS.join("\t"))?;
        for position in 0..PROFILE_LENGTH {
            for reference in BASES {
                for read in BASES {
                    for context in Context::ALL {
                        let count = self.count(position, reference, read, context);
                        if count > 0 {
                            writeln!(
                                writer,
                                "{}\t{}\t{}\t{context}\t{count}",
                                position + 1,
                                reference as char,
                                read as char
                            )?;
                        }
                    }
                }
            }
        }
        writer.flush()?;
        drop(writer);
        staged.persist()
    }

    /// Read the counts another run wrote with [`EndProfile::write`].
    pub fn read(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("failed to open the tagmentation profile: {path:?}"))?;
        let mut lines = BufReader::new(file).lines();
        let header = lines.next().transpose()?.unwrap_or_default();
        if header.split('\t').collect::<Vec<_>>() != PROFILE_COLUMNS {
            bail!(
                "the tagmentation profile {path:?} must start with the columns: {}",
                PROFILE_COLUMNS.join(", ")
            );
        }
        let mut profile = Self::default();
        for (number, line) in lines.enumerate() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            let row = (|| {
                let [position, reference, read, context, count] = fields.as_slice() else {
                    return None;
                };
                let position: usize = position.parse().ok()?;
                let base = |text: &str| match text.as_bytes() {
                    [base] if index(*base).is_some() => Some(*base),
                    _ => None,
                };
                let context = Context::ALL
                    .into_iter()
                    .find(|c| c.to_string() == *context)?;
                let slot = Self::slot(
                    position.checked_sub(1)?,
                    base(reference)?,
                    base(read)?,
                    context,
                )?;
                Some((slot, count.parse::<u64>().ok()?))
            })();
            let Some((slot, count)) = row else {
                bail!(
                    "line {} of the tagmentation profile {path:?} is not a position from 1 to {PROFILE_LENGTH}, two bases, a context, and a count: {line}",
                    number + 2
                );
            };
            profile.counts[slot] += count;
        }
        Ok(profile)
    }
}

/// The columns of a profile file.
const PROFILE_COLUMNS: [&str; 5] = ["position", "ref_base", "read_base", "context", "count"];

/// What a profile says of one class at each position within [`WINDOW`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ClassWeights {
    /// The class's substitutions per base in the interior.
    pub background: f64,
    /// Its substitutions per base at each position.
    pub rate: [f64; WINDOW],
    /// The rate expected at each position without an artifact.
    pub expected: [f64; WINDOW],
    /// How much more often a real allele of the class is found at each
    /// position than the reference base it replaces.
    pub preference: [f64; WINDOW],
    /// The excess rate at each learned position, zero at the others.
    pub excess: [f64; WINDOW],
    /// The substitutions counted at each position.
    pub events: [u64; WINDOW],
}

/// The weights a profile gives every class.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Weights {
    /// Each class's weights.
    pub classes: BTreeMap<EndClass, ClassWeights>,
}

impl Weights {
    /// Each learned class and 1-based position, in class order.
    pub fn learned(&self) -> Vec<(EndClass, usize)> {
        self.classes
            .iter()
            .flat_map(|(class, weights)| {
                (0..WINDOW)
                    .filter(|&p| weights.excess[p] > 0.0)
                    .map(|p| (*class, p + 1))
            })
            .collect()
    }

    /// The learned positions as text, such as `T>G@1,A>G@1`, or `None`
    /// without one.
    pub fn learned_text(&self, keep: impl Fn(EndClass) -> bool) -> Option<String> {
        let learned: Vec<String> = self
            .learned()
            .into_iter()
            .filter(|(class, _)| keep(*class))
            .map(|(class, position)| format!("{class}@{position}"))
            .collect();
        (!learned.is_empty()).then(|| learned.join(","))
    }

    /// The weights a call of `class`, as the forward strand reads it, is
    /// scored with, or `None` when neither end has a learned position.
    pub fn of(&self, class: EndClass) -> Option<CallWeights> {
        let none = ClassWeights::default();
        let left = self.classes.get(&class).unwrap_or(&none);
        let right = self
            .classes
            .get(&class.reverse_complement())
            .unwrap_or(&none);
        let any = |w: &ClassWeights| w.excess.iter().any(|&x| x > 0.0);
        if !any(left) && !any(right) {
            return None;
        }
        let folds = |w: &ClassWeights| w.preference.map(|fold| if fold > 0.0 { fold } else { 1.0 });
        Some(CallWeights {
            excess_left: left.excess,
            excess_right: right.excess,
            fold_left: folds(left),
            fold_right: folds(right),
        })
    }
}

/// The excess rates and preference folds one call's molecules take at each
/// end.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CallWeights {
    excess_left: [f64; WINDOW],
    excess_right: [f64; WINDOW],
    fold_left: [f64; WINDOW],
    fold_right: [f64; WINDOW],
}

impl CallWeights {
    fn at<const N: usize>(table: &[f64; N], distance: Option<usize>, outside: f64) -> f64 {
        distance
            .and_then(|d| table.get(d))
            .copied()
            .unwrap_or(outside)
    }

    /// A molecule's excess rate and preference fold.
    fn molecule(&self, left: Option<usize>, right: Option<usize>) -> (f64, f64) {
        let excess =
            Self::at(&self.excess_left, left, 0.0) + Self::at(&self.excess_right, right, 0.0);
        let fold = Self::at(&self.fold_left, left, 1.0) * Self::at(&self.fold_right, right, 1.0);
        (excess, fold)
    }
}

/// Molecules by how far their site lies from each template end, as far as
/// [`WINDOW`], and their templates' lengths.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ends {
    left: [u32; WINDOW],
    right: [u32; WINDOW],
    molecules: u32,
    measured: u32,
    inverse_lengths: f64,
}

impl Ends {
    fn add(&mut self, molecule: &Molecule) {
        self.molecules += 1;
        if let Some(count) = molecule.left.and_then(|d| self.left.get_mut(d)) {
            *count += 1;
        }
        if let Some(count) = molecule.right.and_then(|d| self.right.get_mut(d)) {
            *count += 1;
        }
        if let Some(length) = molecule.length() {
            self.measured += 1;
            self.inverse_lengths += 1.0 / length.max(WINDOW) as f64;
        }
    }

    /// Add another site's molecules.
    pub fn merge(&mut self, other: &Ends) {
        for (a, b) in self.left.iter_mut().zip(other.left) {
            *a += b;
        }
        for (a, b) in self.right.iter_mut().zip(other.right) {
            *a += b;
        }
        self.molecules += other.molecules;
        self.measured += other.measured;
        self.inverse_lengths += other.inverse_lengths;
    }

    /// The molecules' summed excess and fold, and how many sit at a learned
    /// position.
    fn sums(&self, weights: &CallWeights) -> (f64, f64, u32) {
        let (mut excess, mut fold, mut congruent) = (0.0, f64::from(self.molecules), 0u32);
        let ends = [
            (&self.left, &weights.excess_left, &weights.fold_left),
            (&self.right, &weights.excess_right, &weights.fold_right),
        ];
        for (counts, excesses, folds) in ends {
            for position in 0..WINDOW {
                let n = f64::from(counts[position]);
                excess += n * excesses[position];
                fold += n * (folds[position] - 1.0);
                if excesses[position] > 0.0 {
                    congruent += counts[position];
                }
            }
        }
        (excess, fold.max(0.0), congruent.min(self.molecules))
    }

    /// The mean of one over the templates' lengths, or `None` without one.
    fn inverse_length(&self) -> Option<f64> {
        (self.measured > 0).then(|| self.inverse_lengths / f64::from(self.measured))
    }
}

/// One call's molecules as the filter scores them under one parent base.
#[derive(Clone, Debug, PartialEq)]
pub struct EndEvidence {
    /// The call's class as the forward strand reads it.
    pub class: EndClass,
    /// The reference molecules.
    pub reference: Ends,
    /// The alternate molecules' distances from the leftmost and rightmost
    /// template bases, and their base qualities.
    pub alternate: Vec<(Option<usize>, Option<usize>, u8)>,
    alternate_ends: Ends,
}

impl EndEvidence {
    /// The evidence of a call `ref_base>alt_base` of `class` from its
    /// molecules.
    pub fn of(molecules: &[Molecule], ref_base: u8, alt_base: u8, class: EndClass) -> Self {
        let mut evidence = Self {
            class,
            reference: Ends::default(),
            alternate: Vec::new(),
            alternate_ends: Ends::default(),
        };
        for molecule in molecules {
            if molecule.left.is_none() && molecule.right.is_none() {
                continue;
            }
            if molecule.base == ref_base {
                evidence.reference.add(molecule);
            } else if molecule.base == alt_base {
                evidence.alternate_ends.add(molecule);
                evidence
                    .alternate
                    .push((molecule.left, molecule.right, molecule.quality));
            }
        }
        evidence
    }

    /// The call's score under `weights`, its reference molecules shrunk
    /// toward `pool`, the reference molecules of the calls of its class.
    pub fn score(&self, weights: &CallWeights, pool: Option<&Ends>) -> Score {
        let k = REFERENCE_PSEUDO_MOLECULES;
        let (excess, fold, ref_congruent) = self.reference.sums(weights);
        let inverse_length = pool
            .and_then(Ends::inverse_length)
            .or_else(|| self.reference.inverse_length())
            .or_else(|| self.alternate_ends.inverse_length())
            .unwrap_or(1.0 / NOMINAL_LENGTH);
        let total = |table: &[f64; WINDOW], outside: f64| -> f64 {
            table.iter().map(|value| value - outside).sum()
        };
        let spread_excess =
            (total(&weights.excess_left, 0.0) + total(&weights.excess_right, 0.0)) * inverse_length;
        let spread_fold = (1.0
            + (total(&weights.fold_left, 1.0) + total(&weights.fold_right, 1.0)) * inverse_length)
            .max(MIN_PREFERENCE);
        let (pool_excess, pool_fold, pool_molecules) = pool.map_or((0.0, 0.0, 0.0), |pool| {
            let (excess, fold, _) = pool.sums(weights);
            (excess, fold, f64::from(pool.molecules))
        });
        let pooled_excess = (pool_excess + k * spread_excess) / (pool_molecules + k);
        let pooled_fold = (pool_fold + k * spread_fold) / (pool_molecules + k);
        let n = f64::from(self.reference.molecules);
        let mean_excess = (excess + k * pooled_excess) / (n + k);
        let mean_fold = ((fold + k * pooled_fold) / (n + k)).max(MIN_PREFERENCE);

        let mut alt_congruent = 0;
        let mut ratio = 0.0;
        for &(left, right, quality) in &self.alternate {
            let (x, s) = weights.molecule(left, right);
            if x > 0.0 {
                alt_congruent += 1;
            }
            let e = error_probability(quality);
            let tilt = (x / mean_excess) / (s.max(MIN_PREFERENCE) / mean_fold);
            ratio += ((1.0 - e) * tilt + e).ln();
        }
        Score {
            log_likelihood_ratio: (!self.alternate.is_empty() && mean_excess > 0.0)
                .then_some(ratio),
            alt_molecules: u32::try_from(self.alternate.len()).unwrap_or(u32::MAX),
            alt_congruent,
            ref_molecules: self.reference.molecules,
            ref_congruent,
        }
    }
}

/// The tagmentation filter.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Tagmentation {
    /// A profile to score with in place of the BAM's own, as a small BAM may
    /// borrow that of a larger BAM of its library.
    pub profile: Option<EndProfile>,
    /// Where to write the profile the run scored with.
    pub profile_output: Option<std::path::PathBuf>,
}

impl Tagmentation {
    /// The INFO key for the posterior probability of a true mutation.
    pub const INFO_POSTERIOR: &'static str = "TGAP";
    /// The INFO key for the log10 likelihood ratio, artifact to mutation.
    pub const INFO_RATIO: &'static str = "TGLR";
    /// The INFO key for the alternate molecules at a learned position, and
    /// all alternate molecules measured.
    pub const INFO_ALT: &'static str = "TGAC";
    /// The INFO key for the reference molecules at a learned position, and
    /// all reference molecules measured.
    pub const INFO_REF: &'static str = "TGRC";
    /// The FILTER name.
    pub const FILTER: &'static str = "TagmentationArtifact";

    /// Heterozygous calls whose every called allele is one base.
    pub fn applies_to(gt: &Genotype) -> bool {
        gt.is_het() && gt.calls_are_single_bases()
    }
}

/// A contig's reference bases, read a chunk at a time as the records advance.
struct Chunk {
    reference: Reference,
    contig: String,
    start: usize,
    bases: Vec<u8>,
}

impl Chunk {
    /// Load the 0-based span `start..end` of `contig`, or return `false` for a
    /// contig the reference lacks.
    fn cover(&mut self, contig: &str, start: usize, end: usize) -> Result<bool> {
        if self.reference.length(contig).is_none() {
            return Ok(false);
        }
        let loaded =
            self.contig == contig && self.start <= start && end <= self.start + self.bases.len();
        if !loaded {
            self.bases = self
                .reference
                .bases(contig, start, end.max(start + REFERENCE_CHUNK))?;
            self.contig = contig.to_string();
            self.start = start;
        }
        Ok(true)
    }

    fn base(&self, pos: usize) -> Option<u8> {
        pos.checked_sub(self.start)
            .and_then(|i| self.bases.get(i))
            .copied()
    }
}

/// Counts the 5' ends of a coordinate-sorted BAM's records into an
/// [`EndProfile`].
struct EndScanner {
    options: PileupOptions,
    chunk: Chunk,
    contig: Option<(usize, String)>,
    seen: HashSet<(usize, u8, bool, u8)>,
    prune_at: usize,
    profile: EndProfile,
}

impl EndScanner {
    fn new(options: PileupOptions, reference: Reference) -> Self {
        Self {
            options,
            chunk: Chunk {
                reference,
                contig: String::new(),
                start: 0,
                bases: Vec::new(),
            },
            contig: None,
            seen: HashSet::new(),
            prune_at: 1 << 16,
            profile: EndProfile::default(),
        }
    }

    fn add(&mut self, record: &bam::Record, header: &noodles::sam::Header) -> Result<()> {
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
        let (Some(Ok(contig_id)), Some(Ok(start))) =
            (record.reference_sequence_id(), record.alignment_start())
        else {
            return Ok(());
        };
        let start = usize::from(start) - 1;
        let ops = record
            .cigar()
            .iter()
            .map(|op| op.map(|op| (op.kind(), op.len())))
            .collect::<std::io::Result<Vec<_>>>()?;
        let reverse = flags.is_reverse_complemented();
        let first = if reverse { ops.last() } else { ops.first() };
        let aligned = matches!(
            first,
            Some(&(Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch, len)) if len >= PROFILE_LENGTH
        );
        if !aligned {
            return Ok(());
        }
        let span: usize = ops
            .iter()
            .filter(|(kind, _)| kind.consumes_reference())
            .map(|(_, len)| len)
            .sum();
        if self.contig.as_ref().map(|(id, _)| *id) != Some(contig_id) {
            let name = header
                .reference_sequences()
                .get_index(contig_id)
                .map(|(name, _)| name.to_string())
                .context("a record's contig is not in the BAM header")?;
            self.contig = Some((contig_id, name));
            self.seen.clear();
        }
        let contig = self
            .contig
            .as_ref()
            .map(|(_, name)| name.as_str())
            .unwrap_or_default();
        if !self
            .chunk
            .cover(contig, start.saturating_sub(1), start + span + 1)?
        {
            return Ok(());
        }
        if self.seen.len() >= self.prune_at {
            self.seen.retain(|&(pos, ..)| pos + 1 >= start);
            self.prune_at = (2 * self.seen.len()).max(1 << 16);
        }
        let sequence = record.sequence();
        let qualities = record.quality_scores();
        let qualities = qualities.as_ref();
        let length = sequence.len();
        for position in 0..PROFILE_LENGTH {
            let (query, pos) = if reverse {
                (length - 1 - position, start + span - 1 - position)
            } else {
                (position, start + position)
            };
            let (Some(read), Some(&quality)) = (sequence.get(query), qualities.get(query)) else {
                continue;
            };
            let read = read.to_ascii_uppercase();
            let Some(reference) = self.chunk.base(pos) else {
                continue;
            };
            if quality < self.options.min_base_quality
                || index(read).is_none()
                || index(reference).is_none()
            {
                continue;
            }
            let prev = pos.checked_sub(1).and_then(|p| self.chunk.base(p));
            let context = Context::of(prev, reference, self.chunk.base(pos + 1));
            if read != reference && !self.seen.insert((pos, read, reverse, position as u8)) {
                continue;
            }
            let (reference, read) = if reverse {
                (complement(reference), complement(read))
            } else {
                (reference, read)
            };
            self.profile.add(position, reference, read, context, 1);
        }
        Ok(())
    }
}

/// Profile the 5' ends of a coordinate-sorted BAM's records, or `None` when
/// `stop` is set before the last is read.
pub fn profile_ends(
    bam: &Path,
    reference: &Path,
    options: &PileupOptions,
    stop: &AtomicBool,
) -> Result<Option<EndProfile>> {
    let file = File::open(bam).with_context(|| format!("failed to open BAM: {bam:?}"))?;
    let decoder = bgzf::io::MultithreadedReader::with_worker_count(DECOMPRESSION_THREADS, file);
    let mut reader = bam::io::Reader::from(decoder);
    let header = reader
        .read_header()
        .context("failed to read the BAM header")?;
    let mut scanner = EndScanner::new(*options, Reference::open(reference)?);
    let mut record = bam::Record::default();
    let mut records = 0u64;
    while reader.read_record(&mut record)? != 0 {
        if stop.load(Ordering::Relaxed) {
            return Ok(None);
        }
        records += 1;
        scanner.add(&record, &header)?;
    }
    info!(
        "profiled {} bases at the 5' ends of {records} records for the tagmentation filter",
        scanner.profile.bases()
    );
    Ok(Some(scanner.profile))
}

#[cfg(test)]
mod tests {
    use streampile::testing::{Frag, SamBuilder, Strand};

    use super::*;
    use crate::testing::write_fasta;

    const A: u8 = b'A';
    const C: u8 = b'C';
    const G: u8 = b'G';
    const T: u8 = b'T';

    /// A profile of `bases` reference bases of each kind at every position,
    /// with `rate` of each read as each other base.
    fn flat(bases: u64, rate: f64) -> EndProfile {
        let mut profile = EndProfile::default();
        let events = (bases as f64 * rate).round() as u64;
        for position in 0..PROFILE_LENGTH {
            for reference in BASES {
                for read in BASES {
                    let count = if read == reference {
                        bases - 3 * events
                    } else {
                        events
                    };
                    profile.add(position, reference, read, Context::NonCpG, count);
                }
            }
        }
        profile
    }

    #[test]
    fn test_a_library_without_an_artifact_learns_nothing() {
        let weights = flat(1_000_000, 0.001).learn();
        assert_eq!(weights.learned(), vec![]);
        assert_eq!(weights.learned_text(|_| true), None);
        for class in EndClass::all() {
            assert_eq!(weights.of(class), None, "{class}");
        }
    }

    #[test]
    fn test_an_empty_profile_learns_nothing() {
        let weights = EndProfile::default().learn();
        assert!(weights.classes.is_empty());
        assert_eq!(weights.of(EndClass::of(T, G, Context::NonCpG)), None);
    }

    #[test]
    fn test_a_spike_at_the_first_base_is_learned_with_its_excess() {
        let mut profile = flat(1_000_000, 0.001);
        profile.add(0, T, G, Context::NonCpG, 4_000);
        let weights = profile.learn();
        let class = EndClass::of(T, G, Context::NonCpG);
        assert_eq!(weights.learned(), vec![(class, 1)]);
        assert_eq!(weights.learned_text(|_| true).as_deref(), Some("T>G@1"));
        let learned = weights.classes[&class];
        assert!((learned.background - 0.001).abs() < 1e-9);
        assert!((learned.excess[0] - 0.004).abs() < 1e-4, "{learned:?}");
        assert!(learned.excess[1..].iter().all(|&x| x == 0.0));
        assert!(weights.of(class).is_some());
        assert!(weights.of(class.reverse_complement()).is_some());
        assert_eq!(weights.of(EndClass::of(T, C, Context::NonCpG)), None);
    }

    #[test]
    fn test_too_few_substitutions_are_not_learned() {
        let mut profile = flat(10_000, 0.001);
        profile.add(0, T, G, Context::NonCpG, 30);
        assert_eq!(profile.learn().learned(), vec![]);
    }

    /// A position with four times the interior's share of G and a
    /// quarter of its share of T finds a real T>G allele sixteen times as
    /// often, so sixteen times the interior rate is no excess.
    #[test]
    fn test_a_real_allele_the_preference_enriches_is_not_an_excess() {
        let mut profile = EndProfile::default();
        let rate = 0.001;
        for position in 0..PROFILE_LENGTH {
            for reference in BASES {
                let bases: u64 = match (position, reference) {
                    (0, b'G') => 2_800_000,
                    (0, b'T') => 175_000,
                    (0, _) => 512_500,
                    _ => 1_000_000,
                };
                for read in BASES.into_iter().filter(|&read| read != reference) {
                    let fold = match (position, reference, read) {
                        (0, b'T', b'G') => 16.0,
                        _ => 1.0,
                    };
                    let events = (bases as f64 * rate * fold).round() as u64;
                    profile.add(position, reference, read, Context::NonCpG, events);
                }
                let matches = bases - 3 * (bases as f64 * rate).round() as u64;
                profile.add(position, reference, reference, Context::NonCpG, matches);
            }
        }
        let weights = profile.learn();
        let class = EndClass::of(T, G, Context::NonCpG);
        let learned = weights.classes[&class];
        assert!((learned.preference[0] - 16.0).abs() < 0.5, "{learned:?}");
        assert_eq!(learned.excess[0], 0.0, "{learned:?}");
        assert_eq!(weights.learned(), vec![]);

        profile.add(0, T, G, Context::NonCpG, 175_000 / 20);
        let weights = profile.learn();
        assert_eq!(weights.learned(), vec![(class, 1)]);
        let fold = weights.of(class).unwrap().fold_left[0];
        assert!((fold - 16.0).abs() < 1.0, "{fold}");
    }

    /// A substitution the preference depletes is compared with the interior
    /// rate, never with less.
    #[test]
    fn test_a_depleted_substitution_is_expected_at_the_interior_rate() {
        let mut profile = flat(1_000_000, 0.001);
        profile.add(0, G, G, Context::NonCpG, 3_000_000);
        let weights = profile.learn();
        let learned = weights.classes[&EndClass::of(G, T, Context::NonCpG)];
        assert!(learned.preference[0] < 0.5, "{learned:?}");
        assert!((learned.expected[0] - learned.background).abs() < 1e-12);
    }

    /// A position rich in CpG has more C>T in all, since a CpG's rate is
    /// ten times any other's, but no more within either context.
    #[test]
    fn test_a_position_rich_in_cpg_is_not_an_excess_of_c_to_t() {
        let mut profile = flat(1_000_000, 0.001);
        let cpg_rate = 0.01;
        for position in 0..PROFILE_LENGTH {
            let cpgs: u64 = if position == 2 { 600_000 } else { 50_000 };
            let events = (cpgs as f64 * cpg_rate) as u64;
            profile.add(position, C, T, Context::CpG, events);
            profile.add(position, C, C, Context::CpG, cpgs - events);
        }
        let weights = profile.learn();
        assert_eq!(weights.learned(), vec![]);
        let split: Vec<EndClass> = weights
            .classes
            .keys()
            .copied()
            .filter(|class| class.reference == C && class.read == T)
            .collect();
        assert_eq!(
            split,
            [
                EndClass::of(C, T, Context::CpG),
                EndClass::of(C, T, Context::NonCpG)
            ]
        );
        let pooled_rate = |position: usize| {
            let events: u64 = Context::ALL
                .map(|c| profile.count(position, C, T, c))
                .iter()
                .sum();
            let bases: u64 = BASES
                .iter()
                .flat_map(|&read| Context::ALL.map(|c| profile.count(position, C, read, c)))
                .sum();
            events as f64 / bases as f64
        };
        assert!(pooled_rate(2) > 2.5 * pooled_rate(15));
    }

    #[test]
    fn test_classes_split_only_the_transitions_at_c_and_g_by_context() {
        assert_eq!(EndClass::all().len(), 14);
        assert_eq!(EndClass::of(C, T, Context::CpG).to_string(), "C>T:CpG");
        assert_eq!(EndClass::of(G, A, Context::NonCpG).stratum(), "C>T:non-CpG");
        assert_eq!(EndClass::of(C, A, Context::CpG).to_string(), "C>A");
        assert_eq!(EndClass::of(A, C, Context::NonCpG).stratum(), "T>G");
        assert_eq!(
            EndClass::of(G, A, Context::CpG).reverse_complement(),
            EndClass::of(C, T, Context::CpG)
        );
        assert_eq!(
            strata(),
            ["C>A", "C>G", "C>T:CpG", "C>T:non-CpG", "T>A", "T>C", "T>G"]
        );
    }

    #[test]
    fn test_a_profile_round_trips_through_its_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("profile.tsv");
        let mut profile = flat(1_000, 0.01);
        profile.add(0, T, G, Context::NonCpG, 77);
        profile.add(3, C, T, Context::CpG, 5);
        profile.write(&path).unwrap();
        assert_eq!(EndProfile::read(&path).unwrap(), profile);
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text
            .starts_with("position\tref_base\tread_base\tcontext\tcount\n1\tA\tA\tnon-CpG\t970\n"));
    }

    #[test]
    fn test_a_malformed_profile_file_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("profile.tsv");
        std::fs::write(&path, "position\tcount\n").unwrap();
        let error = EndProfile::read(&path).unwrap_err().to_string();
        assert!(error.contains("must start with the columns"), "{error}");
        let header = PROFILE_COLUMNS.join("\t");
        for row in [
            "0\tA\tC\tnon-CpG\t1",
            "1\tA\tN\tnon-CpG\t1",
            "1\tA\tC\tCpG\tx",
            "1\tA\tC",
        ] {
            std::fs::write(&path, format!("{header}\n{row}\n")).unwrap();
            let error = EndProfile::read(&path).unwrap_err().to_string();
            assert!(error.contains("line 2"), "{error}");
        }
    }

    fn spiked() -> Weights {
        let mut profile = flat(1_000_000, 0.001);
        profile.add(0, T, G, Context::NonCpG, 4_000);
        profile.learn()
    }

    /// A 150 bp molecule whose leftmost base is `left` bases before the site.
    fn at(base: u8, left: usize) -> Molecule {
        Molecule::new(base, 37, left, 149 - left)
    }

    fn spread(base: u8) -> Vec<Molecule> {
        (0..150).map(|left| at(base, left)).collect()
    }

    fn score(molecules: &[Molecule], ref_base: u8, alt_base: u8) -> Score {
        let class = EndClass::of(ref_base, alt_base, Context::NonCpG);
        let weights = spiked().of(class).unwrap();
        EndEvidence::of(molecules, ref_base, alt_base, class).score(&weights, None)
    }

    #[test]
    fn test_alternates_at_the_learned_position_favor_the_artifact() {
        let mut molecules = spread(T);
        molecules.extend([at(G, 0), at(G, 0)]);
        let score = score(&molecules, T, G);
        assert!(score.log_likelihood_ratio.unwrap() > 6.0, "{score:?}");
        assert_eq!((score.alt_congruent, score.alt_molecules), (2, 2));
        assert_eq!((score.ref_congruent, score.ref_molecules), (1, 150));
    }

    #[test]
    fn test_an_alternate_elsewhere_outweighs_those_at_the_learned_position() {
        let mut molecules = spread(T);
        molecules.extend([at(G, 0), at(G, 60)]);
        let score = score(&molecules, T, G);
        assert!(score.log_likelihood_ratio.unwrap() < -2.0, "{score:?}");
        assert_eq!((score.alt_congruent, score.alt_molecules), (1, 2));
    }

    /// The right end reads the call's reverse complement, so an A>C call's
    /// molecules at their rightmost base are the T>G the library makes.
    #[test]
    fn test_the_right_end_reads_the_reverse_complement() {
        let mut molecules = spread(A);
        molecules.push(at(C, 149));
        let right = score(&molecules, A, C);
        assert!(right.log_likelihood_ratio.unwrap() > 3.0, "{right:?}");
        assert_eq!(right.alt_congruent, 1);

        let mut molecules = spread(A);
        molecules.push(at(C, 0));
        let left = score(&molecules, A, C);
        assert!(left.log_likelihood_ratio.unwrap() < -5.0, "{left:?}");
        assert_eq!(left.alt_congruent, 0);
    }

    #[test]
    fn test_a_duplex_and_a_single_strand_molecule_score_alike() {
        let mut duplex = spread(T);
        duplex.push(at(G, 0));
        let single: Vec<Molecule> = duplex
            .iter()
            .map(|m| m.from_strand(crate::classes::Strand::Forward))
            .collect();
        assert_eq!(score(&duplex, T, G), score(&single, T, G));
    }

    #[test]
    fn test_a_call_without_reference_molecules_is_scored_from_the_pool() {
        let class = EndClass::of(T, G, Context::NonCpG);
        let weights = spiked().of(class).unwrap();
        let pool = EndEvidence::of(&spread(T), T, G, class).reference;
        let lone = EndEvidence::of(&[at(G, 0)], T, G, class);
        let pooled = lone.score(&weights, Some(&pool));
        assert!(pooled.log_likelihood_ratio.unwrap() > 3.0, "{pooled:?}");
        assert_eq!(pooled.ref_molecules, 0);
        let alone = lone.score(&weights, None);
        assert!(alone.log_likelihood_ratio.unwrap() > 3.0, "{alone:?}");
    }

    #[test]
    fn test_a_call_without_alternate_molecules_has_no_ratio() {
        let score = score(&spread(T), T, G);
        assert_eq!(score.log_likelihood_ratio, None);
        assert_eq!(score.ref_molecules, 150);
    }

    /// A real allele the preference puts at the first base far more often
    /// than its reference base is weaker evidence there than one it does not.
    #[test]
    fn test_the_preference_fold_weakens_the_evidence_of_an_enriched_allele() {
        let class = EndClass::of(T, G, Context::NonCpG);
        let plain = spiked().of(class).unwrap();
        let mut enriched = plain;
        enriched.fold_left[0] = 16.0;
        let mut molecules = spread(T);
        molecules.push(at(G, 0));
        let evidence = EndEvidence::of(&molecules, T, G, class);
        let (plain, enriched) = (
            evidence.score(&plain, None).log_likelihood_ratio.unwrap(),
            evidence
                .score(&enriched, None)
                .log_likelihood_ratio
                .unwrap(),
        );
        assert!(plain - enriched > 2.0, "{plain} {enriched}");
    }

    fn reference() -> String {
        "ACGTTGCATGCCATGACTGATCGGATTACA".repeat(10)
    }

    #[test]
    fn test_a_bam_is_profiled_from_each_record_s_five_prime_end() {
        let dir = tempfile::tempdir().unwrap();
        let sequence = reference();
        let fasta = write_fasta(dir.path(), "chr1", &sequence);
        let mut reads = SamBuilder::new().read_length(40);
        let bases = &sequence[..40];
        let mut changed = bases.as_bytes().to_vec();
        changed[0] = b'G';
        let changed = String::from_utf8(changed).unwrap();
        for _ in 0..3 {
            reads.add_frag(Frag::at(1).bases(bases.to_string()));
            reads.add_frag(Frag::at(1).bases(changed.clone()));
        }
        reads.add_frag(Frag::at(1).bases(bases.to_string()).strand(Strand::Minus));
        reads.add_frag(Frag::at(1).bases(bases.to_string()).cigar("5S35M"));
        let bam = dir.path().join("reads.bam");
        reads.write_bam(&bam).unwrap();
        let stop = AtomicBool::new(false);
        let profile = profile_ends(&bam, &fasta, &PileupOptions::default(), &stop)
            .unwrap()
            .unwrap();
        assert_eq!(profile.count(0, A, A, Context::NonCpG), 3);
        assert_eq!(profile.count(0, A, G, Context::NonCpG), 1);
        assert_eq!(profile.count(1, C, C, Context::CpG), 6);
        let last = sequence.as_bytes()[39];
        assert_eq!(
            profile.count(0, complement(last), complement(last), Context::NonCpG),
            1
        );
        assert_eq!(profile.bases(), 7 * PROFILE_LENGTH as u64 - 2);

        stop.store(true, Ordering::Relaxed);
        let stopped = profile_ends(&bam, &fasta, &PileupOptions::default(), &stop).unwrap();
        assert_eq!(stopped, None);
    }
}
