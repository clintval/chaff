//! Per-sample metrics: one row per filter and stratum, written as a headered
//! TSV and logged.
//!
//! For each stratum the counts pool every annotated call's molecules. A
//! molecule is congruent when the site lies within the filter's distance of
//! the end its artifact favors, the same end and distance its posterior
//! weighs.
//!
//! The asymmetry test asks whether more alternate molecules are congruent
//! than each call's own reference molecules predict. Call `i` contributes its
//! `a_i` alternate molecules, each congruent under the null with probability
//! `f_i = (ref_congruent_i + 1) / (ref_molecules_i + 2)`, its reference
//! molecules' congruent fraction smoothed by one molecule each way, so a call
//! without reference molecules has NanoSeq's null of one half. The p-value is
//! `P(X >= alt_congruent)` for `X` the sum of the calls' `Binomial(a_i, f_i)`,
//! computed exactly. Because each call is measured against its own reference
//! molecules, skew that a site's alleles share, such as capture or fragment
//! length, stays out of the test, and calls whose alternate molecules follow
//! their own reference molecules cannot add up to an asymmetric stratum.

use std::io::BufWriter;
use std::path::Path;

use anyhow::Result;
use csv::{Terminator, WriterBuilder};
use serde::{Serialize, Serializer};

use crate::io::StagedFile;

/// One metrics row.
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct StratumMetrics {
    /// The sample whose molecules were measured.
    pub sample: String,
    /// The filter: `copied-damage`, `end-repair-fill-in`, or `a-tailing`.
    pub filter: String,
    /// The stratum the artifact fraction is learned in.
    pub stratum: String,
    /// Calls annotated by the filter.
    pub calls: u64,
    /// Calls the filter's FILTER was applied to.
    pub filtered: u64,
    /// The learned artifact fraction of the stratum, empty under fgbio's
    /// per-call prior.
    #[serde(serialize_with = "six_digits_or_empty")]
    pub artifact_fraction: Option<f64>,
    /// The learned artifact fraction of all the filter's calls, which each
    /// stratum's prior shrinks toward; empty under fgbio's per-call prior.
    #[serde(serialize_with = "six_digits_or_empty")]
    pub filter_artifact_fraction: Option<f64>,
    /// The filter's distance in bases, shared by its strata: a decay's scale,
    /// learned or fixed, or a window.
    #[serde(serialize_with = "six_digits")]
    pub distance: f64,
    /// The sum over calls of the posterior probability of an artifact.
    #[serde(serialize_with = "six_digits")]
    pub expected_artifacts: f64,
    /// The sum over calls of the posterior probability of a real mutation,
    /// a call without a posterior counting as one: the stratum's expected
    /// count of real mutations, a burden that weighs each call by how likely
    /// it is real rather than filtering it. With `expected_artifacts` it sums
    /// to `calls`.
    #[serde(serialize_with = "six_digits")]
    pub expected_mutations: f64,
    /// Alternate molecules measured.
    pub alt_molecules: u64,
    /// Alternate molecules congruent with the artifact.
    pub alt_congruent: u64,
    /// Alternate molecules congruent under the asymmetry test's null, the sum
    /// over calls of `a_i f_i`.
    #[serde(serialize_with = "six_digits_or_empty")]
    pub expected_alt_congruent: Option<f64>,
    /// `alt_congruent / alt_molecules`.
    #[serde(serialize_with = "six_digits_or_empty")]
    pub alt_congruent_fraction: Option<f64>,
    /// Reference molecules measured.
    pub ref_molecules: u64,
    /// Reference molecules congruent with the artifact.
    pub ref_congruent: u64,
    /// `ref_congruent / ref_molecules`.
    #[serde(serialize_with = "six_digits_or_empty")]
    pub ref_congruent_fraction: Option<f64>,
    /// The one-sided asymmetry p-value, `P(X >= alt_congruent)`.
    #[serde(serialize_with = "six_digits_or_empty")]
    pub asymmetry_p_value: Option<f64>,
    /// Copied damage only, from the single-strand consensus of a duplex BAM:
    /// duplex changes of the stratum's class per molecule over the library,
    /// leaving out positions with at least 3 changes in 1% of their
    /// molecules as germline or clonal, and positions of 10 molecules or
    /// fewer, where a germline variant would pass for one change.
    #[serde(serialize_with = "six_digits_or_empty")]
    pub change_rate: Option<f64>,
    /// Single-strand changes per molecule with both strands called, each by
    /// at least 2 raw reads, over the same positions.
    #[serde(serialize_with = "six_digits_or_empty")]
    pub single_strand_rate: Option<f64>,
    /// `change_rate / single_strand_rate`.
    #[serde(serialize_with = "six_digits_or_empty")]
    pub conversion_ratio: Option<f64>,
    /// The mean artifact prior of the stratum's calls under the library's
    /// chance model.
    #[serde(serialize_with = "six_digits_or_empty")]
    pub chance_fraction: Option<f64>,
    /// The positions with two or more changes that the library's chance
    /// model expects beyond those observed and three standard deviations of
    /// a Poisson count, at each depth, as a share of those observed: zero
    /// when the model fits.
    #[serde(serialize_with = "six_digits_or_empty")]
    pub chance_excess: Option<f64>,
}

/// Serialize `value` rounded to six significant digits.
fn six_digits<S: Serializer>(value: &f64, serializer: S) -> Result<S::Ok, S::Error> {
    serializer.serialize_f64(format!("{value:.5e}").parse().unwrap_or(*value))
}

/// Serialize `value` rounded to six significant digits, or an empty field for `None`.
fn six_digits_or_empty<S: Serializer>(
    value: &Option<f64>,
    serializer: S,
) -> Result<S::Ok, S::Error> {
    match value {
        Some(value) => six_digits(value, serializer),
        None => serializer.serialize_none(),
    }
}

/// The probability under the asymmetry test's null that one of a call's
/// alternate molecules is congruent: its reference molecules' congruent
/// fraction, smoothed by one molecule each way.
pub fn null_fraction(ref_congruent: u32, ref_molecules: u32) -> f64 {
    (f64::from(ref_congruent) + 1.0) / (f64::from(ref_molecules) + 2.0)
}

/// `P(X >= k)` for `X` the sum of independent `Binomial(n_i, p_i)`, one per
/// `(n_i, p_i)`, by convolving their distributions one trial at a time.
pub fn poisson_binomial_upper_tail(k: u64, trials: &[(u32, f64)]) -> f64 {
    let n: usize = trials.iter().map(|&(n, _)| n as usize).sum();
    if k == 0 {
        return 1.0;
    }
    if k as usize > n {
        return 0.0;
    }
    let mut pmf = vec![0.0; n + 1];
    pmf[0] = 1.0;
    let mut most = 0;
    for &(count, p) in trials {
        for _ in 0..count {
            most += 1;
            for j in (1..=most).rev() {
                pmf[j] = pmf[j] * (1.0 - p) + pmf[j - 1] * p;
            }
            pmf[0] *= 1.0 - p;
        }
    }
    pmf[k as usize..].iter().sum::<f64>().min(1.0)
}

/// A ratio, or `None` over zero.
pub fn fraction(numerator: u64, denominator: u64) -> Option<f64> {
    (denominator > 0).then(|| numerator as f64 / denominator as f64)
}

impl StratumMetrics {
    /// Fill the derived fractions and the asymmetry test from the counts and
    /// each call's alternate molecules and null fraction, as `(a_i, f_i)`.
    pub fn finish(mut self, trials: &[(u32, f64)]) -> Self {
        self.alt_congruent_fraction = fraction(self.alt_congruent, self.alt_molecules);
        self.ref_congruent_fraction = fraction(self.ref_congruent, self.ref_molecules);
        let measured = self.alt_molecules > 0;
        self.expected_alt_congruent =
            measured.then(|| trials.iter().map(|&(n, p)| f64::from(n) * p).sum());
        self.asymmetry_p_value =
            measured.then(|| poisson_binomial_upper_tail(self.alt_congruent, trials));
        self
    }
}

/// The column names of a metrics row.
fn columns() -> Result<csv::StringRecord> {
    let mut row = csv::Writer::from_writer(Vec::new());
    row.serialize(StratumMetrics::default())?;
    let text = row.into_inner()?;
    Ok(csv::Reader::from_reader(text.as_slice()).headers()?.clone())
}

/// Write the rows as a tab-separated file with a header.
pub fn write_metrics(path: &Path, rows: &[StratumMetrics]) -> Result<()> {
    let staged = StagedFile::create(path)?;
    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .terminator(Terminator::Any(b'\n'))
        .from_writer(BufWriter::new(staged.writer()?));
    if rows.is_empty() {
        writer.write_record(&columns()?)?;
    }
    for row in rows {
        writer.serialize(row)?;
    }
    writer.flush()?;
    staged.persist()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row() -> StratumMetrics {
        StratumMetrics {
            sample: "s1".into(),
            filter: "copied-damage".into(),
            stratum: "C>T:CpG".into(),
            calls: 3,
            filtered: 1,
            artifact_fraction: Some(0.4),
            distance: 30.0,
            expected_artifacts: 1.2,
            alt_molecules: 10,
            alt_congruent: 9,
            ref_molecules: 100,
            ref_congruent: 50,
            ..StratumMetrics::default()
        }
    }

    /// One call with all ten alternate molecules at a null of one half.
    const TRIALS: [(u32, f64); 1] = [(10, 0.5)];

    #[test]
    fn test_poisson_binomial_upper_tail() {
        let tail = poisson_binomial_upper_tail;
        assert_eq!(tail(0, &TRIALS), 1.0);
        assert!((tail(10, &TRIALS) - 0.5f64.powi(10)).abs() < 1e-12);
        assert!((tail(9, &TRIALS) - 11.0 / 1024.0).abs() < 1e-12);
        assert_eq!(tail(1, &[(10, 0.0)]), 0.0);
        assert_eq!(tail(3, &[(10, 1.0)]), 1.0);
        assert_eq!(tail(11, &TRIALS), 0.0);
        let mixed = [(1, 0.1), (1, 0.9)];
        assert!((tail(1, &mixed) - 0.91).abs() < 1e-12);
        assert!((tail(2, &mixed) - 0.09).abs() < 1e-12);
    }

    #[test]
    fn test_null_fraction_is_smoothed_by_one_molecule_each_way() {
        assert_eq!(null_fraction(0, 0), 0.5);
        assert_eq!(null_fraction(0, 50), 1.0 / 52.0);
        assert_eq!(null_fraction(9, 10), 10.0 / 12.0);
    }

    #[test]
    fn test_finish_fills_fractions_and_the_asymmetry_test() {
        let row = row().finish(&TRIALS);
        assert_eq!(row.alt_congruent_fraction, Some(0.9));
        assert_eq!(row.ref_congruent_fraction, Some(0.5));
        assert_eq!(row.expected_alt_congruent, Some(5.0));
        assert!((row.asymmetry_p_value.unwrap() - 11.0 / 1024.0).abs() < 1e-12);
        let empty = StratumMetrics {
            alt_molecules: 0,
            alt_congruent: 0,
            ..row
        }
        .finish(&[]);
        assert_eq!(empty.asymmetry_p_value, None);
        assert_eq!(empty.alt_congruent_fraction, None);
    }

    #[test]
    fn test_write_metrics_has_a_header_and_rows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("nested/metrics.tsv");
        write_metrics(&path, &[row().finish(&TRIALS)]).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("sample\tfilter\tstratum\tcalls"));
        assert!(lines[1].starts_with("s1\tcopied-damage\tC>T:CpG\t3\t1\t0.4"));
    }

    #[test]
    fn test_write_metrics_rounds_values_to_six_significant_digits() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("metrics.tsv");
        let learned = StratumMetrics {
            artifact_fraction: Some(1.0 / 3.0),
            expected_artifacts: 2.0 / 3.0,
            expected_mutations: 7.0 / 3.0,
            ..row()
        };
        let fgbio = StratumMetrics {
            artifact_fraction: None,
            ..row()
        };
        write_metrics(&path, &[learned.finish(&TRIALS), fgbio.finish(&TRIALS)]).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let rows: Vec<Vec<&str>> = text
            .lines()
            .skip(1)
            .map(|l| l.split('\t').collect())
            .collect();
        assert_eq!(
            rows[0][5..10],
            ["0.333333", "", "30.0", "0.666667", "2.33333"]
        );
        assert_eq!(rows[0][12], "5.0");
        assert_eq!(rows[0][17], "0.0107422");
        assert_eq!(rows[1][5], "");
    }

    #[test]
    fn test_write_metrics_without_rows_still_writes_the_header() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("metrics.tsv");
        write_metrics(&path, &[]).unwrap();
        let empty = std::fs::read_to_string(&path).unwrap();
        write_metrics(&path, &[row().finish(&TRIALS)]).unwrap();
        let full = std::fs::read_to_string(&path).unwrap();
        assert_eq!(
            empty.lines().collect::<Vec<_>>(),
            [full.lines().next().unwrap()]
        );
        assert!(empty.starts_with("sample\tfilter"));
    }
}
