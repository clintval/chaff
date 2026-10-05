//! Per-sample metrics: one row per filter and stratum, written as a headered
//! TSV and logged.
//!
//! For each stratum the counts pool every annotated call's molecules. A
//! molecule is congruent when it sits where the artifact would put it: nearer
//! the lesion strand's 5' end than its 3' end for the lesion copy filter
//! (NanoSeq's strand assignment by nearest 5' end), within the distance of the
//! relevant template end for the read-end filters.
//!
//! The asymmetry test is a one-sided binomial test of the congruent alternate
//! molecules against the congruent fraction of the reference molecules:
//! `P(X >= alt_congruent)` for `X ~ Binomial(alt_molecules,
//! ref_congruent_fraction)`. When fragment geometry is symmetric the
//! reference fraction is one half, NanoSeq's null; using the reference
//! molecules keeps capture and length skew out of the test.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use anyhow::{Context, Result};
use csv::{Terminator, WriterBuilder};
use serde::Serialize;
use statrs::distribution::{Binomial, DiscreteCDF};

/// One metrics row.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct StratumMetrics {
    /// The sample whose molecules were measured.
    pub sample: String,
    /// The filter: `lesion-copy`, `end-repair-fill-in`, or `a-tailing`.
    pub filter: String,
    /// The stratum the artifact fraction is learned in.
    pub stratum: String,
    /// Calls annotated by the filter.
    pub calls: u64,
    /// Calls the filter's FILTER was applied to.
    pub filtered: u64,
    /// The learned artifact fraction of the stratum.
    pub artifact_fraction: f64,
    /// The sum over calls of the posterior probability of an artifact.
    pub expected_artifacts: f64,
    /// Alternate molecules measured.
    pub alt_molecules: u64,
    /// Alternate molecules congruent with the artifact.
    pub alt_congruent: u64,
    /// `alt_congruent / alt_molecules`.
    pub alt_congruent_fraction: Option<f64>,
    /// Reference molecules measured.
    pub ref_molecules: u64,
    /// Reference molecules congruent with the artifact.
    pub ref_congruent: u64,
    /// `ref_congruent / ref_molecules`.
    pub ref_congruent_fraction: Option<f64>,
    /// The one-sided binomial asymmetry p-value.
    pub asymmetry_p_value: Option<f64>,
}

/// `P(X >= k)` for `X ~ Binomial(n, p)`.
pub fn binomial_upper_tail(k: u64, n: u64, p: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if k > n {
        return 0.0;
    }
    if p <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return 1.0;
    }
    let binomial = Binomial::new(p, n).expect("a probability in (0, 1)");
    binomial.sf(k - 1)
}

/// A ratio, or `None` over zero.
pub fn fraction(numerator: u64, denominator: u64) -> Option<f64> {
    (denominator > 0).then(|| numerator as f64 / denominator as f64)
}

impl StratumMetrics {
    /// Fill the derived fractions and the asymmetry test from the counts.
    pub fn finish(mut self) -> Self {
        self.alt_congruent_fraction = fraction(self.alt_congruent, self.alt_molecules);
        self.ref_congruent_fraction = fraction(self.ref_congruent, self.ref_molecules);
        self.asymmetry_p_value = match (self.alt_molecules, self.ref_congruent_fraction) {
            (n, Some(p)) if n > 0 => Some(binomial_upper_tail(self.alt_congruent, n, p)),
            _ => None,
        };
        self
    }
}

/// Write the rows as a tab-separated file with a header.
pub fn write_metrics(path: &Path, rows: &[StratumMetrics]) -> Result<()> {
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create metrics directory: {parent:?}"))?;
    }
    let file = File::create(path).with_context(|| format!("failed to create metrics: {path:?}"))?;
    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .terminator(Terminator::Any(b'\n'))
        .from_writer(BufWriter::new(file));
    for row in rows {
        writer.serialize(row)?;
    }
    if rows.is_empty() {
        writer.write_record([
            "sample",
            "filter",
            "stratum",
            "calls",
            "filtered",
            "artifact_fraction",
            "expected_artifacts",
            "alt_molecules",
            "alt_congruent",
            "alt_congruent_fraction",
            "ref_molecules",
            "ref_congruent",
            "ref_congruent_fraction",
            "asymmetry_p_value",
        ])?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row() -> StratumMetrics {
        StratumMetrics {
            sample: "s1".into(),
            filter: "lesion-copy".into(),
            stratum: "C>T:CpG".into(),
            calls: 3,
            filtered: 1,
            artifact_fraction: 0.4,
            expected_artifacts: 1.2,
            alt_molecules: 10,
            alt_congruent: 9,
            alt_congruent_fraction: None,
            ref_molecules: 100,
            ref_congruent: 50,
            ref_congruent_fraction: None,
            asymmetry_p_value: None,
        }
    }

    #[test]
    fn test_binomial_upper_tail() {
        assert_eq!(binomial_upper_tail(0, 10, 0.5), 1.0);
        assert!((binomial_upper_tail(10, 10, 0.5) - 0.5f64.powi(10)).abs() < 1e-12);
        assert!((binomial_upper_tail(9, 10, 0.5) - 11.0 / 1024.0).abs() < 1e-12);
        assert_eq!(binomial_upper_tail(1, 10, 0.0), 0.0);
        assert_eq!(binomial_upper_tail(3, 10, 1.0), 1.0);
        assert_eq!(binomial_upper_tail(11, 10, 0.5), 0.0);
    }

    #[test]
    fn test_finish_fills_fractions_and_the_asymmetry_test() {
        let row = row().finish();
        assert_eq!(row.alt_congruent_fraction, Some(0.9));
        assert_eq!(row.ref_congruent_fraction, Some(0.5));
        assert!((row.asymmetry_p_value.unwrap() - 11.0 / 1024.0).abs() < 1e-12);
        let empty = StratumMetrics {
            alt_molecules: 0,
            alt_congruent: 0,
            ..row
        }
        .finish();
        assert_eq!(empty.asymmetry_p_value, None);
        assert_eq!(empty.alt_congruent_fraction, None);
    }

    #[test]
    fn test_write_metrics_has_a_header_and_rows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("nested/metrics.tsv");
        write_metrics(&path, &[row().finish()]).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("sample\tfilter\tstratum\tcalls"));
        assert!(lines[1].starts_with("s1\tlesion-copy\tC>T:CpG\t3\t1\t0.4"));
    }

    #[test]
    fn test_write_metrics_without_rows_still_writes_the_header() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("metrics.tsv");
        write_metrics(&path, &[]).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.starts_with("sample\tfilter"));
    }
}
