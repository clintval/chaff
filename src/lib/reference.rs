//! Reference sequence context from an indexed FASTA.

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use anyhow::{bail, Context as _, Result};
use noodles::core::{Position, Region};
use noodles::fasta;

/// An indexed FASTA (`.fai` alongside it).
pub struct Reference {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    lengths: HashMap<Vec<u8>, usize>,
}

impl Reference {
    /// Open an indexed FASTA.
    pub fn open(path: &Path) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("failed to open indexed FASTA (needs a .fai): {path:?}"))?;
        let lengths = reader
            .index()
            .as_ref()
            .iter()
            .map(|r| (r.name().to_vec(), r.length() as usize))
            .collect();
        Ok(Self { reader, lengths })
    }

    /// The forward-strand bases before, at, and after the 1-based `pos`,
    /// upper-cased; a neighbor past a contig end is `None`.
    pub fn context(&mut self, contig: &str, pos: usize) -> Result<(Option<u8>, u8, Option<u8>)> {
        let length = *self
            .lengths
            .get(contig.as_bytes())
            .with_context(|| format!("contig is not in the reference: {contig}"))?;
        if pos == 0 || pos > length {
            bail!("position is outside the reference contig: {contig}:{pos}");
        }
        let start = pos.saturating_sub(1).max(1);
        let end = (pos + 1).min(length);
        let region = Region::new(
            contig,
            Position::try_from(start)?..=Position::try_from(end)?,
        );
        let record = self
            .reader
            .query(&region)
            .with_context(|| format!("failed to read the reference at {contig}:{pos}"))?;
        let bases: Vec<u8> = record
            .sequence()
            .as_ref()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        let at = |p: usize| bases.get(p - start).copied();
        let base = at(pos).with_context(|| format!("no reference base at {contig}:{pos}"))?;
        let prev = (pos > 1).then(|| at(pos - 1)).flatten();
        let next = (pos < length).then(|| at(pos + 1)).flatten();
        Ok((prev, base, next))
    }

    /// The length of a contig, or `None` for one the FASTA does not hold.
    pub fn length(&self, contig: &str) -> Option<usize> {
        self.lengths.get(contig.as_bytes()).copied()
    }

    /// The upper-cased forward-strand bases of the 0-based, half-open span
    /// `start..end` of a contig, clipped to its end.
    pub fn bases(&mut self, contig: &str, start: usize, end: usize) -> Result<Vec<u8>> {
        let length = self
            .length(contig)
            .with_context(|| format!("contig is not in the reference: {contig}"))?;
        let end = end.min(length);
        if start >= end {
            return Ok(Vec::new());
        }
        let region = Region::new(
            contig,
            Position::try_from(start + 1)?..=Position::try_from(end)?,
        );
        let record = self.reader.query(&region).with_context(|| {
            format!(
                "failed to read the reference at {contig}:{}-{end}",
                start + 1
            )
        })?;
        Ok(record
            .sequence()
            .as_ref()
            .iter()
            .map(u8::to_ascii_uppercase)
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::write_fasta;

    #[test]
    fn test_context_reads_neighbors_and_stops_at_contig_ends() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), "chr1", &"acgT".repeat(40));
        let mut reference = Reference::open(&path).unwrap();
        assert_eq!(
            reference.context("chr1", 1).unwrap(),
            (None, b'A', Some(b'C'))
        );
        assert_eq!(
            reference.context("chr1", 62).unwrap(),
            (Some(b'A'), b'C', Some(b'G'))
        );
        assert_eq!(
            reference.context("chr1", 160).unwrap(),
            (Some(b'G'), b'T', None)
        );
        assert!(reference.context("chr1", 161).is_err());
        assert!(reference.context("chr2", 1).is_err());
    }

    #[test]
    fn test_bases_reads_a_span_clipped_to_the_contig() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), "chr1", &"acgT".repeat(40));
        let mut reference = Reference::open(&path).unwrap();
        assert_eq!(reference.length("chr1"), Some(160));
        assert_eq!(reference.length("chr2"), None);
        assert_eq!(reference.bases("chr1", 0, 6).unwrap(), b"ACGTAC");
        assert_eq!(reference.bases("chr1", 158, 200).unwrap(), b"GT");
        assert!(reference.bases("chr1", 160, 170).unwrap().is_empty());
        assert!(reference.bases("chr2", 0, 1).is_err());
    }
}
