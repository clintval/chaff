//! Molecule evidence at a variant site.
//!
//! The statistics see each template once, as a [`Molecule`]: the base it holds
//! at the site, that base's quality, and the template ends the reads reveal.
//! [`Evidence`] is the seam between the statistics and the reads.
//! [`PileupEvidence`] fills it from any [`PileupSource`], a streaming pileup
//! engine that lists the reads covering a position with the offset of their
//! aligned base there. The engine owns record streaming and CIGAR walking; this
//! module owns the read floors, the template geometry, and the collapse of
//! overlapping mates into one molecule.

use std::collections::HashMap;

use anyhow::{Context as _, Result};
use noodles::core::Position;
use noodles::sam;
use noodles::sam::alignment::Record;

use crate::template::{template_ends, ReadBase, ReadFilter, TemplateEnds};

/// One template's observation at a site.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Molecule {
    /// The upper-cased base on the forward strand, `N` when the mates disagree.
    pub base: u8,
    /// The base quality.
    pub quality: u8,
    /// The template's leftmost base, the forward strand's 5' end, when known.
    pub start: Option<i64>,
    /// The template's rightmost base, the reverse strand's 5' end, when known.
    pub end: Option<i64>,
}

impl Molecule {
    /// A molecule with both ends known.
    pub fn new(base: u8, quality: u8, start: i64, end: i64) -> Self {
        Self {
            base,
            quality,
            start: Some(start),
            end: Some(end),
        }
    }

    /// The 0-based distance of `pos` from the template's leftmost base.
    pub fn from_start(&self, pos: i64) -> Option<i64> {
        self.start.map(|start| pos - start)
    }

    /// The 0-based distance of `pos` from the template's rightmost base.
    pub fn from_end(&self, pos: i64) -> Option<i64> {
        self.end.map(|end| end - pos)
    }

    /// The template length, when both ends are known.
    pub fn length(&self) -> Option<i64> {
        Some(self.end? - self.start? + 1)
    }
}

/// The molecules covering each variant site, asked for in coordinate order.
pub trait Evidence {
    /// One molecule per template covering the 1-based `pos` on `contig`.
    fn molecules(&mut self, contig: &str, pos: Position) -> Result<Vec<Molecule>>;
}

/// A streaming pileup engine over coordinate-sorted records.
pub trait PileupSource {
    /// The alignment record type the engine yields.
    type Record: Record;

    /// The header of the records.
    fn header(&self) -> &sam::Header;

    /// Every read with an aligned base at the 1-based `pos` on `contig`, with
    /// the offset of that base; reads with a deletion or a skip there are left
    /// out. Positions are asked for in coordinate order.
    fn pileup(&mut self, contig: &str, pos: Position) -> Result<Vec<ReadBase<'_, Self::Record>>>;
}

impl<S: streampile::RecordSource> PileupSource for streampile::StreamingPileupBuilder<'_, S> {
    type Record = noodles::bam::Record;

    fn header(&self) -> &sam::Header {
        streampile::StreamingPileupBuilder::header(self)
    }

    fn pileup(&mut self, contig: &str, pos: Position) -> Result<Vec<ReadBase<'_, Self::Record>>> {
        let pileup =
            streampile::StreamingPileupBuilder::pileup(self, contig, usize::from(pos) - 1)?;
        Ok(pileup
            .iter()
            .filter_map(|entry| {
                entry
                    .query_position()
                    .map(|offset| ReadBase::new(entry.record(), offset))
            })
            .collect())
    }
}

/// [`Evidence`] from a [`PileupSource`], after the read floors.
#[derive(Debug)]
pub struct PileupEvidence<P> {
    source: P,
    header: sam::Header,
    filter: ReadFilter,
}

impl<P: PileupSource> PileupEvidence<P> {
    /// Evidence from `source` that keeps the reads and bases `filter` accepts.
    pub fn new(source: P, filter: ReadFilter) -> Self {
        let header = source.header().clone();
        Self {
            source,
            header,
            filter,
        }
    }

    /// The read-level observations at a site, before overlapping mates are
    /// collapsed: one per read passing the floors, with its name.
    pub fn observations(
        &mut self,
        contig: &str,
        pos: Position,
    ) -> Result<Vec<(Option<Vec<u8>>, Molecule)>> {
        let site = usize::from(pos) as i64;
        let entries = self.source.pileup(contig, pos)?;
        let mut observations = Vec::with_capacity(entries.len());
        for entry in entries {
            let record = entry.record;
            if !self.filter.accepts(record)? {
                continue;
            }
            let Some(base) = entry.base() else { continue };
            let quality = entry.quality()?.unwrap_or(u8::MAX);
            if quality < self.filter.min_base_quality {
                continue;
            }
            let ends: TemplateEnds = template_ends(record, &self.header)
                .with_context(|| format!("reading template ends at {contig}:{site}"))?;
            if ends.excludes(site) {
                continue;
            }
            let name = record.name().map(|n| n.to_vec());
            observations.push((
                name,
                Molecule {
                    base,
                    quality,
                    start: ends.start,
                    end: ends.end,
                },
            ));
        }
        Ok(observations)
    }
}

impl<P: PileupSource> Evidence for PileupEvidence<P> {
    fn molecules(&mut self, contig: &str, pos: Position) -> Result<Vec<Molecule>> {
        Ok(collapse_templates(self.observations(contig, pos)?))
    }
}

/// Collapse reads that share a name into one molecule per template.
///
/// Mates that agree on the base become one molecule with the higher of the
/// two qualities; mates that disagree become one molecule holding `N`, which
/// counts as neither allele. Ends one mate cannot see are taken from the other.
/// Unnamed reads are never collapsed. First-seen order is kept.
pub fn collapse_templates(observations: Vec<(Option<Vec<u8>>, Molecule)>) -> Vec<Molecule> {
    let mut molecules: Vec<Molecule> = Vec::with_capacity(observations.len());
    let mut index: HashMap<Vec<u8>, usize> = HashMap::with_capacity(observations.len());
    for (name, molecule) in observations {
        let Some(name) = name else {
            molecules.push(molecule);
            continue;
        };
        match index.get(&name) {
            None => {
                index.insert(name, molecules.len());
                molecules.push(molecule);
            }
            Some(&i) => {
                let kept = &mut molecules[i];
                if kept.base == molecule.base {
                    kept.quality = kept.quality.max(molecule.quality);
                } else {
                    kept.base = b'N';
                    kept.quality = kept.quality.min(molecule.quality);
                }
                kept.start = kept.start.or(molecule.start);
                kept.end = kept.end.or(molecule.end);
            }
        }
    }
    molecules
}

/// [`Evidence`] from a fixed table of molecules per site, for callers that
/// already hold molecule observations.
#[derive(Clone, Debug, Default)]
pub struct MoleculeTable {
    sites: HashMap<(String, usize), Vec<Molecule>>,
}

impl MoleculeTable {
    /// An empty table.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the molecules at a 1-based site.
    pub fn insert(&mut self, contig: &str, pos: usize, molecules: Vec<Molecule>) {
        self.sites.insert((contig.to_string(), pos), molecules);
    }
}

impl Evidence for MoleculeTable {
    fn molecules(&mut self, contig: &str, pos: Position) -> Result<Vec<Molecule>> {
        Ok(self
            .sites
            .get(&(contig.to_string(), usize::from(pos)))
            .cloned()
            .unwrap_or_default())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::{Frag, Pair, SamBuilder};

    fn pos(n: usize) -> Position {
        Position::try_from(n).unwrap()
    }

    fn names(observations: &[(Option<Vec<u8>>, Molecule)]) -> Vec<String> {
        observations
            .iter()
            .map(|(n, _)| String::from_utf8(n.clone().unwrap()).unwrap())
            .collect()
    }

    /// fgbio `PileupBuilderTest`: "filter out reads below the minimum mapping
    /// quality".
    #[test]
    fn test_filter_out_reads_below_the_minimum_mapping_quality() {
        let mut builder = SamBuilder::new().read_length(50);
        for (name, mapq) in [("q1", 9), ("q2", 10), ("q3", 11)] {
            builder.add_frag(Frag {
                name: Some(name.into()),
                start: 101,
                mapq,
                ..Frag::default()
            });
        }
        let filter = ReadFilter {
            min_mapping_quality: 10,
            ..ReadFilter::default()
        };
        let mut evidence = PileupEvidence::new(builder.pileup(), filter);
        let pile = evidence.observations("chr1", pos(105)).unwrap();
        assert_eq!(pile.len(), 2);
        assert!(!names(&pile).contains(&"q1".to_string()));
    }

    /// fgbio `PileupBuilderTest`: "filter out base entries below the minimum
    /// base quality".
    #[test]
    fn test_filter_out_base_entries_below_the_minimum_base_quality() {
        let mut builder = SamBuilder::new().read_length(50).base_quality(19);
        builder.add_frag(Frag {
            name: Some("q1".into()),
            start: 101,
            ..Frag::default()
        });
        for (name, quality) in [("q2", 20), ("q3", 21)] {
            let recs = SamBuilder::new()
                .read_length(50)
                .base_quality(quality)
                .add_frag(Frag {
                    name: Some(name.into()),
                    start: 101,
                    ..Frag::default()
                });
            builder.extend(recs);
        }
        let mut evidence = PileupEvidence::new(builder.pileup(), ReadFilter::default());
        let pile = evidence.observations("chr1", pos(105)).unwrap();
        assert_eq!(pile.len(), 2);
        assert!(!names(&pile).contains(&"q1".to_string()));
    }

    /// fgbio `PileupBuilderTest`: "filter out reads that are not a part of a
    /// mapped pair".
    #[test]
    fn test_filter_out_reads_that_are_not_part_of_a_mapped_pair() {
        let mut builder = SamBuilder::new().read_length(50);
        builder.add_frag(Frag {
            name: Some("q1".into()),
            start: 101,
            ..Frag::default()
        });
        builder.add_pair(Pair {
            name: Some("q2".into()),
            start1: 101,
            start2: 101,
            unmapped2: true,
            ..Pair::default()
        });
        builder.add_pair(Pair {
            name: Some("q3".into()),
            ..Pair::at(101, 300)
        });
        let filter = ReadFilter {
            paired_reads_only: true,
            ..ReadFilter::default()
        };
        let mut evidence = PileupEvidence::new(builder.pileup(), filter);
        let pile = evidence.observations("chr1", pos(105)).unwrap();
        assert_eq!(names(&pile), vec!["q3".to_string()]);
    }

    /// fgbio `PileupBuilderTest`: "remove one half of each overlapping pair".
    #[test]
    fn test_remove_one_half_of_each_overlapping_pair() {
        let mut builder = SamBuilder::new().read_length(50);
        for (name, start1, start2) in [("q1", 100, 110), ("q2", 110, 100), ("q3", 50, 100)] {
            builder.add_pair(Pair {
                name: Some(name.into()),
                ..Pair::at(start1, start2)
            });
        }
        let mut evidence = PileupEvidence::new(builder.pileup(), ReadFilter::default());
        let pile = evidence.observations("chr1", pos(125)).unwrap();
        assert_eq!(pile.len(), 5);
        assert_eq!(collapse_templates(pile).len(), 3);
    }

    /// fgbio `PileupBuilderTest`: "not filter out single-end records when we
    /// are not filter for mapped pairs only and we are not removing records
    /// where the position is outside the insert of FR pairs".
    #[test]
    fn test_keep_single_end_records_at_both_read_ends() {
        let mut builder = SamBuilder::new().read_length(50);
        builder.add_frag(Frag {
            name: Some("q1".into()),
            start: 100,
            ..Frag::default()
        });
        let mut evidence = PileupEvidence::new(builder.pileup(), ReadFilter::default());
        assert_eq!(evidence.observations("chr1", pos(100)).unwrap().len(), 1);
        assert_eq!(evidence.observations("chr1", pos(149)).unwrap().len(), 1);
    }

    /// fgbio `PileupBuilderTest`: "filter out records where a position is
    /// outside the insert for an FR pair".
    #[test]
    fn test_filter_out_positions_outside_the_insert_of_an_fr_pair() {
        let mut builder = SamBuilder::new().read_length(50);
        builder.add_pair(Pair {
            name: Some("q2".into()),
            ..Pair::at(101, 100)
        });
        let mut evidence = PileupEvidence::new(builder.pileup(), ReadFilter::default());
        let mut depth = |p| evidence.observations("chr1", pos(p)).unwrap().len();
        assert_eq!(depth(100), 0);
        assert_eq!(depth(101), 2);
        assert_eq!(depth(149), 2);
        assert_eq!(depth(150), 0);
    }

    /// fgbio `PileupBuilderTest`: "not filter out records where a position is
    /// outside what might look like an 'insert' for a non-FR pair".
    #[test]
    fn test_keep_positions_outside_what_looks_like_an_insert_for_a_non_fr_pair() {
        let mut builder = SamBuilder::new().read_length(50);
        builder.add_pair(Pair {
            name: Some("q2".into()),
            strand1: crate::testing::Strand::Minus,
            strand2: crate::testing::Strand::Plus,
            ..Pair::at(101, 100)
        });
        let mut evidence = PileupEvidence::new(builder.pileup(), ReadFilter::default());
        let mut depth = |p| evidence.observations("chr1", pos(p)).unwrap().len();
        assert_eq!(depth(100), 1);
        assert_eq!(depth(101), 2);
        assert_eq!(depth(149), 2);
        assert_eq!(depth(150), 1);
    }

    #[test]
    fn test_collapse_templates_merges_agreeing_mates_and_masks_disagreeing_ones() {
        let a = |q, start, end| Molecule {
            base: b'A',
            quality: q,
            start,
            end,
        };
        let observations = vec![
            (Some(b"x".to_vec()), a(30, Some(10), None)),
            (Some(b"y".to_vec()), a(30, Some(10), Some(60))),
            (Some(b"x".to_vec()), a(40, None, Some(50))),
            (
                Some(b"y".to_vec()),
                Molecule {
                    base: b'C',
                    ..a(20, Some(10), Some(60))
                },
            ),
            (None, a(10, None, None)),
            (None, a(10, None, None)),
        ];
        let molecules = collapse_templates(observations);
        assert_eq!(molecules.len(), 4);
        assert_eq!(molecules[0], a(40, Some(10), Some(50)));
        assert_eq!(molecules[1].base, b'N');
        assert_eq!(molecules[1].quality, 20);
    }

    #[test]
    fn test_molecule_distances() {
        let m = Molecule::new(b'A', 30, 101, 150);
        assert_eq!(m.from_start(101), Some(0));
        assert_eq!(m.from_end(101), Some(49));
        assert_eq!(m.length(), Some(50));
        let half = Molecule { end: None, ..m };
        assert_eq!(half.from_end(101), None);
        assert_eq!(half.length(), None);
    }

    #[test]
    fn test_molecule_table_returns_inserted_molecules() {
        let mut table = MoleculeTable::new();
        table.insert("chr1", 10, vec![Molecule::new(b'A', 30, 1, 20)]);
        assert_eq!(table.molecules("chr1", pos(10)).unwrap().len(), 1);
        assert!(table.molecules("chr1", pos(11)).unwrap().is_empty());
    }
}
