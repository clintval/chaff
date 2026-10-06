//! Molecule evidence at a variant site.
//!
//! The statistics see each template once, as a [`Molecule`]: the base it holds
//! at the site, that base's quality, and the site's distances from the
//! template ends the reads reveal. [`Evidence`] separates the statistics from
//! the reads. [`PileupEvidence`] fills it from any [`PileupSource`], a
//! streaming pileup engine that lists the reads covering a position with the
//! offset of their aligned base there. The engine owns record streaming and
//! CIGAR walking, streampile counts the distances for any engine's records, and
//! this module owns the read floors and the collapse of overlapping mates into
//! one molecule.

use std::collections::HashMap;

use anyhow::{Context as _, Result};
use noodles::core::Position;
use noodles::sam;
use noodles::sam::alignment::Record;

use crate::template::{ReadBase, ReadFilter};

/// One template's observation at a site.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Molecule {
    /// The upper-cased base on the forward strand, `N` when the mates disagree.
    pub base: u8,
    /// The base quality.
    pub quality: u8,
    /// The template's bases between the site and its leftmost base, the
    /// forward strand's 5' end, when known: 0 at that base.
    pub left: Option<usize>,
    /// The template's bases between the site and its rightmost base, the
    /// reverse strand's 5' end, when known: 0 at that base.
    pub right: Option<usize>,
}

impl Molecule {
    /// A molecule with both distances known.
    pub fn new(base: u8, quality: u8, left: usize, right: usize) -> Self {
        Self {
            base,
            quality,
            left: Some(left),
            right: Some(right),
        }
    }

    /// The template's length in bases, when both distances are known.
    pub fn length(&self) -> Option<usize> {
        Some(self.left? + self.right? + 1)
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
            let Some((left, right)) = entry
                .template_distances(&self.header, pos)
                .with_context(|| format!("reading template ends at {contig}:{pos}"))?
            else {
                continue;
            };
            let name = record.name().map(|n| n.to_vec());
            observations.push((
                name,
                Molecule {
                    base,
                    quality,
                    left,
                    right,
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
/// counts as neither allele. Where the mates of an FR pair overlap, each
/// distance is counted along the read sequenced from that end, so both mates
/// report the same distances. The molecule keeps the first read's distances
/// and takes any it lacks from its mate. Unnamed reads are never collapsed.
/// First-seen order is kept.
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
                kept.left = kept.left.or(molecule.left);
                kept.right = kept.right.or(molecule.right);
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
        let a = |q, left, right| Molecule {
            base: b'A',
            quality: q,
            left,
            right,
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
    fn test_overlapping_mates_report_the_same_template_bases_across_indels() {
        let mut builder = SamBuilder::new().read_length(50);
        builder.add_pair(Pair {
            name: Some("q1".into()),
            cigar1: Some("30M4D20M".into()),
            cigar2: Some("20M2I28M".into()),
            ..Pair::at(101, 121)
        });
        let mut evidence = PileupEvidence::new(builder.pileup(), ReadFilter::default());
        let pile = evidence.observations("chr1", pos(140)).unwrap();
        let distances: Vec<_> = pile.iter().map(|(_, m)| (m.left, m.right)).collect();
        assert_eq!(distances, vec![(Some(35), Some(30)); 2]);
        let molecules = collapse_templates(pile);
        assert_eq!(molecules.len(), 1);
        assert_eq!(
            (molecules[0].left, molecules[0].right),
            (Some(35), Some(30))
        );
    }

    #[test]
    fn test_the_streampile_engine_measures_template_bases_as_the_test_pileup_does() {
        let mut builder = SamBuilder::new().read_length(50).coordinate_sorted();
        for (start1, start2, cigar1, cigar2) in [
            (101, 121, "30M4D20M", "20M2I28M"),
            (96, 131, "5S45M", "40M10S"),
            (111, 141, "5H45M", "45M5H"),
            (121, 120, "50M", "50M"),
        ] {
            let bases = |cigar: &str| "A".repeat(if cigar.contains('H') { 45 } else { 50 });
            builder.add_pair(Pair {
                bases1: Some(bases(cigar1)),
                bases2: Some(bases(cigar2)),
                cigar1: Some(cigar1.into()),
                cigar2: Some(cigar2.into()),
                ..Pair::at(start1, start2)
            });
        }
        let dir = tempfile::tempdir().unwrap();
        let path = builder.write_bam(&dir.path().join("reads.bam"));
        let mut reader = noodles::bam::io::reader::Builder
            .build_from_path(path)
            .unwrap();
        let header = reader.read_header().unwrap();
        let engine = streampile::StreamingPileupBuilder::new(reader, &header).unwrap();
        let mut streamed = PileupEvidence::new(engine, ReadFilter::default());
        let mut expected = PileupEvidence::new(builder.pileup(), ReadFilter::default());
        let mut measured = 0;
        for site in 100..=190 {
            let observations = streamed.observations("chr1", pos(site)).unwrap();
            assert_eq!(
                observations,
                expected.observations("chr1", pos(site)).unwrap()
            );
            measured += observations.len();
        }
        assert!(measured > 300, "{measured}");
    }

    #[test]
    fn test_molecule_length_counts_the_site_and_both_distances() {
        let m = Molecule::new(b'A', 30, 0, 49);
        assert_eq!(m.length(), Some(50));
        let half = Molecule { right: None, ..m };
        assert_eq!(half.length(), None);
    }

    #[test]
    fn test_molecule_table_returns_inserted_molecules() {
        let mut table = MoleculeTable::new();
        table.insert("chr1", 10, vec![Molecule::new(b'A', 30, 9, 10)]);
        assert_eq!(table.molecules("chr1", pos(10)).unwrap().len(), 1);
        assert!(table.molecules("chr1", pos(11)).unwrap().is_empty());
    }
}
