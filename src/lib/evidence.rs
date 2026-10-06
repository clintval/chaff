//! Molecule evidence at a variant site.
//!
//! The statistics see each template once, as a [`Molecule`]: the base it holds
//! at the site, that base's quality, and the site's distances from the
//! template ends its reads reveal. [`Evidence`] separates the statistics from
//! the reads. [`PileupEvidence`] fills it from a streampile pileup builder that
//! leaves out the reads [`PileupOptions`] reject and calls the bases of
//! overlapping mates into one: mates that agree keep the higher quality, and
//! mates that disagree become an `N`, which counts as neither allele.

use std::collections::HashMap;
use std::io;

use anyhow::{Context as _, Result};
use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::record::Flags;
use streampile::{
    AgreementStrategy, AlignmentRecord, DisagreementStrategy, PileupTemplate, RecordSource,
    StreamingPileupBuilder,
};

/// How the quality of agreeing mates is called: the higher of the two.
const AGREEMENT: AgreementStrategy = AgreementStrategy::MaxQual;

/// How the base of disagreeing mates is called: an `N`.
const DISAGREEMENT: DisagreementStrategy = DisagreementStrategy::MaskBoth;

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

/// Which reads and bases count as evidence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PileupOptions {
    /// Reads below this mapping quality are left out.
    pub min_mapping_quality: u8,
    /// Templates whose called base is below this quality are left out, unless
    /// every read's base is at it, as for mates that disagree, called `N`.
    pub min_base_quality: u8,
    /// Keep only paired reads whose mate is also mapped.
    pub paired_reads_only: bool,
}

impl Default for PileupOptions {
    fn default() -> Self {
        Self {
            min_mapping_quality: 20,
            min_base_quality: 20,
            paired_reads_only: false,
        }
    }
}

impl PileupOptions {
    /// A builder that piles up only the reads these options accept: mapped,
    /// primary, not duplicates, at or above the mapping quality floor, and
    /// paired with a mapped mate when only pairs are kept. QC-failed reads are
    /// kept, as fgbio keeps them, and a read without a mapping quality (255)
    /// passes the floor, as in htsjdk.
    pub fn configure<'f, S: RecordSource>(
        &self,
        builder: StreamingPileupBuilder<'f, S>,
    ) -> StreamingPileupBuilder<'f, S> {
        let builder = builder
            .exclude_flags(Flags::SECONDARY | Flags::DUPLICATE | Flags::SUPPLEMENTARY)
            .min_mapping_quality(self.min_mapping_quality)
            .min_base_quality(self.min_base_quality);
        if !self.paired_reads_only {
            return builder;
        }
        builder.read_filter(|record: &S::Record| {
            let flags = record.bam().flags();
            flags.is_segmented() && !flags.is_mate_unmapped()
        })
    }
}

/// [`Evidence`] from a streampile pileup builder over coordinate-sorted reads.
pub struct PileupEvidence<'f, S: RecordSource> {
    builder: StreamingPileupBuilder<'f, S>,
}

impl<'f, S: RecordSource> PileupEvidence<'f, S> {
    /// Evidence from `builder` once `options` configure it.
    pub fn new(builder: StreamingPileupBuilder<'f, S>, options: &PileupOptions) -> Self {
        Self {
            builder: options.configure(builder),
        }
    }
}

impl<S: RecordSource> Evidence for PileupEvidence<'_, S> {
    fn molecules(&mut self, contig: &str, pos: Position) -> Result<Vec<Molecule>> {
        let pileup = self.builder.pileup(contig, usize::from(pos) - 1)?;
        let mut molecules = Vec::new();
        for template in pileup.templates(AGREEMENT, DISAGREEMENT) {
            let molecule = molecule(&template, pileup.min_base_quality())
                .with_context(|| format!("reading template ends at {contig}:{pos}"))?;
            molecules.extend(molecule);
        }
        Ok(molecules)
    }
}

/// The molecule a template shows: its called base and the site's distances
/// from the template's ends, or `None` for a template with no base, a site
/// outside it, or a base under the quality floor. A base is at the floor when
/// its called quality is, or when every read's base here is, so mates that
/// disagree are an `N`, which counts as neither allele.
fn molecule<R: AlignmentRecord>(
    template: &PileupTemplate<'_, R>,
    min_base_quality: u8,
) -> streampile::Result<Option<Molecule>> {
    let (Some(base), Some(quality)) = (template.base(), template.quality()) else {
        return Ok(None);
    };
    if !template.passes(min_base_quality)
        && !template
            .entries()
            .all(|entry| entry.passes(min_base_quality))
    {
        return Ok(None);
    }
    let Some((left, right)) = distances(template)? else {
        return Ok(None);
    };
    Ok(Some(Molecule {
        base,
        quality,
        left,
        right,
    }))
}

/// The site's distances from a template's leftmost and rightmost bases, as
/// `(left, right)`, or `None` for a site outside the template.
///
/// Each distance is counted along the read sequenced from that end where it
/// holds a base here, and otherwise by a read of the other strand, which walks
/// its mate's CIGAR from the `MC` tag. A read whose mate maps to the same contig on the other
/// strand is outside its template wherever a distance is unknown: past its
/// mate's 5' end, or anywhere in a pair whose reads face away from each other.
fn distances<R: AlignmentRecord>(
    template: &PileupTemplate<'_, R>,
) -> streampile::Result<Option<(Option<usize>, Option<usize>)>> {
    let end = |reverse: bool| -> streampile::Result<Option<usize>> {
        let own = template
            .entries()
            .filter(|entry| entry.is_reverse() == reverse)
            .find_map(|entry| entry.five_prime_distance());
        let other = template
            .entries()
            .find(|entry| entry.is_reverse() != reverse);
        match (own, other) {
            (Some(distance), _) => Ok(Some(distance)),
            (None, Some(other)) => other.template_end_distance(),
            (None, None) => Ok(None),
        }
    };
    let (left, right) = (end(false)?, end(true)?);
    let read = template.entries().next().map(|entry| entry.record());
    if (left.is_none() || right.is_none())
        && read.map_or(Ok(false), has_mate_on_the_other_strand)?
    {
        return Ok(None);
    }
    Ok(Some((left, right)))
}

/// Whether a read's mate maps to the same contig on the other strand.
fn has_mate_on_the_other_strand(record: &bam::Record) -> io::Result<bool> {
    let flags = record.flags();
    if !flags.is_segmented()
        || flags.is_mate_unmapped()
        || flags.is_reverse_complemented() == flags.is_mate_reverse_complemented()
    {
        return Ok(false);
    }
    let this = record.reference_sequence_id().transpose()?;
    let mate = record.mate_reference_sequence_id().transpose()?;
    Ok(this == mate)
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
    use noodles::sam::alignment::RecordBuf;
    use rstest::rstest;
    use streampile::testing::{Frag, Pair, SamBuilder, Strand};

    use super::*;

    fn pos(n: usize) -> Position {
        Position::try_from(n).unwrap()
    }

    /// The molecules at a 1-based site of reads piled up under `options`.
    fn molecules_at(reads: &SamBuilder, options: PileupOptions, site: usize) -> Vec<Molecule> {
        let mut evidence = PileupEvidence::new(reads.to_pileup_builder(), &options);
        evidence.molecules("chr1", pos(site)).unwrap()
    }

    /// The site's distances from the template ends of each molecule of these
    /// records at a 1-based site.
    fn distances_at(records: &[RecordBuf], site: usize) -> Vec<(Option<usize>, Option<usize>)> {
        let mut reads = SamBuilder::new();
        reads.extend(records.iter().cloned());
        molecules_at(&reads, PileupOptions::default(), site)
            .iter()
            .map(|m| (m.left, m.right))
            .collect()
    }

    /// The name and number of reads of each template counted as a molecule at
    /// a 1-based site.
    fn templates_at(
        reads: &SamBuilder,
        options: PileupOptions,
        site: usize,
    ) -> Vec<(String, usize)> {
        let mut builder = options.configure(reads.to_pileup_builder());
        let pileup = builder.pileup("chr1", site - 1).unwrap();
        pileup
            .templates(AGREEMENT, DISAGREEMENT)
            .iter()
            .filter(|template| {
                molecule(template, pileup.min_base_quality())
                    .unwrap()
                    .is_some()
            })
            .map(|template| (template.name().to_string(), template.entries().count()))
            .collect()
    }

    fn names(templates: &[(String, usize)]) -> Vec<&str> {
        templates.iter().map(|(name, _)| name.as_str()).collect()
    }

    fn depth(templates: &[(String, usize)]) -> usize {
        templates.iter().map(|(_, reads)| reads).sum()
    }

    /// fgbio `PileupBuilderTest`: "filter out reads below the minimum mapping
    /// quality".
    #[test]
    fn test_filter_out_reads_below_the_minimum_mapping_quality() {
        let mut reads = SamBuilder::new().read_length(50);
        for (name, mapq) in [("q1", 9), ("q2", 10), ("q3", 11)] {
            reads.add_frag(Frag {
                name: Some(name.into()),
                start: 101,
                mapq,
                ..Frag::default()
            });
        }
        let options = PileupOptions {
            min_mapping_quality: 10,
            ..PileupOptions::default()
        };
        assert_eq!(names(&templates_at(&reads, options, 105)), ["q2", "q3"]);
    }

    /// fgbio `PileupBuilderTest`: "filter out base entries below the minimum
    /// base quality".
    #[test]
    fn test_filter_out_base_entries_below_the_minimum_base_quality() {
        let mut reads = SamBuilder::new().read_length(50);
        for (name, quality) in [("q1", 19), ("q2", 20), ("q3", 21)] {
            reads.add_frag(Frag {
                name: Some(name.into()),
                start: 101,
                quals: Some(vec![quality; 50]),
                ..Frag::default()
            });
        }
        let options = PileupOptions::default();
        assert_eq!(names(&templates_at(&reads, options, 105)), ["q2", "q3"]);
    }

    /// fgbio `PileupBuilderTest`: "filter out reads that are not a part of a
    /// mapped pair".
    #[test]
    fn test_filter_out_reads_that_are_not_part_of_a_mapped_pair() {
        let mut reads = SamBuilder::new().read_length(50);
        reads.add_frag(Frag {
            name: Some("q1".into()),
            start: 101,
            ..Frag::default()
        });
        reads.add_pair(Pair {
            name: Some("q2".into()),
            start1: 101,
            start2: 101,
            unmapped2: true,
            ..Pair::default()
        });
        reads.add_pair(Pair {
            name: Some("q3".into()),
            ..Pair::at(101, 300)
        });
        let options = PileupOptions {
            paired_reads_only: true,
            ..PileupOptions::default()
        };
        assert_eq!(names(&templates_at(&reads, options, 105)), ["q3"]);
    }

    /// fgbio `PileupBuilderTest`: "remove one half of each overlapping pair".
    #[test]
    fn test_remove_one_half_of_each_overlapping_pair() {
        let mut reads = SamBuilder::new().read_length(50);
        for (name, start1, start2) in [("q1", 100, 110), ("q2", 110, 100), ("q3", 50, 100)] {
            reads.add_pair(Pair {
                name: Some(name.into()),
                ..Pair::at(start1, start2)
            });
        }
        let templates = templates_at(&reads, PileupOptions::default(), 125);
        assert_eq!(depth(&templates), 5);
        assert_eq!(templates.len(), 3);
    }

    /// fgbio `PileupBuilderTest`: "not filter out single-end records when we
    /// are not filter for mapped pairs only and we are not removing records
    /// where the position is outside the insert of FR pairs".
    #[test]
    fn test_keep_single_end_records_at_both_read_ends() {
        let mut reads = SamBuilder::new().read_length(50);
        reads.add_frag(Frag {
            name: Some("q1".into()),
            start: 100,
            ..Frag::default()
        });
        for site in [100, 149] {
            assert_eq!(
                templates_at(&reads, PileupOptions::default(), site).len(),
                1
            );
        }
    }

    /// fgbio `PileupBuilderTest`: "filter out records where a position is
    /// outside the insert for an FR pair".
    #[test]
    fn test_filter_out_positions_outside_the_insert_of_an_fr_pair() {
        let mut reads = SamBuilder::new().read_length(50);
        reads.add_pair(Pair {
            name: Some("q2".into()),
            ..Pair::at(101, 100)
        });
        let depths: Vec<usize> = [100, 101, 149, 150]
            .map(|site| depth(&templates_at(&reads, PileupOptions::default(), site)))
            .into();
        assert_eq!(depths, [0, 2, 2, 0]);
    }

    /// fgbio `PileupBuilderTest`: "not filter out records where a position is
    /// outside what might look like an 'insert' for a non-FR pair".
    #[test]
    fn test_keep_positions_outside_what_looks_like_an_insert_for_a_non_fr_pair() {
        let mut reads = SamBuilder::new().read_length(50);
        reads.add_pair(Pair {
            name: Some("q2".into()),
            strand1: Strand::Minus,
            strand2: Strand::Plus,
            ..Pair::at(101, 100)
        });
        let depths: Vec<usize> = [100, 101, 149, 150]
            .map(|site| depth(&templates_at(&reads, PileupOptions::default(), site)))
            .into();
        assert_eq!(depths, [1, 2, 2, 1]);
    }

    /// fgbio `PileupTest`: "BaseEntry should report the correct
    /// offsets/positions/bases", and the matching `PileupBuilderTest` case, on
    /// the entries chaff reads. fgbio's 1-based positions in read order are one
    /// more than these 0-based distances from the 5' end.
    #[test]
    fn test_entries_report_offsets_positions_and_bases() {
        let mut reads = SamBuilder::new().read_length(50).base_quality(35);
        reads.add_pair(Pair {
            name: Some("q1".into()),
            start1: 101,
            start2: 201,
            bases1: Some("A".repeat(50)),
            bases2: Some("C".repeat(50)),
            ..Pair::default()
        });
        let mut builder = PileupOptions::default().configure(reads.to_pileup_builder());
        let mut seen = Vec::new();
        for site in [105, 205] {
            let pileup = builder.pileup("chr1", site - 1).unwrap();
            let entry = pileup.get(0).unwrap();
            seen.push((
                entry.base(),
                entry.sequenced_base(),
                entry.quality(),
                entry.query_position(),
                entry.five_prime_distance(),
            ));
        }
        assert_eq!(
            seen,
            [
                (Some(b'A'), Some(b'A'), Some(35), Some(4), Some(4)),
                (Some(b'C'), Some(b'G'), Some(35), Some(4), Some(45)),
            ]
        );
    }

    /// fgbio `BamsTest`: "Bams.insertCoordinates should fail on fragments and
    /// inappropriate pairs". chaff gives such records only their own end where
    /// fgbio raises.
    #[test]
    fn test_template_distances_know_only_their_own_end_on_fragments_and_inappropriate_pairs() {
        let mut reads = SamBuilder::new().read_length(10).base_quality(20);
        let own_only = vec![(Some(4), None)];
        assert_eq!(distances_at(&reads.add_frag(Frag::at(100)), 104), own_only);
        let unmapped_mate = reads.add_pair(Pair {
            start1: 100,
            start2: 100,
            unmapped2: true,
            ..Pair::default()
        });
        assert_eq!(distances_at(&unmapped_mate[..1], 104), own_only);
        let pair = reads.add_pair(Pair::at(100, 200));
        let other_contig = SamBuilder::with_mate_reference_sequence_id(pair[0].clone(), 1);
        assert_eq!(distances_at(&[other_contig], 104), own_only);
    }

    /// fgbio `BamsTest`: "Bams.insertCoordinates should calculate insert
    /// coordinates correctly", as the distances from the template's first and
    /// last bases. chaff reads the far end from `MC`, not `TLEN`.
    #[test]
    fn test_template_distances_of_an_fr_pair() {
        let mut reads = SamBuilder::new().read_length(10).base_quality(20);
        let pair = reads.add_pair(Pair::at(100, 191));
        assert_eq!(distances_at(&pair, 100), [(Some(0), Some(100))]);
        assert_eq!(distances_at(&pair, 200), [(Some(100), Some(0))]);
    }

    /// fgbio `BamsTest`: "Bams.positionFromOtherEndOfTemplate should return
    /// None for anything that's not an FR mapped pair".
    #[test]
    fn test_distance_from_the_other_end_is_none_unless_fr_pair() {
        let mut reads = SamBuilder::new();
        let header = reads.header().clone();
        let check = |records: Vec<RecordBuf>| {
            for r in records.iter().filter(|r| !r.flags().is_unmapped()) {
                for site in [r.alignment_start(), r.alignment_end()] {
                    let position = usize::from(site.unwrap()) - 1;
                    let distance = streampile::template_end_distance(r, &header, position);
                    assert_eq!(distance.unwrap(), None, "{r:?}");
                }
            }
        };
        check(reads.add_frag(Frag::at(100)));
        check(reads.add_pair(Pair {
            start1: 100,
            start2: 200,
            unmapped2: true,
            ..Pair::default()
        }));
        for (strand1, strand2) in [
            (Strand::Plus, Strand::Plus),
            (Strand::Minus, Strand::Minus),
            (Strand::Minus, Strand::Plus),
        ] {
            check(reads.add_pair(Pair {
                start1: 100,
                start2: 200,
                strand1,
                strand2,
                ..Pair::default()
            }));
        }
    }

    /// fgbio `BamsTest`: "Bams.positionFromOtherEndOfTemplate should correctly
    /// calculate the position from the other end of the template for FR pairs",
    /// as the 1-based distance from the far template end.
    #[test]
    fn test_distance_from_the_other_end_of_the_template() {
        let mut reads = SamBuilder::new().read_length(50);
        let pair = reads.add_pair(Pair::at(101, 151));
        let from_other_end = |i: usize, site: usize| {
            streampile::template_end_distance(&pair[i], reads.header(), site - 1)
                .unwrap()
                .map(|d| d + 1)
        };
        assert_eq!(from_other_end(0, 101), Some(100));
        assert_eq!(from_other_end(0, 111), Some(90));
        assert_eq!(from_other_end(0, 151), Some(50));
        assert_eq!(from_other_end(0, 200), Some(1));
        assert_eq!(from_other_end(1, 200), Some(100));
        assert_eq!(from_other_end(1, 190), Some(90));
        assert_eq!(from_other_end(1, 150), Some(50));
        assert_eq!(from_other_end(1, 101), Some(1));
    }

    #[test]
    fn test_template_distances_count_soft_clips_and_not_hard_clips() {
        let mut reads = SamBuilder::new().read_length(50);
        let soft = reads.add_pair(Pair {
            cigar1: Some("5S45M".into()),
            cigar2: Some("40M10S".into()),
            ..Pair::at(101, 151)
        });
        assert_eq!(distances_at(&soft, 101), [(Some(5), Some(99))]);
        assert_eq!(distances_at(&soft, 190), [(Some(94), Some(10))]);
        let hard = reads.add_pair(Pair {
            cigar1: Some("5H45M".into()),
            cigar2: Some("40M10H".into()),
            ..Pair::at(101, 151)
        });
        assert_eq!(distances_at(&hard, 101), [(Some(0), Some(89))]);
        assert_eq!(distances_at(&hard, 190), [(Some(89), Some(0))]);
    }

    /// A deletion between a site and a template end takes its length off the
    /// distance, and an insertion adds its length, whether the read sequenced
    /// from that end holds the site or only its mate does.
    #[test]
    fn test_template_distances_count_indels_by_their_length() {
        let mut reads = SamBuilder::new().read_length(50);
        let pair = reads.add_pair(Pair {
            cigar1: Some("30M4D20M".into()),
            cigar2: Some("20M2I28M".into()),
            ..Pair::at(101, 121)
        });
        let (deleted, inserted) = (4, 2);
        let expected = [(Some(140 - 101 - deleted), Some(168 - 140 + inserted))];
        for records in [&pair[..], &pair[..1], &pair[1..]] {
            assert_eq!(distances_at(records, 140), expected);
        }
    }

    /// A site past its mate's 5' end lies outside the template, and so does
    /// every site of a pair whose reads face away from each other.
    #[test]
    fn test_template_distances_leave_out_sites_past_the_mates_five_prime_end() {
        let mut reads = SamBuilder::new().read_length(50);
        let through = reads.add_pair(Pair::at(101, 100));
        assert_eq!(distances_at(&through, 150), []);
        assert_eq!(distances_at(&through, 100), []);
        assert_eq!(distances_at(&through, 149), [(Some(48), Some(0))]);
        let away = reads.add_pair(Pair {
            strand1: Strand::Minus,
            strand2: Strand::Plus,
            ..Pair::at(100, 200)
        });
        assert_eq!(distances_at(&away, 120), []);
        assert_eq!(distances_at(&away, 220), []);
    }

    #[test]
    fn test_template_distances_of_a_pair_without_a_mate_cigar_are_an_error() {
        let mut reads = SamBuilder::new().read_length(50);
        let pair = reads.add_pair(Pair {
            name: Some("q1".into()),
            ..Pair::at(101, 151)
        });
        let mut stripped = SamBuilder::new();
        stripped.extend(pair.into_iter().map(SamBuilder::without_mate_cigar));
        for site in [111, 161] {
            let mut evidence =
                PileupEvidence::new(stripped.to_pileup_builder(), &PileupOptions::default());
            let error = evidence.molecules("chr1", pos(site)).unwrap_err();
            let message = format!("{error:#}");
            assert!(message.contains("read q1"), "{message}");
            assert!(message.contains("MC"), "{message}");
        }
    }

    #[test]
    fn test_template_distances_ignore_the_insert_size() {
        let mut reads = SamBuilder::new().read_length(50);
        let pair: Vec<RecordBuf> = reads
            .add_pair(Pair::at(101, 151))
            .into_iter()
            .map(|mut r| {
                *r.template_length_mut() = 7;
                r
            })
            .collect();
        assert_eq!(distances_at(&pair, 120), [(Some(19), Some(80))]);
        assert_eq!(distances_at(&pair, 160), [(Some(59), Some(40))]);
    }

    #[rstest]
    #[case(Strand::Plus, "5S45M", (Some(5), None))]
    #[case(Strand::Minus, "45M5S", (None, Some(49)))]
    fn test_template_distances_of_a_fragment_know_only_their_own_end(
        #[case] strand: Strand,
        #[case] cigar: &str,
        #[case] expected: (Option<usize>, Option<usize>),
    ) {
        let mut reads = SamBuilder::new().read_length(50);
        let frag = reads.add_frag(Frag {
            strand,
            cigar: Some(cigar.into()),
            ..Frag::at(100)
        });
        assert_eq!(distances_at(&frag, 100), [expected]);
    }

    #[test]
    fn test_template_distances_of_a_tandem_pair_know_only_their_own_end() {
        let mut reads = SamBuilder::new().read_length(50);
        let pair = reads.add_pair(Pair {
            strand1: Strand::Plus,
            strand2: Strand::Plus,
            ..Pair::at(100, 200)
        });
        assert_eq!(distances_at(&pair, 210), [(Some(10), None)]);
    }

    /// Mates that agree are one molecule at the higher quality, and mates that
    /// disagree are an `N`, which counts as neither allele. A mate under the
    /// quality floor still takes part, so its disagreement leaves no molecule.
    #[test]
    fn test_overlapping_mates_are_called_into_one_molecule() {
        let mut reads = SamBuilder::new().read_length(50);
        for (bases2, quality2) in [('A', 40), ('C', 30), ('C', 10)] {
            reads.add_pair(Pair {
                bases1: Some("A".repeat(50)),
                bases2: Some(bases2.to_string().repeat(50)),
                quals1: Some(vec![30; 50]),
                quals2: Some(vec![quality2; 50]),
                ..Pair::at(101, 121)
            });
        }
        let molecules = molecules_at(&reads, PileupOptions::default(), 130);
        let calls: Vec<(u8, u8)> = molecules.iter().map(|m| (m.base, m.quality)).collect();
        assert_eq!(calls, [(b'A', 40), (b'N', 2)]);
        assert_eq!(molecules[0], Molecule::new(b'A', 40, 29, 40));
    }

    /// Secondary, duplicate, supplementary, and unmapped reads are left out,
    /// QC-failed reads are kept, as fgbio keeps them.
    #[test]
    fn test_reads_are_left_out_by_their_flags() {
        for (flag, kept) in [
            (0x100, false),
            (0x400, false),
            (0x800, false),
            (0x4, false),
            (0x200, true),
        ] {
            let mut built = SamBuilder::new().read_length(50);
            let frag = built.add_frag(Frag::at(101));
            let mut reads = SamBuilder::new();
            reads.extend([SamBuilder::with_flags(frag[0].clone(), flag)]);
            let templates = templates_at(&reads, PileupOptions::default(), 110);
            assert_eq!(templates.len(), usize::from(kept), "flag {flag:#x}");
        }
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
