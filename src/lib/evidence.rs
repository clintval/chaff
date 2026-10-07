//! Molecule evidence at a variant site.
//!
//! The statistics see each template once, as a [`Molecule`]: the base it holds
//! at the site, that base's quality, and the site's distances from the
//! template ends its mates reveal, whether they are reads or consensus.
//! [`Evidence`] separates the statistics from the BAM. [`PileupEvidence`]
//! fills it from a streampile pileup builder that leaves out the records
//! [`PileupOptions`] reject and calls the bases of
//! overlapping mates into one: mates that agree keep the higher quality, and
//! mates that disagree become an `N`, which counts as neither allele. A
//! template that holds a deletion at the site is a [`DELETION`], which counts
//! as neither allele but in the depth, as fgbio counts it. Each molecule also
//! names the strand its mates were copied from, by pair orientation as GATK's
//! `LearnReadOrientationModel` reads it: the first of pair is a copy of the
//! original strand from its 5' end, so an F1R2 template comes from the forward
//! strand and an F2R1 template from the reverse. A duplex consensus, whose
//! records carry fgbio's `aD` and `bD` depths of both strands, comes from
//! both.

use std::collections::HashMap;

use anyhow::{anyhow, Context as _, Result};
use noodles::core::Position;
use noodles::sam::alignment::record::Flags;
use streampile::{
    AgreementStrategy, AlignmentRecord, DisagreementStrategy, EntryKind, PileupEntry,
    PileupTemplate, RecordSource, StreamingPileupBuilder,
};

use crate::classes::Strand;

/// How the quality of agreeing mates is called: the higher of the two.
const AGREEMENT: AgreementStrategy = AgreementStrategy::MaxQual;

/// How the base of disagreeing mates is called: an `N`.
const DISAGREEMENT: DisagreementStrategy = DisagreementStrategy::MaskBoth;

/// The base of a template that holds a deletion at the site, as a VCF writes a
/// spanning deletion.
pub const DELETION: u8 = b'*';

/// One template's observation at a site.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Molecule {
    /// The upper-cased base on the forward strand, `N` when the mates
    /// disagree, or [`DELETION`].
    pub base: u8,
    /// The base quality, 0 for a deletion.
    pub quality: u8,
    /// The template's bases between the site and its leftmost base, the
    /// forward strand's 5' end, when known: 0 at that base.
    pub left: Option<usize>,
    /// The template's bases between the site and its rightmost base, the
    /// reverse strand's 5' end, when known: 0 at that base.
    pub right: Option<usize>,
    /// Whether the mate fgbio keeps for the template, its first mate here at
    /// the quality floor, is reverse: a site equally far from both ends is
    /// nearer that mate's own 5' end.
    pub reverse: bool,
    /// The strand the template's mates were copied from, or `None` for a
    /// duplex consensus, which holds both.
    pub origin: Option<Strand>,
}

impl Molecule {
    /// A duplex molecule with both distances known.
    pub fn new(base: u8, quality: u8, left: usize, right: usize) -> Self {
        Self {
            base,
            quality,
            left: Some(left),
            right: Some(right),
            reverse: false,
            origin: None,
        }
    }

    /// The same molecule copied from one strand.
    pub fn from_strand(self, origin: Strand) -> Self {
        Self {
            origin: Some(origin),
            ..self
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

/// Which records and bases count as evidence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PileupOptions {
    /// Records below this mapping quality are left out.
    pub min_mapping_quality: u8,
    /// Records whose base is below this quality take no part in their
    /// template's base.
    pub min_base_quality: u8,
    /// Keep only paired records whose mate is also mapped.
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
    /// A builder that piles up only the records these options accept: mapped,
    /// primary, not duplicates, at or above the mapping quality floor, and
    /// paired with a mapped mate when only pairs are kept. QC-failed records
    /// are kept, as fgbio keeps them, and a record without a mapping quality
    /// (255) passes the floor, as in htsjdk.
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

/// [`Evidence`] from a streampile pileup builder over a coordinate-sorted BAM.
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
        let pileup = self
            .builder
            .pileup(contig, usize::from(pos) - 1)
            .map_err(|error| match error {
                streampile::Error::Backwards {
                    contig,
                    position,
                    from_contig,
                    from_position,
                } => anyhow!(
                    "the VCF/BCF reaches {contig}:{} after {from_contig}:{}, which the BAM sorts after it, so both must order their contigs alike",
                    position + 1,
                    from_position + 1
                ),
                streampile::Error::UnknownContig(contig) => {
                    anyhow!("the BAM header has no contig {contig}")
                }
                error => error.into(),
            })?;
        let mut molecules = Vec::new();
        for template in pileup.templates(AGREEMENT, DISAGREEMENT) {
            let molecule = molecule(&template, pileup.min_base_quality())
                .map_err(|error| match error {
                    streampile::Error::MissingMateCigar { .. } => {
                        anyhow!("{error}; add mate CIGARs with samtools fixmate")
                    }
                    error => error.into(),
                })
                .with_context(|| format!("reading template ends at {contig}:{pos}"))?;
            molecules.extend(molecule);
        }
        Ok(molecules)
    }
}

/// The molecule a template shows: the base its mates at the quality floor
/// call, or else a deletion any of its mates holds, and the site's distances
/// from the template's ends, or `None` for a template with neither or a site
/// outside it.
fn molecule<R: AlignmentRecord>(
    template: &PileupTemplate<'_, R>,
    min_base_quality: u8,
) -> streampile::Result<Option<Molecule>> {
    let deleted = |e: &PileupEntry<'_, R>| e.is_deletion() || e.is_skip();
    let (base, quality) = match (template.base(), template.quality()) {
        (Some(base), Some(quality)) => (base, quality),
        _ if template.entries().any(|e| deleted(&e)) => (DELETION, 0),
        _ => return Ok(None),
    };
    let Some((left, right)) = distances(template)? else {
        return Ok(None);
    };
    let kept = template
        .entries()
        .find(|e| deleted(e) || (e.kind() == EntryKind::Base && e.passes(min_base_quality)));
    Ok(Some(Molecule {
        base,
        quality,
        left,
        right,
        reverse: kept.is_some_and(|e| e.is_reverse()),
        origin: template.entries().next().and_then(origin),
    }))
}

/// The strand a mate's template was copied from: the first of pair's strand,
/// which for a second of pair is its mate's, or `None` for a duplex consensus,
/// which holds both.
fn origin<R: AlignmentRecord>(entry: PileupEntry<'_, R>) -> Option<Strand> {
    let depth = |tag: &[u8; 2]| {
        entry
            .record()
            .data()
            .get(tag)
            .and_then(Result::ok)
            .and_then(|value| value.as_int())
    };
    if depth(b"aD").is_some_and(|d| d > 0) && depth(b"bD").is_some_and(|d| d > 0) {
        return None;
    }
    let flags = entry.flags();
    let read_one_reverse = if flags.is_segmented() && flags.is_last_segment() {
        flags.is_mate_reverse_complemented()
    } else {
        flags.is_reverse_complemented()
    };
    Some(if read_one_reverse {
        Strand::Reverse
    } else {
        Strand::Forward
    })
}

/// The site's distances from a template's leftmost and rightmost bases, as
/// `(left, right)`, or `None` for a site outside the template.
///
/// Each distance is counted along the mate sequenced from that end where it
/// holds a base here, and otherwise by a mate of the other strand, which walks
/// its mate's CIGAR from the `MC` tag. A mate of an FR pair, as htsjdk 5.0.0
/// classifies it, is outside its template wherever a distance is unknown: past
/// its mate's 5' end. A mate of a pair whose mates face away from each other
/// knows only its own end, as fgbio has it.
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
    let read = template.entries().next();
    if (left.is_none() || right.is_none()) && read.map_or(Ok(false), |r| r.is_fr_pair())? {
        return Ok(None);
    }
    Ok(Some((left, right)))
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

    /// Read 1 is a copy of the strand its template came from: forward for an
    /// F1R2 pair, reverse for F2R1, even where only read 2 covers the site; a
    /// duplex consensus carrying both strands' depths comes from neither one.
    #[test]
    fn test_a_template_s_origin_is_read_one_s_strand_unless_it_is_a_duplex() {
        use crate::classes::Strand as Origin;
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as Field;
        let origin = |pair: Pair, site: usize| {
            let mut reads = SamBuilder::new().read_length(50);
            reads.add_pair(pair);
            molecules_at(&reads, PileupOptions::default(), site)[0].origin
        };
        let f1r2 = || Pair::at(101, 131);
        let f2r1 = || {
            Pair::at(131, 101)
                .strand1(Strand::Minus)
                .strand2(Strand::Plus)
        };
        assert_eq!(origin(f1r2(), 110), Some(Origin::Forward));
        assert_eq!(origin(f1r2(), 175), Some(Origin::Forward));
        assert_eq!(origin(f2r1(), 110), Some(Origin::Reverse));
        assert_eq!(origin(f2r1(), 175), Some(Origin::Reverse));
        let depth = |pair: Pair, a: i32, b: i32| {
            pair.attr(Tag::new(b'a', b'D'), Field::from(a))
                .attr(Tag::new(b'b', b'D'), Field::from(b))
        };
        assert_eq!(origin(depth(f1r2(), 3, 2), 110), None);
        assert_eq!(origin(depth(f2r1(), 3, 0), 110), Some(Origin::Reverse));
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
        let floor = pileup.min_base_quality();
        pileup
            .templates(AGREEMENT, DISAGREEMENT)
            .iter()
            .filter(|template| molecule(template, floor).unwrap().is_some())
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
            reads.add_frag(Frag::at(101).name(name).mapq(mapq));
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
            reads.add_frag(Frag::at(101).name(name).quals(vec![quality; 50]));
        }
        let options = PileupOptions::default();
        assert_eq!(names(&templates_at(&reads, options, 105)), ["q2", "q3"]);
    }

    /// fgbio `PileupBuilderTest`: "filter out reads that are not a part of a
    /// mapped pair".
    #[test]
    fn test_filter_out_reads_that_are_not_part_of_a_mapped_pair() {
        let mut reads = SamBuilder::new().read_length(50);
        reads.add_frag(Frag::at(101).name("q1"));
        reads.add_pair(Pair::at(101, 101).name("q2").unmapped2(true));
        reads.add_pair(Pair::at(101, 300).name("q3"));
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
            reads.add_pair(Pair::at(start1, start2).name(name));
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
        reads.add_frag(Frag::at(100).name("q1"));
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
        reads.add_pair(Pair::at(101, 100).name("q2"));
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
        reads.add_pair(
            Pair::at(101, 100)
                .name("q2")
                .strand1(Strand::Minus)
                .strand2(Strand::Plus),
        );
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
        reads.add_pair(
            Pair::at(101, 201)
                .name("q1")
                .bases1("A".repeat(50))
                .bases2("C".repeat(50)),
        );
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
        let unmapped_mate = reads.add_pair(Pair::at(100, 100).unmapped2(true));
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
        check(reads.add_pair(Pair::at(100, 200).unmapped2(true)));
        for (strand1, strand2) in [
            (Strand::Plus, Strand::Plus),
            (Strand::Minus, Strand::Minus),
            (Strand::Minus, Strand::Plus),
        ] {
            check(reads.add_pair(Pair::at(100, 200).strand1(strand1).strand2(strand2)));
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
        let soft = reads.add_pair(Pair::at(101, 151).cigar1("5S45M").cigar2("40M10S"));
        assert_eq!(distances_at(&soft, 101), [(Some(5), Some(99))]);
        assert_eq!(distances_at(&soft, 190), [(Some(94), Some(10))]);
        let hard = reads.add_pair(Pair::at(101, 151).cigar1("5H45M").cigar2("40M10H"));
        assert_eq!(distances_at(&hard, 101), [(Some(0), Some(89))]);
        assert_eq!(distances_at(&hard, 190), [(Some(89), Some(0))]);
    }

    /// A deletion between a site and a template end takes its length off the
    /// distance, and an insertion adds its length, whether the read sequenced
    /// from that end holds the site or only its mate does.
    #[test]
    fn test_template_distances_count_indels_by_their_length() {
        let mut reads = SamBuilder::new().read_length(50);
        let pair = reads.add_pair(Pair::at(101, 121).cigar1("30M4D20M").cigar2("20M2I28M"));
        let (deleted, inserted) = (4, 2);
        let expected = [(Some(140 - 101 - deleted), Some(168 - 140 + inserted))];
        for records in [&pair[..], &pair[..1], &pair[1..]] {
            assert_eq!(distances_at(records, 140), expected);
        }
    }

    /// A site past its mate's 5' end lies outside the template of an FR pair.
    #[test]
    fn test_template_distances_leave_out_sites_past_the_mates_five_prime_end() {
        let mut reads = SamBuilder::new().read_length(50);
        let through = reads.add_pair(Pair::at(101, 100));
        assert_eq!(distances_at(&through, 150), []);
        assert_eq!(distances_at(&through, 100), []);
        assert_eq!(distances_at(&through, 149), [(Some(48), Some(0))]);
    }

    /// A pair whose aligned 5' ends coincide is FR, as htsjdk 5.0.0 has it, so
    /// only their shared site lies inside its template. A pair whose forward
    /// read starts one base later faces away, and each read knows only its own
    /// end, as fgbio has it.
    #[test]
    fn test_a_five_prime_tie_is_fr_and_a_pair_one_base_apart_faces_away() {
        let mut reads = SamBuilder::new().read_length(50);
        let tie = reads.add_pair(Pair::at(150, 101));
        assert_eq!(distances_at(&tie, 149), []);
        assert_eq!(distances_at(&tie, 150), [(Some(0), Some(0))]);
        assert_eq!(distances_at(&tie, 151), []);
        let away = reads.add_pair(Pair::at(151, 101));
        assert_eq!(distances_at(&away, 150), [(None, Some(0))]);
        assert_eq!(distances_at(&away, 151), [(Some(0), None)]);
        let apart = reads.add_pair(
            Pair::at(100, 200)
                .strand1(Strand::Minus)
                .strand2(Strand::Plus),
        );
        assert_eq!(distances_at(&apart, 120), [(None, Some(29))]);
        assert_eq!(distances_at(&apart, 220), [(Some(20), None)]);
    }

    #[test]
    fn test_template_distances_of_a_pair_without_a_mate_cigar_are_an_error() {
        let mut reads = SamBuilder::new().read_length(50);
        let pair = reads.add_pair(Pair::at(101, 151).name("q1"));
        let mut stripped = SamBuilder::new();
        stripped.extend(pair.into_iter().map(SamBuilder::without_mate_cigar));
        for site in [111, 161] {
            let mut evidence =
                PileupEvidence::new(stripped.to_pileup_builder(), &PileupOptions::default());
            let error = evidence.molecules("chr1", pos(site)).unwrap_err();
            let message = format!("{error:#}");
            assert!(message.contains("read q1"), "{message}");
            assert!(message.contains("MC"), "{message}");
            assert!(
                message.contains("add mate CIGARs with samtools fixmate"),
                "{message}"
            );
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
        let frag = reads.add_frag(Frag::at(100).strand(strand).cigar(cigar));
        assert_eq!(distances_at(&frag, 100), [expected]);
    }

    #[test]
    fn test_template_distances_of_a_tandem_pair_know_only_their_own_end() {
        let mut reads = SamBuilder::new().read_length(50);
        let pair = reads.add_pair(
            Pair::at(100, 200)
                .strand1(Strand::Plus)
                .strand2(Strand::Plus),
        );
        assert_eq!(distances_at(&pair, 210), [(Some(10), None)]);
    }

    /// Mates that agree are one molecule at the higher quality, and mates that
    /// disagree are an `N`, which counts as neither allele. A mate under the
    /// quality floor takes no part, so the other mate's base stands.
    #[test]
    fn test_overlapping_mates_are_called_into_one_molecule() {
        let mut reads = SamBuilder::new().read_length(50);
        for (bases2, quality2) in [('A', 40), ('C', 30), ('C', 10)] {
            reads.add_pair(
                Pair::at(101, 121)
                    .bases1("A".repeat(50))
                    .bases2(bases2.to_string().repeat(50))
                    .quals1(vec![30; 50])
                    .quals2(vec![quality2; 50]),
            );
        }
        let molecules = molecules_at(&reads, PileupOptions::default(), 130);
        let calls: Vec<(u8, u8)> = molecules.iter().map(|m| (m.base, m.quality)).collect();
        assert_eq!(calls, [(b'A', 40), (b'N', 2), (b'A', 30)]);
        assert_eq!(
            molecules[0],
            Molecule::new(b'A', 40, 29, 40).from_strand(crate::classes::Strand::Forward)
        );
    }

    /// A template holding a deletion at the site is a deletion, neither allele.
    #[test]
    fn test_a_template_holding_a_deletion_is_a_deletion_molecule() {
        let mut reads = SamBuilder::new().read_length(50);
        reads.add_pair(Pair::at(101, 151).cigar1("20M5D30M"));
        let molecules = molecules_at(&reads, PileupOptions::default(), 122);
        let bases: Vec<u8> = molecules.iter().map(|m| m.base).collect();
        assert_eq!(bases, [DELETION]);
        assert_eq!(molecules[0].left, Some(20));
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
