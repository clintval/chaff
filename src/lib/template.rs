//! Template geometry from single alignment records.
//!
//! A template is the sequenced molecule. Its leftmost base is the 5' end of
//! its forward strand and its rightmost base is the 5' end of its reverse
//! strand. Every position here is 1-based on the reference.
//!
//! The ends come from the unclipped 5' ends of both mates: the record's own
//! from its position and CIGAR, its mate's from the mate position and the `MC`
//! (mate CIGAR) tag. The insert size (`TLEN`) is never read.

use std::io;

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::alignment::Record;
use noodles::sam::Header;

/// One read's aligned base at a pileup position.
#[derive(Debug)]
pub struct ReadBase<'a, R> {
    /// The read.
    pub record: &'a R,
    /// The 0-based offset of the base in the read's stored sequence.
    pub offset: usize,
}

impl<'a, R> Clone for ReadBase<'a, R> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, R> Copy for ReadBase<'a, R> {}

impl<'a, R: Record> ReadBase<'a, R> {
    /// A read's base at `offset`.
    pub fn new(record: &'a R, offset: usize) -> Self {
        Self { record, offset }
    }

    /// The upper-cased base as aligned to the forward strand.
    pub fn base(&self) -> Option<u8> {
        self.record
            .sequence()
            .get(self.offset)
            .map(|b| b.to_ascii_uppercase())
    }

    /// The base as read by the sequencer: complemented on a reverse-strand read.
    pub fn base_in_read_orientation(&self) -> io::Result<Option<u8>> {
        let reverse = self.record.flags()?.is_reverse_complemented();
        Ok(self.base().map(|b| if reverse { complement(b) } else { b }))
    }

    /// The base quality, or `None` when the read stores no qualities.
    pub fn quality(&self) -> io::Result<Option<u8>> {
        self.record
            .quality_scores()
            .iter()
            .nth(self.offset)
            .transpose()
    }

    /// The 1-based position of the base in the read's stored sequence.
    pub fn position_in_read(&self) -> usize {
        self.offset + 1
    }

    /// The 1-based position of the base counted from the read's 5' end, in
    /// sequencing order; soft-clipped bases count.
    pub fn position_in_read_in_read_order(&self) -> io::Result<usize> {
        if self.record.flags()?.is_reverse_complemented() {
            Ok(self.record.sequence().len() - self.offset)
        } else {
            Ok(self.position_in_read())
        }
    }
}

/// The complement of a base, keeping case-insensitive `N` and others as `N`.
fn complement(base: u8) -> u8 {
    crate::classes::complement(base.to_ascii_uppercase())
}

/// The ends of a template that one record can see.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct TemplateEnds {
    /// The leftmost base, the forward strand's 5' end, when known.
    pub start: Option<i64>,
    /// The rightmost base, the reverse strand's 5' end, when known.
    pub end: Option<i64>,
}

impl TemplateEnds {
    /// Whether `pos` falls outside the template; `false` when an end is unknown.
    pub fn excludes(&self, pos: i64) -> bool {
        matches!(self.start, Some(start) if pos < start)
            || matches!(self.end, Some(end) if pos > end)
    }
}

/// The length of the clipping (soft and hard) at the start and at the end of a
/// list of CIGAR operations.
fn clipping(ops: &[(Kind, usize)]) -> (usize, usize) {
    let is_clip = |kind: &Kind| matches!(kind, Kind::SoftClip | Kind::HardClip);
    let leading = ops
        .iter()
        .take_while(|(kind, _)| is_clip(kind))
        .map(|(_, len)| len)
        .sum();
    let trailing = ops
        .iter()
        .rev()
        .take_while(|(kind, _)| is_clip(kind))
        .map(|(_, len)| len)
        .sum();
    (leading, trailing)
}

/// The number of reference bases a list of CIGAR operations spans.
fn reference_span(ops: &[(Kind, usize)]) -> usize {
    ops.iter()
        .filter(|(kind, _)| kind.consumes_reference())
        .map(|(_, len)| len)
        .sum()
}

/// A record's CIGAR as `(kind, length)` pairs.
fn cigar_ops<R: Record>(record: &R) -> io::Result<Vec<(Kind, usize)>> {
    record
        .cigar()
        .iter()
        .map(|op| op.map(|op| (op.kind(), op.len())))
        .collect()
}

/// Parse a CIGAR string such as `5S40M2D10M` into `(kind, length)` pairs.
pub fn parse_cigar(text: &[u8]) -> io::Result<Vec<(Kind, usize)>> {
    let invalid = || io::Error::new(io::ErrorKind::InvalidData, "invalid CIGAR string");
    let mut ops = Vec::new();
    let mut len: usize = 0;
    let mut has_len = false;
    for &byte in text {
        if byte.is_ascii_digit() {
            len = len
                .checked_mul(10)
                .and_then(|n| n.checked_add(usize::from(byte - b'0')))
                .ok_or_else(invalid)?;
            has_len = true;
            continue;
        }
        let kind = match byte {
            b'M' => Kind::Match,
            b'I' => Kind::Insertion,
            b'D' => Kind::Deletion,
            b'N' => Kind::Skip,
            b'S' => Kind::SoftClip,
            b'H' => Kind::HardClip,
            b'P' => Kind::Pad,
            b'=' => Kind::SequenceMatch,
            b'X' => Kind::SequenceMismatch,
            _ => return Err(invalid()),
        };
        if !has_len {
            return Err(invalid());
        }
        ops.push((kind, len));
        len = 0;
        has_len = false;
    }
    if has_len {
        return Err(invalid());
    }
    Ok(ops)
}

/// The record's unclipped 5' end: the unclipped start of a forward read or the
/// unclipped end of a reverse read, counting soft and hard clips.
pub fn unclipped_five_prime<R: Record>(record: &R) -> io::Result<Option<i64>> {
    let flags = record.flags()?;
    if flags.is_unmapped() {
        return Ok(None);
    }
    let Some(start) = record.alignment_start().transpose()? else {
        return Ok(None);
    };
    let start = usize::from(start) as i64;
    let ops = cigar_ops(record)?;
    let (leading, trailing) = clipping(&ops);
    if flags.is_reverse_complemented() {
        let span = reference_span(&ops).max(1) as i64;
        Ok(Some(start + span - 1 + trailing as i64))
    } else {
        Ok(Some(start - leading as i64))
    }
}

/// The mate's unclipped 5' end from its position, strand, and `MC` tag, or
/// `None` without an `MC` tag.
pub fn mate_unclipped_five_prime<R: Record>(record: &R) -> io::Result<Option<i64>> {
    let flags = record.flags()?;
    let data = record.data();
    let Some(value) = data.get(&Tag::MATE_CIGAR).transpose()? else {
        return Ok(None);
    };
    let Value::String(text) = value else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "the MC tag is not a string",
        ));
    };
    let Some(mate_start) = record.mate_alignment_start().transpose()? else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "a record with an MC tag has no mate position",
        ));
    };
    let mate_start = usize::from(mate_start) as i64;
    let ops = parse_cigar(text)?;
    let (leading, trailing) = clipping(&ops);
    if flags.is_mate_reverse_complemented() {
        let span = reference_span(&ops).max(1) as i64;
        Ok(Some(mate_start + span - 1 + trailing as i64))
    } else {
        Ok(Some(mate_start - leading as i64))
    }
}

/// Whether the record and its mate map to the same reference sequence, both
/// mapped.
fn is_mapped_pair_on_one_contig<R: Record>(record: &R, header: &Header) -> io::Result<bool> {
    let flags = record.flags()?;
    if !flags.is_segmented() || flags.is_unmapped() || flags.is_mate_unmapped() {
        return Ok(false);
    }
    let this = record.reference_sequence_id(header).transpose()?;
    let mate = record.mate_reference_sequence_id(header).transpose()?;
    Ok(this.is_some() && this == mate)
}

/// The template ends a record can see: its own unclipped 5' end always, and
/// the far end from the mate's unclipped 5' end (`MC` tag) when the record is
/// in an FR pair. A record whose mate maps to the same contig on the other
/// strand must carry `MC`.
pub fn template_ends<R: Record>(record: &R, header: &Header) -> io::Result<TemplateEnds> {
    let Some(own) = unclipped_five_prime(record)? else {
        return Ok(TemplateEnds::default());
    };
    let flags = record.flags()?;
    let reverse = flags.is_reverse_complemented();
    let own_only = if reverse {
        TemplateEnds {
            start: None,
            end: Some(own),
        }
    } else {
        TemplateEnds {
            start: Some(own),
            end: None,
        }
    };
    if !is_mapped_pair_on_one_contig(record, header)?
        || reverse == flags.is_mate_reverse_complemented()
    {
        return Ok(own_only);
    }
    let Some(mate) = mate_unclipped_five_prime(record)? else {
        let name = record
            .name()
            .map(|n| String::from_utf8_lossy(n).into_owned())
            .unwrap_or_else(|| "*".to_string());
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("read {name} has a mapped mate but no MC tag to place the far template end"),
        ));
    };
    let (start, end) = if reverse { (mate, own) } else { (own, mate) };
    if start < end {
        Ok(TemplateEnds {
            start: Some(start),
            end: Some(end),
        })
    } else {
        Ok(own_only)
    }
}

/// Which records and bases count as evidence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadFilter {
    /// Records below this mapping quality are left out.
    pub min_mapping_quality: u8,
    /// Bases below this quality are left out.
    pub min_base_quality: u8,
    /// Keep only paired records whose mate is also mapped.
    pub paired_reads_only: bool,
}

impl Default for ReadFilter {
    fn default() -> Self {
        Self {
            min_mapping_quality: 20,
            min_base_quality: 20,
            paired_reads_only: false,
        }
    }
}

impl ReadFilter {
    /// Whether a record passes the static filters: mapped, primary, not a
    /// duplicate, at or above the mapping quality floor, and paired with a
    /// mapped mate when only pairs are kept. A record without a mapping
    /// quality (255) passes the floor, as in htsjdk.
    pub fn accepts<R: Record>(&self, record: &R) -> io::Result<bool> {
        let flags = record.flags()?;
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
        {
            return Ok(false);
        }
        if let Some(mapq) = record.mapping_quality().transpose()? {
            if u8::from(mapq) < self.min_mapping_quality {
                return Ok(false);
            }
        }
        if self.paired_reads_only && (!flags.is_segmented() || flags.is_mate_unmapped()) {
            return Ok(false);
        }
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use noodles::sam::alignment::RecordBuf;

    use super::*;
    use crate::testing::{Frag, Pair, SamBuilder, Strand};

    #[test]
    fn test_parse_cigar() {
        let ops = parse_cigar(b"5S40M2D10M3H").unwrap();
        assert_eq!(
            ops,
            vec![
                (Kind::SoftClip, 5),
                (Kind::Match, 40),
                (Kind::Deletion, 2),
                (Kind::Match, 10),
                (Kind::HardClip, 3)
            ]
        );
        assert_eq!(clipping(&ops), (5, 3));
        assert_eq!(reference_span(&ops), 52);
        assert!(parse_cigar(b"M").is_err());
        assert!(parse_cigar(b"10").is_err());
        assert!(parse_cigar(b"10Q").is_err());
    }

    /// fgbio `BamsTest`: "Bams.insertCoordinates should fail on fragments and
    /// inappropriate pairs". chaff gives such records only their own end where
    /// fgbio raises.
    #[test]
    fn test_template_ends_know_only_their_own_end_on_fragments_and_inappropriate_pairs() {
        let mut builder = SamBuilder::new().read_length(10).base_quality(20);
        let own_only = TemplateEnds {
            start: Some(100),
            end: None,
        };
        for r in builder.add_frag(Frag::at(100)) {
            assert_eq!(template_ends(&r, builder.header()).unwrap(), own_only);
        }
        let recs = builder.add_pair(Pair {
            start1: 100,
            start2: 100,
            unmapped2: true,
            ..Pair::default()
        });
        assert_eq!(template_ends(&recs[0], builder.header()).unwrap(), own_only);
        let recs = builder.add_pair(Pair::at(100, 200));
        let r = SamBuilder::with_mate_reference_sequence_id(recs[0].clone(), 1);
        assert_eq!(template_ends(&r, builder.header()).unwrap(), own_only);
    }

    /// fgbio `BamsTest`: "Bams.insertCoordinates should calculate insert
    /// coordinates correctly". chaff reads the far end from `MC`, not `TLEN`.
    #[test]
    fn test_template_ends_of_an_fr_pair() {
        let mut builder = SamBuilder::new().read_length(10).base_quality(20);
        for r in builder.add_pair(Pair::at(100, 191)) {
            assert_eq!(
                template_ends(&r, builder.header()).unwrap(),
                TemplateEnds {
                    start: Some(100),
                    end: Some(200)
                }
            );
        }
    }

    /// fgbio `BamsTest`: "Bams.positionFromOtherEndOfTemplate should return
    /// None for anything that's not an FR mapped pair".
    #[test]
    fn test_template_ends_have_no_far_end_unless_fr_pair() {
        let mut builder = SamBuilder::new();
        let header = builder.header().clone();
        let check = |recs: Vec<RecordBuf>| {
            for r in recs.iter().filter(|r| !r.flags().is_unmapped()) {
                let ends = template_ends(r, &header).unwrap();
                assert!(ends.start.is_none() || ends.end.is_none(), "{ends:?}");
            }
        };
        check(builder.add_frag(Frag::at(100)));
        check(builder.add_pair(Pair {
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
            check(builder.add_pair(Pair {
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
        let mut builder = SamBuilder::new().read_length(50);
        let recs = builder.add_pair(Pair::at(101, 151));
        let header = builder.header();
        let ends: Vec<TemplateEnds> = recs
            .iter()
            .map(|r| template_ends(r, header).unwrap())
            .collect();
        let from_other_end = |i: usize, pos: i64| match i {
            0 => ends[0].end.unwrap() - pos + 1,
            _ => pos - ends[1].start.unwrap() + 1,
        };
        assert_eq!(from_other_end(0, 101), 100);
        assert_eq!(from_other_end(0, 111), 90);
        assert_eq!(from_other_end(0, 151), 50);
        assert_eq!(from_other_end(0, 200), 1);
        assert_eq!(from_other_end(1, 200), 100);
        assert_eq!(from_other_end(1, 190), 90);
        assert_eq!(from_other_end(1, 150), 50);
        assert_eq!(from_other_end(1, 101), 1);
    }

    /// fgbio `PileupTest`: "BaseEntry should report the correct
    /// offsets/positions/bases", and the matching `PileupBuilderTest` case.
    #[test]
    fn test_read_base_reports_offsets_positions_and_bases() {
        let mut builder = SamBuilder::new().read_length(50).base_quality(35);
        let recs = builder.add_pair(Pair {
            name: Some("q1".into()),
            start1: 101,
            start2: 201,
            bases1: Some("A".repeat(50)),
            bases2: Some("C".repeat(50)),
            ..Pair::default()
        });

        let p1 = ReadBase::new(&recs[0], 4);
        assert_eq!(p1.base(), Some(b'A'));
        assert_eq!(p1.base_in_read_orientation().unwrap(), Some(b'A'));
        assert_eq!(p1.quality().unwrap(), Some(35));
        assert_eq!(p1.offset, 4);
        assert_eq!(p1.position_in_read(), 5);
        assert_eq!(p1.position_in_read_in_read_order().unwrap(), 5);

        let p2 = ReadBase::new(&recs[1], 4);
        assert_eq!(p2.base(), Some(b'C'));
        assert_eq!(p2.base_in_read_orientation().unwrap(), Some(b'G'));
        assert_eq!(p2.quality().unwrap(), Some(35));
        assert_eq!(p2.offset, 4);
        assert_eq!(p2.position_in_read(), 5);
        assert_eq!(p2.position_in_read_in_read_order().unwrap(), 46);
    }

    #[test]
    fn test_template_ends_of_an_fr_pair_from_the_mate_cigar() {
        let mut builder = SamBuilder::new().read_length(50);
        let recs = builder.add_pair(Pair {
            start1: 101,
            start2: 151,
            cigar1: Some("5S45M".into()),
            cigar2: Some("40M10S".into()),
            ..Pair::default()
        });
        let header = builder.header();
        for r in &recs {
            assert_eq!(
                template_ends(r, header).unwrap(),
                TemplateEnds {
                    start: Some(96),
                    end: Some(200)
                }
            );
        }
    }

    #[test]
    fn test_template_ends_of_a_pair_without_a_mate_cigar_are_an_error() {
        let mut builder = SamBuilder::new().read_length(50);
        let recs = builder.add_pair(Pair {
            name: Some("q1".into()),
            ..Pair::at(101, 151)
        });
        for r in recs {
            let r = SamBuilder::without_mate_cigar(r);
            let error = template_ends(&r, builder.header()).unwrap_err();
            assert!(error.to_string().contains("read q1"), "{error}");
            assert!(error.to_string().contains("MC"), "{error}");
        }
    }

    #[test]
    fn test_template_ends_ignore_the_insert_size() {
        let mut builder = SamBuilder::new().read_length(50);
        for mut r in builder.add_pair(Pair::at(101, 151)) {
            *r.template_length_mut() = 7;
            assert_eq!(
                template_ends(&r, builder.header()).unwrap(),
                TemplateEnds {
                    start: Some(101),
                    end: Some(200)
                }
            );
        }
    }

    #[rstest]
    #[case(Strand::Plus, "5S45M", TemplateEnds { start: Some(95), end: None })]
    #[case(Strand::Minus, "45M5S", TemplateEnds { start: None, end: Some(149) })]
    fn test_template_ends_of_a_fragment_know_only_their_own_end(
        #[case] strand: Strand,
        #[case] cigar: &str,
        #[case] expected: TemplateEnds,
    ) {
        let mut builder = SamBuilder::new().read_length(50);
        let recs = builder.add_frag(Frag {
            start: 100,
            strand,
            cigar: Some(cigar.into()),
            ..Frag::default()
        });
        assert_eq!(template_ends(&recs[0], builder.header()).unwrap(), expected);
    }

    #[test]
    fn test_template_ends_of_a_tandem_pair_know_only_their_own_end() {
        let mut builder = SamBuilder::new().read_length(50);
        let recs = builder.add_pair(Pair {
            start1: 100,
            start2: 200,
            strand1: Strand::Plus,
            strand2: Strand::Plus,
            ..Pair::default()
        });
        let ends = template_ends(&recs[1], builder.header()).unwrap();
        assert_eq!(
            ends,
            TemplateEnds {
                start: Some(200),
                end: None
            }
        );
    }

    #[test]
    fn test_template_ends_excludes() {
        let ends = TemplateEnds {
            start: Some(10),
            end: Some(20),
        };
        assert!(ends.excludes(9));
        assert!(!ends.excludes(10));
        assert!(!ends.excludes(20));
        assert!(ends.excludes(21));
        assert!(!TemplateEnds::default().excludes(1));
    }

    #[test]
    fn test_read_filter_flags_and_mapping_quality() {
        let filter = ReadFilter {
            min_mapping_quality: 10,
            ..ReadFilter::default()
        };
        let mut builder = SamBuilder::new();
        let low = builder.add_frag(Frag {
            start: 101,
            mapq: 9,
            ..Frag::default()
        });
        let ok = builder.add_frag(Frag {
            start: 101,
            mapq: 10,
            ..Frag::default()
        });
        assert!(!filter.accepts(&low[0]).unwrap());
        assert!(filter.accepts(&ok[0]).unwrap());
        for flag in [0x100, 0x400, 0x800, 0x4] {
            let r = SamBuilder::with_flags(ok[0].clone(), flag);
            assert!(!filter.accepts(&r).unwrap(), "flag {flag:#x}");
        }
    }
}
