//! Where a read's base sits on its template.
//!
//! A template is the sequenced molecule. Its leftmost base is the 5' end of
//! its forward strand and its rightmost base is the 5' end of its reverse
//! strand. A base's distance from either end counts the template's bases, as
//! streampile counts them: from the read's own 5' end along its sequence, and
//! from its mate's 5' end by walking both reads' CIGARs, the mate's from the
//! `MC` tag, so an indel counts by its length. Soft-clipped bases count and
//! hard-clipped bases do not. The insert size (`TLEN`) is never read.

use std::io;

use noodles::core::Position;
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

    /// The 0-based distance of the base from the read's 5' end, in bases as
    /// sequenced.
    pub fn five_prime_distance(&self) -> io::Result<Option<usize>> {
        streampile::five_prime_distance(self.record, self.offset)
    }

    /// The template's bases between this base, at the 1-based `pos`, and the
    /// template's leftmost and rightmost bases, as `(left, right)`: 0 at an end.
    ///
    /// The read's own 5' end is always known: the leftmost base for a forward
    /// read and the rightmost for a reverse one. The other end is its mate's
    /// 5' end in an FR pair, and unknown otherwise. It is `None` for a base
    /// outside the template: one streampile places past the 5' end of a mate
    /// on the same contig and the other strand, as it places every base of a
    /// pair whose reads face away from each other. A read of an FR pair
    /// without a usable `MC` tag is an error naming the read.
    pub fn template_distances(
        &self,
        header: &Header,
        pos: Position,
    ) -> streampile::Result<Option<(Option<usize>, Option<usize>)>> {
        let own = self.five_prime_distance()?;
        let mate = streampile::template_end_distance(self.record, header, usize::from(pos) - 1)?;
        if mate.is_none() && has_mate_on_the_other_strand(self.record, header)? {
            return Ok(None);
        }
        let reverse = self.record.flags()?.is_reverse_complemented();
        Ok(Some(if reverse { (mate, own) } else { (own, mate) }))
    }
}

/// The complement of a base, keeping case-insensitive `N` and others as `N`.
fn complement(base: u8) -> u8 {
    crate::classes::complement(base.to_ascii_uppercase())
}

/// Whether a read's mate maps to the same contig on the other strand.
fn has_mate_on_the_other_strand<R: Record>(record: &R, header: &Header) -> io::Result<bool> {
    let flags = record.flags()?;
    if !flags.is_segmented()
        || flags.is_mate_unmapped()
        || flags.is_reverse_complemented() == flags.is_mate_reverse_complemented()
    {
        return Ok(false);
    }
    let this = record.reference_sequence_id(header).transpose()?;
    let mate = record.mate_reference_sequence_id(header).transpose()?;
    Ok(this == mate)
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
    use crate::testing::{offset_at, Frag, Pair, SamBuilder, Strand};

    type Distances = Option<(Option<usize>, Option<usize>)>;

    /// The template distances of a record's base at the 1-based `pos`.
    fn distances(builder: &SamBuilder, record: &RecordBuf, pos: usize) -> Distances {
        let offset = offset_at(record, pos).unwrap();
        ReadBase::new(record, offset)
            .template_distances(builder.header(), Position::try_from(pos).unwrap())
            .unwrap()
    }

    /// fgbio `BamsTest`: "Bams.insertCoordinates should fail on fragments and
    /// inappropriate pairs". chaff gives such records only their own end where
    /// fgbio raises.
    #[test]
    fn test_template_distances_know_only_their_own_end_on_fragments_and_inappropriate_pairs() {
        let mut builder = SamBuilder::new().read_length(10).base_quality(20);
        let own_only = Some((Some(4), None));
        for r in builder.add_frag(Frag::at(100)) {
            assert_eq!(distances(&builder, &r, 104), own_only);
        }
        let recs = builder.add_pair(Pair {
            start1: 100,
            start2: 100,
            unmapped2: true,
            ..Pair::default()
        });
        assert_eq!(distances(&builder, &recs[0], 104), own_only);
        let recs = builder.add_pair(Pair::at(100, 200));
        let r = SamBuilder::with_mate_reference_sequence_id(recs[0].clone(), 1);
        assert_eq!(distances(&builder, &r, 104), own_only);
    }

    /// fgbio `BamsTest`: "Bams.insertCoordinates should calculate insert
    /// coordinates correctly", as the distances from the template's first and
    /// last bases. chaff reads the far end from `MC`, not `TLEN`.
    #[test]
    fn test_template_distances_of_an_fr_pair() {
        let mut builder = SamBuilder::new().read_length(10).base_quality(20);
        let recs = builder.add_pair(Pair::at(100, 191));
        assert_eq!(
            distances(&builder, &recs[0], 100),
            Some((Some(0), Some(100)))
        );
        assert_eq!(
            distances(&builder, &recs[1], 200),
            Some((Some(100), Some(0)))
        );
    }

    /// fgbio `BamsTest`: "Bams.positionFromOtherEndOfTemplate should return
    /// None for anything that's not an FR mapped pair".
    #[test]
    fn test_distance_from_the_other_end_is_none_unless_fr_pair() {
        let mut builder = SamBuilder::new();
        let header = builder.header().clone();
        let check = |recs: Vec<RecordBuf>| {
            for r in recs.iter().filter(|r| !r.flags().is_unmapped()) {
                for pos in [r.alignment_start(), r.alignment_end()] {
                    let pos = usize::from(pos.unwrap()) - 1;
                    let distance = streampile::template_end_distance(r, &header, pos);
                    assert_eq!(distance.unwrap(), None, "{r:?}");
                }
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
        let from_other_end = |i: usize, pos: usize| {
            streampile::template_end_distance(&recs[i], builder.header(), pos - 1)
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

    /// fgbio `PileupTest`: "BaseEntry should report the correct
    /// offsets/positions/bases", and the matching `PileupBuilderTest` case.
    /// fgbio's 1-based positions in read order are one more than these
    /// 0-based distances from the 5' end.
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
        assert_eq!(p1.five_prime_distance().unwrap(), Some(4));

        let p2 = ReadBase::new(&recs[1], 4);
        assert_eq!(p2.base(), Some(b'C'));
        assert_eq!(p2.base_in_read_orientation().unwrap(), Some(b'G'));
        assert_eq!(p2.quality().unwrap(), Some(35));
        assert_eq!(p2.offset, 4);
        assert_eq!(p2.five_prime_distance().unwrap(), Some(45));
    }

    #[test]
    fn test_template_distances_count_soft_clips_and_not_hard_clips() {
        let mut builder = SamBuilder::new().read_length(50);
        let soft = builder.add_pair(Pair {
            start1: 101,
            start2: 151,
            cigar1: Some("5S45M".into()),
            cigar2: Some("40M10S".into()),
            ..Pair::default()
        });
        assert_eq!(
            distances(&builder, &soft[0], 101),
            Some((Some(5), Some(99)))
        );
        assert_eq!(
            distances(&builder, &soft[1], 190),
            Some((Some(94), Some(10)))
        );
        let hard = builder.add_pair(Pair {
            start1: 101,
            start2: 151,
            bases1: Some("A".repeat(45)),
            bases2: Some("C".repeat(40)),
            cigar1: Some("5H45M".into()),
            cigar2: Some("40M10H".into()),
            ..Pair::default()
        });
        assert_eq!(
            distances(&builder, &hard[0], 101),
            Some((Some(0), Some(89)))
        );
        assert_eq!(
            distances(&builder, &hard[1], 190),
            Some((Some(89), Some(0)))
        );
    }

    /// A deletion between a base and a template end takes its length off the
    /// distance, and an insertion adds its length, on either read's side.
    #[test]
    fn test_template_distances_count_indels_by_their_length() {
        let mut builder = SamBuilder::new().read_length(50);
        let recs = builder.add_pair(Pair {
            start1: 101,
            start2: 121,
            cigar1: Some("30M4D20M".into()),
            cigar2: Some("20M2I28M".into()),
            ..Pair::default()
        });
        let (deleted, inserted) = (4, 2);
        let (left, right) = (140 - 101 - deleted, 168 - 140 + inserted);
        for r in &recs {
            assert_eq!(distances(&builder, r, 140), Some((Some(left), Some(right))));
        }
    }

    /// A base past its mate's 5' end lies outside the template, and so does
    /// every base of a pair whose reads face away from each other.
    #[test]
    fn test_template_distances_leave_out_bases_past_the_mates_five_prime_end() {
        let mut builder = SamBuilder::new().read_length(50);
        let through = builder.add_pair(Pair::at(101, 100));
        assert_eq!(distances(&builder, &through[0], 150), None);
        assert_eq!(distances(&builder, &through[1], 100), None);
        assert_eq!(
            distances(&builder, &through[0], 149),
            Some((Some(48), Some(0)))
        );
        let away = builder.add_pair(Pair {
            start1: 100,
            start2: 200,
            strand1: Strand::Minus,
            strand2: Strand::Plus,
            ..Pair::default()
        });
        assert_eq!(distances(&builder, &away[0], 120), None);
        assert_eq!(distances(&builder, &away[1], 220), None);
    }

    #[test]
    fn test_template_distances_of_a_pair_without_a_mate_cigar_are_an_error() {
        let mut builder = SamBuilder::new().read_length(50);
        let recs = builder.add_pair(Pair {
            name: Some("q1".into()),
            ..Pair::at(101, 151)
        });
        for r in recs {
            let pos = usize::from(r.alignment_start().unwrap()) + 10;
            let r = SamBuilder::without_mate_cigar(r);
            let error = ReadBase::new(&r, 10)
                .template_distances(builder.header(), Position::try_from(pos).unwrap())
                .unwrap_err();
            assert!(error.to_string().contains("read q1"), "{error}");
            assert!(error.to_string().contains("MC"), "{error}");
        }
    }

    #[test]
    fn test_template_distances_ignore_the_insert_size() {
        let mut builder = SamBuilder::new().read_length(50);
        let recs = builder.add_pair(Pair::at(101, 151));
        let cases = [(&recs[0], 120, (19, 80)), (&recs[1], 160, (59, 40))];
        for (r, pos, (left, right)) in cases {
            let mut r = r.clone();
            *r.template_length_mut() = 7;
            assert_eq!(
                distances(&builder, &r, pos),
                Some((Some(left), Some(right)))
            );
        }
    }

    #[rstest]
    #[case(Strand::Plus, "5S45M", (Some(5), None))]
    #[case(Strand::Minus, "45M5S", (None, Some(49)))]
    fn test_template_distances_of_a_fragment_know_only_their_own_end(
        #[case] strand: Strand,
        #[case] cigar: &str,
        #[case] expected: (Option<usize>, Option<usize>),
    ) {
        let mut builder = SamBuilder::new().read_length(50);
        let recs = builder.add_frag(Frag {
            start: 100,
            strand,
            cigar: Some(cigar.into()),
            ..Frag::default()
        });
        assert_eq!(distances(&builder, &recs[0], 100), Some(expected));
    }

    #[test]
    fn test_template_distances_of_a_tandem_pair_know_only_their_own_end() {
        let mut builder = SamBuilder::new().read_length(50);
        let recs = builder.add_pair(Pair {
            start1: 100,
            start2: 200,
            strand1: Strand::Plus,
            strand2: Strand::Plus,
            ..Pair::default()
        });
        assert_eq!(distances(&builder, &recs[1], 210), Some((Some(10), None)));
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
