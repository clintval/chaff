//! Test builders that mirror fgbio's `SamBuilder` and `VcfBuilder`, and an
//! in-memory pileup over built records for driving the evidence layer. Builder
//! fields carry fgbio's parameter names and defaults.
#![allow(missing_docs)]

use std::collections::BTreeMap;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Result};
use noodles::core::Position;
use noodles::sam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::{Flags, MappingQuality};
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{Cigar, Data, QualityScores, Sequence};
use noodles::sam::alignment::RecordBuf;
use noodles::sam::header::record::value::map::{Map, ReferenceSequence};
use std::num::NonZeroUsize;

use crate::evidence::PileupSource;
use crate::template::{parse_cigar, ReadBase};

/// The strand a built read aligns to.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Strand {
    /// The forward strand.
    #[default]
    Plus,
    /// The reverse strand.
    Minus,
}

/// A read pair to build, with fgbio `SamBuilder.addPair` defaults.
#[derive(Clone, Debug)]
pub struct Pair {
    pub name: Option<String>,
    pub bases1: Option<String>,
    pub bases2: Option<String>,
    pub contig: usize,
    pub contig2: Option<usize>,
    pub start1: usize,
    pub start2: usize,
    pub unmapped1: bool,
    pub unmapped2: bool,
    pub cigar1: Option<String>,
    pub cigar2: Option<String>,
    pub mapq1: u8,
    pub mapq2: u8,
    pub strand1: Strand,
    pub strand2: Strand,
}

impl Default for Pair {
    fn default() -> Self {
        Self {
            name: None,
            bases1: None,
            bases2: None,
            contig: 0,
            contig2: None,
            start1: 0,
            start2: 0,
            unmapped1: false,
            unmapped2: false,
            cigar1: None,
            cigar2: None,
            mapq1: 60,
            mapq2: 60,
            strand1: Strand::Plus,
            strand2: Strand::Minus,
        }
    }
}

impl Pair {
    /// A default FR pair at two starts.
    pub fn at(start1: usize, start2: usize) -> Self {
        Self {
            start1,
            start2,
            ..Self::default()
        }
    }

    /// A default FR pair at two starts with the same base repeated on both reads.
    pub fn filled(start1: usize, start2: usize, base: char, read_length: usize) -> Self {
        Self {
            start1,
            start2,
            bases1: Some(base.to_string().repeat(read_length)),
            bases2: Some(base.to_string().repeat(read_length)),
            ..Self::default()
        }
    }
}

/// A fragment (unpaired read) to build, with fgbio `SamBuilder.addFrag` defaults.
#[derive(Clone, Debug)]
pub struct Frag {
    pub name: Option<String>,
    pub bases: Option<String>,
    pub contig: usize,
    pub start: usize,
    pub unmapped: bool,
    pub cigar: Option<String>,
    pub mapq: u8,
    pub strand: Strand,
}

impl Default for Frag {
    fn default() -> Self {
        Self {
            name: None,
            bases: None,
            contig: 0,
            start: 0,
            unmapped: false,
            cigar: None,
            mapq: 60,
            strand: Strand::Plus,
        }
    }
}

impl Frag {
    /// A default forward fragment at a start.
    pub fn at(start: usize) -> Self {
        Self {
            start,
            ..Self::default()
        }
    }
}

/// Builds alignment records the way fgbio's `SamBuilder` does: contigs `chr1`
/// to `chr22`, `chrX`, `chrY`, `chrM` of 200 Mbp, one read group, sequential
/// names, and full mate information on pairs.
#[derive(Clone, Debug)]
pub struct SamBuilder {
    read_length: usize,
    base_quality: u8,
    header: sam::Header,
    records: Vec<RecordBuf>,
    counter: usize,
    seed: u64,
}

impl Default for SamBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// The contig names of fgbio's default sequence dictionary.
pub fn default_contigs() -> Vec<String> {
    (1..=22)
        .map(|i| i.to_string())
        .chain(["X", "Y", "M"].map(String::from))
        .map(|c| format!("chr{c}"))
        .collect()
}

impl SamBuilder {
    /// A builder with read length 100 and base quality 30.
    pub fn new() -> Self {
        let mut builder = sam::Header::builder();
        for name in default_contigs() {
            builder = builder.add_reference_sequence(
                name,
                Map::<ReferenceSequence>::new(NonZeroUsize::new(200_000_000).unwrap()),
            );
        }
        Self {
            read_length: 100,
            base_quality: 30,
            header: builder.build(),
            records: Vec::new(),
            counter: 0,
            seed: 42,
        }
    }

    /// Set the length of reads with default bases.
    pub fn read_length(mut self, read_length: usize) -> Self {
        self.read_length = read_length;
        self
    }

    /// Set the quality of every base.
    pub fn base_quality(mut self, base_quality: u8) -> Self {
        self.base_quality = base_quality;
        self
    }

    /// Declare the records coordinate sorted (`@HD SO:coordinate`), as fgbio's
    /// `sort = Some(SamOrder.Coordinate)` does.
    pub fn coordinate_sorted(mut self) -> Self {
        use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;
        use noodles::sam::header::record::value::map::Header;
        let hd = Map::<Header>::builder()
            .insert(SORT_ORDER, "coordinate")
            .build()
            .expect("a valid @HD");
        *self.header.header_mut() = Some(hd);
        self
    }

    /// The header of the built records.
    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Every built record, in insertion order.
    pub fn records(&self) -> &[RecordBuf] {
        &self.records
    }

    /// Write the records as a BAM in coordinate order (fgbio's `write`).
    pub fn write_bam(&self, path: &Path) -> PathBuf {
        use noodles::bam;
        use noodles::sam::alignment::io::Write as _;
        let mut writer = bam::io::Writer::new(std::fs::File::create(path).unwrap());
        writer.write_header(&self.header).unwrap();
        for record in &self.pileup().records {
            writer.write_alignment_record(&self.header, record).unwrap();
        }
        writer.try_finish().unwrap();
        path.to_path_buf()
    }

    /// Add records built elsewhere (fgbio's `++=`).
    pub fn extend(&mut self, records: impl IntoIterator<Item = RecordBuf>) {
        self.records.extend(records);
    }

    /// A pileup over the built records, in coordinate order.
    pub fn pileup(&self) -> TestPileup {
        let mut records = self.records.clone();
        records.sort_by_key(|r| (r.reference_sequence_id(), r.alignment_start()));
        TestPileup {
            header: self.header.clone(),
            records,
        }
    }

    fn next_name(&mut self) -> String {
        let name = format!("{:04}", self.counter);
        self.counter += 1;
        name
    }

    fn random_bases(&mut self) -> String {
        (0..self.read_length)
            .map(|_| {
                self.seed = self
                    .seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(self.seed >> 62) as usize] as char
            })
            .collect()
    }

    fn cigar(&self, cigar: &Option<String>) -> Cigar {
        let text = cigar
            .clone()
            .unwrap_or_else(|| format!("{}M", self.read_length));
        parse_cigar(text.as_bytes())
            .unwrap()
            .into_iter()
            .map(|(kind, len)| Op::new(kind, len))
            .collect()
    }

    fn quals(&self, len: usize) -> QualityScores {
        QualityScores::from(vec![self.base_quality; len])
    }

    #[allow(clippy::too_many_arguments)]
    fn record(
        &self,
        name: &str,
        bases: &str,
        contig: usize,
        start: usize,
        unmapped: bool,
        cigar: &Cigar,
        mapq: u8,
        strand: Strand,
        mut flags: Flags,
    ) -> RecordBuf {
        if unmapped {
            flags |= Flags::UNMAPPED;
        }
        if strand == Strand::Minus {
            flags |= Flags::REVERSE_COMPLEMENTED;
        }
        let mut builder = RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_sequence(Sequence::from(bases.as_bytes().to_vec()))
            .set_quality_scores(self.quals(bases.len()))
            .set_data(Data::from_iter([(Tag::READ_GROUP, Value::from("A"))]));
        if !unmapped && start > 0 {
            builder = builder
                .set_reference_sequence_id(contig)
                .set_alignment_start(Position::try_from(start).unwrap())
                .set_cigar(cigar.clone());
            if let Some(mapq) = MappingQuality::new(mapq) {
                builder = builder.set_mapping_quality(mapq);
            }
        }
        builder.build()
    }

    /// Add a pair of reads with full mate information, as htsjdk's
    /// `SamPairUtil.setProperPairAndMateInfo` sets it.
    pub fn add_pair(&mut self, pair: Pair) -> Vec<RecordBuf> {
        let name = pair.name.clone().unwrap_or_else(|| self.next_name());
        let bases1 = pair.bases1.clone().unwrap_or_else(|| self.random_bases());
        let bases2 = pair.bases2.clone().unwrap_or_else(|| self.random_bases());
        let cigar1 = self.cigar(&pair.cigar1);
        let cigar2 = self.cigar(&pair.cigar2);
        let unmapped1 = pair.unmapped1 || pair.start1 == 0;
        let unmapped2 = pair.unmapped2 || pair.start2 == 0;
        let first = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        let second = Flags::SEGMENTED | Flags::LAST_SEGMENT;
        let mut r1 = self.record(
            &name,
            &bases1,
            pair.contig,
            pair.start1,
            unmapped1,
            &cigar1,
            pair.mapq1,
            pair.strand1,
            first,
        );
        let mut r2 = self.record(
            &name,
            &bases2,
            pair.contig2.unwrap_or(pair.contig),
            pair.start2,
            unmapped2,
            &cigar2,
            pair.mapq2,
            pair.strand2,
            second,
        );
        set_mate_info(&mut r1, &mut r2, &self.header);
        let recs = vec![r1, r2];
        self.records.extend(recs.iter().cloned());
        recs
    }

    /// Add an unpaired read.
    pub fn add_frag(&mut self, frag: Frag) -> Vec<RecordBuf> {
        let name = frag.name.clone().unwrap_or_else(|| self.next_name());
        let bases = frag.bases.clone().unwrap_or_else(|| self.random_bases());
        let cigar = self.cigar(&frag.cigar);
        let rec = self.record(
            &name,
            &bases,
            frag.contig,
            frag.start,
            frag.unmapped || frag.start == 0,
            &cigar,
            frag.mapq,
            frag.strand,
            Flags::empty(),
        );
        self.records.push(rec.clone());
        vec![rec]
    }

    /// A copy of `record` claiming its mate maps to another reference sequence.
    pub fn with_mate_reference_sequence_id(mut record: RecordBuf, id: usize) -> RecordBuf {
        *record.mate_reference_sequence_id_mut() = Some(id);
        record
    }

    /// A copy of `record` without its `MC` tag.
    pub fn without_mate_cigar(mut record: RecordBuf) -> RecordBuf {
        record.data_mut().remove(&Tag::MATE_CIGAR);
        record
    }

    /// A copy of `record` with extra flag bits set.
    pub fn with_flags(mut record: RecordBuf, bits: u16) -> RecordBuf {
        let flags = record.flags() | Flags::from_bits_truncate(bits);
        *record.flags_mut() = flags;
        record
    }

    /// A copy of `record` with new bases (fgbio's `rec.bases = ...`).
    pub fn with_bases(mut record: RecordBuf, bases: &str) -> RecordBuf {
        *record.sequence_mut() = Sequence::from(bases.as_bytes().to_vec());
        record
    }
}

/// The aligned 5' end of a record, for the insert size.
fn aligned_five_prime(record: &RecordBuf) -> i64 {
    if record.flags().is_reverse_complemented() {
        usize::from(record.alignment_end().unwrap()) as i64
    } else {
        usize::from(record.alignment_start().unwrap()) as i64
    }
}

/// htsjdk's `SamPairUtil.computeInsertSize`.
fn insert_size(first: &RecordBuf, second: &RecordBuf) -> i32 {
    if first.flags().is_unmapped()
        || second.flags().is_unmapped()
        || first.reference_sequence_id() != second.reference_sequence_id()
    {
        return 0;
    }
    let first5 = aligned_five_prime(first);
    let second5 = aligned_five_prime(second);
    let adjustment = if second5 >= first5 { 1 } else { -1 };
    (second5 - first5 + adjustment) as i32
}

fn cigar_text(record: &RecordBuf) -> String {
    record
        .cigar()
        .as_ref()
        .iter()
        .map(|op| {
            let c = match op.kind() {
                Kind::Match => 'M',
                Kind::Insertion => 'I',
                Kind::Deletion => 'D',
                Kind::Skip => 'N',
                Kind::SoftClip => 'S',
                Kind::HardClip => 'H',
                Kind::Pad => 'P',
                Kind::SequenceMatch => '=',
                Kind::SequenceMismatch => 'X',
            };
            format!("{}{c}", op.len())
        })
        .collect()
}

/// Point each mate at the other, as htsjdk's `SamPairUtil.setMateInfo` does
/// with mate CIGARs on.
fn set_mate_info(r1: &mut RecordBuf, r2: &mut RecordBuf, header: &sam::Header) {
    let _ = header;
    let u1 = r1.flags().is_unmapped();
    let u2 = r2.flags().is_unmapped();
    if u1 && u2 {
        let (rev1, rev2) = (
            r1.flags().is_reverse_complemented(),
            r2.flags().is_reverse_complemented(),
        );
        for (rec, mate_reverse) in [(&mut *r1, rev2), (&mut *r2, rev1)] {
            let mut flags = rec.flags() | Flags::MATE_UNMAPPED;
            flags.set(Flags::MATE_REVERSE_COMPLEMENTED, mate_reverse);
            *rec.flags_mut() = flags;
        }
        return;
    }
    if u1 || u2 {
        let (mapped, unmapped) = if u1 { (r2, r1) } else { (r1, r2) };
        *unmapped.reference_sequence_id_mut() = mapped.reference_sequence_id();
        *unmapped.alignment_start_mut() = mapped.alignment_start();
        *mapped.mate_reference_sequence_id_mut() = unmapped.reference_sequence_id();
        *mapped.mate_alignment_start_mut() = unmapped.alignment_start();
        let mut flags = mapped.flags() | Flags::MATE_UNMAPPED;
        flags.set(
            Flags::MATE_REVERSE_COMPLEMENTED,
            unmapped.flags().is_reverse_complemented(),
        );
        *mapped.flags_mut() = flags;
        *mapped.template_length_mut() = 0;
        *unmapped.mate_reference_sequence_id_mut() = mapped.reference_sequence_id();
        *unmapped.mate_alignment_start_mut() = mapped.alignment_start();
        let mut flags = unmapped.flags();
        flags.set(
            Flags::MATE_REVERSE_COMPLEMENTED,
            mapped.flags().is_reverse_complemented(),
        );
        *unmapped.flags_mut() = flags;
        *unmapped.template_length_mut() = 0;
        let mc = cigar_text(mapped);
        unmapped.data_mut().insert(Tag::MATE_CIGAR, Value::from(mc));
        return;
    }
    let tlen = insert_size(r1, r2);
    let proper = r1.reference_sequence_id() == r2.reference_sequence_id();
    let (mc1, mc2) = (cigar_text(r2), cigar_text(r1));
    let (ref1, ref2) = (r1.reference_sequence_id(), r2.reference_sequence_id());
    let (start1, start2) = (r1.alignment_start(), r2.alignment_start());
    let (rev1, rev2) = (
        r1.flags().is_reverse_complemented(),
        r2.flags().is_reverse_complemented(),
    );
    for (rec, mate_ref, mate_start, mate_rev, mc, tlen) in [
        (&mut *r1, ref2, start2, rev2, mc1, tlen),
        (&mut *r2, ref1, start1, rev1, mc2, -tlen),
    ] {
        *rec.mate_reference_sequence_id_mut() = mate_ref;
        *rec.mate_alignment_start_mut() = mate_start;
        let mut flags = rec.flags();
        flags.set(Flags::MATE_REVERSE_COMPLEMENTED, mate_rev);
        flags.set(Flags::PROPERLY_SEGMENTED, proper);
        *rec.flags_mut() = flags;
        *rec.template_length_mut() = tlen;
        rec.data_mut().insert(Tag::MATE_CIGAR, Value::from(mc));
    }
}

/// The 0-based query offset of the read base aligned to the 1-based `pos`, or
/// `None` when the read does not cover `pos` or has a deletion there.
fn offset_at(record: &RecordBuf, pos: usize) -> Option<usize> {
    let start = usize::from(record.alignment_start()?);
    if pos < start {
        return None;
    }
    let mut ref_pos = start;
    let mut query = 0usize;
    for op in record.cigar().as_ref() {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if pos < ref_pos + len {
                    return Some(query + (pos - ref_pos));
                }
                ref_pos += len;
                query += len;
            }
            Kind::Insertion | Kind::SoftClip => query += len,
            Kind::Deletion | Kind::Skip => {
                if pos < ref_pos + len {
                    return None;
                }
                ref_pos += len;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
    None
}

/// An in-memory pileup over built records, for tests only.
#[derive(Clone, Debug)]
pub struct TestPileup {
    header: sam::Header,
    records: Vec<RecordBuf>,
}

impl TestPileup {
    /// A pileup over the given records.
    pub fn new(header: sam::Header, records: Vec<RecordBuf>) -> Self {
        Self { header, records }
    }
}

impl PileupSource for TestPileup {
    type Record = RecordBuf;

    fn header(&self) -> &sam::Header {
        &self.header
    }

    fn pileup(&mut self, contig: &str, pos: Position) -> Result<Vec<ReadBase<'_, RecordBuf>>> {
        let id = self
            .header
            .reference_sequences()
            .get_index_of(contig.as_bytes())
            .ok_or_else(|| anyhow!("unknown contig: {contig}"))?;
        let pos = usize::from(pos);
        Ok(self
            .records
            .iter()
            .filter(|r| !r.flags().is_unmapped() && r.reference_sequence_id() == Some(id))
            .filter_map(|r| offset_at(r, pos).map(|offset| ReadBase::new(r, offset)))
            .collect())
    }
}

/// A genotype for [`VcfBuilder::add`], by alleles (`C/A`) or indices (`0/1`).
#[derive(Clone, Debug)]
pub struct Gt {
    pub sample: String,
    pub gt: String,
}

/// Shorthand for a [`Gt`].
pub fn gt(sample: &str, gt: &str) -> Gt {
    Gt {
        sample: sample.to_string(),
        gt: gt.to_string(),
    }
}

/// A variant for [`VcfBuilder::add`], with fgbio `VcfBuilder.add` defaults.
#[derive(Clone, Debug, Default)]
pub struct Variant {
    pub chrom: Option<String>,
    pub pos: usize,
    pub alleles: Vec<String>,
    pub info: Vec<(String, String)>,
    pub filters: Vec<String>,
    pub gts: Vec<Gt>,
}

impl Variant {
    /// A variant at `pos` with alleles and genotypes.
    pub fn new(pos: usize, alleles: &[&str], gts: Vec<Gt>) -> Self {
        Self {
            pos,
            alleles: alleles.iter().map(|a| a.to_string()).collect(),
            gts,
            ..Self::default()
        }
    }
}

/// Builds VCFs the way fgbio's `VcfBuilder` does, on its default header.
#[derive(Clone, Debug)]
pub struct VcfBuilder {
    samples: Vec<String>,
    variants: Vec<(usize, Variant)>,
}

impl VcfBuilder {
    /// A builder for these samples.
    pub fn new(samples: &[&str]) -> Self {
        Self {
            samples: samples.iter().map(|s| s.to_string()).collect(),
            variants: Vec::new(),
        }
    }

    /// Add a variant.
    pub fn add(&mut self, variant: Variant) -> &mut Self {
        let chrom = variant.chrom.clone().unwrap_or_else(|| "chr1".into());
        let index = default_contigs().iter().position(|c| *c == chrom).unwrap();
        self.variants.push((index, variant));
        self
    }

    fn header_text(&self) -> String {
        let mut text = String::from("##fileformat=VCFv4.2\n");
        for contig in default_contigs() {
            text.push_str(&format!("##contig=<ID={contig},length=200000000>\n"));
        }
        text.push_str(concat!(
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate allele counts in genotypes.\">\n",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth across all samples.\">\n",
            "##FILTER=<ID=LowQD,Description=\"Low Quality/Depth value\">\n",
            "##FILTER=<ID=LowAB,Description=\"Low/poor allele balance.\">\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype string.\">\n",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depth per allele.\">\n",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality.\">\n",
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred scaled genotype likelihoods.\">\n",
            "##ALT=<ID=NON_REF,Description=\"Represents any non-reference allele.\">\n",
        ));
        text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if !self.samples.is_empty() {
            text.push_str("\tFORMAT");
            for sample in &self.samples {
                text.push('\t');
                text.push_str(sample);
            }
        }
        text.push('\n');
        text
    }

    fn record_text(&self, variant: &Variant) -> String {
        let chrom = variant.chrom.clone().unwrap_or_else(|| "chr1".into());
        let alleles = &variant.alleles;
        let alt = if alleles.len() > 1 {
            alleles[1..].join(",")
        } else {
            ".".to_string()
        };
        let filter = if variant.filters.is_empty() {
            ".".to_string()
        } else {
            variant.filters.join(";")
        };
        let info = if variant.info.is_empty() {
            ".".to_string()
        } else {
            variant
                .info
                .iter()
                .map(|(k, v)| format!("{k}={v}"))
                .collect::<Vec<_>>()
                .join(";")
        };
        let mut line = format!(
            "{chrom}\t{}\t.\t{}\t{alt}\t.\t{filter}\t{info}",
            variant.pos, alleles[0]
        );
        if !self.samples.is_empty() {
            line.push_str("\tGT");
            for sample in &self.samples {
                let call = variant
                    .gts
                    .iter()
                    .find(|g| &g.sample == sample)
                    .map(|g| {
                        let sep = if g.gt.contains('|') { '|' } else { '/' };
                        g.gt.split(sep)
                            .map(|a| {
                                if a == "." || a.chars().all(|c| c.is_ascii_digit()) {
                                    a.to_string()
                                } else {
                                    alleles.iter().position(|x| x == a).unwrap().to_string()
                                }
                            })
                            .collect::<Vec<_>>()
                            .join(&sep.to_string())
                    })
                    .unwrap_or_else(|| "./.".to_string());
                line.push('\t');
                line.push_str(&call);
            }
        }
        line.push('\n');
        line
    }

    /// Write the VCF, records in coordinate order.
    pub fn write(&self, path: &Path) -> PathBuf {
        let mut sorted: BTreeMap<(usize, usize), &Variant> = BTreeMap::new();
        for (index, variant) in &self.variants {
            sorted.insert((*index, variant.pos), variant);
        }
        self.write_records(path, sorted.into_values())
    }

    /// Write the VCF, records in the order they were added.
    pub fn write_unsorted(&self, path: &Path) -> PathBuf {
        self.write_records(path, self.variants.iter().map(|(_, v)| v))
    }

    fn write_records<'a>(
        &self,
        path: &Path,
        variants: impl Iterator<Item = &'a Variant>,
    ) -> PathBuf {
        let mut file = std::fs::File::create(path).unwrap();
        file.write_all(self.header_text().as_bytes()).unwrap();
        for variant in variants {
            file.write_all(self.record_text(variant).as_bytes())
                .unwrap();
        }
        path.to_path_buf()
    }
}

/// Write a FASTA and its index for a single contig.
pub fn write_fasta(dir: &Path, name: &str, sequence: &str) -> PathBuf {
    let path = dir.join("ref.fa");
    let mut fasta = std::fs::File::create(&path).unwrap();
    writeln!(fasta, ">{name}").unwrap();
    let width = 60;
    for chunk in sequence.as_bytes().chunks(width) {
        fasta.write_all(chunk).unwrap();
        fasta.write_all(b"\n").unwrap();
    }
    let offset = name.len() + 2;
    let mut fai = std::fs::File::create(dir.join("ref.fa.fai")).unwrap();
    writeln!(
        fai,
        "{name}\t{}\t{offset}\t{width}\t{}",
        sequence.len(),
        width + 1
    )
    .unwrap();
    path
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_pair_sets_mate_information() {
        let mut builder = SamBuilder::new().read_length(40);
        let recs = builder.add_pair(Pair::at(100, 140));
        let (r1, r2) = (&recs[0], &recs[1]);
        assert_eq!(r1.template_length(), 80);
        assert_eq!(r2.template_length(), -80);
        assert!(r1.flags().is_mate_reverse_complemented());
        assert!(!r2.flags().is_mate_reverse_complemented());
        assert_eq!(r1.mate_alignment_start(), r2.alignment_start());
        assert_eq!(builder.records().len(), 2);
    }

    #[test]
    fn test_offset_at_walks_the_cigar() {
        let mut builder = SamBuilder::new().read_length(10);
        let recs = builder.add_frag(Frag {
            start: 100,
            cigar: Some("2S3M2D3M2I".into()),
            ..Frag::default()
        });
        let r = &recs[0];
        assert_eq!(offset_at(r, 99), None);
        assert_eq!(offset_at(r, 100), Some(2));
        assert_eq!(offset_at(r, 102), Some(4));
        assert_eq!(offset_at(r, 103), None);
        assert_eq!(offset_at(r, 105), Some(5));
        assert_eq!(offset_at(r, 108), None);
    }
}
