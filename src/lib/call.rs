//! The genotype of the sample under test at one VCF record.

use noodles::vcf::variant::record_buf::samples::sample::Value;
use noodles::vcf::variant::RecordBuf;

/// The no-call allele.
pub const NO_CALL: &str = ".";

/// The allele of a deletion that spans the site, as a VCF writes it.
pub const SPANNING_DELETION: &str = "*";

/// Why the filters score no call at a record.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Skip {
    /// An allele of the call is more or less than one base.
    NotSnv,
    /// The genotype calls only the reference allele.
    HomozygousReference,
    /// The genotype calls one alternate allele on two or more copies.
    HomozygousAlternate,
    /// The genotype calls a single copy of an alternate allele.
    HaploidAlternate,
    /// No enabled filter scores the substitution.
    NoFilter,
}

impl Skip {
    /// The reason in words, as a log line gives it.
    pub fn reason(self) -> &'static str {
        match self {
            Skip::NotSnv => "not an SNV",
            Skip::HomozygousReference => "homozygous reference",
            Skip::HomozygousAlternate => "homozygous alternate",
            Skip::HaploidAlternate => "haploid alternate",
            Skip::NoFilter => "no enabled filter scores the substitution",
        }
    }
}

/// One sample's genotype as fgbio models it: the record's alleles (reference
/// first) and the sample's called alleles, `.` for a no-call.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Genotype {
    /// The record's alleles, reference first.
    pub alleles: Vec<String>,
    /// The called alleles, one per ploidy, `.` for a no-call.
    pub calls: Vec<String>,
}

impl Genotype {
    /// A genotype from allele and call strings.
    pub fn new(alleles: &[&str], calls: &[&str]) -> Self {
        Self {
            alleles: alleles.iter().map(|a| a.to_string()).collect(),
            calls: calls.iter().map(|c| c.to_string()).collect(),
        }
    }

    /// The genotype of the sample at `sample_index`, or `None` when the record
    /// has no `GT` for it.
    pub fn from_record(record: &RecordBuf, sample_index: usize) -> Option<Self> {
        let mut alleles = vec![record.reference_bases().to_string()];
        alleles.extend(record.alternate_bases().as_ref().iter().cloned());
        let sample = record.samples().get_index(sample_index)?;
        let Some(Some(Value::Genotype(genotype))) = sample.get("GT") else {
            return None;
        };
        let calls = genotype
            .as_ref()
            .iter()
            .map(|allele| match allele.position() {
                Some(i) => alleles.get(i).cloned().unwrap_or_else(|| NO_CALL.into()),
                None => NO_CALL.into(),
            })
            .collect();
        Some(Self { alleles, calls })
    }

    /// The genotype the filters score at a record, or why they score none.
    ///
    /// A heterozygous SNV is scored as called. A genotype that calls no
    /// allele, `.` or no `GT` at all, as many somatic callers write, is scored
    /// as heterozygous for the reference and the first alternate allele. A
    /// genotype calling only the reference, or only one alternate allele, is
    /// not: the filters compare alternate molecules with the reference
    /// molecules of the sample at the site.
    pub fn scored(record: &RecordBuf, sample_index: usize) -> Result<Self, Skip> {
        let missing = |alleles: Vec<String>| {
            let alt = alleles[1..]
                .iter()
                .find(|a| a.as_str() != SPANNING_DELETION)
                .cloned()
                .ok_or(Skip::NotSnv)?;
            Ok(Self {
                calls: vec![alleles[0].clone(), alt],
                alleles,
            })
        };
        let gt = match Self::from_record(record, sample_index) {
            Some(gt) if gt.called().next().is_some() => gt,
            Some(gt) => missing(gt.alleles)?,
            None => {
                let mut alleles = vec![record.reference_bases().to_string()];
                alleles.extend(record.alternate_bases().as_ref().iter().cloned());
                missing(alleles)?
            }
        };
        if !gt.is_het() {
            let reference = gt.reference().to_string();
            return Err(if gt.called().all(|c| c == reference) {
                Skip::HomozygousReference
            } else if gt.calls.len() == 1 {
                Skip::HaploidAlternate
            } else {
                Skip::HomozygousAlternate
            });
        }
        if !gt.calls_are_single_bases() {
            return Err(Skip::NotSnv);
        }
        Ok(gt)
    }

    /// The reference allele.
    pub fn reference(&self) -> &str {
        &self.alleles[0]
    }

    /// The called alleles, without no-calls or spanning deletions, as fgbio's
    /// `calledAlleles`.
    pub fn called(&self) -> impl Iterator<Item = &str> {
        self.calls
            .iter()
            .map(String::as_str)
            .filter(|c| *c != NO_CALL && *c != SPANNING_DELETION)
    }

    /// At least two called alleles differ.
    pub fn is_het(&self) -> bool {
        let mut called = self.called();
        match called.next() {
            Some(first) => called.any(|c| c != first),
            None => false,
        }
    }

    /// Heterozygous with no reference allele called.
    pub fn is_het_non_ref(&self) -> bool {
        self.is_het() && !self.calls.iter().any(|c| c == self.reference())
    }

    /// Every call, no-calls included, is one base long, as fgbio checks it.
    pub fn calls_are_single_bases(&self) -> bool {
        self.calls.iter().all(|c| c.len() == 1)
    }

    /// The first called allele that is not the reference.
    pub fn first_alt(&self) -> Option<&str> {
        self.called().find(|c| *c != self.reference())
    }

    /// The called alternate alleles.
    pub fn alts(&self) -> impl Iterator<Item = &str> {
        let reference = self.reference().to_string();
        self.called().filter(move |c| *c != reference)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_het_and_hom() {
        assert!(Genotype::new(&["A", "T"], &["A", "T"]).is_het());
        assert!(!Genotype::new(&["A"], &["A", "A"]).is_het());
        assert!(!Genotype::new(&["A", "T"], &["T", "T"]).is_het());
        assert!(!Genotype::new(&["A", "T"], &[".", "T"]).is_het());
        assert!(!Genotype::new(&["A", "T"], &[".", "."]).is_het());
        assert!(!Genotype::new(&["A", "T", "*"], &["A", "*"]).is_het());
        assert!(!Genotype::new(&["A", "T", "*"], &["T", "*"]).is_het());
        assert!(Genotype::new(&["A", "T", "*"], &["A", "T", "*"]).is_het());
        assert!(Genotype::new(&["A", "C", "T"], &["C", "T"]).is_het_non_ref());
        assert!(!Genotype::new(&["A", "C"], &["A", "C"]).is_het_non_ref());
    }

    /// Records of one sample, `tumor`, from VCF body lines.
    fn records(lines: &[&str]) -> Vec<RecordBuf> {
        use noodles::vcf;
        let text = format!(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n{}\n",
            lines.join("\n")
        );
        let mut reader = vcf::io::Reader::new(text.as_bytes());
        let header = reader.read_header().unwrap();
        reader
            .record_bufs(&header)
            .map(|record| record.unwrap())
            .collect()
    }

    /// Heterozygous SNVs are scored as called, and SNVs whose genotype calls
    /// nothing, or that have none, as heterozygous for the first alternate
    /// allele; homozygous, haploid alternate, and non-SNV calls are not.
    #[test]
    fn test_which_genotypes_are_scored() {
        let lines = [
            "chr1\t10\t.\tC\tT\t.\t.\t.\tGT\t0/1",
            "chr1\t11\t.\tC\tT\t.\t.\t.\tGT\t.",
            "chr1\t12\t.\tC\tT,A\t.\t.\t.\tGT\t./.",
            "chr1\t13\t.\tC\tT\t.\t.\t.\tAD\t5,2",
            "chr1\t14\t.\tC\tT\t.\t.\t.\tGT\t1/1",
            "chr1\t15\t.\tC\tT\t.\t.\t.\tGT\t0/0",
            "chr1\t16\t.\tC\tT\t.\t.\t.\tGT\t1",
            "chr1\t17\t.\tC\tT\t.\t.\t.\tGT\t0",
            "chr1\t18\t.\tCT\tC\t.\t.\t.\tGT\t0/1",
            "chr1\t19\t.\tCT\tC\t.\t.\t.\tGT\t.",
            "chr1\t20\t.\tC\t*\t.\t.\t.\tGT\t.",
        ];
        let scored: Vec<Result<Option<String>, Skip>> = records(&lines)
            .iter()
            .map(|r| Genotype::scored(r, 0).map(|gt| gt.first_alt().map(String::from)))
            .collect();
        let alt = |a: &str| Ok(Some(a.to_string()));
        assert_eq!(
            scored,
            [
                alt("T"),
                alt("T"),
                alt("T"),
                alt("T"),
                Err(Skip::HomozygousAlternate),
                Err(Skip::HomozygousReference),
                Err(Skip::HaploidAlternate),
                Err(Skip::HomozygousReference),
                Err(Skip::NotSnv),
                Err(Skip::NotSnv),
                Err(Skip::NotSnv),
            ]
        );
    }

    #[test]
    fn test_first_alt_and_single_bases() {
        let gt = Genotype::new(&["A", "G"], &["A", "G", "TT"]);
        assert_eq!(gt.first_alt(), Some("G"));
        assert!(!gt.calls_are_single_bases());
        assert!(Genotype::new(&["A", "G"], &["A", "."]).calls_are_single_bases());
    }
}
