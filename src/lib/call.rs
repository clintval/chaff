//! The genotype of the sample under test at one VCF record.

use noodles::vcf::variant::record_buf::samples::sample::Value;
use noodles::vcf::variant::RecordBuf;

/// The no-call allele.
pub const NO_CALL: &str = ".";

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

    /// The reference allele.
    pub fn reference(&self) -> &str {
        &self.alleles[0]
    }

    /// The called alleles, without no-calls.
    pub fn called(&self) -> impl Iterator<Item = &str> {
        self.calls
            .iter()
            .map(String::as_str)
            .filter(|c| *c != NO_CALL)
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
        assert!(Genotype::new(&["A", "C", "T"], &["C", "T"]).is_het_non_ref());
        assert!(!Genotype::new(&["A", "C"], &["A", "C"]).is_het_non_ref());
    }

    #[test]
    fn test_first_alt_and_single_bases() {
        let gt = Genotype::new(&["A", "G"], &["A", "G", "TT"]);
        assert_eq!(gt.first_alt(), Some("G"));
        assert!(!gt.calls_are_single_bases());
        assert!(Genotype::new(&["A", "G"], &["A", "."]).calls_are_single_bases());
    }
}
