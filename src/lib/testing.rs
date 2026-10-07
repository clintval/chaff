//! Test builders that mirror fgbio's `VcfBuilder` on its default contigs, and
//! a reference FASTA writer. Builder fields carry fgbio's parameter names and
//! defaults; reads are built with `streampile::testing`.
#![allow(missing_docs)]

use std::collections::BTreeMap;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use crate::simplex::{negative_binomial, StratumProfile, MAX_CHANGES};

/// The contig names of fgbio's default sequence dictionary.
fn default_contigs() -> Vec<String> {
    (1..=22)
        .map(|i| i.to_string())
        .chain(["X", "Y", "M"].map(String::from))
        .map(|c| format!("chr{c}"))
        .collect()
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

/// A stratum whose duplex and single-strand changes fall as Poisson damage at
/// `rate` per molecule over `positions` positions of each depth in `depths`,
/// its duplex changes joined at each depth by `real` positions with 2.
pub fn poisson_stratum(depths: &[u32], positions: f64, rate: f64, real: u64) -> StratumProfile {
    let mut stratum = StratumProfile::default();
    for &n in depths {
        for k in 0..=MAX_CHANGES {
            let mean = f64::from(n) * rate;
            let count = (positions * negative_binomial(k, mean, f64::INFINITY)).round() as u64;
            let extra = if k == 2 { real } else { 0 };
            stratum.positions.insert((n, k), count + extra);
            stratum.strand_positions.insert((n, k), count);
        }
    }
    stratum
}
