//! VCF and BCF reading and writing, header lines, and value formatting.

use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use noodles::bcf;
use noodles::bgzf;
use noodles::vcf;
use noodles::vcf::header::record::value::map::info::{Number, Type};
use noodles::vcf::header::record::value::map::{Filter, Info, Map};
use noodles::vcf::variant::io::Write as _;
use noodles::vcf::variant::RecordBuf;

/// Whether a path names a BCF file.
fn is_bcf(path: &Path) -> bool {
    path.extension().is_some_and(|ext| ext == "bcf")
}

/// A reader over VCF (plain or BGZF) or BCF records.
pub enum VariantReader {
    /// A VCF reader.
    Vcf(vcf::io::Reader<Box<dyn BufRead>>),
    /// A BCF reader.
    Bcf(bcf::io::Reader<bgzf::io::Reader<File>>),
}

impl VariantReader {
    /// Open a VCF or BCF by its extension.
    pub fn open(path: &Path) -> Result<Self> {
        if is_bcf(path) {
            let file = File::open(path).with_context(|| format!("failed to open BCF: {path:?}"))?;
            Ok(Self::Bcf(bcf::io::Reader::new(file)))
        } else {
            let reader = vcf::io::reader::Builder::default()
                .build_from_path(path)
                .with_context(|| format!("failed to open VCF: {path:?}"))?;
            Ok(Self::Vcf(reader))
        }
    }

    /// Read the header.
    pub fn read_header(&mut self) -> io::Result<vcf::Header> {
        match self {
            Self::Vcf(r) => r.read_header(),
            Self::Bcf(r) => r.read_header(),
        }
    }

    /// Read the next record into `record`, returning zero at the end.
    pub fn read_record(
        &mut self,
        header: &vcf::Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        match self {
            Self::Vcf(r) => r.read_record_buf(header, record),
            Self::Bcf(r) => r.read_record_buf(header, record),
        }
    }
}

/// A writer of VCF (plain or BGZF) or BCF records.
pub enum VariantWriter {
    /// A VCF writer.
    Vcf(vcf::io::Writer<Box<dyn Write>>),
    /// A BCF writer.
    Bcf(Box<bcf::io::Writer<bgzf::io::Writer<BufWriter<File>>>>),
}

impl VariantWriter {
    /// Create a VCF or BCF by its extension; `-` writes VCF to standard output.
    pub fn create(path: &Path) -> Result<Self> {
        if path == Path::new("-") {
            let stdout: Box<dyn Write> = Box::new(BufWriter::new(io::stdout()));
            return Ok(Self::Vcf(vcf::io::Writer::new(stdout)));
        }
        if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("failed to create output directory: {parent:?}"))?;
        }
        if is_bcf(path) {
            let file =
                File::create(path).with_context(|| format!("failed to create BCF: {path:?}"))?;
            let writer = bcf::io::Writer::new(BufWriter::new(file));
            Ok(Self::Bcf(Box::new(writer)))
        } else {
            let writer = vcf::io::writer::Builder::default()
                .build_from_path(path)
                .with_context(|| format!("failed to create VCF: {path:?}"))?;
            Ok(Self::Vcf(writer))
        }
    }

    /// Write the header.
    pub fn write_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        match self {
            Self::Vcf(w) => w.write_header(header),
            Self::Bcf(w) => w.write_header(header),
        }
    }

    /// Write one record.
    pub fn write_record(&mut self, header: &vcf::Header, record: &RecordBuf) -> io::Result<()> {
        match self {
            Self::Vcf(w) => w.write_variant_record(header, record),
            Self::Bcf(w) => w.write_variant_record(header, record),
        }
    }

    /// Flush and close the stream, finishing any BGZF blocks.
    pub fn finish(self) -> io::Result<()> {
        match self {
            Self::Vcf(mut w) => {
                w.get_mut().flush()?;
                drop(w);
                Ok(())
            }
            Self::Bcf(w) => {
                let mut inner = w.into_inner().finish()?;
                inner.flush()
            }
        }
    }
}

/// Add an INFO header line, replacing any with the same ID.
pub fn add_info(header: &mut vcf::Header, id: &str, number: Number, ty: Type, description: &str) {
    header
        .infos_mut()
        .insert(id.to_string(), Map::<Info>::new(number, ty, description));
}

/// Add a FILTER header line, replacing any with the same ID.
pub fn add_filter(header: &mut vcf::Header, id: &str, description: &str) {
    header
        .filters_mut()
        .insert(id.to_string(), Map::<Filter>::new(description));
}

/// Round a value the way htsjdk's `VCFEncoder.formatVCFDouble` prints it (two
/// decimals at or above one, three decimals down to 0.01, four significant
/// digits below that, and zero under 1e-20), so a value read back matches
/// fgbio's output.
pub fn vcf_float(value: f64) -> f32 {
    if !value.is_finite() {
        return value as f32;
    }
    let rounded = if value >= 1.0 {
        (value * 100.0).round() / 100.0
    } else if value >= 0.01 {
        (value * 1000.0).round() / 1000.0
    } else if value.abs() >= 1e-20 {
        let exponent = value.abs().log10().floor();
        let scale = 10f64.powf(3.0 - exponent);
        (value * scale).round() / scale
    } else {
        0.0
    };
    rounded as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vcf_float_rounds_like_htsjdk() {
        assert_eq!(vcf_float(0.99999), 1.0);
        assert_eq!(vcf_float(1.234567), 1.23);
        assert_eq!(vcf_float(0.123456), 0.123);
        assert_eq!(vcf_float(0.0012345678), 0.001_235);
        assert_eq!(vcf_float(1.5e-25), 0.0);
        assert_eq!(vcf_float(-12.3456), -12.35);
        assert_eq!(vcf_float(-0.0012345), -0.001_235);
    }

    #[test]
    fn test_header_lines_are_added() {
        let mut header = vcf::Header::default();
        add_info(&mut header, "XX", Number::Count(1), Type::Float, "a test");
        add_filter(&mut header, "Bad", "a filter");
        assert!(header.infos().contains_key("XX"));
        assert!(header.filters().contains_key("Bad"));
    }
}
