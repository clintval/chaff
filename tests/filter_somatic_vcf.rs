//! fgbio's `FilterSomaticVcfTest`, ported.
//!
//! Each test runs the whole command (VCF in, merge-join with molecules from a
//! coordinate-sorted pileup over built reads, priors, VCF out) and asserts what
//! fgbio's test asserts. fgbio runs every case under both of its BAM access
//! patterns; chaff only streams, so each case runs once. Cases run under
//! `--prior fgbio` assert fgbio's values exactly; the expected values were
//! produced by fgbio 4.1.1 on the same reads. Cases where the learned prior
//! changes the outcome on purpose say so in their names.

use std::path::{Path, PathBuf};

use chaff::filter::{run_filter_on, FilterArgs, FilterKind, FilterOptions};
use chaff::io::VariantReader;
use chaff::prior::PriorMode;
use chaff::read_end::{ATailing, EndRepairFillIn};
use chaff::template::ReadFilter;
use chaff::testing::{gt, Pair, SamBuilder, Variant, VcfBuilder};
use noodles::vcf::variant::record_buf::info::field::Value;
use noodles::vcf::variant::RecordBuf;
use tempfile::TempDir;

const RLEN: usize = 40;

fn read_vcf(path: &Path) -> (noodles::vcf::Header, Vec<RecordBuf>) {
    let mut reader = VariantReader::open(path).unwrap();
    let header = reader.read_header().unwrap();
    let mut records = Vec::new();
    let mut record = RecordBuf::default();
    while reader.read_record(&header, &mut record).unwrap() != 0 {
        records.push(record.clone());
    }
    (header, records)
}

fn float(record: &RecordBuf, key: &str) -> Option<f32> {
    match record.info().get(key) {
        Some(Some(Value::Float(f))) => Some(*f),
        _ => None,
    }
}

fn tumor_vcfs(dir: &Path) -> (PathBuf, PathBuf) {
    let build = |samples: &[&str], name: &str| {
        let normal = samples.contains(&"normal");
        let gts = |t: &str, n: &str| {
            let mut gts = vec![gt("tumor", t)];
            if normal {
                gts.push(gt("normal", n));
            }
            gts
        };
        let mut builder = VcfBuilder::new(samples);
        builder.add(Variant::new(100, &["C", "A"], gts("C/A", "C/C")));
        builder.add(Variant::new(200, &["G", "A"], gts("G/A", "G/G")));
        builder.add(Variant::new(300, &["AAA", "A"], gts("AAA/A", "AAA/AAA")));
        builder.add(Variant::new(400, &["A", "T"], gts("A/T", "A/A")));
        builder.add(Variant::new(500, &["C", "G"], gts("C/G", "C/C")));
        builder.write(&dir.join(name))
    };
    (
        build(&["tumor"], "tumor.vcf"),
        build(&["tumor", "normal"], "tumor_normal.vcf"),
    )
}

/// The reads of fgbio's shared BAM: artifact signal at 100, 400, and 500, and
/// a low-fraction G>A with alternate molecules spread evenly at 200.
fn tumor_bam() -> SamBuilder {
    let mut b = SamBuilder::new()
        .read_length(RLEN)
        .base_quality(40)
        .coordinate_sorted();
    let pair = |b: &mut SamBuilder, s1: usize, s2: usize, base: char| {
        b.add_pair(Pair::filled(s1, s2, base, RLEN));
    };
    for pos in 61..=100 {
        for _ in 1..=3 {
            pair(&mut b, pos, pos + RLEN, 'C');
            pair(&mut b, pos - RLEN, pos, 'C');
        }
    }
    pair(&mut b, 20, 62, 'A');
    pair(&mut b, 21, 61, 'A');
    pair(&mut b, 16, 63, 'A');
    add_even_g_to_a(&mut b);
    for pos in 361..=400 {
        for _ in 1..=3 {
            pair(&mut b, pos, pos + RLEN, 'A');
            pair(&mut b, pos - RLEN, pos, 'A');
        }
    }
    for (s1, s2) in [(398, 430), (398, 440), (399, 437), (400, 449), (400, 441)] {
        pair(&mut b, s1, s2, 'T');
    }
    for pos in 461..=500 {
        for _ in 1..=3 {
            pair(&mut b, pos, pos + RLEN, 'C');
            pair(&mut b, pos - RLEN, pos, 'C');
        }
    }
    for (s1, s2) in [(498, 530), (498, 540), (499, 537), (500, 549), (500, 541)] {
        pair(&mut b, s1, s2, 'G');
    }
    b
}

/// A G>A at 200 with a low alternate fraction but alternate molecules spread
/// evenly through the reads.
fn add_even_g_to_a(b: &mut SamBuilder) {
    for pos in 161..=200 {
        for i in 1..=3 {
            b.add_pair(Pair::filled(pos, pos + RLEN, 'G', RLEN));
            b.add_pair(Pair::filled(pos - RLEN, pos, 'G', RLEN));
            if i == 1 && pos % 4 == 0 {
                b.add_pair(Pair::filled(pos, pos + RLEN, 'A', RLEN));
                b.add_pair(Pair::filled(pos - RLEN, pos, 'A', RLEN));
            }
        }
    }
}

fn fgbio_options() -> FilterOptions {
    FilterOptions {
        filters: vec![FilterKind::EndRepairFillIn, FilterKind::ATailing],
        ..FilterOptions::default()
    }
}

fn run(
    dir: &TempDir,
    input: &Path,
    reads: &SamBuilder,
    options: FilterOptions,
) -> anyhow::Result<Vec<RecordBuf>> {
    let output = dir.path().join("filtered.vcf");
    let args = FilterArgs {
        input: input.to_path_buf(),
        output: output.clone(),
        bam: PathBuf::from("reads.bam"),
        reference: None,
        metrics: None,
        read_filter: ReadFilter::default(),
        options,
    };
    run_filter_on(&args, reads.pileup())?;
    Ok(read_vcf(&output).1)
}

fn options(sample: Option<&str>, prior: PriorMode) -> FilterOptions {
    FilterOptions {
        sample: sample.map(String::from),
        prior,
        ..fgbio_options()
    }
}

fn has_filter(record: &RecordBuf, name: &str) -> bool {
    record.filters().as_ref().contains(name)
}

/// The annotations fgbio's tests check on the five calls of the shared VCF:
/// ATAP on the SNVs with an A or T alternate allele, ERFAP on every SNV, and
/// nothing on the deletion.
fn assert_annotated_as_fgbio(records: &[RecordBuf]) {
    assert_eq!(records.len(), 5);
    let atap: Vec<bool> = records
        .iter()
        .map(|r| float(r, ATailing::INFO).is_some())
        .collect();
    let erfap: Vec<bool> = records
        .iter()
        .map(|r| float(r, EndRepairFillIn::INFO).is_some())
        .collect();
    assert_eq!(atap, vec![true, true, false, true, false]);
    assert_eq!(erfap, vec![true, true, false, true, true]);
}

/// fgbio: "work on an empty VCF".
#[test]
fn test_work_on_an_empty_vcf() {
    let dir = TempDir::new().unwrap();
    let input = VcfBuilder::new(&["tumor"]).write(&dir.path().join("empty.vcf"));
    let records = run(&dir, &input, &tumor_bam(), fgbio_options()).unwrap();
    assert!(records.is_empty());
    let (header, _) = read_vcf(&dir.path().join("filtered.vcf"));
    assert!(header.infos().contains_key(ATailing::INFO));
    assert!(header.filters().contains_key(ATailing::FILTER));
    assert!(header.infos().contains_key(EndRepairFillIn::INFO));
    assert!(header.filters().contains_key(EndRepairFillIn::FILTER));
}

/// fgbio: "raise an exception when stream BAM is true and variant records are
/// not coordinate ordered".
#[test]
fn test_raise_an_error_when_variant_records_are_not_coordinate_ordered() {
    let dir = TempDir::new().unwrap();
    let mut vcf = VcfBuilder::new(&["sample1"]);
    vcf.add(Variant::new(2, &["C", "A"], vec![gt("sample1", "C/A")]));
    vcf.add(Variant::new(1, &["C", "A"], vec![gt("sample1", "C/A")]));
    let input = vcf.write_unsorted(&dir.path().join("unsorted.vcf"));
    let mut reads = SamBuilder::new().coordinate_sorted();
    reads.add_frag(chaff::testing::Frag::at(1));
    let error = run(&dir, &input, &reads, fgbio_options()).unwrap_err();
    assert!(
        error.to_string().contains("not coordinate sorted"),
        "{error}"
    );
}

/// The single call of fgbio's INFO and FILTER preservation tests: a G>A at 200
/// with a low alternate fraction but alternate molecules spread evenly.
fn even_g_to_a(
    dir: &TempDir,
    info: Vec<(String, String)>,
    filters: Vec<String>,
) -> (PathBuf, SamBuilder) {
    let mut vcf = VcfBuilder::new(&["sample1"]);
    vcf.add(Variant {
        info,
        filters,
        ..Variant::new(200, &["G", "A"], vec![gt("sample1", "G/A")])
    });
    let input = vcf.write(&dir.path().join("in.vcf"));
    let mut reads = SamBuilder::new().read_length(RLEN).coordinate_sorted();
    add_even_g_to_a(&mut reads);
    (input, reads)
}

/// fgbio: "not change existing INFO data when adding new INFO data".
#[test]
fn test_not_change_existing_info_data_when_adding_new_info_data() {
    for prior in [PriorMode::Fgbio, PriorMode::Learned] {
        let dir = TempDir::new().unwrap();
        let (input, reads) = even_g_to_a(&dir, vec![("DP".into(), "2".into())], vec![]);
        let records = run(&dir, &input, &reads, options(None, prior)).unwrap();
        assert_eq!(records.len(), 1);
        let annotated = &records[0];
        assert_eq!(usize::from(annotated.variant_start().unwrap()), 200);
        assert_eq!(annotated.info().get("DP"), Some(Some(&Value::Integer(2))));
        assert_eq!(float(annotated, ATailing::INFO), Some(1.0));
        assert_eq!(float(annotated, EndRepairFillIn::INFO), Some(1.0));
        assert!(annotated.filters().as_ref().is_empty());
    }
}

/// fgbio: "not change existing FILTER data when adding new FILTER data".
#[test]
fn test_not_change_existing_filter_data_when_adding_new_filter_data() {
    for prior in [PriorMode::Fgbio, PriorMode::Learned] {
        let dir = TempDir::new().unwrap();
        let (input, reads) = even_g_to_a(&dir, vec![], vec!["LowQD".into()]);
        let options = FilterOptions {
            a_tailing_threshold: Some(1.0),
            end_repair_fill_in_threshold: Some(1.0),
            ..options(None, prior)
        };
        let records = run(&dir, &input, &reads, options).unwrap();
        let annotated = &records[0];
        assert_eq!(float(annotated, ATailing::INFO), Some(1.0));
        assert_eq!(float(annotated, EndRepairFillIn::INFO), Some(1.0));
        let mut filters: Vec<&str> = annotated
            .filters()
            .as_ref()
            .iter()
            .map(String::as_str)
            .collect();
        filters.sort_unstable();
        assert_eq!(
            filters,
            vec![ATailing::FILTER, EndRepairFillIn::FILTER, "LowQD"]
        );
    }
}

/// fgbio: "not remove PASS from the FILTER field if no new filters are added".
#[test]
fn test_not_remove_pass_if_no_new_filters_are_added() {
    for prior in [PriorMode::Fgbio, PriorMode::Learned] {
        let dir = TempDir::new().unwrap();
        let (input, reads) = even_g_to_a(&dir, vec![], vec!["PASS".into()]);
        let records = run(&dir, &input, &reads, options(None, prior)).unwrap();
        let annotated = &records[0];
        assert_eq!(float(annotated, ATailing::INFO), Some(1.0));
        assert_eq!(float(annotated, EndRepairFillIn::INFO), Some(1.0));
        assert!(annotated.filters().is_pass());
    }
}

/// fgbio: "remove PASS from the FILTER field if a new filter is added".
#[test]
fn test_remove_pass_if_a_new_filter_is_added() {
    for prior in [PriorMode::Fgbio, PriorMode::Learned] {
        let dir = TempDir::new().unwrap();
        let (input, reads) = even_g_to_a(&dir, vec![], vec!["PASS".into()]);
        let options = FilterOptions {
            a_tailing_threshold: Some(1.0),
            end_repair_fill_in_threshold: Some(1.0),
            ..options(None, prior)
        };
        let records = run(&dir, &input, &reads, options).unwrap();
        let annotated = &records[0];
        let mut filters: Vec<&str> = annotated
            .filters()
            .as_ref()
            .iter()
            .map(String::as_str)
            .collect();
        filters.sort_unstable();
        assert_eq!(filters, vec![ATailing::FILTER, EndRepairFillIn::FILTER]);
    }
}

/// fgbio: "work on a single sample VCF when BAM access is: Streaming" (and
/// RandomAccess). The values are fgbio 4.1.1's.
#[test]
fn test_work_on_a_single_sample_vcf() {
    let dir = TempDir::new().unwrap();
    let (tumor, _) = tumor_vcfs(dir.path());
    let records = run(&dir, &tumor, &tumor_bam(), options(None, PriorMode::Fgbio)).unwrap();
    assert_annotated_as_fgbio(&records);
    assert!(!records.iter().any(|r| has_filter(r, ATailing::FILTER)));
    assert!(!records
        .iter()
        .any(|r| has_filter(r, EndRepairFillIn::FILTER)));
    let atap: Vec<Option<f32>> = records.iter().map(|r| float(r, ATailing::INFO)).collect();
    let erfap: Vec<Option<f32>> = records
        .iter()
        .map(|r| float(r, EndRepairFillIn::INFO))
        .collect();
    assert_eq!(
        atap,
        vec![Some(3.732e-3), Some(1.0), None, Some(0.715), None]
    );
    assert_eq!(
        erfap,
        vec![
            Some(3.218e-5),
            Some(1.0),
            None,
            Some(1.239e-5),
            Some(1.239e-5)
        ]
    );
}

/// fgbio: "fail on a single-sample VCF if an invalid sample name is provided".
#[test]
fn test_fail_on_a_single_sample_vcf_with_an_invalid_sample_name() {
    let dir = TempDir::new().unwrap();
    let (tumor, _) = tumor_vcfs(dir.path());
    let error = run(
        &dir,
        &tumor,
        &tumor_bam(),
        options(Some("WhoDis"), PriorMode::Fgbio),
    )
    .unwrap_err();
    assert!(error.to_string().contains("WhoDis"), "{error}");
}

/// fgbio: "fail on a multi-sample VCF if no sample name is provided".
#[test]
fn test_fail_on_a_multi_sample_vcf_without_a_sample_name() {
    let dir = TempDir::new().unwrap();
    let (_, tumor_normal) = tumor_vcfs(dir.path());
    let error = run(
        &dir,
        &tumor_normal,
        &tumor_bam(),
        options(None, PriorMode::Fgbio),
    )
    .unwrap_err();
    assert!(error.to_string().contains("--sample"), "{error}");
}

/// fgbio: "fail on a multi-sample VCF if an invalid sample name is provided".
#[test]
fn test_fail_on_a_multi_sample_vcf_with_an_invalid_sample_name() {
    let dir = TempDir::new().unwrap();
    let (_, tumor_normal) = tumor_vcfs(dir.path());
    let error = run(
        &dir,
        &tumor_normal,
        &tumor_bam(),
        options(Some("WhoDis"), PriorMode::Fgbio),
    )
    .unwrap_err();
    assert!(error.to_string().contains("WhoDis"), "{error}");
}

/// fgbio: "work on a multi-sample VCF if a sample name is given".
#[test]
fn test_work_on_a_multi_sample_vcf_with_a_sample_name() {
    let dir = TempDir::new().unwrap();
    let (_, tumor_normal) = tumor_vcfs(dir.path());
    let records = run(
        &dir,
        &tumor_normal,
        &tumor_bam(),
        options(Some("tumor"), PriorMode::Fgbio),
    )
    .unwrap();
    assert_annotated_as_fgbio(&records);
    assert!(!records.iter().any(|r| has_filter(r, ATailing::FILTER)));
    assert!(!records
        .iter()
        .any(|r| has_filter(r, EndRepairFillIn::FILTER)));
}

fn thresholded(prior: PriorMode) -> FilterOptions {
    FilterOptions {
        a_tailing: ATailing { distance: 4 },
        a_tailing_threshold: Some(0.001),
        end_repair_fill_in_threshold: Some(0.001),
        ..options(Some("tumor"), prior)
    }
}

/// fgbio: "apply filters if filter-specific p-value thresholds are supplied",
/// under fgbio's prior. The values are fgbio 4.1.1's.
#[test]
fn test_apply_filters_with_thresholds_under_the_fgbio_prior() {
    let dir = TempDir::new().unwrap();
    let (tumor, _) = tumor_vcfs(dir.path());
    let records = run(&dir, &tumor, &tumor_bam(), thresholded(PriorMode::Fgbio)).unwrap();
    assert_annotated_as_fgbio(&records);
    let atap: Vec<bool> = records
        .iter()
        .map(|r| has_filter(r, ATailing::FILTER))
        .collect();
    let erfap: Vec<bool> = records
        .iter()
        .map(|r| has_filter(r, EndRepairFillIn::FILTER))
        .collect();
    assert_eq!(atap, vec![true, false, false, true, false]);
    assert_eq!(erfap, vec![true, false, false, true, true]);
    let atap: Vec<Option<f32>> = records.iter().map(|r| float(r, ATailing::INFO)).collect();
    assert_eq!(
        atap,
        vec![Some(7.669e-8), Some(1.0), None, Some(5.265e-10), None]
    );
}

/// fgbio: "apply filters if filter-specific p-value thresholds are supplied".
/// Intended difference: under the learned prior, 3 to 5 alternate molecules
/// within 15 bp of an end, where 37.5% of the reference molecules also are,
/// leave ERFAP at 0.027 and 0.0037, above the 0.001 threshold; fgbio's
/// `(2 * maf)^2` prior alone drives them under it. A-tailing still filters 100
/// and 400, whose alternate molecules all sit where 5% of reference ones do.
#[test]
fn test_apply_filters_with_thresholds_under_the_learned_prior_spares_weak_end_repair_evidence() {
    let dir = TempDir::new().unwrap();
    let (tumor, _) = tumor_vcfs(dir.path());
    let records = run(&dir, &tumor, &tumor_bam(), thresholded(PriorMode::Learned)).unwrap();
    assert_annotated_as_fgbio(&records);
    let atap: Vec<bool> = records
        .iter()
        .map(|r| has_filter(r, ATailing::FILTER))
        .collect();
    let erfap: Vec<bool> = records
        .iter()
        .map(|r| has_filter(r, EndRepairFillIn::FILTER))
        .collect();
    assert_eq!(atap, vec![true, false, false, true, false]);
    assert_eq!(erfap, vec![false, false, false, false, false]);
    let erfap: Vec<Option<f32>> = records
        .iter()
        .map(|r| float(r, EndRepairFillIn::INFO))
        .collect();
    assert_eq!(
        erfap,
        vec![Some(0.027), Some(1.0), None, Some(3.718e-3), Some(3.718e-3)]
    );
}

/// The metrics rows of the thresholded run: one per filter and substitution
/// class, pooled over the calls.
#[test]
fn test_metrics_rows_of_the_shared_vcf() {
    let dir = TempDir::new().unwrap();
    let (tumor, _) = tumor_vcfs(dir.path());
    let output = dir.path().join("filtered.vcf");
    let metrics = dir.path().join("metrics.tsv");
    let args = FilterArgs {
        input: tumor,
        output,
        bam: PathBuf::from("reads.bam"),
        reference: None,
        metrics: Some(metrics.clone()),
        read_filter: ReadFilter::default(),
        options: thresholded(PriorMode::Learned),
    };
    run_filter_on(&args, tumor_bam().pileup()).unwrap();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&metrics)
        .unwrap();
    let rows: Vec<Vec<String>> = reader
        .records()
        .map(|r| r.unwrap().iter().map(String::from).collect())
        .collect();
    let keys: Vec<(String, String)> = rows.iter().map(|r| (r[1].clone(), r[2].clone())).collect();
    assert_eq!(
        keys,
        vec![
            ("a-tailing".into(), "C>A".into()),
            ("a-tailing".into(), "C>T".into()),
            ("a-tailing".into(), "T>A".into()),
            ("end-repair-fill-in".into(), "C>A".into()),
            ("end-repair-fill-in".into(), "C>G".into()),
            ("end-repair-fill-in".into(), "C>T".into()),
            ("end-repair-fill-in".into(), "T>A".into()),
        ]
    );
}

/// fgbio `PileupBuilderTest`: "raise an exception if the input SAM source is
/// not coordinate sorted when BAM access is `Streaming`".
#[test]
fn test_raise_an_error_if_the_reads_are_not_coordinate_sorted() {
    let dir = TempDir::new().unwrap();
    let (tumor, _) = tumor_vcfs(dir.path());
    let reads = SamBuilder::new();
    let error = run(
        &dir,
        &tumor,
        &reads,
        options(Some("tumor"), PriorMode::Fgbio),
    )
    .unwrap_err();
    assert!(error.to_string().contains("coordinate sorted"), "{error}");
}
