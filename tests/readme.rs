//! The README's examples, run on `tests/data` by the real binary.

use std::fs;
use std::path::Path;

use assert_cmd::Command;
use tempfile::TempDir;

const README: &str = include_str!("../README.md");

fn chaff(args: &[&str]) {
    Command::cargo_bin("chaff")
        .unwrap()
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .env("NO_COLOR", "1")
        .args(args)
        .assert()
        .success();
}

/// Whether the README shows a line with these whitespace-separated fields.
fn shows(fields: &[&str]) -> bool {
    README
        .lines()
        .any(|line| line.split_whitespace().eq(fields.iter().copied()))
}

/// Columns 1, 2, 4, 5, 7, and 8 of each record, as the README cuts them.
fn records(path: &Path) -> Vec<Vec<String>> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            [0, 1, 3, 4, 6, 7].map(|i| fields[i].to_string()).to_vec()
        })
        .collect()
}

fn assert_shown(rows: &[Vec<String>]) {
    assert!(!rows.is_empty());
    for row in rows {
        let fields: Vec<&str> = row.iter().map(String::as_str).collect();
        assert!(
            shows(&fields),
            "the README does not show: {}",
            fields.join(" ")
        );
    }
}

#[test]
fn test_the_readme_shows_the_scored_calls_header_and_metrics() {
    let dir = TempDir::new().unwrap();
    let output = dir.path().join("calls.chaff.vcf");
    let metrics = dir.path().join("tumor.chaff.tsv");
    chaff(&[
        "--input",
        "tests/data/calls.vcf",
        "--bam",
        "tests/data/tumor.bam",
        "--ref",
        "tests/data/ref.fa",
        "--sample",
        "tumor",
        "--output",
        output.to_str().unwrap(),
        "--metrics",
        metrics.to_str().unwrap(),
        "--copied-damage-threshold",
        "0.05",
    ]);
    assert_shown(&records(&output));

    let filter = fs::read_to_string(&output)
        .unwrap()
        .lines()
        .find(|line| line.starts_with("##FILTER=<ID=CopiedDamage"))
        .unwrap()
        .to_string();
    assert!(README.lines().any(|line| line == filter), "{filter}");

    let rows: Vec<Vec<String>> = fs::read_to_string(&metrics)
        .unwrap()
        .lines()
        .take(3)
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            [1, 2, 3, 4, 5, 6, 13]
                .map(|i| fields[i].to_string())
                .to_vec()
        })
        .collect();
    assert_shown(&rows);
}

#[test]
fn test_the_readme_shows_the_fgbio_prior_calls() {
    let dir = TempDir::new().unwrap();
    let output = dir.path().join("fgbio.vcf");
    chaff(&[
        "--input",
        "tests/data/calls.vcf",
        "--bam",
        "tests/data/tumor.bam",
        "--sample",
        "tumor",
        "--output",
        output.to_str().unwrap(),
        "--filters",
        "end-repair-fill-in,a-tailing",
        "--end-repair-fill-in-threshold",
        "0.001",
        "--prior",
        "fgbio",
    ]);
    assert_shown(&records(&output));
}

#[test]
fn test_the_readme_choosing_filters_example_parses() {
    let section = README.split("## Choosing Filters").nth(1).unwrap();
    let block = section.split("```console\n").nth(1).unwrap();
    let command = block.split("```").next().unwrap().replace("\\\n", " ");
    let args: Vec<&str> = command.split_whitespace().skip(1).collect();
    assert!(args.contains(&"--copied-damage-classes"), "{args:?}");
    let dir = TempDir::new().unwrap();
    let output = Command::cargo_bin("chaff")
        .unwrap()
        .current_dir(dir.path())
        .env("NO_COLOR", "1")
        .args(&args)
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(stderr.contains("failed to open BAM"), "{stderr}");
}
