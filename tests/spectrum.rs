//! The spectrum PDF, written by the real binary on `tests/data`.

use std::fs;
use std::path::Path;
use std::process::Command;

use tempfile::TempDir;

/// Each run is its own process, so a writer that follows a randomly seeded
/// hash map's order writes different bytes from one run to the next.
#[test]
fn test_the_spectrum_is_byte_identical_across_runs() {
    let dir = TempDir::new().unwrap();
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let spectra: Vec<Vec<u8>> = (0..8)
        .map(|run| {
            let spectrum = dir.path().join(format!("{run}.spectrum.pdf"));
            let output = Command::new(env!("CARGO_BIN_EXE_chaff"))
                .arg("--input")
                .arg(data.join("calls.vcf"))
                .arg("--output")
                .arg(dir.path().join(format!("{run}.vcf")))
                .arg("--bam")
                .arg(data.join("tumor.bam"))
                .arg("--ref")
                .arg(data.join("ref.fa"))
                .args(["--sample", "tumor", "--filters", "a-tailing"])
                .arg("--spectrum")
                .arg(&spectrum)
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            fs::read(spectrum).unwrap()
        })
        .collect();
    for (run, spectrum) in spectra.iter().enumerate().skip(1) {
        assert!(
            *spectrum == spectra[0],
            "run {run} wrote a spectrum that differs from run 0's"
        );
    }
}
