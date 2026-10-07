//! The README's examples, run on `tests/data` by the real binary.
//!
//! Every `console` block but the installation runs in order, under bash, in
//! one working directory holding a copy of `tests/data`, with `chaff` on the
//! `PATH`. The `text` block right after a `console` block is what that block
//! prints, compared line by line and field by field in both directions, and a
//! `console` block without one must print nothing.

use std::fs;
use std::path::Path;
use std::process::Command;

use tempfile::TempDir;

const README: &str = include_str!("../README.md");

/// The README's fenced blocks, as their language and their lines.
fn blocks() -> Vec<(&'static str, &'static str)> {
    let mut blocks = Vec::new();
    let mut rest = README;
    while let Some(start) = rest.find("```") {
        let (language, body) = rest[start + 3..].split_once('\n').unwrap();
        let end = body.find("```").unwrap();
        blocks.push((language, &body[..end]));
        rest = &body[end + 3..];
    }
    blocks
}

/// The whitespace-separated fields of each line that has any.
fn fields(text: &str) -> Vec<Vec<&str>> {
    text.lines()
        .map(|line| line.split_whitespace().collect::<Vec<_>>())
        .filter(|fields| !fields.is_empty())
        .collect()
}

/// A working directory holding a copy of `tests/data`.
fn work_dir() -> TempDir {
    let dir = TempDir::new().unwrap();
    let data = dir.path().join("tests/data");
    fs::create_dir_all(&data).unwrap();
    let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    for entry in fs::read_dir(source).unwrap() {
        let entry = entry.unwrap();
        fs::copy(entry.path(), data.join(entry.file_name())).unwrap();
    }
    dir
}

#[test]
fn test_every_readme_example_prints_what_the_readme_shows() {
    let dir = work_dir();
    let binary = Path::new(env!("CARGO_BIN_EXE_chaff")).parent().unwrap();
    let path = format!("{}:{}", binary.display(), std::env::var("PATH").unwrap());
    let blocks = blocks();
    let mut shown = 0;
    for (i, (language, body)) in blocks.iter().enumerate() {
        if *language != "console" || body.starts_with("cargo install") {
            continue;
        }
        let output = Command::new("bash")
            .args(["-e", "-c", &body.replace("| column -t", "")])
            .current_dir(dir.path())
            .env("PATH", &path)
            .env("NO_COLOR", "1")
            .output()
            .unwrap();
        let stdout = String::from_utf8(output.stdout).unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{body}\n{stderr}");
        match blocks.get(i + 1) {
            Some(("text", expected)) => {
                assert_eq!(fields(&stdout), fields(expected), "{body}");
                shown += 1;
            }
            _ => assert_eq!(stdout, "", "the README shows no output of:\n{body}"),
        }
    }
    assert!(shown >= 3, "only {shown} outputs were checked");
}
