//! Cucumber behavior specs for chaff.
//!
//! Every feature file under `tests/features/` is an executable acceptance spec:
//! it drives the real `chaff` binary on small fixtures and asserts on its exit
//! code and stdout/stderr. The statistics and VCF logic are validated by the
//! inline unit tests in `src/lib` and the ported fgbio suite in
//! `tests/filter_somatic_vcf.rs`.
//!
//! Each scenario runs in its own temporary working directory: `Given a file
//! "..." containing:` writes a fixture there, `Given a BAM "..." sorted by
//! "..."` writes an empty BAM with that `@HD SO`, and the `I run` step executes
//! the binary with that directory as the cwd. A `\t` in a fixture or an
//! expected stdout is a tab.

use std::fs;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::Output;

use assert_cmd::Command;
use cucumber::gherkin::Step;
use cucumber::{given, then, when, World};
use noodles::bam;
use noodles::sam;
use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;
use noodles::sam::header::record::value::map::{Header, Map, ReferenceSequence};
use tempfile::TempDir;

/// Shared state across steps: the scenario's working directory and the most
/// recent invocation.
#[derive(Debug, Default, World)]
struct ChaffWorld {
    /// Per-scenario working directory, created lazily on first use.
    dir: Option<TempDir>,
    /// Arguments of the most recent invocation.
    args: Vec<String>,
    /// Exit code of the most recent invocation.
    code: Option<i32>,
    /// Captured stdout (lossy UTF-8).
    stdout: String,
    /// Captured stderr (lossy UTF-8).
    stderr: String,
}

impl ChaffWorld {
    /// The scenario's working directory, created on first use.
    fn work_dir(&mut self) -> PathBuf {
        self.dir
            .get_or_insert_with(|| TempDir::new().expect("create a temp working directory"))
            .path()
            .to_path_buf()
    }
}

fn run_chaff(dir: &Path, args: &[String]) -> Output {
    Command::cargo_bin("chaff")
        .expect("binary `chaff` builds")
        .current_dir(dir)
        .env("NO_COLOR", "1")
        .args(args)
        .output()
        .expect("chaff runs")
}

#[given(regex = r#"^a file "([^"]+)" containing:$"#)]
async fn given_file(world: &mut ChaffWorld, name: String, step: &Step) {
    let dir = world.work_dir();
    let content = format!("{}\n", step.docstring.clone().unwrap_or_default().trim());
    let content = content.replace("\\t", "\t");
    fs::write(dir.join(&name), content).expect("write fixture file");
}

#[given(regex = r#"^a BAM "([^"]+)" sorted by "([^"]+)"$"#)]
async fn given_bam(world: &mut ChaffWorld, name: String, sort_order: String) {
    let dir = world.work_dir();
    let hd = Map::<Header>::builder()
        .insert(SORT_ORDER, sort_order)
        .build()
        .expect("build @HD");
    let header = sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence(
            "chr1",
            Map::<ReferenceSequence>::new(NonZeroUsize::new(1_000).unwrap()),
        )
        .build();
    let file = fs::File::create(dir.join(&name)).expect("create BAM");
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).expect("write BAM header");
    writer.try_finish().expect("finish BAM");
}

#[when(regex = r"^I run `chaff ?(.*)`$")]
async fn i_run(world: &mut ChaffWorld, arg_line: String) {
    let dir = world.work_dir();
    let args: Vec<String> = arg_line.split_whitespace().map(str::to_owned).collect();
    let out = run_chaff(&dir, &args);
    world.code = out.status.code();
    world.stdout = String::from_utf8_lossy(&out.stdout).into_owned();
    world.stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    world.args = args;
}

#[then(regex = r"^the exit code is (\d+)$")]
async fn exit_code_is(world: &mut ChaffWorld, expected: i32) {
    assert_eq!(
        world.code,
        Some(expected),
        "args={:?}\nstderr:\n{}",
        world.args,
        world.stderr
    );
}

#[then(regex = r#"^stdout contains "(.*)"$"#)]
async fn stdout_contains(world: &mut ChaffWorld, needle: String) {
    let needle = needle.replace("\\t", "\t");
    assert!(
        world.stdout.contains(&needle),
        "stdout did not contain {needle:?}\nstdout:\n{}",
        world.stdout
    );
}

#[then(regex = r#"^stderr contains "(.*)"$"#)]
async fn stderr_contains(world: &mut ChaffWorld, needle: String) {
    assert!(
        world.stderr.contains(&needle),
        "stderr did not contain {needle:?}\nstderr:\n{}",
        world.stderr
    );
}

fn main() {
    futures::executor::block_on(ChaffWorld::cucumber().run_and_exit("tests/features"));
}
