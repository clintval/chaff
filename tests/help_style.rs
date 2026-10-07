//! The magenta one-line description at the top of `--help` and the license
//! footer at the bottom. These are styled on the `clap::Command` at startup
//! (`src/main.rs`), so the tests drive the real binary and inspect the raw ANSI
//! in its output.

use assert_cmd::Command;
use clap::builder::styling::{AnsiColor, Style};

/// Mirrors `TITLE` in `src/main.rs`.
const TITLE: Style = AnsiColor::Magenta.on_default();
/// Mirrors `FOOTER` in `src/main.rs`.
const FOOTER: Style = AnsiColor::Magenta.on_default();

/// The tool's one-line description; the first line of the help text.
const DESCRIPTION: &str = "Flag somatic variant calls that are library-preparation artifacts.";

/// Run the binary with `args` and return its captured stdout.
fn help_stdout(args: &[&str], no_color: bool) -> String {
    let mut cmd = Command::cargo_bin("chaff").unwrap();
    if no_color {
        cmd.env("NO_COLOR", "1");
    } else {
        cmd.env_remove("NO_COLOR");
    }
    let output = cmd
        .args(args)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    String::from_utf8(output).unwrap()
}

#[test]
fn long_help_opens_with_the_magenta_description_line() {
    let stdout = help_stdout(&["--help"], false);
    let title = format!("{}{DESCRIPTION}{}\n", TITLE.render(), TITLE.render_reset());
    assert!(stdout.starts_with(&title), "got:\n{stdout}");
}

#[test]
fn help_ends_with_the_colored_license_footer() {
    let stdout = help_stdout(&["--help"], false);
    let footer = format!(
        "{}MIT License 2026 · Clint Valentine{}",
        FOOTER.render(),
        FOOTER.render_reset()
    );
    assert!(stdout.trim_end().ends_with(&footer), "got:\n{stdout}");
}

#[test]
fn no_color_strips_every_escape() {
    let stdout = help_stdout(&["--help"], true);
    assert!(!stdout.contains('\u{1b}'), "got:\n{stdout}");
    assert!(stdout.starts_with(DESCRIPTION));
}

#[test]
fn backtick_terms_are_stripped_of_their_backticks() {
    assert!(!help_stdout(&["--help"], true).contains('`'));
}
