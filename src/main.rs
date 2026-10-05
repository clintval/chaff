//! Separate somatic variant calls from library-preparation damage artifacts.
use std::process;

use anyhow::{Error, Result};
use clap::builder::styling::{AnsiColor, Effects, Style, Styles};
use clap::{CommandFactory, FromArgMatches, Parser};
use env_logger::Env;
use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

pub(crate) const HEADER: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
pub(crate) const USAGE: Style = AnsiColor::Green.on_default().effects(Effects::BOLD);
pub(crate) const LITERAL: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
pub(crate) const PLACEHOLDER: Style = AnsiColor::Cyan.on_default();
pub(crate) const ERROR: Style = AnsiColor::Red.on_default().effects(Effects::BOLD);
pub(crate) const VALID: Style = AnsiColor::Cyan.on_default().effects(Effects::BOLD);
pub(crate) const INVALID: Style = AnsiColor::Yellow.on_default().effects(Effects::BOLD);

/// The tool's one-line description at the top of the help: magenta, the same
/// color code as the [`FOOTER`], so the help opens and closes in one color.
pub(crate) const TITLE: Style = AnsiColor::Magenta.on_default();
/// Indented code examples and backtick-wrapped terms in help text: tertiary
/// color (yellow), a named ANSI slot that follows the user's terminal palette.
pub(crate) const CODE: Style = AnsiColor::Yellow.on_default();
/// The license/attribution footer at the bottom of the help: magenta, shared
/// with the [`TITLE`].
pub(crate) const FOOTER: Style = AnsiColor::Magenta.on_default();
/// The right-hand description column of an indented two-column table: a faded
/// gray (bright-black), so the left-hand code term (in [`CODE`]) stands out.
pub(crate) const FADED: Style = AnsiColor::BrightBlack.on_default();

/// Cargo's color style.
/// [source](https://github.com/crate-ci/clap-cargo/blob/master/src/style.rs)
pub(crate) const CARGO_STYLING: Styles = Styles::styled()
    .header(HEADER)
    .usage(USAGE)
    .literal(LITERAL)
    .placeholder(PLACEHOLDER)
    .error(ERROR)
    .valid(VALID)
    .invalid(INVALID);

/// Separate somatic variant calls from library-preparation damage artifacts.
///
/// Library preparation can turn DNA damage into a base change that both
/// strands of a duplex agree on. This tool scores each somatic call by where
/// its alternate allele sits on the molecules that carry it, compared with the
/// molecules that carry the reference allele at the same site.
#[derive(Debug, Parser)]
#[command(author, version, color = clap::ColorChoice::Always, verbatim_doc_comment, arg_required_else_help = true)]
#[clap(styles = CARGO_STYLING)]
struct Cli {}

/// The ANSI escape that starts `style`, or an empty string when `color` is off
/// (honoring `NO_COLOR`). Its matching reset comes from [`esc_reset`].
fn esc(style: Style, color: bool) -> String {
    if color {
        style.render().to_string()
    } else {
        String::new()
    }
}

/// The reset escape for `style` (a full SGR reset), or empty when `color` is off.
fn esc_reset(style: Style, color: bool) -> String {
    if color {
        style.render_reset().to_string()
    } else {
        String::new()
    }
}

/// Paint the first line of `text` in the [`TITLE`] color (the tool's one-line
/// description), leaving the rest of the text untouched.
fn title_first_line(text: &str, color: bool) -> String {
    let title = esc(TITLE, color);
    let reset = esc_reset(TITLE, color);
    match text.split_once('\n') {
        Some((first, rest)) => format!("{title}{first}{reset}\n{rest}"),
        None => format!("{title}{text}{reset}"),
    }
}

/// Drop the backticks around terms (`` `like this` ``), painting each term with
/// the `code` escape. `after` is the escape that restores the surrounding text's
/// color once a term ends.
fn paint_backtick_terms(text: &str, code: &str, after: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(open) = rest.find('`') {
        out.push_str(&rest[..open]);
        let tail = &rest[open + 1..];
        match tail.find('`') {
            Some(close) => {
                out.push_str(code);
                out.push_str(&tail[..close]);
                out.push_str(after);
                rest = &tail[close + 1..];
            }
            None => {
                out.push('`');
                rest = tail;
                break;
            }
        }
    }
    out.push_str(rest);
    out
}

/// Split an indented two-column row into (left, gap, right) at the first run of
/// three or more spaces that follows the left-hand term.
fn split_two_column(content: &str) -> Option<(&str, &str, &str)> {
    let bytes = content.as_bytes();
    let mut seen_term = false;
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b' ' {
            let start = i;
            while i < bytes.len() && bytes[i] == b' ' {
                i += 1;
            }
            if seen_term && i - start >= 3 {
                return Some((&content[..start], &content[start..i], &content[i..]));
            }
        } else {
            seen_term = true;
            i += 1;
        }
    }
    None
}

/// Style one help string, dropping backticks throughout:
///
/// - indented two-column rows: the left code term in [`CODE`], the right
///   description column in [`FADED`];
/// - indented description continuations (deeper indent, no term): all [`FADED`];
/// - indented single-column examples (a bare command): all [`CODE`];
/// - prose: default color, with backtick terms painted in [`CODE`].
fn style_help_text(text: &str, color: bool) -> String {
    let code = esc(CODE, color);
    let faded = esc(FADED, color);
    let reset = esc_reset(CODE, color);
    text.split_inclusive('\n')
        .map(|line| {
            let (content, newline) = match line.strip_suffix('\n') {
                Some(content) => (content, "\n"),
                None => (line, ""),
            };
            if !content.starts_with("  ") || content.trim().is_empty() {
                return format!("{}{newline}", paint_backtick_terms(content, &code, &reset));
            }
            if let Some((left, gap, right)) = split_two_column(content) {
                let left = paint_backtick_terms(left, &code, &code);
                let right = paint_backtick_terms(right, &code, &faded);
                format!("{code}{left}{reset}{gap}{faded}{right}{reset}{newline}")
            } else if content.len() - content.trim_start().len() > 2 {
                let painted = paint_backtick_terms(content, &code, &faded);
                format!("{faded}{painted}{reset}{newline}")
            } else {
                let painted = paint_backtick_terms(content, &code, &code);
                format!("{code}{painted}{reset}{newline}")
            }
        })
        .collect()
}

/// Style a command's help: paint the first line of the about text in the
/// [`TITLE`] color, apply [`style_help_text`] to the abouts and every option's
/// help, and add the license footer. Subcommands are styled the same way.
fn decorate_help(cmd: clap::Command, color: bool) -> clap::Command {
    let about = cmd.get_about().map(ToString::to_string);
    let long_about = cmd
        .get_long_about()
        .map(ToString::to_string)
        .or_else(|| about.clone());

    let footer = format!(
        "{}MIT License 2026 · Clint Valentine{}",
        esc(FOOTER, color),
        esc_reset(FOOTER, color),
    );

    let choice = if color {
        clap::ColorChoice::Always
    } else {
        clap::ColorChoice::Never
    };
    let mut cmd = cmd.color(choice);
    if let Some(about) = about {
        cmd = cmd.about(title_first_line(&style_help_text(&about, color), color));
    }
    if let Some(long_about) = long_about {
        cmd = cmd.long_about(title_first_line(
            &style_help_text(&long_about, color),
            color,
        ));
    }
    cmd = cmd
        .next_line_help(true)
        .after_help(footer.clone())
        .after_long_help(footer);
    cmd.mut_args(|mut arg| {
        if let Some(help) = arg.get_help().map(ToString::to_string) {
            arg = arg.help(style_help_text(&help, color));
        }
        if let Some(long_help) = arg.get_long_help().map(ToString::to_string) {
            arg = arg.long_help(style_help_text(&long_help, color));
        }
        arg
    })
    .mut_subcommands(|sub| decorate_help(sub, color))
}

/// Main binary entrypoint.
#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), Error> {
    let color = std::env::var_os("NO_COLOR").is_none();

    let env = Env::default().default_filter_or("info");
    let write_style = if color {
        env_logger::WriteStyle::Auto
    } else {
        env_logger::WriteStyle::Never
    };
    env_logger::Builder::from_env(env)
        .write_style(write_style)
        .init();

    // Doc comments are hand-wrapped to 80 columns and `decorate_help` injects
    // ANSI color into them, which clap's wrapping would count as visible width.
    let cmd = decorate_help(Cli::command().term_width(usize::MAX), color);
    let matches = cmd.get_matches();
    let _cli = Cli::from_arg_matches(&matches).unwrap_or_else(|e| e.exit());
    process::exit(0);
}
