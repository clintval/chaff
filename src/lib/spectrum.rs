//! The trinucleotide spectrum of the SNVs a run scores, before and after its
//! filters, written as a PDF.
//!
//! Each SNV falls in one of 96 channels: its substitution read from the
//! pyrimidine of the base pair (C>A, C>G, C>T, T>A, T>C, or T>G) and the
//! bases 5' and 3' of it on that strand. The spectrum counts every
//! heterozygous SNV of the sample, whatever its FILTER, an SNV without a
//! genotype counting as heterozygous, before filtering. After filtering it
//! counts the SNVs that no filter flagged when any filter has a threshold,
//! and otherwise weighs each SNV by the product of the posteriors of the
//! enabled filters, its expected count as a real mutation.

use std::io::Write;
use std::path::Path;

use anyhow::{anyhow, Context as _, Result};
use kuva::plot::BarPlot;
use kuva::render::annotations::TextAnnotation;
use kuva::render::figure::Figure;
use kuva::render::layout::Layout;
use kuva::render::plots::Plot;
use kuva::render::render::Primitive;

use crate::classes::complement;
use crate::io::StagedFile;

/// The six substitution classes, read from the pyrimidine.
pub const CLASSES: [&str; 6] = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"];

/// The colours of the six classes, as mutational signature plots draw them.
pub const CLASS_COLORS: [&str; 6] = [
    "#1ebff0", "#050708", "#e62725", "#cbcacb", "#a1cf64", "#edc8c5",
];

/// The number of channels.
pub const CHANNELS: usize = 96;

const BASES: [u8; 4] = *b"ACGT";

/// The colour of the CpG C>T labels.
const CPG_COLOR: &str = "#12875c";

/// The channel of a forward-strand SNV and its flanking bases, read from the
/// pyrimidine: `16 class + 4 five_prime + three_prime`, with the classes in
/// [`CLASSES`] order and the flanking bases in `ACGT` order. `None` without
/// both flanking bases or for a base that is not one of `ACGT`.
pub fn channel(prev: Option<u8>, ref_base: u8, alt_base: u8, next: Option<u8>) -> Option<usize> {
    let upper = |b: u8| b.to_ascii_uppercase();
    let (mut five, mut reference, mut alternate, mut three) =
        (upper(prev?), upper(ref_base), upper(alt_base), upper(next?));
    if matches!(reference, b'G' | b'A') {
        (five, reference, alternate, three) = (
            complement(three),
            complement(reference),
            complement(alternate),
            complement(five),
        );
    }
    let index = |b: u8| BASES.iter().position(|&x| x == b);
    let class = CLASSES
        .iter()
        .position(|c| c.as_bytes() == [reference, b'>', alternate])?;
    Some(16 * class + 4 * index(five)? + index(three)?)
}

/// A channel's trinucleotide, its pyrimidine between its flanking bases, such
/// as `ACG`.
pub fn context(channel: usize) -> String {
    let class = CLASSES[channel / 16].as_bytes()[0];
    let bases = [BASES[(channel % 16) / 4], class, BASES[channel % 4]];
    String::from_utf8_lossy(&bases).into_owned()
}

/// Whether a channel is a C>T at a CpG.
pub fn is_cpg_c_to_t(channel: usize) -> bool {
    channel / 16 == 2 && channel % 4 == 2
}

/// What the spectrum after filtering counts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum After {
    /// The SNVs that no filter flagged, when any filter has a threshold.
    Passing,
    /// Each SNV weighed by the product of its posteriors, without thresholds.
    Weighted,
}

/// The spectrum of one sample's scored SNVs, before and after filtering.
#[derive(Clone, Debug, PartialEq)]
pub struct Spectrum {
    /// The scored SNVs in each channel.
    pub before: [f64; CHANNELS],
    /// The SNVs in each channel after filtering, as [`Spectrum::after_kind`]
    /// counts them.
    pub after: [f64; CHANNELS],
    /// What [`Spectrum::after`] counts.
    pub after_kind: After,
}

impl Spectrum {
    /// An empty spectrum.
    pub fn new(after_kind: After) -> Self {
        Self {
            before: [0.0; CHANNELS],
            after: [0.0; CHANNELS],
            after_kind,
        }
    }

    /// Count one scored SNV in `channel`, flagged by some filter or not, with
    /// the product of its posteriors.
    pub fn add(&mut self, channel: usize, flagged: bool, weight: f64) {
        self.before[channel] += 1.0;
        self.after[channel] += match self.after_kind {
            After::Passing => f64::from(u8::from(!flagged)),
            After::Weighted => weight,
        };
    }

    /// The figure: two panels of 96 bars on one axis, each under a strip of
    /// the class colours, with the CpG C>T contexts named after filtering.
    fn figure(&self, sample: &str) -> Figure {
        let top = self.before.iter().copied().fold(1.0, f64::max);
        let after_title = match self.after_kind {
            After::Passing => "After chaff: Calls Passing Every Threshold",
            After::Weighted => "After chaff: Expected Real Calls, Each Weighted by Its Posteriors",
        };
        let (before, before_layout) =
            panel(&self.before, "Before chaff: Every Heterozygous SNV", top);
        let (after, mut after_layout) = panel(&self.after, after_title, top);
        for i in (0..CHANNELS).filter(|&i| is_cpg_c_to_t(i)) {
            let label = TextAnnotation::new(context(i), i as f64 + 1.0, self.after[i] + top * 0.03)
                .with_color(CPG_COLOR)
                .with_font_size(10);
            after_layout = after_layout.with_annotation(label);
        }
        Figure::new(2, 1)
            .with_title(title(sample))
            .with_plots(vec![before, after])
            .with_layouts(vec![before_layout, after_layout])
            .with_cell_size(1100.0, 360.0)
    }

    /// Write the spectrum of `sample` as a PDF to `path`, its title bold.
    pub fn write_pdf(&self, path: &Path, sample: &str) -> Result<()> {
        let title = title(sample);
        let mut scene = self.figure(sample).render();
        for element in &mut scene.elements {
            if let Primitive::Text { content, bold, .. } = element {
                if *content == title {
                    *bold = true;
                }
            }
        }
        let bytes = kuva::backend::pdf::PdfBackend
            .render_scene(&scene)
            .map_err(|e| anyhow!("failed to render the spectrum for {path:?}: {e}"))?;
        let write = || -> Result<()> {
            let staged = StagedFile::create(path)?;
            staged.writer()?.write_all(&bytes)?;
            staged.persist()
        };
        write().with_context(|| format!("failed to write the spectrum: {path:?}"))
    }
}

/// The figure's title.
fn title(sample: &str) -> String {
    format!("Trinucleotide Spectrum of {sample}")
}

/// A round step that splits `top` into about four ticks.
fn tick_step(top: f64) -> f64 {
    let raw = top / 4.0;
    let magnitude = 10f64.powf(raw.log10().floor());
    [1.0, 2.0, 2.5, 5.0, 10.0]
        .into_iter()
        .map(|m| m * magnitude)
        .find(|step| *step >= raw)
        .unwrap_or(10.0 * magnitude)
}

/// One panel: a bar per channel in its class's colour, under the class names
/// and a strip of their colours. The axis ticks stop at the first tick at or
/// above `top`, the tallest bar before filtering, and the strip and names sit
/// between it and the next tick, which the axis never reaches. The strip is a
/// stacked bar per channel on an unpainted base.
fn panel(values: &[f64; CHANNELS], title: &str, top: f64) -> (Vec<Plot>, Layout) {
    let step = tick_step(top);
    let ceiling = (top / step).ceil() * step;
    let mut bars = BarPlot::new();
    let mut strip = BarPlot::new().with_width(1.0).with_stacked();
    for (i, value) in values.iter().enumerate() {
        let color = CLASS_COLORS[i / 16];
        bars = bars.with_colored_bar(context(i), *value, color);
        strip = strip.with_group(
            context(i),
            [(ceiling + step * 0.25, "none"), (step * 0.15, color)],
        );
    }
    let plots = vec![Plot::Bar(strip), Plot::Bar(bars)];
    let mut layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_y_label("SNVs")
        .with_y_axis_min(0.0)
        .with_y_axis_max(ceiling + step * 0.95)
        .with_y_tick_step(step)
        .with_show_grid(false)
        .with_x_tick_rotate(90.0)
        .with_tick_size(8);
    for (class, name) in CLASSES.iter().enumerate() {
        let label = TextAnnotation::new(*name, 16.0 * class as f64 + 8.5, ceiling + step * 0.7)
            .with_color("black")
            .with_font_size(11);
        layout = layout.with_annotation(label);
    }
    (plots, layout)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_channel_reads_the_pyrimidine_strand() {
        assert_eq!(channel(Some(b'A'), b'C', b'A', Some(b'A')), Some(0));
        assert_eq!(
            channel(Some(b'A'), b'C', b'T', Some(b'G')),
            Some(16 * 2 + 2)
        );
        assert_eq!(
            channel(Some(b'C'), b'G', b'A', Some(b'T')),
            Some(16 * 2 + 2)
        );
        assert_eq!(
            channel(Some(b't'), b't', b'g', Some(b't')),
            Some(16 * 5 + 15)
        );
        assert_eq!(
            channel(Some(b'A'), b'A', b'C', Some(b'A')),
            Some(16 * 5 + 15)
        );
        assert_eq!(channel(None, b'C', b'T', Some(b'G')), None);
        assert_eq!(channel(Some(b'N'), b'C', b'T', Some(b'G')), None);
        assert_eq!(channel(Some(b'A'), b'C', b'C', Some(b'G')), None);
    }

    #[test]
    fn test_context_and_cpg_name_every_channel() {
        assert_eq!(context(0), "ACA");
        assert_eq!(context(16 * 2 + 2), "ACG");
        assert_eq!(context(95), "TTT");
        let cpg: Vec<String> = (0..CHANNELS)
            .filter(|&c| is_cpg_c_to_t(c))
            .map(context)
            .collect();
        assert_eq!(cpg, ["ACG", "CCG", "GCG", "TCG"]);
        for c in 0..CHANNELS {
            let bases = context(c).into_bytes();
            let class = CLASSES[c / 16].as_bytes();
            assert_eq!(
                channel(Some(bases[0]), class[0], class[2], Some(bases[2])),
                Some(c)
            );
        }
    }

    #[test]
    fn test_spectrum_counts_before_passing_and_weighted() {
        let acg = 16 * 2 + 2;
        let mut passing = Spectrum::new(After::Passing);
        passing.add(acg, true, 0.01);
        passing.add(acg, false, 0.9);
        passing.add(0, false, 1.0);
        assert_eq!((passing.before[acg], passing.after[acg]), (2.0, 1.0));
        assert_eq!((passing.before[0], passing.after[0]), (1.0, 1.0));

        let mut weighted = Spectrum::new(After::Weighted);
        weighted.add(acg, true, 0.01);
        weighted.add(acg, false, 0.9);
        assert_eq!(weighted.before[acg], 2.0);
        assert!((weighted.after[acg] - 0.91).abs() < 1e-12);
        assert_eq!(weighted.before.iter().sum::<f64>(), 2.0);
    }

    #[test]
    fn test_write_pdf_writes_a_pdf() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("nested/tumor.spectrum.pdf");
        let mut spectrum = Spectrum::new(After::Weighted);
        for c in 0..CHANNELS {
            spectrum.add(c, c % 3 == 0, 0.5);
        }
        spectrum.write_pdf(&path, "tumor").unwrap();
        let bytes = std::fs::read(&path).unwrap();
        assert!(bytes.starts_with(b"%PDF-"), "{:?}", &bytes[..8]);
        let blocked = path.join("tumor.spectrum.pdf");
        let error = spectrum.write_pdf(&blocked, "tumor").unwrap_err();
        assert!(
            format!("{error}").contains("tumor.spectrum.pdf/tumor.spectrum.pdf"),
            "{error}"
        );
    }

    /// The axis steps split the tallest bar into about four round ticks.
    #[test]
    fn test_tick_steps_are_round() {
        for (top, step) in [
            (810.0, 250.0),
            (1.0, 0.25),
            (3.0, 1.0),
            (47.0, 20.0),
            (9.5, 2.5),
        ] {
            assert_eq!(tick_step(top), step, "{top}");
        }
    }
}
