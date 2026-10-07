//! The trinucleotide spectrum of the SNVs a run scores, before and after its
//! filters, written as a PDF.
//!
//! Each SNV falls in one of 96 channels: its substitution read from the
//! pyrimidine of the base pair (C>A, C>G, C>T, T>A, T>C, or T>G) and the
//! bases 5' and 3' of it on that strand. One PDF page stacks its views on one
//! scale, so their heights compare:
//!
//! - **Before filtering:** every heterozygous SNV of the sample, whatever its
//!   FILTER, an SNV without a genotype counting as heterozygous.
//! - **Expected real:** each SNV weighed by the product of the posteriors of
//!   the enabled filters, its expected count as a real mutation, a call
//!   without a posterior weighing one.
//! - **Passing:** the SNVs no filter flagged, drawn only when some filter has
//!   a threshold.
//!
//! The views share a page rather than taking a page each because `kuva`
//! renders one scene to one single-page PDF.

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

/// The spectrum of one sample's scored SNVs, before and after filtering.
#[derive(Clone, Debug, PartialEq)]
pub struct Spectrum {
    /// The heterozygous SNVs in each channel.
    pub before: [f64; CHANNELS],
    /// The expected real SNVs in each channel: each weighed by the product
    /// of its posteriors.
    pub weighted: [f64; CHANNELS],
    /// The SNVs in each channel that no filter flagged, when some filter has
    /// a threshold.
    pub passing: Option<[f64; CHANNELS]>,
}

impl Spectrum {
    /// An empty spectrum, counting passing SNVs when `thresholded`.
    pub fn new(thresholded: bool) -> Self {
        Self {
            before: [0.0; CHANNELS],
            weighted: [0.0; CHANNELS],
            passing: thresholded.then_some([0.0; CHANNELS]),
        }
    }

    /// Count one scored SNV in `channel`, flagged by some filter or not, with
    /// the product of its posteriors.
    pub fn add(&mut self, channel: usize, flagged: bool, weight: f64) {
        self.before[channel] += 1.0;
        self.weighted[channel] += weight;
        if let Some(passing) = &mut self.passing {
            passing[channel] += f64::from(u8::from(!flagged));
        }
    }

    /// The figure: a panel of 96 bars per view, stacked on one axis, each
    /// under a strip of the class colours, with the CpG C>T contexts named
    /// after filtering.
    fn figure(&self, sample: &str) -> Figure {
        let top = self.before.iter().copied().fold(1.0, f64::max);
        let mut views = vec![
            (
                &self.before,
                "Before chaff: Every Heterozygous SNV",
                "SNVs",
                false,
            ),
            (
                &self.weighted,
                "After chaff: Expected Real SNVs, Each SNV Weighted by Its Posteriors",
                "Expected Real SNVs",
                true,
            ),
        ];
        if let Some(passing) = &self.passing {
            views.push((
                passing,
                "After chaff: SNVs Passing Every Threshold",
                "Passing SNVs",
                true,
            ));
        }
        let (mut plots, mut layouts) = (Vec::new(), Vec::new());
        for (values, title, y_label, named) in views {
            let (panel_plots, mut layout) = panel(values, title, y_label, top);
            if named {
                for i in (0..CHANNELS).filter(|&i| is_cpg_c_to_t(i)) {
                    let label =
                        TextAnnotation::new(context(i), i as f64 + 1.0, values[i] + top * 0.03)
                            .with_color(CPG_COLOR)
                            .with_font_size(10);
                    layout = layout.with_annotation(label);
                }
            }
            plots.push(panel_plots);
            layouts.push(layout);
        }
        Figure::new(plots.len(), 1)
            .with_title(title(sample))
            .with_plots(plots)
            .with_layouts(layouts)
            .with_cell_size(1100.0, 380.0)
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

/// The share of a panel's axis above its last tick that holds the class
/// strip and names.
const STRIP_SHARE: f64 = 0.24;

/// A panel's y axis for bars up to `top`: the tick step, the last tick at or
/// above `top`, and the axis maximum. The step is the smallest round one,
/// from a fifth of `top` up, whose next tick lies past the maximum, so the
/// band above the last tick, [`STRIP_SHARE`] of the axis, holds no tick.
fn axis(top: f64) -> (f64, f64, f64) {
    let magnitude = 10f64.powf((top / 5.0).log10().floor());
    let steps = (0..4).flat_map(|e| [1.0, 2.0, 2.5, 5.0].map(|m| m * magnitude * 10f64.powi(e)));
    for step in steps.filter(|s| *s >= top / 5.0) {
        let ceiling = (top / step).ceil() * step;
        let max = ceiling / (1.0 - STRIP_SHARE);
        if max - ceiling < step {
            return (step, ceiling, max);
        }
    }
    (top, top, top / (1.0 - STRIP_SHARE))
}

/// One panel: a bar per channel in its class's colour, under the class names
/// and a strip of their colours. The axis ticks stop at the first tick at or
/// above `top`, the tallest bar before filtering, and the strip and names sit
/// between it and the next tick, which the axis never reaches, with room above
/// the names. The strip is a stacked bar per channel on an unpainted base.
fn panel(values: &[f64; CHANNELS], title: &str, y_label: &str, top: f64) -> (Vec<Plot>, Layout) {
    let (step, ceiling, max) = axis(top);
    let band = max - ceiling;
    let mut bars = BarPlot::new();
    let mut strip = BarPlot::new().with_width(1.03).with_stacked();
    for (i, value) in values.iter().enumerate() {
        let color = CLASS_COLORS[i / 16];
        bars = bars.with_colored_bar(context(i), *value, color);
        strip = strip.with_group(
            context(i),
            [(ceiling + band * 0.12, "none"), (band * 0.15, color)],
        );
    }
    let plots = vec![Plot::Bar(strip), Plot::Bar(bars)];
    let mut layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_y_label(y_label)
        .with_y_axis_min(0.0)
        .with_y_axis_max(max)
        .with_y_tick_step(step)
        .with_show_grid(false)
        .with_x_tick_rotate(90.0)
        .with_tick_size(8);
    for (class, name) in CLASSES.iter().enumerate() {
        let label = TextAnnotation::new(*name, 16.0 * class as f64 + 8.5, ceiling + band * 0.48)
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
    fn test_spectrum_counts_before_weighted_and_passing() {
        let acg = 16 * 2 + 2;
        let mut spectrum = Spectrum::new(true);
        spectrum.add(acg, true, 0.01);
        spectrum.add(acg, false, 0.9);
        spectrum.add(0, false, 1.0);
        assert_eq!((spectrum.before[acg], spectrum.before[0]), (2.0, 1.0));
        assert!((spectrum.weighted[acg] - 0.91).abs() < 1e-12);
        assert_eq!(spectrum.weighted[0], 1.0);
        let passing = spectrum.passing.unwrap();
        assert_eq!((passing[acg], passing[0]), (1.0, 1.0));

        let mut unthresholded = Spectrum::new(false);
        unthresholded.add(acg, true, 0.01);
        assert_eq!(unthresholded.passing, None);
        assert_eq!(unthresholded.weighted[acg], 0.01);
    }

    /// The figure holds a panel per view, each with its own y label.
    #[test]
    fn test_the_figure_draws_one_panel_per_view() {
        let labels = |spectrum: &Spectrum| -> Vec<String> {
            let scene = spectrum.figure("tumor").render();
            ["SNVs", "Expected Real SNVs", "Passing SNVs"]
                .into_iter()
                .filter(|label| {
                    scene
                        .elements
                        .iter()
                        .any(|e| matches!(e, Primitive::Text { content, .. } if content == label))
                })
                .map(String::from)
                .collect()
        };
        assert_eq!(
            labels(&Spectrum::new(false)),
            ["SNVs", "Expected Real SNVs"]
        );
        assert_eq!(
            labels(&Spectrum::new(true)),
            ["SNVs", "Expected Real SNVs", "Passing SNVs"]
        );
    }

    #[test]
    fn test_write_pdf_writes_a_pdf() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("nested/tumor.spectrum.pdf");
        let mut spectrum = Spectrum::new(true);
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

    /// The axis steps round, ends on a tick at or above the tallest bar, and
    /// leaves the band above it for the strip with no tick in it.
    #[test]
    fn test_the_axis_keeps_ticks_out_of_the_strip() {
        for (top, step, ceiling) in [
            (810.0, 500.0, 1000.0),
            (300.0, 100.0, 300.0),
            (1.0, 0.5, 1.0),
            (47.0, 20.0, 60.0),
        ] {
            let (s, c, max) = axis(top);
            assert_eq!((s, c), (step, ceiling), "{top}");
            assert!(c >= top && c + s > max, "{top}");
            assert!(((max - c) / max - STRIP_SHARE).abs() < 1e-12, "{top}");
        }
    }
}
