#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = ["matplotlib>=3.10", "numpy>=2.2", "pysam>=0.23"]
# ///
"""Redraw the README figures: simulate a damaged duplex sample, run chaff on it, and plot."""

import array
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

REPO = Path(__file__).resolve().parents[1]
WORK, OUT = REPO / "target" / "figures", REPO / ".github" / "img"
SEED, REAL, ARTIFACTS, ALT_MOLECULES = 11, 6000, 2000, (2, 3, 5, 10)
DEPTH, DISPERSION, MEDIAN, SIGMA, FILL_IN, CLOCK = 400, 20, 200, 0.35, 30, 0.12
SPACING, READ_LENGTH, THRESHOLD = 1200, 150, 0.05
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
CLASS_COLOR = {"C>A": "#1ebff0", "C>G": "#050708", "C>T": "#e62725", "T>A": "#cbcacb", "T>C": "#a1cf64", "T>G": "#edc8c5"}
CHANNELS = [f"{a}[{c}]{b}" for c in CLASSES for a in "ACGT" for b in "ACGT"]
CPG_CT = [ch for ch in CHANNELS if ch[2:5] == "C>T" and ch[-1] == "G"]
ALT_COLOR, REF_COLOR, REAL_COLOR, GRAY, GREEN = "#e34948", "#8c8b86", "#2a78d6", "0.45", "#12875c"
RAMP = ["#b7aee8", "#8a7cd3", "#5f4fbb", "#33267f"]
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
NOTE = (
    f"Simulated duplex sample: {REAL:,} real mutations and {ARTIFACTS:,} copied-damage calls at CpG C>T, "
    f"a quarter of each with {', '.join(map(str, ALT_MOLECULES[:-1]))}, or {ALT_MOLECULES[-1]} alternate molecules;\n"
    f"about {DEPTH} duplex molecules per site and a median fragment of {MEDIAN} bp;\n"
    f"fill-in copies a lesion onto the partner strand from the lesion strand's 5′ end, over an exponential mean of {FILL_IN} bp."
)


def spectrum(rng):
    weights = {"C>A": 0.12, "C>G": 0.10, "C>T": 0.28, "T>A": 0.12, "T>C": 0.25, "T>G": 0.13}
    broad = np.concatenate([weights[c] * rng.dirichlet(np.full(16, 8.0)) for c in CLASSES])
    clock = np.array([ch in CPG_CT for ch in CHANNELS], dtype=float)
    p = CLOCK * clock / clock.sum() + (1 - CLOCK) * broad / broad.sum()
    return p / p.sum()


def fragment(rng, site):
    while True:
        length = int(np.clip(round(rng.lognormal(np.log(MEDIAN), SIGMA)), 60, 600))
        offset, read_length = int(rng.integers(0, length)), min(READ_LENGTH, length)
        if offset < read_length or offset >= length - read_length:
            return site - offset, length, read_length


def pair(header, name, reference, start, length, read_length, site, base, rng):
    end, first_forward = start + length - 1, bool(rng.integers(0, 2))
    a, b = (int(d) for d in rng.integers(1, 6, size=2))
    for forward in (True, False):
        read_start = start if forward else end - read_length + 1
        sequence = list(reference[read_start:read_start + read_length])
        qualities = rng.integers(30, 61, size=read_length).astype(np.uint8)
        if read_start <= site < read_start + read_length:
            sequence[site - read_start] = base
            qualities[site - read_start] = 2 if base == "N" else qualities[site - read_start]
        read = pysam.AlignedSegment(header)
        read.query_name, read.query_sequence = name, "".join(sequence)
        read.flag = 0x3 | (0x20 if forward else 0x10) | (0x40 if forward == first_forward else 0x80)
        read.reference_id, read.reference_start, read.mapping_quality = 0, read_start, 60
        read.cigartuples = [(0, read_length)]
        read.next_reference_id = 0
        read.next_reference_start = end - read_length + 1 if forward else start
        read.template_length = length if forward else -length
        read.query_qualities = array.array("B", qualities.tobytes())
        read.set_tags([("MC", f"{read_length}M"), ("RG", "tumor"), ("aD", a), ("bD", b), ("cD", a + b)])
        yield read


def simulate():
    """Write a reference, reads, and calls of real mutations and copied damage, and return the truth by 1-based position."""
    rng = np.random.default_rng(SEED)
    p = spectrum(rng)
    sites = [{"kind": "real", "channel": str(rng.choice(CHANNELS, p=p)), "n": ALT_MOLECULES[i % 4]} for i in range(REAL)]
    sites += [{"kind": "artifact", "channel": str(rng.choice(CPG_CT)), "n": ALT_MOLECULES[i % 4]} for i in range(ARTIFACTS)]
    rng.shuffle(sites)
    reference = np.array(list("ACGT"))[rng.choice(4, size=SPACING * len(sites), p=[0.295, 0.205, 0.205, 0.295])]
    for i, s in enumerate(sites):
        s["site"], s["forward"] = i * SPACING + SPACING // 2, bool(rng.integers(0, 2))
        context, s["ref"], s["alt"] = s["channel"][0] + s["channel"][2] + s["channel"][6], s["channel"][2], s["channel"][4]
        if not s["forward"]:
            context, s["ref"], s["alt"] = (x.translate(COMPLEMENT) for x in (context[::-1], s["ref"], s["alt"]))
        reference[s["site"] - 1:s["site"] + 2] = list(context)
    reference = "".join(reference)
    WORK.mkdir(parents=True, exist_ok=True)
    (WORK / "ref.fa").write_text(">chr1\n" + "".join(reference[i:i + 80] + "\n" for i in range(0, len(reference), 80)))
    pysam.faidx(str(WORK / "ref.fa"))
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(reference)}], "RG": [{"ID": "tumor", "SM": "tumor"}]})
    with pysam.AlignmentFile(str(WORK / "reads.bam"), "wb", header=header) as bam:
        for i, s in enumerate(sites):
            depth, molecules, alt = int(rng.negative_binomial(DISPERSION, DISPERSION / (DISPERSION + DEPTH))), [], 0
            while alt < s["n"]:
                start, length, read_length = fragment(rng, s["site"])
                five_prime = s["site"] - start if s["forward"] else start + length - 1 - s["site"]
                copied = s["kind"] == "real" or rng.random() < np.exp(-five_prime / FILL_IN)
                molecules.append((start, length, read_length, s["alt"] if copied else "N"))
                alt += copied
            molecules += [(*fragment(rng, s["site"]), s["ref"]) for _ in range(max(depth, len(molecules)) - len(molecules))]
            reads = [r for j, m in enumerate(molecules) for r in pair(header, f"s{i}m{j}", reference, *m[:3], s["site"], m[3], rng)]
            for read in sorted(reads, key=lambda r: r.reference_start):
                bam.write(read)
    pysam.index(str(WORK / "reads.bam"))
    lines = ["##fileformat=VCFv4.2", f"##contig=<ID=chr1,length={len(reference)}>", '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor"]
    lines += [f"chr1\t{s['site'] + 1}\t.\t{s['ref']}\t{s['alt']}\t.\t.\t.\tGT\t0/1" for s in sites]
    (WORK / "calls.vcf").write_text("\n".join(lines) + "\n")
    return {s["site"] + 1: s for s in sites}


def run_chaff():
    subprocess.run(["cargo", "build", "--release", "--manifest-path", str(REPO / "Cargo.toml")], check=True)
    for model in ("chaff", "fgbio"):
        subprocess.run([str(REPO / "target" / "release" / "chaff"), "--input", WORK / "calls.vcf", "--bam", WORK / "reads.bam", "--ref", WORK / "ref.fa",
                        "--sample", "tumor", "--filters", "copied-damage", "--copied-damage-threshold", str(THRESHOLD), "--model", model,
                        "--output", WORK / f"{model}.vcf", "--metrics", WORK / f"{model}.tsv"], check=True)


def calls(model):
    with pysam.VariantFile(str(WORK / f"{model}.vcf")) as vcf:
        return {r.pos: {"filtered": "CopiedDamageArtifact" in r.filter.keys(), "cdlr": r.info.get("CDLR"),
                        "cdac": r.info.get("CDAC"), "cdrc": r.info.get("CDRC")} for r in vcf}


def distances(truth):
    """Each molecule's distance from the lesion strand's 5' and 3' ends, by kind, at the C>T and G>T sites."""
    out = defaultdict(lambda: ([], []))
    with pysam.AlignmentFile(str(WORK / "reads.bam")) as bam:
        for t in truth.values():
            forward = (t["ref"], t["alt"]) in {("C", "T"), ("G", "T")}
            if not forward and (t["ref"], t["alt"]) not in {("G", "A"), ("C", "A")}:
                continue
            templates = defaultdict(list)
            for read in bam.fetch("chr1", t["site"], t["site"] + 1):
                templates[read.query_name].append(read)
            for reads in templates.values():
                r = reads[0]
                start, end = (r.next_reference_start, r.reference_end - 1) if r.is_reverse else (r.reference_start, r.next_reference_start + int(r.get_tag("MC")[:-1]) - 1)
                bases = {x.query_sequence[t["site"] - x.reference_start] for x in reads if x.reference_start <= t["site"] < x.reference_end}
                base = bases.pop() if len(bases) == 1 else "N"
                kind = t["kind"] if base == t["alt"] else "reference" if base == t["ref"] else None
                if kind:
                    left, right = t["site"] - start, end - t["site"]
                    out[kind][0].append(left if forward else right)
                    out[kind][1].append(right if forward else left)
    return out


def label(ax, x, y, parts, **kwargs):
    """Text in parts of their own styles, each after the last, from axes coordinates `x` and `y`."""
    text = ax.text(x, y, parts[0][0], transform=ax.transAxes, **parts[0][1], **kwargs)
    for words, style in parts[1:]:
        text = ax.annotate(words, xy=(1, 0), xycoords=text, **style, **kwargs)


def share(called, truth, kind, field):
    pairs = [called[p][field] for p, t in truth.items() if called[p][field] is not None and kind in (None, t["kind"])]
    return 100 * sum(a for a, _ in pairs) / sum(n for _, n in pairs)


def finish(fig, name, title, note):
    box = fig.get_tightbbox(fig.canvas.get_renderer())
    width, height = fig.get_figwidth(), fig.get_figheight()
    fig.text(box.x0 / width, (box.y0 - 0.15) / height, note, fontsize=8.5, color=GRAY, va="top", linespacing=1.4)
    box = fig.get_tightbbox(fig.canvas.get_renderer())
    fig.text(0.5 * (box.x0 + box.x1) / width, (box.y1 + 0.12) / height, title, ha="center", va="bottom", fontsize=11.5, fontweight="bold")
    fig.savefig(OUT / name, dpi=200, bbox_inches="tight", pad_inches=0.15, facecolor="white")
    plt.close(fig)


def ends_figure(truth, called):
    scale = float(next(line.split("\t")[7] for line in (WORK / "chaff.tsv").read_text().splitlines()[1:]))
    groups = [("artifact", "Copied damage, alternate molecules", ALT_COLOR, "cdac"), ("real", "Real mutations, alternate molecules", REAL_COLOR, "cdac"),
              ("reference", "Reference molecules", REF_COLOR, "cdrc")]
    measured, bins = distances(truth), np.arange(0, 401, 10)
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.2), sharey=True)
    for ax, end, end_label in zip(axes, (0, 1), ("5′", "3′")):
        for kind, _, color, _ in groups:
            counts, _ = np.histogram(measured[kind][end], bins=bins)
            counts = 100 * counts / len(measured[kind][end])
            if kind == "reference":
                ax.fill_between(bins[:-1], counts, step="post", color=color, alpha=0.3, lw=0)
            ax.step(bins[:-1], counts, where="post", color=color, lw=1.2 if kind == "reference" else 1.8)
        ax.set_xlim(0, 400)
        ax.set_xlabel(f"Distance from the lesion strand's {end_label} end (bp)")
        bold = {"fontweight": "bold"}
        label(ax, 0, 1.04, [("From the lesion strand's ", bold), (f"{end_label} end", {**bold, "color": GREEN})], fontsize=10, va="bottom")
    axes[0].set_ylabel("Molecules per 10 bp bin (%)")
    axes[0].set_ylim(0, None)
    axes[0].axvline(scale, color=GRAY, lw=0.9, ls="--", zorder=1)
    axes[0].text(scale + 8, axes[0].get_ylim()[1] * 0.92, f"Learned scale, {scale:.1f} bp", fontsize=8, color=GRAY, va="top")
    handles = []
    for kind, name, color, field in groups:
        text = f"{name}: {share(called, truth, None if kind == 'reference' else kind, field):.0f}%" + (f" within {scale:.0f} bp" if kind == "artifact" else "")
        handles.append(Patch(facecolor=color, alpha=0.3, edgecolor=color, label=text) if kind == "reference" else Line2D([], [], color=color, lw=1.8, label=text))
    axes[1].legend(handles=handles, loc="upper right", fontsize=8.5, handlelength=1.6, bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout(w_pad=2.0)
    finish(fig, "copied-damage-ends.png", "Copied Damage Crowds the Lesion Strand's 5′ End; Real Mutations Follow the Reference",
           NOTE + "\nThe dashed line is the decay scale chaff learned; the legend's shares are its CDAC and CDRC counts.")


def spectrum_row(ax, truth, keep, title, detail):
    real, artifact = np.zeros(len(CHANNELS)), np.zeros(len(CHANNELS))
    for pos, t in truth.items():
        if keep(pos):
            (real if t["kind"] == "real" else artifact)[CHANNELS.index(t["channel"])] += 1
    colors = [CLASS_COLOR[ch[2:5]] for ch in CHANNELS]
    ax.bar(range(len(CHANNELS)), real, width=0.78, color=colors, linewidth=0)
    ax.bar(range(len(CHANNELS)), artifact, bottom=real, width=0.78, color="white", edgecolor=colors, hatch="//////", linewidth=0.5)
    ax.set_xlim(-0.7, len(CHANNELS) - 0.3)
    ax.set_xticks([])
    ax.set_ylabel("Calls")
    label(ax, 0.01, 0.9, [(title, {"fontweight": "bold"}), (detail, {})], fontsize=9.5, va="bottom")
    return (real + artifact).max()


def roc(called, truth, n):
    def scores(kind):
        return np.array([-np.inf if called[p]["cdlr"] is None else called[p]["cdlr"] for p, t in truth.items()
                         if t["kind"] == kind and t["n"] == n and (kind == "artifact" or t["channel"] in CPG_CT)])
    art, real = scores("artifact"), scores("real")
    cuts = np.unique(np.concatenate([art, real]))[::-1]
    return [0.0] + [100 * np.mean(real >= c) for c in cuts], [0.0] + [100 * np.mean(art >= c) for c in cuts]


def operating_point(called, truth, n):
    def rate(kind):
        return 100 * np.mean([called[p]["filtered"] for p, t in truth.items() if t["kind"] == kind and t["n"] == n and (kind == "artifact" or t["channel"] in CPG_CT)])
    return rate("real"), rate("artifact")


def outcome_figure(truth, chaff, fgbio):
    fig = plt.figure(figsize=(10.0, 4.5))
    grid = fig.add_gridspec(2, 2, width_ratios=[2.4, 1], hspace=0.12, wspace=0.16)
    before = fig.add_subplot(grid[0, 0])
    after = fig.add_subplot(grid[1, 0], sharey=before)
    top = spectrum_row(before, truth, lambda p: True, "Before chaff", ": all calls")
    spectrum_row(after, truth, lambda p: not chaff[p]["filtered"], "After chaff", f": the calls it passes at a threshold of {THRESHOLD}")
    before.set_ylim(0, top * 1.02)
    for i, cls in enumerate(CLASSES):
        before.add_patch(plt.Rectangle((i * 16 - 0.45, top * 1.04), 15.9, top * 0.05, color=CLASS_COLOR[cls], clip_on=False, lw=0))
        before.text(i * 16 + 7.5, top * 1.12, cls, ha="center", va="bottom", fontsize=8.5)
    after.set_xticks(range(len(CHANNELS)), [ch[0] + ch[2] + ch[6] for ch in CHANNELS], rotation=90, fontsize=5, fontfamily=["Menlo", "DejaVu Sans Mono"])
    after.tick_params(axis="x", length=0, pad=2)
    for label, ch in zip(after.get_xticklabels(), CHANNELS):
        label.set_fontweight("bold" if ch in CPG_CT else "normal")
    handles = [Patch(facecolor=CLASS_COLOR["C>T"], label="Real mutations"), Patch(facecolor="white", edgecolor=CLASS_COLOR["C>T"], hatch="//////", label="Copied damage")]
    after.legend(handles=handles, loc="upper right", fontsize=8.5, ncol=2, handlelength=1.4, bbox_to_anchor=(1.0, 1.0))
    ax = fig.add_subplot(grid[:, 1])
    for n, color in zip(ALT_MOLECULES, RAMP):
        ax.plot(*roc(chaff, truth, n), color=color, lw=1.6, drawstyle="steps-post", zorder=2)
        for called, face in ((chaff, color), (fgbio, "white")):
            ax.scatter(*operating_point(called, truth, n), s=34, facecolors=face, edgecolors=color, linewidths=1.4, zorder=4)
    ax.set_xscale("symlog", linthresh=1, linscale=0.6)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 102)
    ax.set_xticks([0, 1, 10, 100], ["0", "1", "10", "100"])
    ax.set_xlabel("Real C>T at CpG filtered (%)")
    ax.set_ylabel("Copied damage filtered (%)")
    handles = [Line2D([], [], color=c, lw=1.6, label=f"{n} alternate molecules") for n, c in zip(ALT_MOLECULES, RAMP)]
    handles += [Line2D([], [], ls="", marker="o", ms=6, mfc=face, mec="black", mew=1.4 if face == "white" else 1.0, label=f"--model {m} at {THRESHOLD}")
                for m, face in (("chaff", "black"), ("fgbio", "white"))]
    ax.legend(handles=handles, loc="lower right", fontsize=8.5, handlelength=1.4)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.12, top=0.9)
    finish(fig, "copied-damage-filtering.png", "Chaff Removes Most Copied Damage and Keeps Real Mutations; 2 Molecules Are Its Limit",
           NOTE + "\nCurves sweep the threshold under the chaff model, and dots mark a threshold of 0.05. Costs count the real C>T at CpG, the stratum copied damage shares.")


if __name__ == "__main__":
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8, "legend.frameon": False,
    })
    truth = simulate()
    run_chaff()
    chaff = calls("chaff")
    ends_figure(truth, chaff)
    outcome_figure(truth, chaff, calls("fgbio"))
