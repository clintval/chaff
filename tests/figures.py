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
DEPTH, DISPERSION, MEDIAN, SIGMA, FILL_IN, NICK, CLOCK = 400, 20, 200, 0.35, 30, 0.2, 0.12
SPACING, READ_LENGTH, THRESHOLD = 1200, 150, 0.05
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
CLASS_COLOR = {"C>A": "#1ebff0", "C>G": "#050708", "C>T": "#e62725", "T>A": "#cbcacb", "T>C": "#a1cf64", "T>G": "#edc8c5"}
CHANNELS = [f"{a}[{c}]{b}" for c in CLASSES for a in "ACGT" for b in "ACGT"]
CPG_CT = [ch for ch in CHANNELS if ch[2:5] == "C>T" and ch[-1] == "G"]
ALT_COLOR, REF_COLOR, REAL_COLOR, GRAY, GREEN = "#e34948", "#8c8b86", "#2a78d6", "0.45", "#12875c"
RAMP = ["#b7aee8", "#8a7cd3", "#5f4fbb", "#33267f"]
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
LEARNING_TITLE = "Where Copied Damage Is Heavy, the Learned Prior Understates It"
BURDEN_TITLE = "Weighting Calls by CDAP Removes Most of the Burden Inflation; a Threshold Removes Little"


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


def simulate(work=WORK, seed=SEED, real=REAL, artifacts=ARTIFACTS, fill_in=FILL_IN, cpg_only=False, counts=ALT_MOLECULES, nick=NICK):
    """Write a reference, reads, and calls of real mutations and copied damage, and return the truth by 1-based position."""
    rng = np.random.default_rng(seed)
    p = np.array([ch in CPG_CT for ch in CHANNELS], dtype=float) / 4 if cpg_only else spectrum(rng)
    sites = [{"kind": "real", "channel": str(rng.choice(CHANNELS, p=p)), "n": counts[i % len(counts)]} for i in range(real)]
    sites += [{"kind": "artifact", "channel": str(rng.choice(CPG_CT)), "n": counts[i % len(counts)]} for i in range(artifacts)]
    rng.shuffle(sites)
    reference = np.array(list("ACGT"))[rng.choice(4, size=SPACING * len(sites), p=[0.295, 0.205, 0.205, 0.295])]
    for i, s in enumerate(sites):
        s["site"], s["forward"] = i * SPACING + SPACING // 2, bool(rng.integers(0, 2))
        context, s["ref"], s["alt"] = s["channel"][0] + s["channel"][2] + s["channel"][6], s["channel"][2], s["channel"][4]
        if not s["forward"]:
            context, s["ref"], s["alt"] = (x.translate(COMPLEMENT) for x in (context[::-1], s["ref"], s["alt"]))
        reference[s["site"] - 1:s["site"] + 2] = list(context)
    reference = "".join(reference)
    work.mkdir(parents=True, exist_ok=True)
    (work / "ref.fa").write_text(">chr1\n" + "".join(reference[i:i + 80] + "\n" for i in range(0, len(reference), 80)))
    pysam.faidx(str(work / "ref.fa"))
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(reference)}], "RG": [{"ID": "tumor", "SM": "tumor"}]})
    with pysam.AlignmentFile(str(work / "reads.bam"), "wb", header=header) as bam:
        for i, s in enumerate(sites):
            depth, molecules, alt = int(rng.negative_binomial(DISPERSION, DISPERSION / (DISPERSION + DEPTH))), [], 0
            while alt < s["n"]:
                start, length, read_length = fragment(rng, s["site"])
                if s["kind"] == "artifact" and rng.random() >= nick:
                    while rng.random() >= np.exp(-(s["site"] - start if s["forward"] else start + length - 1 - s["site"]) / fill_in):
                        molecules.append((start, length, read_length, "N"))
                        start, length, read_length = fragment(rng, s["site"])
                molecules.append((start, length, read_length, s["alt"]))
                alt += 1
            molecules += [(*fragment(rng, s["site"]), s["ref"]) for _ in range(max(depth, len(molecules)) - len(molecules))]
            reads = [r for j, m in enumerate(molecules) for r in pair(header, f"s{i}m{j}", reference, *m[:3], s["site"], m[3], rng)]
            for read in sorted(reads, key=lambda r: r.reference_start):
                bam.write(read)
    pysam.index(str(work / "reads.bam"))
    lines = ["##fileformat=VCFv4.2", f"##contig=<ID=chr1,length={len(reference)}>", '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor"]
    lines += [f"chr1\t{s['site'] + 1}\t.\t{s['ref']}\t{s['alt']}\t.\t.\t.\tGT\t0/1" for s in sites]
    (work / "calls.vcf").write_text("\n".join(lines) + "\n")
    return {s["site"] + 1: s for s in sites}


def run_chaff(work=WORK):
    for model in ("chaff", "fgbio"):
        subprocess.run([str(REPO / "target" / "release" / "chaff"), "--input", work / "calls.vcf", "--bam", work / "reads.bam", "--ref", work / "ref.fa",
                        "--sample", "tumor", "--filters", "copied-damage", "--copied-damage-threshold", str(THRESHOLD), "--model", model,
                        "--output", work / f"{model}.vcf", "--metrics", work / f"{model}.tsv"], check=True, stderr=subprocess.DEVNULL)


def calls(model, work=WORK):
    with pysam.VariantFile(str(work / f"{model}.vcf")) as vcf:
        return {r.pos: {"filtered": "CopiedDamageArtifact" in r.filter.keys(), "cdlr": r.info.get("CDLR"), "cdap": r.info.get("CDAP"),
                        "cdac": r.info.get("CDAC"), "cdrc": r.info.get("CDRC")} for r in vcf}


def metrics(work, column):
    lines = [line.split("\t") for line in (work / "chaff.tsv").read_text().splitlines()]
    return float(dict(zip(lines[0], lines[1]))[column])


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


def finish(fig, name, title):
    box = fig.get_tightbbox(fig.canvas.get_renderer())
    width, height = fig.get_figwidth(), fig.get_figheight()
    fig.text(0.5 * (box.x0 + box.x1) / width, (box.y1 + 0.12) / height, title, ha="center", va="bottom", fontsize=11.5, fontweight="bold")
    fig.savefig(OUT / name, dpi=200, bbox_inches="tight", pad_inches=0.15, facecolor="white")
    plt.close(fig)


def ends_figure(truth, called):
    scale = metrics(WORK, "distance")
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
        label(ax, 0, 1.04, [("From the Lesion Strand's ", bold), (f"{end_label} End", {**bold, "color": GREEN})], fontsize=10, va="bottom")
    axes[0].set_ylabel("Molecules per 10 bp bin (%)")
    axes[0].set_ylim(0, None)
    axes[0].axvline(scale, color=GRAY, lw=0.9, ls="--", zorder=1)
    axes[0].text(scale + 8, axes[0].get_ylim()[1] * 0.92, f"Learned scale, {scale:.1f} bp", fontsize=8, color=GRAY, va="top")
    handles = []
    for kind, name, color, field in groups:
        text = f"{name}: {share(called, truth, None if kind == 'reference' else kind, field):.0f}%" + (f" within {scale:.0f} bp" if kind == "artifact" else "")
        handles.append(Patch(facecolor=color, alpha=0.3, edgecolor=color, label=text) if kind == "reference" else Line2D([], [], color=color, lw=1.8, label=text))
    axes[0].legend(handles=handles, loc="upper right", fontsize=8.5, handlelength=1.6, bbox_to_anchor=(1.0, 0.8))
    fig.tight_layout(w_pad=2.0)
    finish(fig, "copied-damage-ends.png", "Copied Damage Crowds the Lesion Strand's 5′ End; Real Mutations Follow the Reference")


def spectrum_row(ax, truth, keep, title, detail):
    real, artifact = np.zeros(len(CHANNELS)), np.zeros(len(CHANNELS))
    for pos, t in truth.items():
        if keep(pos):
            (real if t["kind"] == "real" else artifact)[CHANNELS.index(t["channel"])] += 1
    colors = [CLASS_COLOR[ch[2:5]] for ch in CHANNELS]
    ax.bar(range(len(CHANNELS)), real, width=0.68, color=colors, edgecolor=colors, linewidth=0.5)
    ax.bar(range(len(CHANNELS)), artifact, bottom=real, width=0.68, color="white", edgecolor=colors, hatch="//////", linewidth=0.5)
    ax.set_xlim(-0.7, len(CHANNELS) - 0.3)
    ax.set_xticks([])
    ax.set_ylabel("Calls")
    label(ax, 0.01, 0.9, [(title, {"fontweight": "bold"}), (detail, {})], fontsize=9.5, va="bottom")
    return real + artifact


def roc(called, truth, n):
    def scores(kind):
        return np.array([-np.inf if called[p]["cdlr"] is None else called[p]["cdlr"] for p, t in truth.items()
                         if t["kind"] == kind and t["n"] == n and (kind == "artifact" or t["channel"] in CPG_CT)])
    art, real = scores("artifact"), scores("real")
    cuts = np.unique(np.concatenate([art, real]))[::-1]
    return [0.0] + [100 * np.mean(real >= c) for c in cuts], [0.0] + [100 * np.mean(art >= c) for c in cuts]


def auc(called, truth, n):
    xs, ys = (np.array(v) / 100 for v in roc(called, truth, n))
    return float(np.sum(np.diff(xs) * (ys[1:] + ys[:-1]) / 2) + (1 - xs[-1]) * ys[-1])


def operating_point(called, truth, n):
    def rate(kind):
        return 100 * np.mean([called[p]["filtered"] for p, t in truth.items() if t["kind"] == kind and t["n"] == n and (kind == "artifact" or t["channel"] in CPG_CT)])
    return rate("real"), rate("artifact")


def outcome_figure(truth, chaff, fgbio):
    fig = plt.figure(figsize=(10.0, 4.5))
    grid = fig.add_gridspec(2, 2, width_ratios=[2.4, 1], hspace=0.12, wspace=0.16)
    before = fig.add_subplot(grid[0, 0])
    after = fig.add_subplot(grid[1, 0], sharey=before)
    top = spectrum_row(before, truth, lambda p: True, "Before chaff", ": All Calls").max()
    heights = spectrum_row(after, truth, lambda p: not chaff[p]["filtered"], "After chaff", f": Calls Passing a Threshold of {THRESHOLD}")
    for ch in CPG_CT:
        i = CHANNELS.index(ch)
        after.text(i, heights[i] + top * 0.05, ch[0] + "CG", rotation=90, ha="center", va="bottom", fontsize=7.5, color="#12875c",
                   fontfamily=["Menlo", "DejaVu Sans Mono"], bbox={"boxstyle": "round,pad=0.25,rounding_size=0.4", "facecolor": "#1baf7a", "alpha": 0.18, "edgecolor": "none"})
    before.set_ylim(0, top * 1.02)
    for i, cls in enumerate(CLASSES):
        before.add_patch(plt.Rectangle((i * 16 - 0.45, top * 1.04), 15.9, top * 0.05, color=CLASS_COLOR[cls], clip_on=False, lw=0))
        before.text(i * 16 + 7.5, top * 1.12, cls, ha="center", va="bottom", fontsize=8.5)
    after.set_xticks(range(len(CHANNELS)), [ch[0] + ch[2] + ch[6] for ch in CHANNELS], rotation=90, fontsize=5, fontfamily=["Menlo", "DejaVu Sans Mono"])
    after.tick_params(axis="x", length=0, pad=2)
    for tick, ch in zip(after.get_xticklabels(), CHANNELS):
        tick.set_fontweight("bold" if ch in CPG_CT else "normal")
    handles = [Patch(facecolor=GRAY, edgecolor=GRAY, label="Real mutations"), Patch(facecolor="white", edgecolor=GRAY, hatch="//////", label="Copied damage")]
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
    handles = [Line2D([], [], color=c, lw=1.6, label=f"{n} alternate molecules (AUC {auc(chaff, truth, n):.2f})") for n, c in zip(ALT_MOLECULES, RAMP)]
    handles += [Line2D([], [], ls="", marker="o", ms=6, mfc=face, mec="black", mew=1.4 if face == "white" else 1.0, label=f"--model {m} at {THRESHOLD}")
                for m, face in (("chaff", "black"), ("fgbio", "white"))]
    ax.legend(handles=handles, loc="lower right", fontsize=8.5, handlelength=1.4)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.12, top=0.9)
    finish(fig, "copied-damage-filtering.png", "A Threshold Spares Real Mutations but Misses Most Copied Damage at 2 or 3 Molecules")


LIBRARIES = [(f, 30) for f in (0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)] + [(0.4, 15), (0.4, 60)]
LIBRARY_CALLS, LIBRARY_COUNTS = 400, (2, 2, 2, 2, 3, 3, 3, 4, 4, 5)
BINS = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0001]


def libraries():
    """Simulate small libraries of C>T calls at CpG with known copied-damage fractions and scales, and run chaff on each."""
    out = []
    for i, (fraction, scale) in enumerate(LIBRARIES):
        work = WORK / "libraries" / f"{fraction}-{scale}"
        artifacts = round(LIBRARY_CALLS * fraction)
        truth = simulate(work, SEED + i + 1, LIBRARY_CALLS - artifacts, artifacts, scale, cpg_only=True, counts=LIBRARY_COUNTS)
        run_chaff(work)
        out.append({"fraction": fraction, "scale": scale, "truth": truth, "learned": metrics(work, "artifact_fraction"),
                    "learned_scale": metrics(work, "distance"), "chaff": calls("chaff", work), "fgbio": calls("fgbio", work)})
    return out


def calibration(libs, model):
    pairs = [(lib[model][p]["cdap"], t["kind"] == "real") for lib in libs for p, t in lib["truth"].items() if lib[model][p]["cdap"] is not None]
    posterior, real = (np.array(v, dtype=float) for v in zip(*pairs))
    bins = np.digitize(posterior, BINS) - 1
    return [(posterior[bins == b].mean(), real[bins == b].mean(), int(np.sum(bins == b))) for b in range(len(BINS) - 1) if np.sum(bins == b) >= 15]


def learning_figure(libs):
    fig, (left, right) = plt.subplots(1, 2, figsize=(8.6, 4.2))
    for ax in (left, right):
        ax.set_box_aspect(1)
    left.plot([0, 70], [0, 70], color=GRAY, lw=0.9, ls="--", zorder=1)
    for lib in libs:
        marker = {15: "s", 30: "o", 60: "D"}[lib["scale"]]
        left.scatter(100 * lib["fraction"], 100 * lib["learned"], s=36, marker=marker, color=ALT_COLOR, zorder=3)
    left.set_xlim(0, 70)
    left.set_ylim(0, 70)
    left.set_xlabel("True copied damage (% of calls)")
    left.set_ylabel("Learned artifact fraction (%)")
    left.set_title("chaff Understates Heavy Damage", fontsize=10, fontweight="bold", loc="left")
    inset = left.inset_axes([0.6, 0.1, 0.36, 0.36])
    inset.plot([0, 80], [0, 80], color=GRAY, lw=0.8, ls="--", zorder=1)
    for lib in libs:
        if lib["fraction"] > 0:
            marker = {15: "s", 30: "o", 60: "D"}[lib["scale"]]
            inset.scatter(lib["scale"], lib["learned_scale"], s=16, marker=marker, color=ALT_COLOR, zorder=3)
    inset.set_xlim(0, 80)
    inset.set_ylim(0, 80)
    inset.set_xticks([0, 30, 60])
    inset.set_yticks([0, 30, 60])
    inset.tick_params(labelsize=7)
    inset.set_xlabel("True scale (bp)", fontsize=7.5)
    inset.set_ylabel("Learned (bp)", fontsize=7.5)
    right.plot([0, 1], [0, 1], color=GRAY, lw=0.9, ls="--", zorder=1)
    for model, face, style in (("chaff", REAL_COLOR, "-"), ("fgbio", "white", "--")):
        x, y, _ = zip(*calibration(libs, model))
        right.plot(x, y, color=REAL_COLOR, lw=1.2, ls=style, zorder=2)
        right.scatter(x, y, s=30, facecolors=face, edgecolors=REAL_COLOR, linewidths=1.3, zorder=3, label=f"--model {model}")
    right.set_xlim(0, 1)
    right.set_ylim(0, 1)
    right.set_xlabel("CDAP, the posterior that a call is real")
    right.set_ylabel("Calls that are real")
    right.set_title("Its Posteriors Lean Toward Real", fontsize=10, fontweight="bold", loc="left")
    right.legend(loc="lower right", fontsize=8.5, handlelength=1.4)
    shapes = [Line2D([], [], ls="", marker=m, ms=6, color=ALT_COLOR, label=f"{scale} bp fill-in") for m, scale in (("s", 15), ("o", 30), ("D", 60))]
    left.legend(handles=shapes, loc="upper left", fontsize=8.5, handlelength=1.0)
    fig.tight_layout(w_pad=3.0)
    finish(fig, "copied-damage-learning.png", LEARNING_TITLE)


def burden(lib, estimator):
    """A library's estimate of its real calls, as a multiple of the true count."""
    called, truth = lib["chaff"], lib["truth"]
    real = sum(t["kind"] == "real" for t in truth.values())
    if estimator == "raw":
        estimate = len(truth)
    elif estimator == "threshold":
        estimate = sum(not called[p]["filtered"] for p in truth)
    else:
        estimate = sum(1.0 if called[p]["cdap"] is None else called[p]["cdap"] for p in truth)
    return estimate / real


def burden_figure(libs):
    libs = [lib for lib in libs if lib["scale"] == FILL_IN]
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.axhline(1, color=GRAY, lw=0.9, ls="--", zorder=1)
    x = [100 * lib["fraction"] for lib in libs]
    for estimator, name, color in (("raw", "Every call", REF_COLOR), ("threshold", f"Calls passing a CDAP threshold of {THRESHOLD}", ALT_COLOR),
                                   ("weighted", "Calls weighted by CDAP", REAL_COLOR)):
        y = [burden(lib, estimator) for lib in libs]
        ax.plot(x, y, color=color, lw=1.6, marker="o", ms=5, label=name, zorder=3)
    ax.set_xlim(0, 62)
    ax.set_ylim(0.8, None)
    ax.set_xlabel("True copied damage (% of calls)")
    ax.set_ylabel("Estimated real calls / true real calls")
    ax.legend(loc="upper left", fontsize=8.5, handlelength=1.6)
    fig.tight_layout()
    finish(fig, "copied-damage-burden.png", BURDEN_TITLE)


if __name__ == "__main__":
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8, "legend.frameon": False,
    })
    subprocess.run(["cargo", "build", "--release", "--manifest-path", str(REPO / "Cargo.toml")], check=True)
    truth = simulate()
    run_chaff()
    chaff = calls("chaff")
    ends_figure(truth, chaff)
    outcome_figure(truth, chaff, calls("fgbio"))
    libs = libraries()
    learning_figure(libs)
    burden_figure(libs)
