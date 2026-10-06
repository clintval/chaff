# chaff

[![Build Status](https://github.com/clintval/chaff/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/clintval/chaff/actions/workflows/ci.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/clintval/chaff/badge.svg?branch=main)](https://coveralls.io/github/clintval/chaff?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-rust-dea588.svg)](https://www.rust-lang.org/)

Separate somatic variant calls from library-preparation damage artifacts.

## Installation

```console
cargo install --git https://github.com/clintval/chaff
```

## Quickstart

`chaff` scores each somatic call by where its alternate molecules sit on their templates, against the reference molecules at the same site.
The VCF and the BAM must be coordinate sorted, and neither needs an index.
Each template counts once, and a read whose mate maps to the same contig needs the `MC` tag.
The examples run on fgbio's `FilterSomaticVcf` test data in [`tests/data`](tests/data): five tumor/normal calls on `chr1` and the tumor's reads.

### Scoring Calls

```console
chaff \
    --input tests/data/calls.vcf \
    --bam tests/data/tumor.bam \
    --ref tests/data/ref.fa \
    --sample tumor \
    --output calls.chaff.vcf \
    --metrics tumor.chaff.tsv \
    --copied-damage-threshold 0.05
```

Each filter writes the posterior probability that the call is a true mutation (`CDAP`, `ATAP`, `ERFAP`) and applies its FILTER at or below its threshold:

```console
grep -v '^#' calls.chaff.vcf | cut -f 1,2,4,5,7,8 | column -t
```

```text
chr1  100  C    A  CopiedDamageArtifact  CDAP=0.042;CDLR=1.307;CDAC=3,3;CDRC=120,240;ATAP=0.963;ERFAP=0.023
chr1  200  G    A  .                     CDAP=0.999;CDLR=-2.869;CDAC=10,20;CDRC=120,240;ATAP=1;ERFAP=1
chr1  300  AAA  A  .                     .
chr1  400  A    T  .                     ATAP=1;ERFAP=0.0033
chr1  500  C    G  .                     ERFAP=0.0033
```

All 3 alternate molecules at 100 sit nearer the lesion strand's 5' end (`CDAC=3,3`), where half of the reference molecules sit (`CDRC=120,240`).
The 20 alternate molecules at 200 split like the reference molecules (`CDAC=10,20`), so that call is likely a true mutation.
The deletion at 300 is not scored.

The header records each filter's model, prior, and threshold:

```console
grep '^##FILTER=<ID=CopiedDamage' calls.chaff.vcf
```

```text
##FILTER=<ID=CopiedDamageArtifact,Description="Call is likely damage copied onto both strands, with damage classes C>T,G>T and a 30 bp copy scale, at or below a posterior of 0.05.">
```

### Reading the Metrics

The `--metrics` file has one row per filter and stratum, with the stratum's and the filter's learned artifact fractions and a test of whether more alternate molecules sit where the artifact puts them than each call's own reference molecules predict:

```console
cut -f 2,3,6,7,16 tumor.chaff.tsv | head -3 | column -t
```

```text
filter         stratum      artifact_fraction  filter_artifact_fraction  asymmetry_p_value
copied-damage  C>T:non-CpG  0.443744           0.488011                  0.588099
copied-damage  G>T:CpG      0.530754           0.488011                  0.125
```

### Comparing to fgbio

With `--prior fgbio`, chaff writes fgbio 4.1.1's values and FILTERs wherever overlapping mates agree in base and quality, no read has an indel or soft clip between the call and its mate's 5' end, and every base is Q2 or better, as on this data:

```console
chaff \
    --input tests/data/calls.vcf \
    --bam tests/data/tumor.bam \
    --sample tumor \
    --output fgbio.vcf \
    --filters end-repair-fill-in,a-tailing \
    --end-repair-fill-in-threshold 0.001 \
    --prior fgbio
grep -v '^#' fgbio.vcf | cut -f 1,2,4,5,7,8 | column -t
```

```text
chr1  100  C    A  EndRepairFillInArtifact  ATAP=0.003732;ERFAP=0.00003218
chr1  200  G    A  .                        ATAP=1;ERFAP=1
chr1  300  AAA  A  .                        .
chr1  400  A    T  EndRepairFillInArtifact  ATAP=0.715;ERFAP=0.00001239
chr1  500  C    G  EndRepairFillInArtifact  ERFAP=0.00001239
```

The default learned prior leaves 100, 400, and 500 at `ERFAP=0.023`, `0.0033`, and `0.0033`, above 0.001: their 3 to 5 alternate molecules sit within 15 bp of a template end, but so do 37.5% of the reference molecules.

## Choosing Filters

Each filter models one library-preparation step, so turn on the ones your preparation has.

| Filter | When the preparation | Tune |
| --- | --- | --- |
| `end-repair-fill-in` | blunts fragment ends with end repair after shearing or enzymatic fragmentation | `--end-repair-fill-in-distance` |
| `a-tailing` | A-tails ends for T-overhang adapter ligation | `--a-tailing-distance` |
| `copied-damage` | has a polymerase work on double-stranded DNA before the strands are tagged: end-repair fill-in, nick translation, or gap filling | `--copied-damage-classes` (`C>T` for deamination from heat, storage, or formalin; `G>T` for oxidation from shearing or heat) and `--copied-damage-scale` |

A sheared or enzymatically fragmented, end-repaired, A-tailed duplex library with suspected deamination:

```console
chaff \
    --input tests/data/calls.vcf \
    --bam tests/data/tumor.bam \
    --ref tests/data/ref.fa \
    --sample tumor \
    --output calls.chaff.vcf.gz \
    --metrics tumor.chaff.tsv \
    --filters copied-damage,end-repair-fill-in,a-tailing \
    --copied-damage-classes 'C>T' \
    --copied-damage-threshold 0.05 \
    --end-repair-fill-in-threshold 0.001 \
    --a-tailing-threshold 0.001
```

The distances, scale, read floors, and prior are left at their defaults.
An option of a filter that `--filters` leaves out is a usage error, and so is `--ref` without `copied-damage`.
A VCF that already declares an enabled filter's INFO or FILTER, from an earlier run of chaff or fgbio, is refused, so a FILTER never outlives the run that applied it; remove them first, as with `bcftools annotate -x`.

## Filters

| Filter | INFO | FILTER | Scores |
| --- | --- | --- | --- |
| `copied-damage` | `CDAP`, `CDLR`, `CDAC`, `CDRC` | `CopiedDamageArtifact` | heterozygous SNVs in a damage class |
| `a-tailing` | `ATAP` | `ATailingArtifact` | heterozygous SNVs to `A` or `T` |
| `end-repair-fill-in` | `ERFAP` | `EndRepairFillInArtifact` | heterozygous SNVs |

- `copied-damage`: a copy reaches distance `d` from the lesion strand's 5' end with probability `w(d) = exp(-d / s)`, so `LLR = Σ ln((1 - e) w(d) / W + e)` over the alternate molecules, with `W` the mean `w(d)` of the reference molecules and `e` the base error. `CDLR` is the LLR in log10 units, and `CDAC` and `CDRC` count the molecules nearer the lesion strand's 5' end, out of all measured.
- `a-tailing` and `end-repair-fill-in`: fgbio's windowed likelihoods. `--end-repair-fill-in-scale` replaces the window with the decay above, from the nearest template end, in the posterior.

## Priors

- `learned` (default): EM learns each filter's artifact fraction `π_f` per sample from all its calls, `π_f = (Σ r_i + 1) / (n + 2)` with `r_i = σ(LLR_i + logit π_f)`, then each stratum's `π` from its own calls and 10 pseudo-calls at `π_f`, `π = (Σ r_i + 10 π_f) / (n + 10)`, so a stratum of one or two calls mostly inherits `π_f`. Strata are the damage class and CpG context for `copied-damage`, and the six pyrimidine substitution classes for the others.
- `fgbio`: fgbio's per-call mutation prior, `min((2 * maf)^2, 0.9999)`.

## Differences From fgbio

chaff matches fgbio where fgbio's choices are arbitrary: a deletion at the site counts in the depth of its prior, a spanning deletion `*` is no called allele, and an A-tailing site equally far from both template ends is nearer the kept read's own end.
It differs on purpose here:

- The learned prior: fgbio's mutation prior is near zero at duplex allele fractions, so alternate molecules inside the window make a call an artifact however many reference molecules sit there too.
- End repair extends a recessed 3' end across a 5' overhang; fgbio's docs describe filling a 3' overhang.
- Distances from both template ends count template bases: chaff walks both reads' CIGARs, the mate's from its `MC` tag, as [fgbio #1172](https://github.com/fulcrumgenomics/fgbio/pull/1172) does for clipping, so an indel counts by its length. Soft clips count and hard clips do not. fgbio measures the far end by insert size, and chaff never reads `TLEN`.
- Overlapping mates are called into one base: mates that agree keep the higher quality, and mates that disagree count as neither allele. fgbio keeps the first read of each name, so values differ where overlapping mates differ in base or quality.
- A base's error probability is capped at 0.75, a random base's, so a Q0 or Q1 base cannot zero a likelihood.
- A call with alternate but no reference molecules gets no INFO value; fgbio writes `NaN`.
- The BAM is always streamed, never queried by index.
- Values keep htsjdk's rounding but are written in decimal: `0.00003218` for fgbio's `3.218e-05`.

chaff builds on [Briggs et al. 2007](https://doi.org/10.1073/pnas.0704665104), [NanoSeq](https://doi.org/10.1038/s41586-021-03477-4), [Duplex-Repair](https://doi.org/10.1093/nar/gkab855), [fgbio](https://github.com/fulcrumgenomics/fgbio), and GATK's [`LearnReadOrientationModel`](https://gatk.broadinstitute.org/hc/en-us/articles/360037593911-LearnReadOrientationModel).

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
