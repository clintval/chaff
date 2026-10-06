# chaff

[![Build Status](https://github.com/clintval/chaff/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/clintval/chaff/actions/workflows/ci.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/clintval/chaff/badge.svg?branch=main)](https://coveralls.io/github/clintval/chaff?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-rust-dea588.svg)](https://www.rust-lang.org/)

Separate somatic variant calls from library-preparation damage artifacts.

## Introduction

Duplex and UMI sequencing suppress most sequencing and amplification errors, but not damage that library preparation copies onto both strands of a molecule before the strands are tagged.
The tool `chaff` scores each somatic call by where its alternate allele sits on the molecules that carry it, compared with the molecules that carry the reference allele at the same site.
It reimplements fgbio's [`FilterSomaticVcf`](https://fulcrumgenomics.github.io/fgbio/tools/latest/FilterSomaticVcf.html) without the JVM and adds a strand-aware filter for damage that polymerase copied onto both strands.

<details>
<summary>
Its filters build on prior work.
</summary>

<br>

- [Briggs et al. 2007](https://doi.org/10.1073/pnas.0704665104): damage in single-stranded overhangs is copied by end repair, so C>T sits at 5' ends and G>A at 3' ends
- [Abascal et al. 2021](https://doi.org/10.1038/s41586-021-03477-4) (NanoSeq): strand assignment by the nearest 5' end, and a binomial test of the asymmetry
- [Xiong et al. 2022](https://doi.org/10.1093/nar/gkab855) (Duplex-Repair): resynthesis from 3' ends copies damage into the complementary strand
- [Jiang et al. 2020](https://doi.org/10.1101/gr.261396.120): jagged single-stranded ends of plasma DNA, and how far end-repair fill-in reaches into them
- [Chen et al. 2017](https://doi.org/10.1126/science.aai8690): oxidative damage (8-oxoguanine) as a widespread source of G>T artifacts
- [fgbio `FilterSomaticVcf`](https://github.com/fulcrumgenomics/fgbio): the end repair fill-in (`ERFAP`) and A-tailing (`ATAP`) filters
- [GATK `LearnReadOrientationModel`](https://gatk.broadinstitute.org/hc/en-us/articles/360037593911-LearnReadOrientationModel): learning artifact priors from all calls with EM

</details>

Install from source:

```bash
cargo install --git https://github.com/clintval/chaff
```

## Quick Start

Annotate the calls of one sample, apply the lesion copy FILTER at a posterior of 0.05 or below, and write the per-sample metrics:

```bash
chaff filter \
    --input "calls.vcf.gz" \
    --bam "tumor.bam" \
    --ref "reference.fa" \
    --sample "tumor" \
    --output "calls.chaff.vcf.gz" \
    --metrics "tumor.chaff.tsv" \
    --lesion-copy-threshold 0.05
```

The VCF and the BAM must be coordinate sorted; neither needs an index, because `chaff` merge-joins them in one stream.
Every template counts once: overlapping mates become one molecule, and duplicate, secondary, and supplementary reads are left out.

## Filters

Each filter compares the alternate molecules of a call with its reference molecules.
Under a true mutation both alleles sit on the molecules the same way, so the reference molecules at a site calibrate the null and absorb capture and fragment-length skew.
Under an artifact the alternate molecules crowd the template end that made the artifact.
Each call gets a likelihood ratio of artifact to mutation, and the INFO field reports the posterior probability that the call is a true mutation: lower values mean a likely artifact, as in fgbio.

| Filter | INFO | FILTER | Applies to |
| --- | --- | --- | --- |
| `lesion-copy` | `LCAP`, `LCLR`, `LCAC`, `LCRC` | `LesionCopyArtifact` | heterozygous SNVs in a lesion class |
| `a-tailing` | `ATAP` | `ATailingArtifact` | heterozygous SNVs to `A` or `T` |
| `end-repair-fill-in` | `ERFAP` | `EndRepairFillInArtifact` | heterozygous SNVs |

A FILTER is applied only with a threshold (`--lesion-copy-threshold`, `--a-tailing-threshold`, `--end-repair-fill-in-threshold`), at or below it.

###### Lesion Copy

A lesion on one strand, such as a deaminated cytosine or an 8-oxoguanine, is templated into the other strand when polymerase resynthesizes it by end-repair fill-in, nick translation, or gap filling before the strands are tagged.
Both strands then carry the change, and duplex consensus agrees on it.
Resynthesis runs 5' to 3' along the new strand, so copies sit near the 5' end of the lesion strand and are depleted near its 3' end, over tens to more than a hundred bases.

The lesion strand comes from the substitution class (`--lesion-copy-classes`, default `C>T,G>T`): `C>T` puts the lesion on the strand carrying the reference `C`, so a forward-strand `C>T` and a reverse-strand `G>A` are the same class; `G>T` puts it on the strand carrying the reference `G`.
For each molecule with both template ends known, `d` is the distance of the site from the lesion strand's 5' end.
A lesion at distance `d` is copied with probability `w(d) = exp(-d / s)`, an exponential resynthesis length with mean `s` (`--lesion-copy-scale`, default 30 bp).
The artifact's alternate distances follow the reference distances tilted by `w`, so with `W` the mean of `w` over the reference molecules and `e` an alternate base's error probability, the log likelihood ratio is

```
LLR = sum over alternate molecules of ln((1 - e) * w(d) / W + e)
```

`LCLR` reports it in log10 units; `LCAC` and `LCRC` count the alternate and reference molecules nearer the lesion strand's 5' end than its 3' end, out of all measured.
Calls are stratified by class and by CpG context from the reference.

###### End Repair Fill-in

End repair blunts a fragment: polymerase extends a recessed 3' end across the opposite strand's 5' overhang, and an exonuclease trims a 3' overhang.
Damage in the single-stranded 5' overhang is copied into the extended strand, so after amplification both strands carry the change near a template end.
A molecule is congruent when the site is within `--end-repair-fill-in-distance` (default 15) of its nearest template end, and the likelihoods are fgbio's.
With `--end-repair-fill-in-scale` the window becomes the decay model above, measured from the nearest end.

###### A-tailing

End repair can over-digest a 3' end and leave it recessed; A-tailing then fills it with adenines.
On the forward strand that is a `T` within `--a-tailing-distance` (default 2) of the leftmost template end or an `A` within it of the rightmost one.
The likelihoods are fgbio's.

## Priors

By default the prior is learned.
Within each sample and stratum (lesion class and CpG context for `lesion-copy`, the six pyrimidine substitution classes for the others), the calls form a two-component mixture with an unknown artifact fraction `pi`, estimated by expectation-maximization as GATK's `LearnReadOrientationModel` learns its priors:

```
E-step: r_i = 1 / (1 + exp(-(LLR_i + logit(pi))))
M-step: pi  = (sum r_i + 1) / (n + 2)
```

The `+1` and `+2` are a Beta(2, 2) prior that keeps `pi` inside (0, 1) when a stratum holds few calls.
The posterior probability of a true mutation is then `1 / (1 + exp(LLR_i + logit(pi)))`.

`--prior fgbio` uses fgbio's per-call mutation prior, `min((2 * maf)^2, 0.9999)`, and reproduces fgbio's `ERFAP` and `ATAP` values.

## Differences From fgbio

- **The prior.** fgbio's `(2 * maf)^2` prior is near zero at duplex allele fractions, so any call whose alternate molecules all sit inside the window becomes an artifact however often the reference molecules sit there too. The learned prior needs the molecules themselves to carry the evidence. `--prior fgbio` restores fgbio's.
- **The end repair mechanism.** fgbio describes fill-in of single-stranded 3' overhangs. Polymerase cannot extend a 3' overhang; it extends a recessed 3' end opposite a 5' overhang, and the damage copied is in that 5' overhang.
- **Template ends.** fgbio measures from the read's own 5' end and from the far end by insert size. `chaff` uses the unclipped 5' ends of both mates, the mate's from the `MC` tag, which a read with a mapped mate must carry; it never reads `TLEN`.
- **Overlapping mates.** fgbio keeps the first read of each name. `chaff` keeps one molecule per template: mates that agree keep the higher quality, and mates that disagree count as neither allele.
- **Missing evidence.** fgbio writes `NaN` when a call has alternate molecules but no reference molecules, and a posterior from `1 / depth` when it has neither. `chaff` leaves the INFO field out in the first case and reports the prior in the second.
- **Access.** fgbio can query an indexed BAM; `chaff` always streams.
- **Number formatting.** Values match fgbio's to the four significant digits fgbio prints, but are written in decimal rather than scientific notation.

## Metrics

`--metrics` writes one row per filter and stratum: calls, filtered calls, the learned artifact fraction, the expected number of artifact calls, and the alternate and reference molecules congruent with the artifact (for `lesion-copy`, nearer the lesion strand's 5' end).
The `asymmetry_p_value` is a one-sided binomial test of the congruent alternate molecules against the congruent fraction of the reference molecules, NanoSeq's test with the reference molecules in place of a fixed one half.

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
