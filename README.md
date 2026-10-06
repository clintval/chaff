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

## How Damage Becomes a Call

A DNA fragment is two complementary strands, each running 5′ to 3′.
Its *template ends* are its outermost bases, the 5′ ends of its two strands, and its *template* is the fragment as a read pair sees it.

```text
5′-CTGCCAGGATCC-3′   top strand
3′-GACGGTCCTAGG-5′   bottom strand
```

A *lesion* is a damaged base on one strand that a polymerase copies as another base:

- 5-methylcytosine deaminates to thymine, so a methylated CpG reads C>T.
- Cytosine deaminates to uracil, so any C can read C>T, but proofreading polymerases of the Pfu family stall at template uracil and often leave it uncopied.
- Guanine oxidizes to 8-oxoguanine, which pairs with A, so a G reads G>T, a C>A on the other strand.

Duplex sequencing tags both strands of a fragment with a UMI and keeps a base only where the two strands agree, so a lesion on one strand is outvoted by its partner.
End repair can defeat that.
Fragmentation leaves single-stranded overhangs, and end repair's polymerase extends each recessed 3′ end across the 5′ overhang opposite it, and from any nick, copying a lesion into the partner before the UMIs and adapters go on:

```text
      a C damaged to read as T
      v
5′-CTGTCAGGATCC-3′   top strand, overhanging at its 5′ end
3′-      CCTAGG-5′   bottom strand, recessed at its 3′ end
5′-CTGTCAGGATCC-3′
3′-GACAGTCCTAGG-5′   end repair fills the bottom strand in, copying the T as an A
      ^
```

Both strands now read T at the lesion, so the duplex consensus agrees on a C>T that looks real.
The polymerase copies from the partner's recessed end toward the lesion strand's 5′ end, so copied lesions sit near the 5′ end of the strand that carries them and rarely near its 3′ end, while a true mutation's molecules sit wherever the reference molecules at the site do.
This is *copied damage*.

A-tailing then adds a non-templated A to each 3′ end, for adapters with a T overhang to ligate to.
Where end repair over-digested a 3′ end, that A stands in for a lost base, so copies of the strand begin with a T where another base belongs: a T near the template's left end or, from the other strand, an A near its right end.

```text
5′-CTGCCAGGATCC-3′
3′- ACGGTCCTAGG-5′   end repair over-digests the bottom strand's 3′ end
3′-AACGGTCCTAGG-5′   A-tailing fills it with an A where a G belongs
```

Each chaff filter models one of these steps and writes its posterior probability that a call is a true mutation, applying its FILTER at or below a threshold you set:

| Filter | Models | INFO | FILTER |
| --- | --- | --- | --- |
| `copied-damage` | lesions copied onto the partner strand | `CDAP`, `CDLR`, `CDAC`, `CDRC` | `CopiedDamageArtifact` |
| `end-repair-fill-in` | errors on the strand end repair extends, near its 3′ end | `ERFAP` | `EndRepairFillInArtifact` |
| `a-tailing` | an A added to an over-digested 3′ end | `ATAP` | `ATailingArtifact` |

Each filter scores heterozygous SNVs: the copied damage filter those in its damage classes, `C>T` and `G>T` by default, fgbio's A-tailing filter those to A or T, and fgbio's end repair fill-in filter all of them.

## Which Filters to Use

All three filters run by default, and none applies its FILTER without a threshold.
Keep the ones your library preparation has:

| Filter | Keep it when the library |
| --- | --- |
| `end-repair-fill-in` | was end-repaired by a polymerase before adapter ligation, as most ligation preps after mechanical or enzymatic fragmentation are |
| `a-tailing` | was A-tailed for T-overhang adapters, unlike blunt-end ligation or transposase (tagmentation) preps |
| `copied-damage` | carries UMIs or duplex tags added after any polymerase fills ends, nicks, or gaps; check your prep's order of steps |

The copied damage filter needs the reference FASTA, `--ref`, for CpG context, and paired reads, whose template ends it measures.

The data say to keep a filter when somatic calls show an excess of C>T or C>A with few supporting molecules, when their alternate bases crowd one fragment end, or when the filter's metrics, below, learn an artifact fraction well above zero.
Suspect copied damage in old, stored, or degraded specimens, and when C>T calls crowd CpGs.

## Setting and Tuning

Most runs need only the filters and a threshold, since the prior is learned per sample.

### 1. Look

Profile the raw reads, before consensus, with the `error` tool of [Riker](https://github.com/fulcrumgenomics/riker), as `riker error -i raw.bam -r ref.fa -o raw`, whose default strata report the mismatch rate by cycle, by read number, and by 3 bp context.
A C>T or G>A excess rising toward read starts points to end repair fill-in, a C>A excess on one read number and not the other to 8-oxoguanine, and an A or T excess at read ends to A-tailing.
How far from read starts an excess reaches is a check on the distances chaff learns for end repair fill-in and copied damage.
Damage copied onto both strands before the UMIs went on reads as a real base to every read-level metric, so only the next step can see it.

### 2. Measure

Run the filters without thresholds, so chaff annotates calls without filtering them, and write the metrics:

```console
chaff \
    --input tests/data/calls.vcf \
    --bam tests/data/tumor.bam \
    --ref tests/data/ref.fa \
    --sample tumor \
    --output calls.annotated.vcf.gz \
    --metrics tumor.chaff.tsv
```

The metrics have one row per filter and *stratum*, the subset of calls a prior is learned in: the damage class and CpG context for `copied-damage`, and the substitution for the others.
Each row's *artifact fraction* is the share of its calls that are artifacts, learned from the calls themselves, its *distance* is the filter's decay scale or window in bases, learned for the decays, and its asymmetry p-value tests whether more alternate molecules sit within that distance of the artifact's end than each call's own reference molecules predict:

```console
cut -f 2,3,6,8,17 tumor.chaff.tsv | column -t
```

```text
filter              stratum      artifact_fraction  distance  asymmetry_p_value
copied-damage       C>T:non-CpG  0.454448           4.25027   0.745398
copied-damage       G>T:CpG      0.545324           4.25027   0.000289011
a-tailing           C>A          0.19272            2.0       0.00246167
a-tailing           C>T          0.189318           2.0       1.0
a-tailing           T>A          0.189332           2.0       0.000231639
end-repair-fill-in  C>A          0.393831           3.10398   0.000155019
end-repair-fill-in  C>G          0.302956           3.10398   1.0
end-repair-fill-in  C>T          0.302956           3.10398   0.668562
end-repair-fill-in  T>A          0.302956           3.10398   1.0
```

Every filter learns a fraction well above zero on these calls, which were built to carry A-tailing and end repair artifacts, so all three stay on; across many calls, a fraction near zero says a library lacks that artifact.
The counts behind a p-value are in the row:

```console
grep -e ^sample -e a-tailing tumor.chaff.tsv | cut -f 3,10,11,12,17 | column -t
```

```text
stratum  alt_molecules  alt_congruent  expected_alt_congruent  asymmetry_p_value
C>A      3              2              0.0867769               0.00246167
C>T      20             0              0.578512                1.0
T>A      5              3              0.144628                0.000231639
```

Three of the 5 alternate molecules of the A>T at position 400 sit where A-tailing puts them, against 0.14 expected, yet the other 2 cannot be A-tailing, so the call is spared: the p-value asks whether a library has the artifact, and the posterior whether one call is it.

### 3. Set and Check

A posterior is the probability that a call is a true mutation, given its molecules and the learned prior, so a threshold of 0.05 filters the calls with at most a 5% chance of being real.
Here each filter applies its FILTER at 0.05, and the output is BGZF-compressed by its extension:

```console
chaff \
    --input tests/data/calls.vcf \
    --bam tests/data/tumor.bam \
    --ref tests/data/ref.fa \
    --sample tumor \
    --output calls.chaff.vcf.gz \
    --copied-damage-threshold 0.05 \
    --end-repair-fill-in-threshold 0.05 \
    --a-tailing-threshold 0.05
gzip -dc calls.chaff.vcf.gz | grep -v '^#' | cut -f 2,4,5,7 | column -t
```

```text
100  C    A  CopiedDamageArtifact;EndRepairFillInArtifact
200  G    A  .
300  AAA  A  .
400  A    T  .
500  C    G  .
```

The posteriors behind those FILTERs are in the INFO of either run:

```console
gzip -dc calls.chaff.vcf.gz | grep -v '^#' | cut -f 2,4,5,8 | column -t
```

```text
100  C    A    CDAP=0.0003579;CDLR=3.367;CDAC=3,3;CDRC=15,240;ATAP=0.963;ERFAP=0.0003782
200  G    A    CDAP=1;CDLR=-47.16;CDAC=1,20;CDRC=15,240;ATAP=1;ERFAP=1
300  AAA  A    .
400  A    T    ATAP=1;ERFAP=1
500  C    G    ERFAP=1
```

The `CDAP`, `ATAP`, and `ERFAP` fields are the posteriors, `CDLR` is the log10 likelihood ratio of copied damage to a true mutation, and `CDAC` and `CDRC` count the alternate and reference molecules within the learned distance of the lesion strand's 5′ end, out of all measured.
All 3 alternate molecules of the C>A at position 100 sit within it, where 15 of the 240 reference molecules do, while 1 of the 20 alternate molecules of the G>A at position 200 does.
The deletion at position 300 is not scored.

To check a threshold, run chaff on germline heterozygous calls from the same reads: they are real, so the share it filters estimates how often it filters real somatic calls.
Where a matched normal or a replicate library exists, the somatic calls it shares are a second check.

## Models

The prior is the chance that a call is an artifact before its molecules are seen, and each filter's distance says how near its end an artifact sits.

- The chaff model, the default: each filter first learns its artifact fraction `π_f` from all its calls, `π_f = (Σ r_i + 1) / (n + 2)` with `r_i = σ(LLR_i + logit π_f)`, then each stratum learns `π = (Σ r_i + 10 π_f) / (n + 10)` from its own calls and 10 pseudo-calls at `π_f`, so a stratum of one or two calls mostly inherits `π_f`. End repair fill-in and copied damage decay with distance, and their scales are learned with `π_f`. In a clean library the fraction is near zero and real low-fraction calls are spared; where a library or stratum is damaged it is high, and the filter grows stricter there.
- The fgbio model, `--model fgbio`: a mutation prior of `min((2 * maf)^2, 0.9999)` from the call's alternate molecule fraction alone, so every low-fraction call is presumed an artifact whatever the library, and fgbio's windows from either template end. With duplex base qualities the posterior then nearly becomes a rule, an artifact whenever every alternate molecule sits inside the window. Use it to reproduce fgbio's values, or to compare with a pipeline built on fgbio.

The same calls under each model, fgbio's in the first three columns and chaff's in the last two:

```console
for model in fgbio chaff; do
    chaff \
        --input tests/data/calls.vcf \
        --bam tests/data/tumor.bam \
        --sample tumor \
        --output $model.vcf \
        --filters end-repair-fill-in,a-tailing \
        --end-repair-fill-in-threshold 0.001 \
        --model $model
done
paste fgbio.vcf chaff.vcf | grep -v '^#' | cut -f 2,7,8,18,19 | column -t
```

```text
100  EndRepairFillInArtifact  ATAP=0.003732;ERFAP=0.00003218  EndRepairFillInArtifact  ATAP=0.963;ERFAP=0.0003782
200  .                        ATAP=1;ERFAP=1                  .                        ATAP=1;ERFAP=1
300  .                        .                               .                        .
400  EndRepairFillInArtifact  ATAP=0.715;ERFAP=0.00001239     .                        ATAP=1;ERFAP=1
500  EndRepairFillInArtifact  ERFAP=0.00001239                .                        ERFAP=1
```

The 3 to 5 alternate molecules of the calls at positions 100, 400, and 500 sit within 3 bp of a template end, and fgbio's window filters all three; its prior presumes each low-fraction call an artifact.
The templates here are F1R2, copied from the forward strand, whose 3′ end, the one end repair extends, is the rightmost: the chaff model filters the call at position 100, whose alternate molecules sit there, and spares those at positions 400 and 500, whose alternate molecules sit at the leftmost end.
Because the chaff model's prior is fitted to each library, a threshold weighs a call against that library's own artifact rate and means much the same across samples; under fgbio's prior it moves with each call's allele fraction.

With `--model fgbio`, chaff writes fgbio 4.1.1's values and FILTERs wherever overlapping mates agree in base and quality, no read has an indel or soft clip between the call and its mate's 5′ end, and every base is Q2 or better, as on this data.

## Options

| Option | Default | Sets |
| --- | --- | --- |
| `--filters` | all three | the filters to run |
| `--model` | `chaff` | the model, `chaff` or `fgbio` |
| `--copied-damage-threshold`, `--end-repair-fill-in-threshold`, `--a-tailing-threshold` | none | the posterior at or below which a filter applies its FILTER |
| `--copied-damage-classes` | `C>T,G>T` | the damage classes, damaged base `>` read base: `C>T` for deamination from heat, storage, or formalin; `G>T` for oxidation from shearing or heat |
| `--copied-damage-distance` | `learned` | the decay scale in bases from the lesion strand's 5′ end, the mean length a polymerase copies a lesion strand over |
| `--end-repair-fill-in-distance` | `learned`, or 15 under `--model fgbio` | the decay scale in bases from the 3′ end of the strand each template was copied from, or under `--model fgbio` the window from either template end |
| `--a-tailing-distance` | 2 | the window from the template end, in bases |
| `--min-mapping-quality`, `--min-base-quality` | 20 | the read and base floors |
| `--paired-reads-only` | off | keep only reads whose mate is mapped |

The VCF/BCF and the BAM must be coordinate sorted, and neither needs an index.
Each template counts once, and a read whose mate maps to the same contig needs the mate's CIGAR in its `MC` tag.
Both ends of a template are measured for an FR pair, whose forward read starts at or before its reverse read's 5′ end; a read of any other pair knows only its own end.
An option of a filter that `--filters` leaves out is a usage error, and so is `--ref` without `copied-damage`.
A VCF that already declares an enabled filter's INFO or FILTER, from an earlier run of chaff or fgbio, is refused, so a FILTER never outlives the run that applied it; remove them first, as with `bcftools annotate -x`.

## Likelihoods

- Copied damage, and end repair fill-in under the chaff model: a copy reaches distance `d` from its end with probability `w(d) = exp(-d / s)`, so `LLR = Σ ln((1 - e) w(d) / W + e)` over the alternate molecules, with `W` the mean `w(d)` of the reference molecules and `e` the base error. The scale `s` maximizes the filter's marginal likelihood with `π_f` solved exactly at each scale, under a log-normal prior centered on 30 bp for copied damage and 15 bp for end repair fill-in. A call without both a measured reference and a measured alternate molecule gets no posterior.
- A-tailing, and end repair fill-in under the fgbio model: fgbio's windowed likelihoods, which compare the alternate molecules inside the window with the share of reference molecules there.

## Differences From fgbio

chaff matches fgbio where fgbio's choices are arbitrary: a deletion at the site counts in the depth of its prior, a spanning deletion `*` is no called allele, and an A-tailing site equally far from both template ends is nearer the kept read's own end.
It differs on purpose here:

- The chaff model: fgbio's mutation prior is near zero at duplex allele fractions, so alternate molecules inside the window make a call an artifact however many reference molecules sit there too, and its window counts either template end, where a fill-in error sits only near the 3′ end of the strand a template was copied from.
- End repair extends a recessed 3′ end across a 5′ overhang; fgbio's docs describe filling a 3′ overhang.
- Distances from both template ends count template bases: chaff walks both reads' CIGARs, the mate's from its `MC` tag, as [fgbio #1172](https://github.com/fulcrumgenomics/fgbio/pull/1172) does for clipping, so an indel counts by its length. Soft clips count and hard clips do not. fgbio measures the far end by insert size, and chaff never reads `TLEN`.
- Overlapping mates are called into one base: mates that agree keep the higher quality, and mates that disagree count as neither allele. fgbio keeps the first read of each name, so values differ where overlapping mates differ in base or quality.
- A base's error probability is capped at 0.75, a random base's, so a Q0 or Q1 base cannot zero a likelihood.
- A call with alternate but no reference molecules gets no INFO value; fgbio writes `NaN`.
- The BAM is always streamed, never queried by index.
- Values keep htsjdk's rounding but are written in decimal: `0.00003218` for fgbio's `3.218e-05`.

The examples run on fgbio's `FilterSomaticVcf` test data in [`tests/data`](tests/data): five tumor/normal calls on `chr1` at positions 100 to 500 in `calls.vcf`, the tumor's reads in `tumor.bam`, with artifact signal at positions 100, 400, and 500, and the reference in `ref.fa`.

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.

chaff builds on [Briggs et al. 2007](https://doi.org/10.1073/pnas.0704665104), [NanoSeq](https://doi.org/10.1038/s41586-021-03477-4), [Duplex-Repair](https://doi.org/10.1093/nar/gkab855), GATK's [`LearnReadOrientationModel`](https://gatk.broadinstitute.org/hc/en-us/articles/360037593911-LearnReadOrientationModel), and the `FilterSomaticVcf` of [fgbio](https://github.com/fulcrumgenomics/fgbio), whose filters, likelihoods, and tests it ports.
