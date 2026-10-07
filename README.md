# chaff

[![Install with bioconda](https://img.shields.io/badge/Install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/chaff/README.html)
[![Anaconda Version](https://anaconda.org/bioconda/chaff/badges/version.svg)](http://bioconda.github.io/recipes/chaff/README.html)
[![Build Status](https://github.com/clintval/chaff/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/clintval/chaff/actions/workflows/ci.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/clintval/chaff/badge.svg?branch=main)](https://coveralls.io/github/clintval/chaff?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-rust-dea588.svg)](https://www.rust-lang.org/)

Separate somatic variant calls from library-preparation damage artifacts.

Install with mamba, conda, or run directly with pixi:

```bash
pixi exec \
    -c conda-forge -c bioconda \
    chaff --help
```

## How DNA Damage Becomes a Variant Call

A DNA fragment is two complementary strands, each running 5′ to 3′.
Its *template ends* are its outermost bases, the 5′ ends of its two strands, and its *template* is the fragment as a read pair sees it.

A *lesion* is a damaged base on one strand that a polymerase copies as another base:

- **5-methylcytosine deaminates to thymine**, so a methylated CpG reads C>T.
- **Cytosine deaminates to uracil**, so any C can read C>T, but proofreading polymerases of the Pfu family stall at template uracil and often leave it uncopied.
- **Guanine oxidizes to 8-oxoguanine**, which pairs with A, so a G reads G>T, a C>A on the other strand.

Duplex Sequencing tags both strands of a fragment with a UMI and keeps a base only where the two strands agree.
A duplex consensus therefore removes an error on one strand, such as a lesion left uncopied, but not damage that library preparation copied onto both strands before the UMI-bearing adapters were ligated.

## The Filters

Each filter models one library-preparation step that leaves an artifact near a known end of the fragment, shown here on one fragment with the lesion and its copy, a misincorporated base, and an added A in red, and the bases end repair filled in in blue:

![One fragment with the end each artifact sits near: copied damage, a lesion marked on its base, near the lesion strand's 5′ end, end repair fill-in errors in the new bases near the extended strand's 3′ end, and an added A at the last base of a 3′ end.](.github/img/reference-points.svg)

The length a polymerase fills in varies from fragment to fragment, so the evidence for copied damage and end repair fill-in fades with distance from the end, with no cliff at any one distance, while A-tailing changes only the last base or two of a 3′ end and so is scored within a window.
All three filters run by default, but their thresholds default to none, so a filter annotates calls and applies no FILTER until you give it a threshold.
Each filter scores heterozygous SNVs: copied damage those in its damage classes, C>T and G>T by default, A-tailing those to A or T, and end repair fill-in all of them.

### Copied Damage

Fragmentation leaves single-stranded overhangs, and end repair's polymerase extends each recessed 3′ end across the 5′ overhang opposite it, and from any nick, copying a lesion onto the partner strand before the UMI-bearing adapters are ligated.
Both strands then read the damaged base, so the duplex consensus agrees on a change that looks real.
This is the filter for Duplex Sequencing, since a duplex consensus cannot remove what both strands carry.

![A lesion, a C deaminated to uracil that reads as T, sits near one strand's 5′ end; end repair fill-in copies it onto the partner strand as an A; UMI-bearing adapters are ligated after the copy; and the duplex consensus of both strands agrees on a C>T.](.github/img/copied-damage.svg)

- **Measured from:** the 5′ end of the lesion strand, the strand that carries the reference C of a C>T or the reference G of a G>T.
- **Scored with:** a decay whose scale is learned per library.
- **Writes:** `CDAP`, the posterior probability that the call is a real mutation rather than damage copied onto both strands; `CDLR`, the log10 likelihood ratio of copied damage to a real mutation; `CDAC` and `CDRC`, the alternate and reference molecules within the decay scale of the lesion strand's 5′ end, and all of each measured; and the FILTER `CopiedDamageArtifact`.
- **Use it when:** the UMI-bearing adapters are ligated after any polymerase fills in ends, nicks, or gaps, as in Duplex Sequencing; it needs the reference FASTA, `--ref`, for CpG context.

The copy runs from the partner's recessed end toward the lesion strand's 5′ end, so copied lesions sit near that end, while a real mutation's molecules sit wherever the reference molecules at the site do.
On a simulated duplex sample of 6,000 real mutations and 2,000 copied-damage calls at CpG C>T, with about 400 duplex molecules per site, fragments of median length 200 bp, and fill-in that copies a lesion onto the partner strand from the lesion strand's 5′ end, over an exponential length with a mean of 30 bp, `chaff` learns a scale of 30.3 bp, and 65% of the copied damage's alternate molecules sit within it, against 17% of the real mutations' and of the reference molecules:

![Distances of alternate and reference molecules from the lesion strand's 5′ and 3′ ends: copied damage piles up near the 5′ end and avoids the 3′ end, while real mutations follow the reference molecules.](.github/img/copied-damage-ends.png)

Copied damage has one blind spot.
A lesion copied onto the partner strand from an internal nick, by nick translation or strand displacement, or across the gap an abasic site leaves, can sit anywhere in the template, so its alternate molecules carry no signal of an end: no per-call score can separate them from a real mutation's, and they show only as an excess of the damage class across a library.

### End Repair Fill-In

![End repair fill-in makes an error: a polymerase fills in a recessed 3′ end and misincorporates a C opposite a T, the UMI-bearing adapters are ligated, and only the filled-in strand carries the error, so the two strands disagree and a duplex consensus masks it.](.github/img/end-repair-fill-in.svg)

End repair's polymerase can misincorporate a base, or copy a damaged base from the overhang it fills, so the error sits on the strand it extended, near that strand's 3′ end.
Only that strand carries it, so a duplex consensus mostly removes it.

- **Measured from:** the 3′ end of the strand each template was copied from, the strand of its read 1: forward for an F1R2 pair and reverse for F2R1, as GATK's `LearnReadOrientationModel` reads them, or the nearer end of a duplex consensus, whose reads carry the `aD` and `bD` depths of both strands.
- **Scored with:** a decay whose scale is learned per library.
- **Writes:** `ERFAP`, the posterior probability that the call is a real mutation rather than an end repair fill-in artifact, and the FILTER `EndRepairFillInArtifact`.
- **Use it when:** a polymerase end-repaired the library before adapter ligation, as in most ligation preps after mechanical or enzymatic fragmentation, and its reads are not a duplex consensus.

### A-Tailing

![A-tailing makes an error: end repair trims a 3′ end one base too far, A-tailing adds a non-templated A where a C belongs, the UMI-bearing adapters are ligated, and only that strand reads A at the last base of its 3′ end, so the two strands disagree.](.github/img/a-tailing.svg)

A-tailing adds a non-templated A to each 3′ end, for adapters with a T overhang to ligate to.
Where end repair over-digested a 3′ end, that A stands in for a lost base, so copies of the strand begin with a T where another base belongs: a T near the template's left end or, from the other strand, an A near its right end.
Only one strand carries it, so a duplex consensus mostly removes it.

- **Measured from:** the template end where the added A reads, the left end for a T and the right end for an A.
- **Scored with:** a 2 bp window, since the artifact changes only the last base or two of a 3′ end.
- **Writes:** `ATAP`, the posterior probability that the call is a real mutation rather than an A-tailing artifact, and the FILTER `ATailingArtifact`.
- **Use it when:** the library was A-tailed for T-overhang adapters, unlike blunt-end ligation or transposase (tagmentation) preps, and its reads are not a duplex consensus.

### Choosing Filters for Your Library

Start from each filter's **Use it when** line, which follows from how your library was prepared, then let the data confirm the choice: keep a filter when somatic calls show an excess of C>T or C>A with few supporting molecules, when their alternate bases crowd one fragment end, or when the filter's metrics, below, learn an artifact fraction well above zero.
Suspect copied damage in old, stored, or degraded specimens, and when C>T calls crowd CpGs.

## What `chaff` Writes

Each filter writes its posterior into the INFO of every call it scores and applies its FILTER where the posterior is at or below its threshold, and copied damage also writes the molecule counts behind its posterior.
The examples run on the five calls at positions 100 to 500 in [`tests/data`](tests/data), with the tumor's reads and the reference beside them.
Here copied damage applies its FILTER at a posterior of 0.05, and the call at position 100 is split into its FILTER and its INFO values, one per line:

```console
chaff \
    --input tests/data/calls.vcf \
    --bam tests/data/tumor.bam \
    --ref tests/data/ref.fa \
    --sample tumor \
    --output annotated.vcf \
    --copied-damage-threshold 0.05
grep -v '^#' annotated.vcf | awk '$2 == 100' | cut -f 7,8 | tr '\t;' '\n\n'
```

```text
CopiedDamageArtifact
CDAP=0.022
CDLR=1.591
CDAC=3,3
CDRC=69,240
ATAP=0.963
ERFAP=0.007738
```

Its `CDAP` of 0.022, at or below the threshold of 0.05, puts the copied damage FILTER on the call, and its `CDLR` of 1.591 favors copied damage.
Its `CDAC` of 3,3 and `CDRC` of 69,240 say that all 3 alternate molecules sit within the learned scale of the lesion strand's 5′ end, against 69 of the 240 reference molecules.
Without thresholds of their own, the `ATAP` and `ERFAP` posteriors of A-tailing and end repair fill-in annotate the call without filtering it.

Each posterior is the probability that the call is a real mutation, so a lower value means a call more likely to be an artifact.
The tool learns per library how common each artifact is and how far it reaches, so most runs need only the filters and a threshold.
The model behind the posteriors, and where it differs from fgbio on purpose, are described in the crate documentation, in [`src/lib/mod.rs`](src/lib/mod.rs).

## Setting and Tuning

Tuning a run comes down to three steps: look at the raw reads, measure the calls, then set and check a threshold.

### 1. Look

Profile the raw reads, before consensus, with the `error` tool of [Riker](https://github.com/fulcrumgenomics/riker), as `riker error -i raw.bam -r ref.fa -o raw`, whose default strata report the mismatch rate by cycle, by read number, and by 3 bp context.
A C>T or G>A excess rising toward read starts points to end repair fill-in, a C>A excess on one read number and not the other to 8-oxoguanine, and an A or T excess at read ends to A-tailing.
How far from read starts an excess reaches is a check on the decay scales `chaff` learns.
Damage copied onto both strands before the UMI-bearing adapters were ligated reads as a real base to every read-level metric, so only the next step can see it.

### 2. Measure

Run the filters without thresholds, so `chaff` annotates calls without filtering them, and write the metrics:

```console
chaff \
    --input tests/data/calls.vcf \
    --bam tests/data/tumor.bam \
    --ref tests/data/ref.fa \
    --sample tumor \
    --output calls.annotated.vcf.gz \
    --metrics tumor.chaff.tsv
```

The metrics have one row per filter and *stratum*, a group of calls that share an artifact rate: the damage class and CpG context for copied damage, and the substitution for the others.
In each row, the artifact fraction is the share of the stratum's calls that look like artifacts, the distance is how far from the artifact's end the filter looks, in bases, and the asymmetry p-value is small when more alternate molecules sit within that distance than the reference molecules at the same sites predict:

```console
cut -f 2,3,6,8,17 tumor.chaff.tsv | column -t
```

```text
filter              stratum      artifact_fraction  distance  asymmetry_p_value
copied-damage       C>T:non-CpG  0.448724           22.7678   0.72908
copied-damage       G>T:CpG      0.53767            22.7678   0.0242018
a-tailing           C>A          0.19272            2.0       0.00246167
a-tailing           C>T          0.189318           2.0       1.0
a-tailing           T>A          0.189332           2.0       0.000231639
end-repair-fill-in  C>A          0.391718           12.1245   0.00451579
end-repair-fill-in  C>G          0.301512           12.1245   1.0
end-repair-fill-in  C>T          0.301512           12.1245   0.665409
end-repair-fill-in  T>A          0.301512           12.1245   1.0
```

Every filter learns a fraction well above zero on these calls, which were built with alternate molecules near template ends, so all three stay on; across many calls, a fraction near zero says a library lacks that artifact.
The decays learn scales of 22.8 and 12.1 bp: the alternate molecules here sit within 3 bp of an end, but 2 and 4 calls move a scale only part of the way from its default.
A small p-value says that a library has the artifact, while a posterior says whether one call is it.

### 3. Set and Check

A threshold of 0.05 filters the calls with at most a 5% chance of being real.
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

The call at position 100 is filtered as copied damage and end repair fill-in, the other calls pass, and the deletion at position 300 is not scored.

To check a threshold, run `chaff` on germline heterozygous calls from the same reads: they are real, so the share it filters estimates how often it filters real somatic calls.
Where a matched normal or a replicate library exists, the somatic calls it shares are a second check.

On the simulated sample of the copied damage section, with a quarter of each kind of call at 2, 3, 5, or 10 alternate molecules, a threshold of 0.05 filters 1,427 of the 2,000 copied-damage calls and 23 of the 1,065 real C>T at CpG, and none of the 4,935 calls in other channels.
That threshold filters 95% of the copied damage with 10 alternate molecules and 86% with 5, but only 40% with 2, where it also filters 6 of 253 real C>T at CpG:

![The SBS96 spectrum of the simulated sample before and after `chaff`, where the copied damage at CpG C>T mostly leaves and the other channels stay, and the share of copied damage filtered against the share of real C>T at CpG filtered, by alternate molecules per call and model.](.github/img/copied-damage-filtering.png)

## Options

| Option&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Sets |
| --- | --- |
| `--filters` | The filters to run (default all three). |
| `--model` | The model, `chaff`, or `fgbio` for reproducing fgbio's values (default `chaff`). |
| `--copied-damage-threshold` | The posterior at or below which copied damage applies its FILTER (default none). |
| `--end-repair-fill-in-threshold` | The posterior at or below which end repair fill-in applies its FILTER (default none). |
| `--a-tailing-threshold` | The posterior at or below which A-tailing applies its FILTER (default none). |
| `--copied-damage-classes` | The damage classes, damaged base `>` read base: `C>T` for deamination from heat, storage, or formalin, and `G>T` for oxidation from shearing or heat (default `C>T,G>T`). |
| `--copied-damage-distance` | The decay scale in bases from the lesion strand's 5′ end, the mean length over which a polymerase copies a lesion strand onto its partner, or `learned` (default `learned`). |
| `--end-repair-fill-in-distance` | The decay scale in bases from the 3′ end of the strand each template was copied from, or `learned` (default `learned`). |
| `--a-tailing-distance` | The window from the template end, in bases (default 2). |
| `--min-mapping-quality` | The mapping quality floor of a read (default 20). |
| `--min-base-quality` | The base quality floor at the call (default 20). |
| `--paired-reads-only` | Keep only reads whose mate is mapped (default off). |

The VCF/BCF and the BAM must be coordinate sorted, and neither needs an index.
Each template counts once, and a read whose mate maps to the same contig needs the mate's CIGAR in its `MC` tag.
Both ends of a template are measured for an FR pair, whose forward read starts at or before its reverse read's 5′ end; a read of any other pair knows only its own end.
An option of a filter that `--filters` leaves out is a usage error, and so is `--ref` without `copied-damage`.
A VCF that already declares an enabled filter's INFO or FILTER, from an earlier run, is refused, so a FILTER never outlives the run that applied it; remove them first, as with `bcftools annotate -x`.

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
For compatibility with fgbio's `FilterSomaticVcf`, end repair fill-in and A-tailing write its INFO keys and FILTER names, `--model fgbio` reproduces its values, and the README examples run on its test data.

## References

The tool `chaff` ports the filters, likelihoods, and tests of fgbio's `FilterSomaticVcf`, and builds on these papers and tools:

1. Briggs AW, et al. 2007. Patterns of damage in genomic DNA sequences from a Neandertal. *Proceedings of the National Academy of Sciences* 104(37):14616–14621. [https://doi.org/10.1073/pnas.0704665104](https://doi.org/10.1073/pnas.0704665104)
2. Schmitt MW, et al. 2012. Detection of ultra-rare mutations by next-generation sequencing. *Proceedings of the National Academy of Sciences* 109(36):14508–14513. [https://doi.org/10.1073/pnas.1208715109](https://doi.org/10.1073/pnas.1208715109)
3. Abascal F, et al. 2021. Somatic mutation landscapes at single-molecule resolution. *Nature* 593(7859):405–410. [https://doi.org/10.1038/s41586-021-03477-4](https://doi.org/10.1038/s41586-021-03477-4)
4. Xiong K, et al. 2022. Duplex-Repair enables highly accurate sequencing, despite DNA damage. *Nucleic Acids Research* 50(1):e1. [https://doi.org/10.1093/nar/gkab855](https://doi.org/10.1093/nar/gkab855)
5. GATK's `LearnReadOrientationModel`: [https://gatk.broadinstitute.org/hc/en-us/articles/360057439111-LearnReadOrientationModel](https://gatk.broadinstitute.org/hc/en-us/articles/360057439111-LearnReadOrientationModel)
6. fgbio's `FilterSomaticVcf`: [https://fulcrumgenomics.github.io/fgbio/tools/latest/FilterSomaticVcf.html](https://fulcrumgenomics.github.io/fgbio/tools/latest/FilterSomaticVcf.html), from [https://github.com/fulcrumgenomics/fgbio](https://github.com/fulcrumgenomics/fgbio)
