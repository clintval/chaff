Feature: CLI contract
  `chaff` scores somatic calls against library-preparation artifacts. These
  scenarios hold at the CLI contract layer.

  Background:
    Given a file "calls.vcf" containing:
      """
      ##fileformat=VCFv4.2
      ##contig=<ID=chr1,length=1000>
      ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\tnormal
      chr1\t100\t.\tC\tA\t.\tPASS\t.\tGT\t0/1\t0/0
      """
    And a BAM "sorted.bam" sorted by "coordinate"
    And a BAM "unsorted.bam" sorted by "unsorted"

  Scenario: reports its version
    When I run `chaff --version`
    Then the exit code is 0
    And stdout contains "chaff"

  Scenario: documents every filter and the prior
    When I run `chaff --help`
    Then the exit code is 0
    And stdout contains "copied-damage"
    And stdout contains "end-repair-fill-in"
    And stdout contains "a-tailing"
    And stdout contains "--prior"

  Scenario: the input, output, and BAM are required
    When I run `chaff`
    Then the exit code is 2

  Scenario: the BAM is required
    When I run `chaff -i calls.vcf -o out.vcf`
    Then the exit code is 2
    And stderr contains "--bam <BAM>"

  Scenario: a damage class and its reverse complement cannot both be given
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --copied-damage-classes C>T,G>A`
    Then the exit code is 2
    And stderr contains "error: damage classes C>T and G>A describe the same change on opposite strands"

  Scenario: thresholds are probabilities
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --a-tailing-threshold 2`
    Then the exit code is 2
    And stderr contains "expected a probability from 0 to 1"

  Scenario: fgbio's p-value flag names are accepted
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --filters a-tailing --a-tailing-p-value 0.01 --sample tumor`
    Then the exit code is 0

  Scenario: a multi-sample VCF needs --sample
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --filters a-tailing`
    Then the exit code is 1
    And stderr contains "--sample must name the one whose reads are in the BAM"

  Scenario: an unknown sample is an error
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --filters a-tailing --sample WhoDis`
    Then the exit code is 1
    And stderr contains "WhoDis"

  Scenario: the copied damage filter needs a reference
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --sample tumor`
    Then the exit code is 1
    And stderr contains "needs a reference FASTA"

  Scenario: the BAM must be coordinate sorted
    When I run `chaff -i calls.vcf -o out.vcf -b unsorted.bam --filters a-tailing --sample tumor`
    Then the exit code is 1
    And stderr contains "must be coordinate sorted"

  Scenario: the input VCF cannot be standard input
    When I run `chaff -i - -o out.vcf -b sorted.bam --filters a-tailing`
    Then the exit code is 1
    And stderr contains "must be a file"

  Scenario: an option of a filter left out of --filters is an error
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --sample tumor --filters a-tailing --copied-damage-threshold 0.05`
    Then the exit code is 2
    And stderr contains "error: the argument '--copied-damage-threshold <P>' applies only to the copied-damage filter, which '--filters' leaves out"

  Scenario: the end repair fill-in scale replaces the window only in the posterior
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --sample tumor --filters end-repair-fill-in --end-repair-fill-in-distance 10 --end-repair-fill-in-scale 15`
    Then the exit code is 0
    When I run `chaff --help`
    Then stdout contains "The distance window still sets the congruent molecules in the metrics."
