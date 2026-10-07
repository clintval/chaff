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

  Scenario: documents every filter and the model
    When I run `chaff --help`
    Then the exit code is 0
    And stdout contains "copied-damage"
    And stdout contains "end-repair-fill-in"
    And stdout contains "a-tailing"
    And stdout contains "--model"

  Scenario: the input, output, and BAM are required
    When I run `chaff`
    Then the exit code is 2

  Scenario: the BAM is required
    When I run `chaff -i calls.vcf -o out.vcf`
    Then the exit code is 2
    And stderr contains "--bam <BAM>"

  Scenario: a damage class and its reverse complement cannot both be given
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam -r ref.fa --copied-damage-classes C>T,G>A`
    Then the exit code is 2
    And stderr contains "error: damage classes C>T and G>A describe the same change on opposite strands"

  Scenario: thresholds are probabilities
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --a-tailing-threshold 2`
    Then the exit code is 2
    And stderr contains "expected a probability from 0 to 1"

  Scenario: fgbio's threshold option names are accepted as hidden aliases
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --filters a-tailing --a-tailing-p-value 0.01 --sample tumor`
    Then the exit code is 0

  Scenario: an output of - writes VCF to standard output
    When I run `chaff -i calls.vcf -o - -b sorted.bam --sample tumor --filters a-tailing`
    Then the exit code is 0
    And stdout contains "##INFO=<ID=ATAP,"
    And stdout contains "chr1\t100\t.\tC\tA\t.\tPASS\t.\tGT\t0/1\t0/0"

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
    Then the exit code is 2
    And stderr contains "error: the copied-damage filter needs a reference FASTA: '--ref <FASTA>'"

  Scenario: a REF that differs from the reference FASTA is an error
    Given a file "ref.fa" containing:
      """
      >chr1
      GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
      """
    And a file "ref.fa.fai" containing:
      """
      chr1\t120\t6\t120\t121
      """
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam -r ref.fa --sample tumor`
    Then the exit code is 1
    And stderr contains "the call at chr1:100 has REF C, but the reference FASTA has G there"

  Scenario: the BAM must be coordinate sorted
    When I run `chaff -i calls.vcf -o out.vcf -b unsorted.bam --filters a-tailing --sample tumor`
    Then the exit code is 1
    And stderr contains "must be coordinate sorted"

  Scenario: the input VCF cannot be standard input
    When I run `chaff -i - -o out.vcf -b sorted.bam --filters a-tailing`
    Then the exit code is 2
    And stderr contains "error: invalid value '-' for '--input <VCF>': the input is read twice, so it must be a file, not standard input"

  Scenario: an output that names the input is an error
    When I run `chaff -i calls.vcf -o ./calls.vcf -b sorted.bam --sample tumor --filters a-tailing`
    Then the exit code is 2
    And stderr contains "error: '--input' and '--output' name the same file"

  Scenario: an option of a filter left out of --filters is an error
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --sample tumor --filters a-tailing --copied-damage-threshold 0.05`
    Then the exit code is 2
    And stderr contains "error: the argument '--copied-damage-threshold <P>' applies only to the copied-damage filter, which '--filters' leaves out"

  Scenario: a reference without the copied damage filter is an error
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --sample tumor --filters a-tailing --ref ref.fa`
    Then the exit code is 2
    And stderr contains "error: the argument '--ref <FASTA>' applies only to the copied-damage filter, which '--filters' leaves out"

  Scenario: the spectrum needs a reference
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --sample tumor --filters a-tailing --spectrum spectrum.pdf`
    Then the exit code is 2
    And stderr contains "error: '--spectrum' needs a reference FASTA: '--ref <FASTA>'"

  Scenario: the spectrum is written as a PDF
    Given a file "ref.fa" containing:
      """
      >chr1
      AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGAAAAAAAAAAAAAAAAAAA
      """
    And a file "ref.fa.fai" containing:
      """
      chr1\t120\t6\t120\t121
      """
    When I run `chaff -i calls.vcf -o out.vcf -b sorted.bam --sample tumor --filters a-tailing --ref ref.fa --spectrum spectrum.pdf`
    Then the exit code is 0
    And the file "spectrum.pdf" is a PDF
