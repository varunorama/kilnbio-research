# Targeted methylome and exome capture

Use `kilnbio capture-panel` to measure depth over a capture design or to ask which measured sites or called variants lie inside it. Supply the **covered target BED for the kit and revision you used**. A methylome kit and an exome kit each have their own target file. Probes can tile or extend beyond their intended targets, so their optional BED has a separate checksum and size.

| Input | Result |
|---|---|
| Target BED and FASTA index or chromosome sizes | Merged, unpadded target intervals and number of unique target bases |
| Coordinate-sorted DNA BAM plus targets | Mean depth, zero-depth bases and fraction meeting the chosen minimum over **all** target bases. Requires SAMtools (tested with 1.24). |
| Called small-variant VCF | Record-level eligibility of PASS biallelic A/C/G/T alleles using the entire reference span. A missing variant is unknown. This does not assign genotypes or count false negatives. |
| Strand-collapsed KilnBio CpG counts | Whether both bases of each covered CpG lie within the targets. Counts are retained as modified/valid. |
| BED or bedGraph | Full-containment flag per deposited interval. A percentage without a read count is suitable for interval eligibility only. |

### Choose the target design

For the named hg38 presets, download the exact covered-target file from Twist. An updated vendor file or a panel with spike-ins can be supplied as `--preset custom` with its own explicit design name, matching reference and complete target BED. The software records file SHA-256 hashes in each report. It performs no liftover or contig renaming.

* `twist-human-methylome-hg38` uses `covered_targets_Twist_Methylome_hg38_annotated_collapsed.bed`. [Twist Human Methylome target BED](https://www.twistbioscience.com/resources/bed-files/twist-human-methylome-panel-target-bed-file)
* `twist-exome-2-hg38` uses `hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed`. [Twist Exome 2.0 BED files](https://www.twistbioscience.com/resources/bed-files/twist-exome-20-bed-files)

### Run the audit

```bash
kilnbio capture-panel --preset twist-exome-2-hg38 --name Twist-Exome-v2.0.2 \
  --assembly GRCh38 --data-assembly GRCh38 --sizes genome.fa.fai \
  --targets hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed \
  --bam sample.sorted.bam --min-depth 20 --measurements sample.vcf.gz \
  --format vcf --out exome-qc

kilnbio capture-panel --preset twist-human-methylome-hg38 --name Twist-Methylome \
  --assembly GRCh38 --data-assembly GRCh38 --sizes genome.fa.fai \
  --targets covered_targets_Twist_Methylome_hg38_annotated_collapsed.bed \
  --measurements cpg-counts.tsv --format cpg --out methylome-qc

kilnbio capture-panel --name lab-panel-v3 --assembly GRCh38 --sizes genome.fa.fai \
  --targets custom.targets.bed --baits custom.baits.bed --out design-qc
```

The output directory must be new. In `features.tsv`, coordinates are zero-based and half-open. Each row refers to one input feature. VCF insertions use the VCF anchor. CpG count input has the tab-separated header `chr pos N X` (1-based forward C, valid calls, modified calls). BED features keep their original interval. Alignments use mapping and base quality 20 by default, suppress overlapping mates with `samtools depth -s` and use samtools' default flag exclusions. DNA depth and valid modification depth measure different things. [SAMtools depth](https://www.htslib.org/doc/samtools-depth.html)

For variant accuracy, compare a callset from an actual capture specimen to matching truth inside the intersection of the design and confident regions. A design overlap or depth threshold alone cannot measure variant sensitivity. [NIST Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle)


### Public worked input

For a small exome exercise, use paired FASTQ from run DRR635838, the standard HG002 exome in PRJDB20167. The study used Twist Exome 2.0 plus a Comprehensive Exome spike-in. Download the two FASTQ files listed on the run record, retain the same first 10,000 complete records in each mate file and check paired identifiers. Align against the full GRCh38 no-alt assembly. The following commands take those two matched files and a Bowtie 2 genome index. Bowtie 2 v2.4.4 and SAMtools 1.24 were used for the example.

```bash
set -euo pipefail
bowtie2 --very-sensitive --seed 31 -p 4 -x "$GRCH38_INDEX" \
    -1 mate1.fq -2 mate2.fq | samtools sort -o subset.sorted.bam -
kilnbio capture-panel --preset twist-exome-2-hg38 --name HG002-core-exome \
    --assembly GRCh38 --data-assembly GRCh38 --sizes "$GRCH38_FAI" \
    --targets "$TWIST_EXOME_TARGETS" --bam subset.sorted.bam --out subset-core-qc
```

Here the denominator is the core design only. A shallow prefix cannot measure whole-exome sequencing performance. [ENA run DRR635838](https://www.ebi.ac.uk/ena/browser/view/DRR635838) · [Study design and spike-in](https://doi.org/10.1038/s10038-025-01403-4)
