# Raw-read and spike-in validation protocol

Frozen before alignment/results, 2026-11-01. These are reproducible recovery
benchmarks, not newly discovered biology or general false-discovery calibration.
All downloads, read subsets, alignments, numerical results and figures stay in
the caller's study directory outside Git. Full libraries are retained.

## Human bulk RNA-seq

GSE52778: four paired primary human airway smooth muscle donors, dexamethasone
versus vehicle. Runs SRR1039508/09, SRR1039512/13, SRR1039516/17 and SRR1039520/21.
Use **five million randomly selected complete read pairs per sample**, seed
74129, paired sampling with identical seeds and verified matching identifiers.
Use the complete GRCh37 primary assembly and Ensembl release 75 GTF, with the
manufacturer's ERCC92 sequences/annotation. The older assembly is deliberately
matched to the original airway comparison. It is not a recommended universal
default for new studies. Record assembly and annotation hashes.

The archive contains orphan FASTQs for four runs. Preserve and inventory them;
the primary RSEM/featureCounts comparison uses the same complete pairs and
explicitly excludes orphans. It does not represent every read in the library.
The original paper reports 12-base clipping of 75-base reads. Confirm deposited
lengths before trimming; do not clip another 12 bases if already 63 bases long.

Run the pinned RNA adapter with STAR/RSEM, adapter/quality QC, inferred strand,
and real genomic BAMs through featureCounts. For low-memory machines the STAR
suffix-array density can be reduced without truncating the genome; record the
index settings and actual resource use. Analyze RSEM with tximport/DESeq2 and
featureCounts integers with both DESeq2 and edgeR, adjusting for paired donor.
GEO and the paper identify **three** mix-1-positive libraries (SRR1039513,
SRR1039520, SRR1039521); the other five have `ercc_mix: -`. Audit detected
signal separately in both sets, without requiring positive counts in negatives.
ERCCs support mapping/QC; they are excluded from endogenous differential tests
and are not assumed to establish absolute transcripts per cell.

Primary response: **CRISPLD2 (ENSG00000103196) increases**, paired log2 fold change
positive and BH-adjusted p <0.05; report all four donor directions. Compare
matched-gene effects between counting engines and against the released-count
example, disclosing unequal depths and mapping methods. Mapping/QC failures or
discordant controls are retained, not removed based on the desired response.

## ChIP-seq perturbation and fly spike-ins

GSE64243, PC9 cells, DMSO versus 1 micromolar GSK126 for five days:

| Run | Mark | Condition |
|---|---|---|
| SRR2532668 | H3K27me3 | DMSO |
| SRR2532669 | H3K27me3 | EZH2 inhibitor |
| SRR2532670 | H3K4me3 | DMSO |
| SRR2532667 | H3K4me3 | EZH2 inhibitor |

Select **five million reads per library**, seed 74129, before preprocessing.
These four runs are campaign R2. **Metadata amendment before alignment:** the
older four GEO titles omit the spike-in label, but their supplementary filenames
match identifiers in both tabs of the paper's Supplementary Table S4. R1 uses
SRR1721586/87 (H3K27me3 vehicle/inhibitor, G3/G4) and SRR1721588/89 (H3K4me3
vehicle/inhibitor, G87/G89). R2 IDs are G90/G93 and G76/G77 respectively.
Verify each raw library has usable fly reads. Evaluate the same fixed-depth
endpoints separately in each campaign, then combine the two campaigns with a
campaign-blocked design per mark if both pass input QC. No significance tests
are made from a single campaign and no matched input library is invented.

Map competitively to the complete hg19 + dm3 reference with disjoint species
prefixes. Use BWA aln/samse for the 50-nt single-end libraries, at most two
alignment mismatches, MAPQ >=30, exclude unmapped/secondary/supplementary/QC-fail
records, and mark/remove positional duplicates for primary counts. Report a
duplicate-retaining sensitivity. Preserve all read losses and both species'
counts; do not align to human only before estimating the fly fraction.

The stated mixture is 30 micrograms human chromatin and 750 ng OSS fly
chromatin with H2Av antibody. Relative exposure is
`fly_counts * human_input_amount / fly_input_amount`. Scale human region counts
by this exposure; do not subsequently impose TMM or median-ratio normalization
that can erase the global shift. Units are per input chromatin mass, not per cell.
Cross-reactive antibodies/reference recovery remain potential biases. Whole-fly
scaling is the primary available benchmark; H2Av-only scaling requires a
documented independent reference-region set and is not silently substituted.

Prespecified controls: H3K27me3 inhibitor/control ratio <0.5; H3K4me3 between
0.5 and 2; H3K27me3/H3K4me3 relative response <0.5. Report continuous ratios
and sensitivity, including failures. Inspect MYT1 and CCND2 promoter signal
without selecting them from peak rankings. In [Egan et al., Fig. 2](https://doi.org/10.1371/journal.pone.0166438.g002),
ChIP-qPCR before library amplification supported MYT1 loss; CCND2 appears
in a different panel measuring amplified libraries, where the treatment
difference was no longer apparent. Neither promoter window is the original
amplicon. These are exploratory loci, not two independently validated positives.
Compare reference scaling with ordinary human library-size scaling. The latter
fixes the genome-wide sum by definition; its inability to show a global shift
is a normalization constraint, not a new empirical discovery.

For the replicated exploratory regional analysis use fixed 10 kb autosomal
tiles, exclude assembly gaps and hg19 blacklist overlaps before testing, and
count primary nonduplicate leftmost read starts once. Compare edgeR QL with
ordinary TMM and experimental spike exposures using the same integer matrix,
contrast and `~ campaign + treatment` design. Report within-campaign effect
directions and continuous regional effects, not a false-discovery truth claim;
the blocked model has only one residual degree of freedom. Count-table
reproduction from S4 is an additional, distinct evidence tier, not raw alignment.

## Evidence and limits

- [Himes et al., 2014](https://doi.org/10.1371/journal.pone.0099625), PMID 24926665.
- [Egan et al., 2016](https://doi.org/10.1371/journal.pone.0166438), PMID 27875550;
  supplementary tables distinguish whole-fly and H2Av-restricted factors.
- [Patel et al., 2024](https://doi.org/10.1038/s41587-024-02377-y), PMID 39271835.
- [STAR parameters](https://github.com/alexdobin/STAR/blob/master/source/parametersDefault);
  [seqtk paired sampling](https://github.com/lh3/seqtk).
- Bibliographic/notice records: [PubMed MCP evidence](raw_citations.json).

The companion synthetic test verifies recovery of a known global loss and
direct DESeq2/edgeR equivalence with external factors. It does not validate
biological variance or global FDR. These fixed-depth examples do not establish
full-library throughput, independent clinical generalization, or a peak-caller
default. Their starting point and limits must accompany any status claim.

## Run the examples on another machine

Use a Linux study path without whitespace, outside this repository, and the
RNA/DE environments in the [user guide](../../docs/USER_GUIDE.md). The
[ENA manifest](public_fastq_manifest.tsv) contains accession, URL, size and MD5
for every archived file in the comparison; it contains no analysis outputs.
The [RNA design](airway_design.csv) reproduces the donor, condition and ERCC
labels in [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778).
Sampling uses seqtk 1.5 and pigz. All alignments, checksums and analysis outputs
are written under the analyst's chosen study path.

```bash
set -euo pipefail
K=$PWD                       # KilnBio checkout
B=$K/benchmarks/biological
D=/data/kilnbio-validation   # new directory outside Git
mkdir -p "$D/sources"
python3 "$B/fetch_fastq.py" "$B/public_fastq_manifest.tsv" "$D/fastq" --range-workers 4
python3 "$B/raw_inputs.py" sample "$D/fastq" "$D" "$(command -v seqtk)"

# Save the original reference downloads linked below in $D/sources first.
S=$D/sources
python3 "$B/airway_reference.py" "$S/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz" \
  "$S/Homo_sapiens.GRCh37.75.gtf.gz" "$S/ERCC92.zip" "$D/rna-ref"
python3 "$B/raw_indexes.py" airway "$D/rna-ref"
python3 "$B/raw_inputs.py" rna-plan "$D" "$D/rna-ref"
bash "$K/run.sh" rnaseq "$D/rna-inputs/samples.csv" "$D/rna-inputs/params.json" \
  "$D/rna-run" --cpus 4 --memory-gb 11
python3 "$B/airway_raw.py" prepare "$D/rna-run" "$D/rna-ref/genome.gtf" \
  "$B/airway_design.csv" "$D/rna-analysis"
# AIRWAY_RELEASED points to the completed gene-count comparator described in README.md.
python3 "$B/airway_raw.py" fit "$D/rna-analysis" "$AIRWAY_RELEASED"
# Secondary within-positive-library dose ordering; no absolute per-cell claims.
python3 "$B/ercc_dose.py" "$S/ERCC_Controls_Analysis.txt" "$D/rna-run/results/star_rsem" \
  "$D/rna-analysis" "$D/ercc-dose"

python3 "$B/raw_indexes.py" pc9 "$S/hg19.fa.gz" "$S/chromFa.tar.gz" "$D/chip-ref"
python3 "$B/chip_regions.py" freeze "$D/rna-ref/genome.gtf" "$D/hg19-regions"
for c in R1 R2; do
  python3 "$B/chip_spikein.py" "$D/sampled-5m" "$D/chip-ref" "$D/chip-$c" --campaign "$c"
done
python3 "$B/chip_read_qc.py" "$D/chip-R1" "$D/chip-R2" "$D/chip-read-qc"
python3 "$B/chip_regions.py" assemble "$D/hg19-regions" "$D/chip-R1" "$D/chip-R2" "$D/chip-combined"
python3 "$B/chip_regions.py" fit "$D/chip-combined"
```

Reference sources: [Ensembl 75 human primary assembly and CHECKSUMS](https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/),
[Ensembl 75 annotation and CHECKSUMS](https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/),
[Thermo Fisher ERCC92](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip),
[manufacturer's ERCC concentration table](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt) (save as `ERCC_Controls_Analysis.txt`),
[UCSC hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/),
and [UCSC dm3](https://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/).
Use `chromFa.tar.gz` for dm3. The ChIP recipe verifies both original archives
against UCSC MD5s before adding species prefixes. The RNA recipe records the
assembly/annotation hashes, validates exon coordinates and compares the STAR and
RSEM transcript ID sets before alignment. The ChIP mapping and mask methods
are described in [BWA](https://doi.org/10.1093/bioinformatics/btp324),
[SAMtools](https://doi.org/10.1093/gigascience/giab008) and the
[ENCODE Blacklist](https://doi.org/10.1038/s41598-019-45839-z).
