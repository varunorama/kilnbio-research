# Human biological validation

Protocol frozen 2026-11-01 before these new analyses. Each example reports its
actual starting point. A count-matrix or BAM-start experiment does not validate
FASTQ mapping. Recovery of a published response is a positive control, not a
new biological discovery or an estimate of the false-discovery rate.

| Route | Public experiment / endpoint | Unit and scope | Additional validation required |
|---|---|---|---|
| Bulk RNA | GSE52778 / airway: dexamethasone response, CRISPLD2 induction | Four paired primary human airway smooth muscle donors; full released integer counts → donor-adjusted DESeq2 and edgeR | Human FASTQ → mapping/RSEM; real BAM → featureCounts is assessed separately |
| 10x gene expression | GSE96583 / Kang: IFN-beta response in PBMCs | Eight paired donors; released cell-called raw counts → QC → donor/cell-type pseudobulk → paired DE | Raw human FASTQ, ambient RNA, doublet sensitivity, and V3/V4 chemistry validation |
| ATAC | ENCSR095QNB / ENCSR868FGK, GM12878 versus K562; cell-matched DNase | Two culture replicates per cell line; independently defined regions and biological specificity | Population inference, differential error calibration, full raw-read route and global shifts |
| Histone ChIP | ENCSR000AKD / ENCSR000AKQ H3K27me3; ENCSR000DRW / ENCSR000DWB H3K36me3 | Two culture replicates per line, matched input controls; independent reference regions | Whole-genome differential biology, externally normalized global changes, deeper libraries |
| WGBS / RRBS / DML / DMR | Select protocol-matched public independent biological replicates and orthogonal controls | Conversion and coverage checks separate from differential inference | Existing synthetic inference and tiny upstream tests do not close this route |
| ONT modBAM / POD5 | Select real MM/ML-tagged human reads and matched methylation or control material | Real-data import/pileup separate from raw-signal basecalling | Dorado/GPU route, channel calibration, independent human truth |
| biomodal | Official HG001 duet +modC and evoC v1.5.0 demo, including conversion controls | Vendor reports → explicit modC / mC / hmC channels; single reference donor | Public bucket requires authentication; vendor raw processing remains a separate licensed route |
| Proteomics / metabolomics | Existing Quartet reference-material, held-out-laboratory recipes | Released abundance matrices → reference adjustment → held-out material classification | Raw mass spectrometry is not implemented; clinical generalization requires new cohorts |
| Aging and cell-state models | Existing model-specific interfaces and erythroid donor-held-out examples | Preserve each model's required feature space and preprocessing | Independent cohort calibration for non-DNAm clocks; lineage and tissue transfer |

## Prespecified endpoints

Airway: CRISPLD2 (ENSG00000103196) must increase with dexamethasone in the
paired analysis. Report its effect, FDR and each donor's normalized direction
under both engines. Other established glucocorticoid-response genes are
secondary checks; do not choose endpoints from the new result ranking. Compare
both engines over the same tested genes, and show per-donor behavior rather
than relying only on a p-value. Four pairs are not an empirical FDR calibration.

Kang: IFIT1, ISG15 and MX1 induction in author-annotated CD14 monocytes are
prespecified positive controls. Use donor as a blocking factor. Report the
counts of retained cells and donor-condition groups, and preserve author cell
labels. A pooled condition library can confound treatment with capture batch;
this example is a known-response recovery test, not proof that batch was removed.
Keep author's singlet/donor assignments explicit; do not relabel every cell as
an independent replicate. Report dropout/QC exclusions before differential fits.

ATAC/ChIP: define regions independently of the contrast being evaluated. Use
released matched controls and reference assembly consistently. Do not declare
cell-line cultures to be independent human donors. No genome-wide peak-caller
default change follows from a regional or single-study result alone.

Every run retains source identifiers, download checksums, software versions,
configuration, count/filter diagnostics, complete tested features, and failures.
Do not silently replace unavailable measurements with zero or force published
expected effects by selecting thresholds after seeing the results.

Additional methylation protocol frozen before fitting: use all chr21 loci from
the unsmoothed `bsseqData::BS.cancer.ex` counts, three paired colon tumor/normal
patients, DSS with patient blocking and no smoothing. Report median within-pair
fraction changes on common covered loci and the direction of tested loci.
The published broad hypomethylation is the positive control; this is not a
base-resolution truth set. Check the native reference coordinates before any
annotation transfer. The current dmrseq route rejects paired designs, so this
experiment cannot validate its region q-values. Source: [Hansen et al. 2011](https://doi.org/10.1038/ng.865),
PMID 21706001; [bsseqData documentation](https://bioconductor.org/packages/bsseqData/).

## Primary sources

- Himes et al. 2014, [RNA-seq and experimental CRISPLD2 validation](https://doi.org/10.1371/journal.pone.0099625), PMID 24926665;
  [airway data package](https://bioconductor.org/packages/airway/).
- Kang et al. 2018, [genetic multiplexing and the IFN-beta PBMC experiment](https://doi.org/10.1038/nbt.4042), PMID 29227470;
  [muscData](https://bioconductor.org/packages/muscData/).
- [ENCODE ATAC standards and resources](https://www.encodeproject.org/atac-seq/);
  exact files and original study records are in [the chromatin manifest](../chromatin/files.tsv).
- [Biomodal official demo data](https://biomodal.com/kb_articles/demo-data-instructions/).
  v1.5.0 +modC/evoC controls differ from the newer mosaic controls.
- Patel et al. 2024, [The Wild West of Spike-in Normalization](https://doi.org/10.1038/s41587-024-02377-y), PMID 39271835.

Numerical analyses and manuscript outputs are written to a caller-supplied
directory outside the repository. This protocol does not claim that pending
routes have passed biological validation.
