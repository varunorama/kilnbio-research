# Public-method exome processing and clonal hematopoiesis

Research use only; not validated for clinical diagnosis or patient management. See [intended use](RESEARCH_USE.md) and the [cross-omic integration path](INTEGRATING_OMICS.md).

KilnBio now provides a pinned nf-core/sarek adapter for paired-end human
autosomal capture sequencing, alongside separate downstream CH interfaces.
These are research routes. CH means **clonal hematopoiesis**, not ChIP-seq.

| Route | Input → output | Established component / KilnBio contribution |
|---|---|---|
| `exome` | Paired FASTQ + reference + capture BED → alignment, QC and germline small variants | Sarek 3.12.0; KilnBio checks targets, samples, resource choices and resume provenance |
| `ch-somatic` | Blood FASTQ + capture BED + somatic resources → Mutect2 candidates/filtering | Sarek/Mutect2; a restricted research configuration, not a validated low-VAF CHIP assay |
| `ch-screen` | Biallelic annotated VCF/TSV + curated driver alleles → candidate and exclusion audit | Transparent public-method thresholds; no automatic clinical diagnosis |
| `ch-import` | Existing provider calls + mutation table + explicit specimen map → audited classifications and genomic events | Preserves source classifications, unknowns, VAF and batch metadata |
| `genomic-panel`, `genomic-associate` | Called variants or curated alteration states + assay matrix → genotype/trait associations | Existing KilnBio linkage and denominator checks; see [genomic linkage commands](../CLI.txt) |

The preprocessing backbone is [nf-core/sarek 3.12.0](https://github.com/nf-core/sarek/releases/tag/3.12.0),
released 11 November 2026, with its schema and license retained in the repository.
The workflow is published in [Hanssen et al. 2024](https://doi.org/10.1093/nargab/lqae031),
PMID 38666213. KilnBio does not claim to have invented its aligners or callers.

## Germline exome or a custom capture panel

```bash
kilnbio exome /study/exome.csv /study/exome.json /study/run01 --dry-run
kilnbio exome /study/exome.csv /study/exome.json /study/run01 --resume
```

Start from [examples/exome.csv](../examples/exome.csv) and
[examples/exome.json](../examples/exome.json). Replace every placeholder with a
local file. Linux/WSL2, Nextflow 25.10.4 and a working Docker/Apptainer/Singularity
runtime are required. Use paths without spaces for inputs and run directories.
Supply the native Sarek header:

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
DONOR1,XX,0,SAMPLE1,L1,/data/R1.fastq.gz,/data/R2.fastq.gz
```

`sex` is XX, XY or NA; germline status is 0. Multiple lanes must have consistent
metadata and distinct files. References, a matching FASTA index and a 0-based
half-open capture BED are explicit. There is no inferred genome or capture kit.
Twist Exome 2.0 and custom panels use the same interface; see
[capture-panel QC](CAPTURE_PANELS.md).

This initial adapter requires **autosomal targets**. Sex-chromosome ploidy/PAR
handling, mitochondrial calling and full-genome WGS are outside this supported
configuration. It fails on sex-chromosome targets rather than silently using
diploid assumptions. Make an explicitly documented autosomal subset of the
original panel and retain the original BED for provenance.

Choose `bwa-mem`, `bwa-mem2` or `dragmap` explicitly. Small-variant choices are
`haplotypecaller`, `deepvariant` and `strelka`, including comma-separated
comparisons. At least one small-variant caller is required. The example enables
HaplotypeCaller and DeepVariant; this is a comparison option, not evidence that
either is universally best. WES-compatible Manta and VEP are optional upstream
tools; their inclusion is not an independent SV/annotation benchmark.

BQSR requires declared dbSNP and known-indel resources with indexes. An explicit
`skip_tools: "baserecalibrator"` is allowed for a justified protocol; it must not
be an undocumented default. Mapping, duplicate marking, recalibration, calling
and upstream QC are handled by the pinned workflow. Review capture coverage,
contamination, identity, missing genotypes and filtering before association.
Run plans retain original/resolved inputs, resource limits and execution traces.

**Known empty-callset edge case:** Sarek's HaplotypeCaller CNN/tranche filtering
fails on a VCF without variants. The tiny smoke therefore explicitly sets
`skip_tools: "haplotypecaller_filter"`. That omits post-calling filtering and is
recorded in the plan; it is not the general study default. Any nonempty output
from such a run still needs an independently specified filtering procedure.
Do not claim post-filter validation or accuracy from the empty-callset smoke.

VEP requires a local versioned cache and declared version, species, assembly and
arguments. Use public Ensembl resources; optional predictors have their own
licenses. This preset does not enable SpliceAI, dbNSFP, Mastermind, Condel or
phenotype plugins. No proprietary model, company threshold, private cohort
reference or polygenic score was imported from the supplied exome context.

Gene-level CNV burden/gCNV, ancestry, relatedness, PRS and clinical pathogenicity
classification are **not** implemented by this adapter. A free executable does
not make every annotation resource or trained model freely redistributable.
Download public resources directly, record versions and preserve their terms.

## CH somatic calling is a different problem

```bash
kilnbio ch-somatic /study/blood.csv /study/ch-somatic.json /study/ch01 --dry-run
```

The [somatic examples](../examples/ch-somatic.json) require an assay-appropriate
panel of normals, germline AF resource and common-SNP contamination resource,
with indexes. Status 1 is Sarek's somatic-calling input convention; it does not
label the donor as having cancer. The adapter starts from FASTQ, not an arbitrary
germline VCF. Somatic origin, low-VAF sensitivity, artifact rejection and clinical
context still require a validated protocol. A typical germline exome callset
cannot be assumed to detect small blood clones adequately.

## Screen already-called variants

```bash
kilnbio ch-screen --samples blood_samples.tsv --variants annotated.vcf.gz \
  --drivers curated_drivers.tsv --assembly GRCh38 --out new_ch_screen
```

The sample roster has `sample,donor,specimen,tissue` (tab-separated, tissue `blood`),
one independent specimen per donor. Driver TSV columns are
`assembly,chrom,pos,ref,alt,gene,evidence`: exact normalized biallelic alleles,
1-based positions, a public curation citation and an assembly matching the VCF.
KilnBio does not guess that every variant in a CH-associated gene is a driver.
Build the curated list from a documented public source and review consequence,
transcript and indel representation. The list and input files are hashed.

VCF input must have sample columns, PASS filtering, FORMAT DP/AD, either AF or
enough information to derive VAF, and alternate support in both read strands
(ADF/ADR) or both pair orientations (F1R2/F2R1). These are different measurements
and the audit labels which was supplied. The population-frequency INFO field
defaults to `gnomAD_AF`; change `--population-af-field` for the actual annotation.
Mutect2 output alone need not contain that annotation. **Missing population
frequency is unknown, not zero.** Multiallelics must be split/normalized against
the matching reference first. Do not silently change genome builds or contig names.

An equivalent TSV needs `sample,assembly,chrom,pos,ref,alt,filter,dp,alt_count,vaf,
alt_forward,alt_reverse,population_af,orientation_kind`. Orientation kind is
`read_strand` or `pair_orientation`. Threshold defaults are depth≥20,
alternate reads≥3, VAF≥0.02, alternate support≥1 in each supplied direction and
population AF≤0.001. They are a transparent screening subset of public methods,
not a reconstruction of every hotspot exception, rescue filter or manual review
in [Uddin et al. 2024](https://doi.org/10.1038/s41467-024-52302-9), PMID 39251642.
Thresholds are CLI parameters and are written to `summary.json`.

`variants.tsv` gives each reason for exclusion. `samples.tsv` reports candidate
count and maximum VAF (not the sum of clone sizes). No candidate is promoted to
a confirmed somatic driver or clinical CHIP diagnosis. No qualifying record
means `no_qualifying_variant_in_supplied_records`, **not a negative call**.
Coverage/callability, tissue, germline alternatives, sequencing artifacts,
cytopenias and hematologic malignancy must be reviewed. The hg38 U2AF1 reference
duplication issue needs a separately validated method; there is no automatic
mapping-quality relaxation here.

## Import existing provider classifications

```bash
kilnbio ch-import --calls provider_calls.csv --variants provider_variants.csv \
  --samples private_linkage.tsv --out new_provider_import
```

Supported calls schema: `Biobank_Subject_ID,hasCHIP,nCHIP,Gene,VAF,Batch_CHIP_call`.
Mutation schema: `Biobank_Subject_ID,varID,Gene,VAF,DP,AD_Alt,ADF_Alt,ADR_Alt`.
The tab-separated linkage supplies `sample,donor,specimen,source_subject`.
Aliases are chosen locally; original subject identifiers need not enter output.
Counts/flags must agree with returned variant records. A subject absent from the
provider table remains unknown. Provider negatives are labeled as reported
negatives, not independently verified absence.

`genomic_events.tsv` exposes reported altered/reference/unknown states to
`genomic-associate`; candidate-only screens deliberately emit unknown states so
they cannot manufacture reference controls. For covariate-adjusted multi-layer
analysis, construct a reviewed donor trait from verified calls and use
[multiomic-layers](MULTIOMIC_AND_OLINK.md). Match donor **and tissue/specimen** as
appropriate; a blood CH event cannot be treated as a tumor-cell mutation.

No participant records, source-return archive or private provider documentation
are distributed with the code. The import interface operates on the user's own
authorized data, while public teaching/testing examples remain synthetic.

## Public-data examples and manuscript scope

1. **Exome engineering and accuracy are separate tests.**
   `python3 tests/public_exome_smoke.py /absolute/new_directory` stages a tiny,
   pinned public nf-core fixture and prints its KilnBio command. It includes
   mapping/BQSR/HaplotypeCaller/QC; it is not a full human exome benchmark. For
   accuracy, the existing [Twist guide](CAPTURE_PANELS.md) identifies public HG002
   exome data (DRR635838/PRJDB20167). Compare against GIAB truth only inside the
   intersection of the exact capture BED and confident truth regions, stratified
   by variant type/difficult context. Follow
   [GA4GH benchmarking recommendations](https://doi.org/10.1038/s41587-019-0054-x).
   A full accuracy run is not claimed by this release.
2. **CH + Olink + chromatin/RNA/methylation is a strong candidate demonstration.**
   [Ismail et al., Blood Advances](https://doi.org/10.1182/bloodadvances.2024014467)
   (PMID 39631069; 2025 issue) describes GSE210435 and supplies Olink data in its
   supplements. The paper reports genotype, single-cell RNA/ATAC, methylation
   and protein measurements with public [analysis code](https://github.com/LabFunEpi/CC_multiomics).
   These are targeted CH measurements, not an exome-calling truth set. Availability
   is documented in the publication; a runnable KilnBio cohort recipe still needs
   file checks, cross-assay donor mapping, covariates and an independent endpoint.
   Published IL32 regulation is a reproduction target, not a new KilnBio discovery.
3. **Exome-derived mutations + RNA/protein can use the existing bladder recipe.**
   It starts from open processed somatic calls and molecular data. TCGA raw
   germline-sensitive sequencing requires controlled access; do not describe the
   recipe as an unrestricted FASTQ-to-exome demonstration. Adding genotype-aware
   cohort/target sensitivity can test what a genetic layer changes without
   redistributing restricted raw data.

Suitable manuscript additions now are the support table, methods, interfaces and
explicit software validation. New biological conclusions need an executed,
frozen public reanalysis with the same donors in the competing arms, independent
validation, sample-exclusion accounting and comparison with the source paper.
Do not label synthetic fixtures or imported provider calls as novel biology.
The existing manuscript results and versioned archives are not silently replaced.

## Verification and attribution

`tests/integration_platform_extensions.py` exercises NPX integrity, missing LOD,
confounding, CH filters, VCF orientation and dependence-aware partial nulls.
`tests/integration_multiomic_layers.R` executes limma/voom and paired models.
`tests/public_olink_smoke.R` exports the vendor's explicitly synthetic
[npx_data1 fixture](https://github.com/Olink-Proteomics/OlinkRPackage/blob/7385f477050b65e9470e709d30b1b23c6fcb5702/OlinkAnalyze/man/npx_data1.Rd).
Other callers, full WES accuracy and de novo CH sensitivity are separate validation
claims, not implied by adapter availability. Logs/results belong in the local
analysis archive; only reusable code and public recipes belong in GitHub.
