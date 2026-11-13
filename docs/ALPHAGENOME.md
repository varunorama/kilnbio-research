# Optional AlphaGenome annotations

Research use only; not validated for clinical diagnosis or patient management. See [intended use](RESEARCH_USE.md) and the [cross-omic integration path](INTEGRATING_OMICS.md).

AlphaGenome is a scientific DNA-sequence prediction model, not a general-purpose
LLM. This optional command adds predicted regulatory effects to an explicitly
selected variant list. Ordinary KilnBio preprocessing and statistics do not use it.
It is also different from AlphaFold2, which predicts protein structure.

## What is implemented

`kilnbio alphagenome CONFIG OUT` validates and writes a local plan. Only adding
`--execute` makes a request to the external AlphaGenome service. The adapter pins
the Python SDK to **2.0.0** and selects an explicit AlphaGenome 2 checkpoint.
It supports human GRCh38 variant scoring for reference tissues/cell types, or
conditioning on a deliberately selected set of human scRNA-seq cells. The provider's
supported scorer names include RNA, ATAC, histone/TF ChIP, methylation and splicing.
This is not raw-read processing, a generic upload of every omic, or AlphaFold access.

The SDK's [v2.0.0 release](https://github.com/google-deepmind/alphagenome/releases/tag/v2.0.0)
was published on 11 November 2026. The integration follows its
[getting-started example](https://www.alphagenomedocs.com/getting_started.html)
and [scoring interface](https://www.alphagenomedocs.com/api/generated/alphagenome.models.v2.dna_client.DnaClient.score_variant.html).
AlphaGenome 1's [Nature paper](https://doi.org/10.1038/s41586-025-10014-0)
is not evidence that this new v2 adapter has been biologically validated.

## Try the supplied public variant

Python 3.11 or later is needed for the optional SDK. Install it in its own environment:

```bash
python3 -m venv ../kilnbio-alphagenome-env
../kilnbio-alphagenome-env/bin/python -m pip install -r environments/alphagenome-requirements.txt
export KILNBIO_PYTHON="$PWD/../kilnbio-alphagenome-env/bin/python"
bash run.sh alphagenome examples/alphagenome/reference.json glb1l
```

The last command is offline and creates `glb1l.plan.json`. The variant and blood
monocyte condition come from the provider's tutorial. This demonstrates the
interface; it is not a new biological discovery. Check assembly, alleles, tissue,
scorers and request count in the plan. Obtain your own API key through the
[provider](https://alphagenome.google/api) and set `ALPHA_GENOME_API_KEY` locally.
Keep it out of configuration files, command histories and Git.

```bash
bash run.sh alphagenome examples/alphagenome/reference.json glb1l --execute
# Resume only with the same input bytes, settings and adapter code:
bash run.sh alphagenome examples/alphagenome/reference.json glb1l --execute --resume
```

`scores.tsv` contains model predictions and provider metadata. Per-variant files,
hashes, `plan.json` and `receipt.json` retain the request context and completion
status. An interrupted run is incomplete, returns a failure, and does not produce
a complete combined table. Completed per-variant requests can be reused on resume.
Pinning the SDK/checkpoint name does not guarantee the hosted service never changes;
retain outputs and timestamps for reproducibility.

## Bring variants from your study

Supply a CSV/TSV with exactly `variant_id,chromosome,position,ref,alt`. Positions
are **one-based**, as in VCF; prediction windows in the plan are zero-based half-open.
Use normalized, biallelic variants with `chr`-prefixed primary GRCh38 chromosomes.
Do not include donor identifiers, genotypes or clinical columns. The adapter does
not perform liftOver, VCF normalization or local reference-allele checking.
Validate the reference and allele representation upstream. Structural/symbolic
variants and windows crossing chromosome boundaries are rejected.

An exome/somatic variant can be scored here and then joined to measured RNA,
protein or regulatory evidence through its stable variant/gene mapping. For
ATAC peaks or DMRs, first identify the variants of interest overlapping those
regions; this command does not send peak heights, beta values or NPX as sequence.
Predicted effect units differ among scorers; inspect the provider's metadata.
Never treat these values as p-values, clinical classifications, or additional
independent assays in a multi-assay significance test. Model training overlap
also matters when designing an accuracy or held-out benchmark.

## Optional single-cell expression conditioning

Instead of `condition`, set `cell_counts` to a local `.h5ad` and
`cell_condition_name` to a generic label. Provide 1–256 explicitly selected cells,
raw integer UMI counts in `X`, and unique human Ensembl gene IDs in `var_names`.
Normalized/fractional counts, empty cells and duplicate IDs after version removal
are rejected. KilnBio performs no hidden subsampling or conversion from log counts.
The adapter reconstructs the object with only counts and gene IDs, removing donor
metadata. The expression values still leave your computer: execution additionally
requires `--allow-cell-upload`. Use only data you may send to the service.
`SC_GENE_EXPRESSION_DENOISED` requires this mode. See the provider's
[conditioning tutorial](https://www.alphagenomedocs.com/colabs/v2/conditioning_on_single_cell_data.html).

## Availability, terms and validation

The SDK is Apache-2.0, but the hosted service and outputs have separate
[provider terms](https://alphagenome.google/terms). The provider's
[README](https://github.com/google-deepmind/alphagenome/blob/v2.0.0/README.md#terms)
describes its free API as noncommercial and restricts model-training use, with
specific exceptions and a separate commercial route. KilnBio's own license does
not override those terms. These outputs are not automatically exported for the
foundation venture or used as training labels.

Offline validation covers coordinates, explicit execution, SDK signatures,
count validation, metadata removal, completion accounting and resume integrity.
A credentialed public-variant run and biological accuracy benchmark remain separate
validation steps. No live API success is claimed without a successful receipt.

For coding variants, the separate [variant annotation interfaces](VARIANT_ANNOTATIONS.md)
support VEP consequences, ClinVar assertions, AlphaMissense scores and AlphaFold
Database metadata. They retain their own evidence types and licenses; protein
structure confidence alone does not establish a variant's pathogenicity. These
are separate commands, not outputs of an AlphaGenome query.

The original AlphaGenome model is described by [Avsec et al., Nature 2026](https://doi.org/10.1038/s41586-025-10014-0)
([PubMed](https://pubmed.ncbi.nlm.nih.gov/41606153/)). That paper is background for
regulatory prediction; it does not validate this adapter or the newer version 2
API. Version 2 request behavior is documented by the pinned provider SDK above.
