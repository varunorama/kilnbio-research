# Run an integrated study with one command

`kilnbio study` orchestrates declared KilnBio commands from one study JSON. It runs prerequisites before dependent analyses, stops on failure, retains logs and checks expected output files. A single invocation can include processing, statistical analyses, cross-assay integration and optional annotations **when their inputs and handoffs are explicitly configured**. It does not infer a study design from arbitrary files.

```bash
kilnbio study /study/study.json /results/study01
```

This runner executes stages sequentially to bound resource use. Existing Nextflow preprocessing workflows retain their own task parallelism and container settings. KilnBio's older `workflows/analysis.nf` remains a dispatcher for independent statistical jobs; it is not the dependency runner described here.

## A beginner can run a supplied recipe

Install the [direct-analysis dependencies](MULTIOMIC_AND_OLINK.md#install-the-direct-analysis-dependencies), including NumPy/pandas/SciPy and R limma/edgeR/jsonlite. In the source checkout:

```bash
python3 examples/multiomic/make_study_demo.py ../my-demo
kilnbio study ../my-demo/study.json ../my-demo-run
```

This entirely fictional example imports delivered-style Olink NPX, fits a matched-donor contrast across RNA, protein and methylation, then tests evidence across layers. It requires no API key, LLM or external study download. It is a teaching and software-integration example, not a biological validation study. The methylation orientation is prescribed by the simulation, not a universal biological relationship.

Start by reading `study-state.json` for completion and logs, `layers/layers_design.tsv` for donor/feature denominators, and `joint/genes.tsv` and `joint/summary.json` for the evidence and its interpretation. The table preserves untestable features. Selection does not establish a clinical biomarker or a causal mechanism.

## Preview and resume

```bash
kilnbio study study.json study01 --dry-run
kilnbio study study.json study01 --resume
```

The preview writes the resolved plan and checks the listed inputs without running analysis stages. Resume skips completed stages only when their declared output hashes still match. Changed study configuration, listed inputs, runtime selection or analysis source requires a new run. Failed/interrupted stage directories are moved into `attempts/` before retry; nothing automatically deletes raw reads. For a failed Nextflow stage, this initial study runner retries that stage from a new directory rather than claiming task-level cache recovery. Individual preprocessing commands expose their own `--resume` facility.

A `.study.lock` prevents two runners writing the same study. After an abrupt machine shutdown, first confirm no runner is active before removing a stale lock. Outputs that were not declared are not hash-checked by the study runner; retain the stage's own provenance as well.

## Configure a new study without writing analysis code

The JSON has `schema_version: 1`, `research_use: true`, `inputs`, optional `configs`, and `steps`. List **every external study input** in `inputs`; each file is hashed for resume checks. Large FASTQs and references take time to hash. JSON templates named in `configs` are hashed and copied into `RUN/configs/NAME.json` after resolving `{{study}}` and `{{run}}` placeholders. Paths in templates should use these placeholders explicitly.

Each step has:

| Field | Meaning |
| --- | --- |
| `id` | Unique simple name; its output directory is `RUN/id` |
| `command` | An allowed KilnBio module, such as `rnaseq`, `de`, `multiomic-layers`, `multiomic-meta` or `variant-annotations` |
| `args` | Argument list using `{{study}}` for the recipe directory and `{{run}}` for the output directory |
| `depends_on` | IDs of prerequisite stages; cycles and unknown IDs fail before execution |
| `outputs` | Expected files beneath the stage ID, e.g. `layers/manifest.tsv`; missing files prevent completion |

The generated demonstration supplies a complete configuration. Intermediate users can change contrasts, covariates, methods, mappings and eligible layers using the underlying modules' documented options. Raw-read outputs often require study-specific count manifests or reviewed feature maps before downstream integration; this runner does not yet generate every possible handoff automatically. Use an existing end-to-end recipe or prepare those explicit interfaces with a bioinformatician.

Remote AlphaGenome/AlphaFold stages require both the module's `--execute` flag and `allow_remote_scientific_services: true`. Credentials remain in the local environment. General-purpose LLMs are not invoked. The runner uses argument lists and never evaluates arbitrary shell fragments from the study configuration.

## Reuse modules in advanced workflows

Advanced users can call each R/Python module or the CLI directly, use the existing Nextflow workflows on supported infrastructure, and consume TSV/JSON outputs in a larger workflow. The study runner is an optional orchestration layer over the same functions, not a separate statistical implementation. Complex designs, new assay formats, custom transformations and clinical workflows are not inferred or validated by adding a step to a JSON file. See [integration principles](INTEGRATING_OMICS.md) and [research-use limits](RESEARCH_USE.md).
