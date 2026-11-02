# Worked analyses from public inputs

Each recipe downloads or verifies public measurements, computes new results locally and refuses to overwrite its output directory. These are retrospective research examples at the deposited data level, not FASTQ-to-publication pipelines or validated clinical predictors. Cite the original studies and follow each source's access and use terms. Run commands from the top of the unpacked release.

| Recipe | Question and comparison | Main outputs |
|---|---|---|
| Glioma | With an equal follow-up budget, how many candidates nominated from joint, RNA or protein discovery summaries agree in independent IDH-labelled glioma cohorts? Test directions and track targets that could not be evaluated. | Frozen panels, per-gene eligibility and confirmation, paired differences |
| Quartet | Does adjustment reduce same-material technical error without increasing error in known biological contrasts in held-out labs? RNA, protein and metabolite are evaluated separately. | Per-lab errors and paired comparisons |
| K562 | Across held-out CRISPR enhancer–gene measurements, does CpG methylation improve predictions beyond chromatin and coverage features on the same set of elements? | Coverage table, average precision, shortlists and conditional resampling interval |
| Lung adenocarcinoma | Separate the effect of restricting the RNA cohort from the effect of adding protein; compare an exact-ID join with audited pairing and deliberate label conflicts. Freeze rankings before external RNA evaluation. | Pairing audit, frozen targets, external tests and target-by-target changes |

Use Python 3.12 and a separate virtual environment per recipe to avoid conflicts. Shell examples below write to a new directory under your home directory. Change the output name for a repeat run.

```bash
WORK="$HOME/kilnbio-public-recipes"
mkdir -p "$WORK"
python3.12 -m venv "$WORK/glioma-env"
"$WORK/glioma-env/bin/python" -m pip install -r recipes/glioma-quartet/requirements.txt
"$WORK/glioma-env/bin/python" recipes/glioma-quartet/reproduce.py \
  --fetch --cache "$WORK/glioma-cache" --out "$WORK/glioma-quartet-01" --rscript Rscript

python3.12 -m venv "$WORK/k562-env"
"$WORK/k562-env/bin/python" -m pip install -r recipes/k562/requirements.txt
"$WORK/k562-env/bin/python" recipes/k562/reproduce.py \
  --cache "$WORK/k562-cache" --out "$WORK/k562-01"

# Lung uses numpy and pandas from the first environment and R with limma.
"$WORK/glioma-env/bin/python" recipes/lung/reproduce.py \
  --rscript Rscript --r-library "$HOME/.local/share/kilnbio-research/R-library" \
  --out "$WORK/lung-01"
```

**Sources.** The [glioma and Quartet manifest](glioma-quartet/sources.tsv) pins discovery summaries, Erasmus RNA ([GSE16011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16011)), [Huashan proteomics](https://doi.org/10.1038/s41467-024-54352-5), Quartet reference inputs and an HGNC archive by checksum. The [K562 manifest](k562/public_sources.json) pins published enhancer perturbations and ENCODE assay files (about 1.3 GB of compressed public inputs). The lung recipe pins [TCGA/cBioPortal discovery inputs](lung/provenance/discovery_sources.json) and [GSE31210 external inputs](lung/provenance/external_sources.json). The deliberate lung label conflict tests error detection; it is not evidence of a naturally occurring specimen swap.

For offline verification, populate the glioma or K562 cache with the named files; omit `--fetch` for glioma. Lung accepts `--cache PATH` pointing to a directory with `discovery/` and `external/` subdirectories matching the manifests. It verifies each source before that stage and keeps the external evaluation after the ranking freeze. Reproduction here describes these declared populations, endpoints and inputs; it does not imply that adding a given assay generally improves discovery.
