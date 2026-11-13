# Olink and a shared question across molecular assays

These commands add delivered Olink NPX to KilnBio and fit the same declared
contrast across RNA, protein and other eligible assay matrices. They do not use
LLM calls. Alignment, R statistics and Python integration remain explicit code.

| Command | Input | Output |
|---|---|---|
| `olink-import` | Delivered long NPX export and explicit sample roster | NPX matrix, eligibility, measurement QC, technical-design audit and input hashes |
| `multiomic-layers` | Matrices, donor/specimen maps and a model specification | limma/voom feature effects and two-sided p-values, design audits and session information |
| `multiomic-meta` | Feature results, prespecified gene maps and directions | Directional partial-conjunction tests, BY-adjusted values and nominations |

The statistical methods are established. KilnBio contributes the interfaces,
design checks, missingness policy and traceable connection between the steps.
This is neither a new statistical test nor proof that an integrated target is
better than a single-assay target. That requires a separate held-out endpoint.

## Install the direct-analysis dependencies

Use a dedicated Python environment with NumPy, pandas and SciPy; add pyarrow for
Parquet. R needs limma, edgeR and jsonlite. The complete preprocessing installation
does not install these direct-analysis packages automatically. Point
`KILNBIO_PYTHON` and `KILNBIO_RSCRIPT` at the chosen environments; an optional
`KILNBIO_R_LIBRARY` selects a private R library.

```bash
python3 -m venv ../kilnbio-stats-env
../kilnbio-stats-env/bin/python -m pip install -r environments/multiomic-requirements.txt
export KILNBIO_PYTHON="$PWD/../kilnbio-stats-env/bin/python"
Rscript -e 'if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager",repos="https://cloud.r-project.org"); BiocManager::install(c("limma","edgeR","jsonlite"),ask=FALSE,update=FALSE)'
```

Record package versions for a study. The R command installs into the selected R
library; use a private library if managing several projects.

## Olink import

```bash
kilnbio olink-import delivered_NPX.csv --samples samples.tsv --out new_olink
# Optional, prespecified detection fraction; missing LOD then makes an assay ineligible:
kilnbio olink-import delivered_NPX.parquet --samples samples.tsv --min-detect 0.5 --out new_filtered
```

The roster starts with `sample` matching `SampleID` exactly and may include
design columns. The long CSV/TSV or Parquet needs `SampleID`, `OlinkID`, `NPX`.
Only requested samples enter the output. Known sample-type fields exclude
controls. If sample type is absent, the user must supply a biological roster.
Older `Sample_Type`, `QC_Warning` and `Assay_Warning` names are recognized.
Two measurements of the same SampleID/OlinkID are rejected rather than averaged.
A sample can occur on different panels/plates when OlinkIDs are distinct;
the measurement audit retains those identities and the roster records their union.

The default retains delivered NPX, including its existing vendor normalization.
NPX is a relative log2 scale, not concentration; a one-unit contrast is a twofold
change on that relative scale. This does not make different Olink products,
projects or proteins interchangeable. Cross-project bridging requires a defined
reference/bridge design and the relevant Olink procedure, not column-wise z-scores.
See the [official OlinkAnalyze documentation](https://olink-proteomics.github.io/OlinkRPackage/).

`SampleQC=FAIL` excludes the sample. `AssayQC=FAIL` creates a missing measurement;
warnings remain visible. `npx.tsv` preserves missing values; `eligible_npx.tsv`
contains complete assays in the retained roster. `features.tsv` gives reasons.
Missing QC is marked unreported, not PASS. No values are imputed.

LOD is taken from the delivered NPX-scale `LOD` column, or an explicit TSV with
unique `OlinkID` and one of `LOD_NPX`/`LOD_ExtNPX`. ExtNPX-scale LOD needs delivered
ExtNPX for conversion. No LOD is invented from a few controls. Without LOD and
without `--min-detect`, detection is reported unknown, rather than used to discard
measurements. Olink detection rules depend on product and normalization; see the
[official detection-limit vignette](https://cran.r-project.org/web/packages/OlinkAnalyzeVignettes/vignettes/lod.html).

`--renormalize study` is an optional single-plate centering sensitivity analysis,
not the default or a validated batch correction. It requires complete ExtNPX and
a constant delivered ExtNPX–NPX offset per assay; LOD moves by the same offset.
It cannot resolve a confounded study design.

`aliasing.tsv` distinguishes categorical association (Cramér's V) from exact
rank confounding. V≥0.5 alone does not establish non-identifiability. This compact
audit is not a comprehensive batch model. `measurement_qc.tsv` retains plate/well
IDs and row-level QC. Review failed samples and eligibility before modeling.

Raw Olink barcode counts/BCL/FASTQ decoding, proprietary NPX Map resources and
cross-product normalization are outside this importer. The public interface
starts at an already delivered NPX export.

## Fit the same contrast in each assay

`kilnbio multiomic-layers config.json NEW_OUTPUT` accepts:

* `donors`: TSV with unique first column `donor`, the trait and covariates.
* `estimand`: `association` or `paired`.
* Association: numeric `trait`, numeric `covariates` (explicitly encode categories),
  one independent specimen per donor per layer.
* Paired: `condition` with `column`, `reference`, `test`; exactly one specimen per
  donor/condition, donor fixed effects. Additional covariates are not implemented.
* `cohort_policy`: `matched` by default, requiring identical donor sets across
  layers. `available` permits different sets, explicitly recorded in the audit;
  this changes populations and needs a separate matched-cohort sensitivity check.
* Each `layers` entry: `name`, `matrix`, `specimens`, `scale`, `feature_map`,
  `orientation`. Matrices start with `feature`; specimen maps have `specimen,donor`
  and, for paired data, the declared condition column. All columns must map exactly.

`counts` uses integer counts, edgeR expression filtering/TMM and limma-voom.
`log2` and `mvalue` use limma with a mean–variance trend. `beta` converts array-style
beta values to M values (clip at 10^-6 and 1−10^-6), then fits limma. This beta
route does **not** model WGBS/RRBS read depth. Use DSS or an appropriate count model
first and supply reviewed results to `multiomic-meta`. Similarly, global chromatin
shifts, spike-in offsets, repeated tissues and scRNA cells need their dedicated
models or donor pseudobulk, not automatic treatment as ordinary independent rows.

All declared features remain in the output; untestable features have missing
p-values. A rank-deficient model fails. Results contain real coefficients, not
signed z-scores relabeled as fold changes. `*_design.tsv` and `layers_design.tsv`
record which donors and features were used. Study-specific mixed models and
exchangeability-aware permutations are not implemented here.

## Test cross-assay evidence

```bash
kilnbio multiomic-meta layers/manifest.tsv --min-layers 2 --out new_joint
```

Manifest TSV columns: `layer,results,feature_map,orientation`, with optional
positive `weight` used only for exploratory ACAT. Results have unique
`feature,effect,p`; p-values must be two-sided from a valid symmetric directional
test. Paths are relative to the manifest. Maps have `feature,gene`, or use
`identity` when the feature already is the intended gene identifier. Olink assay
labels and UniProt IDs are not silently converted to gene symbols.

For each gene and direction, Bonferroni combines the prespecified mapped
features within an assay. A layer with missing mapped features is untestable for
that gene. For m measured assays, the partial-conjunction p-value in one
direction is `min(1,(m-r+1)*p_(r))`. The smaller directional result is doubled
for testing both directions. BY adjustment covers **all genes in the declared
maps**, including untestable genes as p=1 for the adjustment. This conservative
construction does not require assay independence, but does require valid
directional feature tests and justified measurement/eligibility selection.

The null is that fewer than r layers support the specified directional
association. It is not the global null of no association anywhere. In contrast,
ACAT is recorded as an exploratory global-null diagnostic and cannot be used
as evidence of convergence or as the module's discovery selection rule.

Freeze gene maps, r, orientations, assay selection and eligibility before looking
at the effects. Methylation is not universally inverse to expression: an
orientation of −1 needs a context-specific hypothesis, not a convenient sign
flip after analysis. Same-donor measurements are correlated and do not constitute
independent replication. Partial conjunction also does not prove mediation,
causality or a shared molecular mechanism.

`genes.tsv` records per-assay directional and combined evidence, number of
measured layers, q-values and classification. `nominations.tsv` records compatible
directions for follow-up; it does not claim every listed assay independently
passes its own FDR threshold. Missing genes remain visible as untestable.

## Run a complete teaching example

From a source checkout with the dependencies above:

```bash
"$KILNBIO_PYTHON" examples/multiomic/make_synthetic.py ../multiomic-demo
bash run.sh olink-import ../multiomic-demo/olink_export.tsv --samples ../multiomic-demo/olink_samples.tsv --out ../multiomic-demo/olink
bash run.sh multiomic-layers ../multiomic-demo/config.json ../multiomic-demo/layers
bash run.sh multiomic-meta ../multiomic-demo/layers/manifest.tsv --min-layers 2 --out ../multiomic-demo/joint
```

This simulates an exome-derived carrier trait linked to RNA, Olink and methylation.
It teaches integration; it contains no real mutations and validates neither
variant calling nor biological target discovery. For public-data candidates see
[the exome and CH guide](EXOME_AND_CH.md#public-data-examples-and-manuscript-scope).

## Methods and review notes

* Ritchie et al. (2015), limma: [10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007).
* Law et al. (2014), voom: [10.1186/gb-2014-15-2-r29](https://doi.org/10.1186/gb-2014-15-2-r29).
* Robinson and Oshlack (2010), TMM: [10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25).
* Benjamini and Heller (2008), partial conjunction: [10.1111/j.1541-0420.2007.00984.x](https://doi.org/10.1111/j.1541-0420.2007.00984.x), PMID 18261164.
* Benjamini and Yekutieli (2001), dependent multiple testing: [10.1214/aos/1013699998](https://doi.org/10.1214/aos/1013699998).
* Liu and Xie (2020), ACAT: [10.1080/01621459.2018.1554485](https://doi.org/10.1080/01621459.2018.1554485), PMID 33012899.

The supplied Olink prototype was reviewed and revised: vendor normalization is
the default; an approximate control-based LOD is not inferred; categorical
association is not called exact confounding; effects remain model coefficients;
and observed cross-assay correlation does not select an independence-dependent
test. The prototype's adaptive-filtering power claims and covariate-ignoring
joint trait permutations are not release claims. These are substantive method
changes, not merely renamed implementation details.
