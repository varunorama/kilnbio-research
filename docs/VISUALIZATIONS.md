# Automatic analysis reports

KilnBio writes figures as part of `multiomic-meta`, `enrich` and `pathway-compare`.
No separate plotting command, notebook or assistant is required. Reports run
locally and make no LLM calls or external network requests.

Install the analysis dependencies with `bash setup.sh --profile analysis`.
Rerun setup from the updated checkout to install this version; an older installed
bundle will retain its original behavior. Existing runs and archived releases
are not rewritten.

## What is generated automatically

| Command | Figures | What to open |
|---|---|---|
| `multiomic-meta` | Gene-by-assay directional evidence; joint ranking; pairwise concordance counts; measurement patterns; leave-one-assay-out rank sensitivity | `NEW_OUTPUT/report.html` |
| `enrich` | Pathway NES ranking; leading-edge gene membership; up to six fgsea running enrichment curves | `NEW_OUTPUT/report.html` and `enrichment-curves.pdf` |
| `pathway-compare` | Pathway-by-assay evidence; measured pathway coverage; shared versus assay-specific contributing genes; pathway overlap network | `NEW_OUTPUT/report.html` |

Each report also includes a multipage `visualizations.pdf`, individual SVG and
PNG files in `figures/`, data in `plotted_data/`, and `visualization.json` with
captions, display choices, software versions, source hash and artifact hashes.
The HTML embeds its images and works offline. Keep the output folder together
to retain its downloadable TSV/PDF/SVG links. Nothing is uploaded by the report.

Existing `de` output includes `results/QC.pdf`; `scrna` includes its exploratory
UMAP. FASTQ workflows have their upstream QC reports. The [visualization
review](VISUALIZATION_REVIEW.md) distinguishes those outputs from additional
visualizations that require future adapters or fitted models.

## A complete synthetic example

From the source checkout, using the Python environment installed for analysis:

```bash
python3 examples/pathways/make_demo.py ../pathway-demo
kilnbio study ../pathway-demo/study.json ../pathway-demo-run
# Open ../pathway-demo-run/pathways/report.html in a browser.
```

The generator uses only the Python standard library. It creates fictional
pathways and p-values for teaching, including missing measurements and opposing
directions. It is a software example, not biological validation. The study
runner records the report outputs and can also be executed through the configured
KilnBio assistant connector.

For matrices and fitted statistics, use the [RNA/Olink/methylation
example](MULTIOMIC_AND_OLINK.md#run-a-complete-teaching-example). The single-command study
generator is `examples/multiomic/make_study_demo.py`; its final joint stage now
generates and records the report automatically.

## Gene-level figures: interpretation

The evidence heatmap shows the sign of the declared orientation times
`-log10(layer p)`. These are gene-level p-values after the existing feature
aggregation, not raw fold changes. Gray indicates missing or untestable evidence.
Prespecify orientations and feature maps; do not automatically invert all
methylation or assume every accessible region activates its nearest gene.

The ranking plots at most 30 genes, ordered by joint raw p, with BY q on the
horizontal axis. Display scores are capped at 16 to keep figures legible;
statistical tables retain the original numbers. All declared genes remain in the
testing family, including untestable genes assigned p=1 for multiplicity control.

Concordance counts cover all declared genes and show same direction, opposite
direction, tied direction, and absence of paired evidence separately. Agreement
between assays from the same donors is not independent replication.

Rank sensitivity removes one assay at a time and holds the required number of
layers, r, fixed. Both rankings use exactly the common testable genes and average
ranks for ties. Positive rank change means a gene ranks better with the removed
assay restored. Lost testability is recorded separately in
`plotted_data/rank_eligibility.tsv`. With two assays and r=2 this comparison is
not estimable; the report says so. These are evidence-rank sensitivity plots,
not estimates of held-out predictive performance or causal assay benefit.

## Compare pathways across assays

First perform the appropriate within-assay analysis using the same scientific
contrast. Choose a frozen gene-set collection and identifier system. RNA and
protein can use full signed gene rankings with `enrich`. Methylation and genomic
regions need bias-aware upstream models and reviewed mappings; see the review.

For two completed `enrich` runs using the identical GMT:

```json
{
  "gmt": "sets.gmt",
  "collection": "My frozen collection and version",
  "species": "human",
  "gene_id_type": "Ensembl gene ID without version",
  "contrast": "treated versus control, same prespecified donor comparison",
  "annotation_reviewed": true,
  "min_assays": 2,
  "alpha": 0.05,
  "joint_test": false,
  "layers": [
    {"layer": "RNA", "enrichment_dir": "rna_enrichment"},
    {"layer": "protein", "enrichment_dir": "protein_enrichment"}
  ]
}
```

```bash
kilnbio pathway-compare pathway_config.json NEW_OUTPUT
```

Paths are relative to the configuration. The adapter checks completed-run status,
GMT identity, species and identifiers, reads all pathways, restores size-filtered
sets as untestable, and obtains each assay's background from its complete ranking.
The original rank and GMT files must remain available. It uses the original
raw fgsea p-values and recomputes within-assay BY over the full declared GMT,
so its `q_within` can differ from the BH `padj` in the original fgsea report.

Other upstream methods can supply the same interface:

```json
{
  "layer": "methylation",
  "results": "methylation_pathways.tsv",
  "background": "methylation_measured_genes.tsv",
  "method": "Name, version and settings of the reviewed upstream pathway test",
  "score_kind": "unsigned",
  "complete_family_reviewed": true
}
```

The result TSV has exactly `pathway`, `p`, `score`, `contributing_genes`. Include
every pathway from the GMT once; use `NA` p/score for untestable pathways.
`p` is the valid raw pathway-test p-value, never an adjusted q-value.
`score` is a descriptive upstream score; `score_kind` is `signed` or `unsigned`.
Contributing genes are semicolon-separated hits or leading-edge genes, or empty
if the upstream test does not identify them. They must occur in both the GMT
pathway and the assay background. The background TSV has one `gene_id` column
and lists the full set of eligible measured genes, not just significant genes.

This is an import interface for correctly analyzed GOmeth/GREAT/other results;
it does not run those programs or repair a biased upstream test. Use the same
pathway identifiers and gene-set definitions for all assays. Verify sample
matching, feature mappings, contrasts, covariates and upstream test validity.

### Optional joint pathway evidence

Default `joint_test=false` produces descriptive comparisons. Setting it to true
explicitly requests the established **nondirectional Bonferroni
partial-conjunction test**, followed by Benjamini–Yekutieli correction over all
declared pathways. For m declared assays and r required assays:

`p_joint = min(1, (m-r+1) * p_(r))`

Here `p_(r)` is the rth smallest raw pathway p-value, with untestable assays
assigned p=1. Pathways with fewer than r tested assays remain untestable in the
output but stay in the multiplicity family. Fixed m is deliberately conservative.
This combination tolerates arbitrary dependence if the input p-values are valid
under the respective upstream nulls. [Benjamini and Heller, 2008](https://pubmed.ncbi.nlm.nih.gov/18261164/);
[Benjamini and Yekutieli, 2001](https://doi.org/10.1214/aos/1013699998).

The conclusion concerns pathway evidence in at least r assays. It does not
require the same genes or the same direction across assays. The figures expose
those differences. A signed NES is not a comparable effect size across modalities,
and positive unsigned enrichment is not evidence of pathway activation.
This is not fixed/random-effects meta-analysis across studies. A forest plot of
pooled effects would require another model and compatible effect estimates.

### Which figures answer which questions?

* **Pathway evidence:** which assays contribute evidence for each pathway?
* **Coverage:** could the assay measure enough of that pathway to test it?
* **Shared gene support:** is agreement driven by the same genes or different ones?
* **Pathway overlap:** are apparently separate findings largely overlapping
  annotation sets? Edges are GMT Jaccard overlap, not inferred molecular interactions.

Automatic analysis remains explicit Python and R code. These reports do not
generate biological narratives with language models.
