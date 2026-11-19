# Automatic-report validation: November 18, 2026

This is software verification using synthetic inputs. It is not new biological
validation, a model-accuracy benchmark or a manuscript result.

## Implemented behavior

The public launcher generates typed offline reports for completed analysis
commands. Existing integrated-evidence/pathway reports remain transactional.
Other reports stage their output before publication, preserve completed analysis
tables on failure and return a failed command if reporting fails. Studies hash
the generated report assets as part of verified resume. Every installation
profile includes the shared Python plotting requirements.

The public gallery contains 119 rendered panels in 24 command families and
three native fitted R diagnostic examples. Shared plot types are reused across
commands, as identified in the guide. Fixtures, plotted values and captions are
public; no private manuscript data or observations are included.

## Checks performed

* Python regression suite: 128 tests run, 126 passed, two optional checks skipped.
  The skipped tests require optional R/samtools configuration in that unit-test run;
  they are not counted as passed.
* Private working-branch regression suite: 122 tests run, 120 passed and the same
  two optional checks skipped.
* Portable installation into a clean prefix completed, and the installed command
  generated a synthetic assay-audit report with its plotting dependencies.
* Actual DESeq2 and edgeR fits: 1,200 synthetic genes, eight independent samples,
  two groups. Both exported model-transformed matrices, native model diagnostics
  and completed automatic reports.
* Actual DSS fit: 400 synthetic covered CpGs, eight independent samples, two
  groups. Locus statistics, fractions, depths and automatic reports completed.
* Complete synthetic Olink/RNA/methylation study: 60 donors, 200 features per layer; import, three fitted layers
  and joint gene evidence completed with automatic reports and study receipts.
* Actual Scanpy processing: 180 synthetic cells from six libraries/donors and
  120 genes completed QC, fitted embedding and reports. Reviewed synthetic cell
  labels produced six donor-level pseudobulks and an automatic pseudobulk report.
* Native fgsea ran on a complete 200-gene synthetic signed ranking and three
  fictional sets. Ranking, leading-edge and native running curves were rendered.
* Report tests verify unchanged source statistics, missing values, bounded
  reads/columns, custom method-agreement alpha, failed-report preservation,
  command-line dispatch, and detection of a modified figure on study resume.
* All initial 118 gallery PDF pages were rendered with Poppler and visually
  inspected. The revised separated AlphaGenome score panels and the three native
  R diagnostics were inspected individually. PNG/TSV links and artifact hashes
  were checked. Different regulatory scoring strata are not pooled on one scale.

Reproduce with the configured analysis/singlecell dependencies:

```bash
python -m unittest discover -s tests -p 'test_*.py'
python tests/integration_auto_reports.py REPO_ROOT NEW_MODEL_CHECK
python tests/integration_scrna_reports.py REPO_ROOT NEW_SINGLECELL_CHECK
python examples/visualizations/make_gallery.py NEW_GALLERY --with-enrichment
```

The integration scripts fit real software models to generated synthetic inputs;
the gallery primarily tests renderers using fictional output fixtures. These
are different levels of verification. Raw FASTQ workflows, remote scientific
model APIs and every optional clock were not rerun for this reporting change.
Their native QC/previous validation scope remains documented separately.

## Scientific boundaries

Figures retain model-specific units, directions, measured backgrounds and
missingness. Plotting does not add inferential tests. Color thresholds are
labeled; statistical output tables remain authoritative. Bounded views declare
their limits and do not claim whole-genome or whole-cohort coverage when capped.
Method agreement does not become a new FDR procedure. Region-only DMR results
do not acquire invented single-CpG tests. Imported pathway evidence still needs
a valid assay-specific upstream test; this change does not execute GOmeth/GREAT.
Unsupported analyses such as 16S, velocity and spatial inference are not implied
by the expanded gallery. See the cited [methods review](../docs/VISUALIZATION_REVIEW.md).
