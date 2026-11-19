# Scientific report design verification — November 18, 2026

This update changes presentation and interpretation aids. Statistical analyses,
FDR procedures, candidate-selection rules and manuscript results are unchanged.

## Delivered design

- Shared typography, white space, direction colors and export styling across
  typed reports; original upstream R diagnostic plots keep their native design.
- Aligned gene-by-assay evidence, measured-layer coverage and joint BY evidence.
- Volcano plots with the effect distribution, displayed denominators, direction
  groups and bounded feature annotations; no FDR line invented on a raw-p axis.
- Exclusive method-call combinations with counts, and explicit assay coverage
  combinations. Missing measurements remain visually distinct from zero calls.
- Detection-fraction legends and displayed cluster counts for single-cell dots.
- Offline HTML navigation, responsive cards, direct downloads and explanations
  of why each panel matters and how to read it. Interpretation text is fixed in
  code; it does not invoke an LLM or infer a finding from a picture.
- Refreshed fictional gallery fixtures, including 24 genes across four assays.
  These illustrate varied coverage and direction. They are not new biological
  findings or replacements for the manuscript examples.

## Verification

The full research regression suite ran 129 tests: 127 passed and two optional
R/samtools checks were skipped. Tests compare statistical output bytes with and
without reports, exercise missing-only inputs and check report artifacts.
A new regression check preserves the heatmap's scientific color-bar label when
cell-value annotations are added.

All 119 gallery panels were exported, and all PDF pages were rendered for visual
review. Review identified and corrected a color-bar labeling collision and two
legend placements. Corrected panels were checked at full resolution. The gallery
contains 24 report families and three previously fitted native R diagnostics.
Every panel has an explanation and downloadable plotted values; 431 report
artifact hashes and the gallery/documentation links were verified.

The offline report was inspected in a headless browser at 1,440- and 600-pixel
widths. A clean portable installation from the private checkout generated the
new report through the installed CLI. Native R diagnostic images use the
previously verified synthetic fits; model fitting is not repeated by styling.

The current gallery generator and regression tests are included in the source.
Original run folders and versioned manuscript releases remain intact.
