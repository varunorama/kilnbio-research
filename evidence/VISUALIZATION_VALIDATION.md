# Automatic report validation — November 18, 2026

Scope: reusable software and synthetic inputs. These checks establish report
execution and selected statistical invariants, not biological discoveries.

* Public checkout regression suite: 121 tests run; 119 passed and two
  environment-dependent tests skipped. New report-specific tests ran (six test cases).
* Synthetic RNA/Olink/methylation study: import, real limma/voom layer models,
  directional gene integration and automatic report completed through `study`.
* Synthetic pathway study: configuration, declared input hashes, execution and
  report artifacts completed through the same study runner.
* Two signed synthetic gene rankings: real R `fgseaMultilevel` execution,
  running enrichment curves, Python reports and direct `enrichment_dir`
  cross-assay comparison completed.
* Figures were exported to SVG/PNG/PDF; the generated PDFs were rendered and
  visually inspected, including the running enrichment curves.

The gene-report test compares `genes.tsv` and `nominations.tsv` byte-for-byte
against execution with rendering disabled by a test mock. It checks missing
measurements, opposing directions, full-family concordance denominators,
HTML-safe labels, artifact hashes and failure cleanup. A two-assay/r=2 case
must report that leave-one-out ranking is not estimable.

The pathway tests check the fixed-family partial-conjunction calculation,
untestable rows, complete-family rejection, incompatible contributing genes,
order invariance and a boundary-null example with perfectly dependent null
p-values. Opposing signed pathway scores can support a nondirectional joint
test; the report preserves their signs. This is deliberately not called
directional pathway convergence.

Reproduce the reusable checks in an analysis environment:

```bash
python -m unittest discover -s tests -p 'test_*.py'
python examples/pathways/make_demo.py ../pathway-demo
kilnbio study ../pathway-demo/study.json ../pathway-run
python examples/multiomic/make_study_demo.py ../matrix-demo
kilnbio study ../matrix-demo/study.json ../matrix-run
```

Python reporting was exercised with Python 3.12, Matplotlib 3.11.2 and pandas
3.1.1. R integration used fgsea 1.38.0 and the existing managed Bioconductor
environment. Each actual report writes its own rendering versions and hashes.
MCP uses the study command interface; this change does not require a separate
assistant tool for each new report. An assistant-client UI run was not repeated
as part of these report-specific checks.
