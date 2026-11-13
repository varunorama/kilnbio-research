# Olink, public exome methods and CH: evidence register

Reviewed 12 November 2026. [Live PubMed metadata](platform_extensions_pubmed.json)
were fetched through the local PubMed MCP adapter. Each record includes its
implementation role, exact DOI and indexed notice links. Metadata verification
does not establish correct implementation or novelty.

| Source / evidence | Implementation decision | Boundary still requiring evidence |
|---|---|---|
| Sarek 3.12.0 release, parameter/input schemas and published Sarek paper | Versioned autosomal capture adapter; open aligners/callers; explicit reference and resource provenance | Other caller branches, full exome accuracy, sex chromosomes and CH low-VAF detection |
| Olink public export examples and detection documentation | Preserve vendor NPX, product-specific LOD supplied rather than invented, explicit roster/QC | Cross-product bridging and batch correction need a separate design and validation |
| limma, voom and TMM | Real coefficients for shared donor-level contrasts; integer counts retain a count-aware mean–variance model | Study-specific repeated measures, confounding, compositional/global shifts and prospective eligibility |
| Benjamini–Heller partial conjunction; Benjamini–Yekutieli | Directional Bonferroni feature aggregation and partial conjunction, BY across the declared gene family | Power may be conservative; dependent-null simulation is a check, not proof of every input model's calibration |
| Liu–Xie Cauchy combination | ACAT retained as exploratory global-null diagnostic only | A combined global-null p-value is not evidence that multiple assays replicate |
| Uddin et al. 2024 CH methods | Explicit depth/read/VAF/population-frequency screening with curated alleles | Not all published exceptions/rescue steps; no clinical classification or somatic-origin validation |
| Ismail et al. 2025 issue, GSE210435 and Olink supplement | Candidate genotype/RNA/chromatin/methylation/protein public example | Download and match the individual files before freezing analysis; existing IL32 results are prior findings |
| GA4GH benchmarking recommendations and GIAB | Restrict future WES accuracy benchmarks to capture × confident-truth regions | No precision/recall claim from a tiny or empty callset |

The integration deliberately changes several assumptions of the supplied Olink
prototype. The [user guide](../docs/MULTIOMIC_AND_OLINK.md#methods-and-review-notes)
records the method changes; the prototype's adaptive-filtering performance is
not evidence for this implementation. [Exome/CH scope](../docs/EXOME_AND_CH.md)
separates public methods from private return files. The latter are not evidence
of an openly reproducible benchmark and are not included in this repository.

Reusable software checks are `tests/test_exome.py`,
`tests/integration_platform_extensions.py`, `tests/integration_multiomic_layers.R`,
`tests/public_olink_smoke.R` and `tests/public_exome_smoke.py`. The full logs,
failed attempts, numerical observations and manuscript notes are retained in
the private local archive. Existing scientific release gates remain in force;
adding citations does not automatically satisfy them.
