# Build biological context across assays

KilnBio's purpose is to help researchers integrate molecular evidence without writing a new processing and statistical workflow for every study. An additional assay is useful when it helps resolve a biological question: whether a genomic alteration accompanies an expression change, which regulatory regions are plausible drivers, whether an RNA signal reaches the protein level, or which candidate merits an experiment.

Integration does not mean putting all measurements onto one numerical scale. Variant genotypes, methylated-read fractions, accessibility counts, RNA counts and Olink NPX describe different quantities. KilnBio processes and models each appropriately, then connects them through declared biological units, genomic coordinates and feature maps.

## Follow one question through the framework

1. **Define the comparison.** Specify the phenotype, treatment, genomic event or cell state. Record donor, specimen, tissue, time point and assay identifiers. Match within the biological unit appropriate to the question, and retain the denominator lost when an assay is missing. Use the [study audits](../CLI.txt).
2. **Process each assay and retain its measurement meaning.** Use the [assay routes](../README.md), [Olink importer](MULTIOMIC_AND_OLINK.md), [modification imports](../README.md) or [genomic interfaces](EXOME_AND_CH.md). Coverage, QC failures and unmeasured features remain visible.
3. **Connect features explicitly.** Use gene identifiers for RNA/protein associations, declared CpG/region-to-gene maps for regulatory hypotheses, and reviewed transcript/protein mappings for variants. A nearest gene is a mapping assumption, not proof of regulation. Do not silently treat two different biopsies, isoforms or genome builds as the same entity.
4. **Estimate effects under the same scientific contrast.** `multiomic-layers` fits supported contrasts with assay-specific models; depth-aware methylation results can be supplied from DSS or dmrseq. `genomic-associate` provides a separate starter carrier-versus-reference comparison. Retain covariates, direction, scale, uncertainty and missingness rather than replacing everything with a z-score.
5. **Assess convergent and conflicting evidence.** `multiomic-meta` tests a prespecified requirement for directional support across layers using established partial-conjunction methods. It preserves the declared gene family, including untestable genes. Same-donor assays are not independent replication, and agreement does not establish a mechanism. Inspect disagreements as potential biology, mapping problems or measurement limitations.
6. **Add external context and test the hypothesis.** Join exact variant, transcript and protein identifiers to [VEP, ClinVar, AlphaMissense and AlphaFold Database annotations](VARIANT_ANNOTATIONS.md); optionally obtain [AlphaGenome predictions](ALPHAGENOME.md). Keep those annotations distinct from measured effects. Freeze candidate nominations before independent validation with `candidate-evidence`, or compare assay combinations on the same held-out donors.

The result is a traceable set of linked measurements, effects, annotations and validation records that can support a biological hypothesis. It is not an automatic causal explanation or a clinical interpretation. There is no universal biological direction linking methylation, accessibility, expression and protein abundance; declare any expected relationship for the relevant locus and context before testing it.

## A practical RNA and Olink entry point

```bash
kilnbio olink-import delivered_NPX.csv --samples samples.tsv --out olink
kilnbio multiomic-layers shared-contrast.json layers
kilnbio multiomic-meta layers/manifest.tsv --min-layers 2 --out joint-evidence
```

The [shared-contrast guide](MULTIOMIC_AND_OLINK.md) defines the sample maps, matrices, models and feature maps required by `shared-contrast.json`. Add other supported layers through that declared design, rather than concatenating incompatible raw matrices. The [genomic linkage recipe](EXOME_AND_CH.md) shows how called genomic states can define a comparison against RNA and methylation. The [public research recipes](../recipes/README.md) demonstrate regulatory-link ranking, candidate follow-up and held-out model evaluation.

A beginner can run a supplied recipe after preparing the required inputs. A bioinformatician is still needed for a new study's mapping assumptions, complex repeated-measure designs, confounding, model choice and interpretation. KilnBio reduces repeated implementation work; it does not infer those scientific decisions from a folder of files. Research-use limits apply to every route: [intended use](RESEARCH_USE.md).
