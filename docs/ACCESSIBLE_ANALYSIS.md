# Start with a recipe, expand when needed

Research use only; not validated for clinical diagnosis or patient management. See [intended use](RESEARCH_USE.md) and the [cross-omic integration path](INTEGRATING_OMICS.md).

KilnBio's purpose is to let a researcher run a documented analysis without writing
their own aligner orchestration, statistical scripts or cross-assay matching code.
More supported methods should expand the available questions, not make every user
configure every method. This remains a command-line toolkit, not a point-and-click
application or an autonomous analyst.

| Starting point | What a supported recipe helps you do | Where expertise still matters |
|---|---|---|
| Wet-lab researcher with processed counts or abundances | Fill in a sample sheet, run an example, inspect QC/exclusions and export results | Biological unit, contrast, confounding, normalization provenance and interpretation |
| New computational biologist | Process a supported assay, compare methods on the same features, connect matched assays | Reference/chemistry choices, resource sizing, unusual data and statistical assumptions |
| Experienced bioinformatician | Reuse the input contracts and checks, change a documented method, construct reproducible comparisons | Validate new adapters, study-specific models and claims on independent endpoints |

For a first run, choose one documented recipe and use its supplied example before
substituting your data. A small processed-data example is usually easier to inspect
than a full FASTQ workflow. Complete preprocessing also needs Linux/WSL, a container
runtime, references and sufficient storage. The direct-analysis routes need their
listed R/Python packages; installing the launcher alone does not install all of them.

The [Olink/shared-contrast guide](MULTIOMIC_AND_OLINK.md) includes a synthetic
RNA/protein/methylation example. It demonstrates how to connect valid feature
results to cross-assay evidence. It does not establish a discovery in real people.
The [exome/CH guide](EXOME_AND_CH.md) explains the stricter reference, assay and
callability requirements for those routes. Optional [AlphaGenome annotations](ALPHAGENOME.md)
are a separate prediction step and are not needed for these analyses.

Users should receive effects, uncertainty, QC, reasons for exclusions and provenance,
not only a significant-feature list. Agreement across methods can show robustness
to an analysis choice; agreement across correlated assays is not independent
replication. A recipe can enforce useful checks but cannot recover information lost
when biological condition and technical batch are completely confounded.

## Runtime policy

The pipeline does not invoke general-purpose LLMs to write code, choose results,
interpret data or execute an analysis. R, Python and command-line programs implement
the workflows. Explicitly selected scientific models are permitted: for example,
DeepVariant and the GATK CNN filter within their configured routes, or the optional
AlphaGenome service. These are not required across all workflows. No hidden LLM
service is a dependency of a recipe.

`python3 ci/check_runtime_policy.py` checks shipped source for known LLM imports,
endpoints and CLI invocation patterns. This is a regression check, not a proof about
all third-party dependencies or a network firewall. Documentation about using a
coding assistant externally does not make that assistant part of pipeline runtime.
