# KilnBio · research release

<img src="docs/brand/kilnbio-icon.png" alt="KilnBio kiln" width="100">

KilnBio helps researchers process multiple assays and measure whether an assay or analysis choice earns its place in a study. For example, if adding proteomics changes the targets you would follow up, KilnBio lets you separate the effect of adding the measurement from the effect of excluding donors who lack it. It records which specimens can be paired, which targets remain measurable and which directions agree in an independent cohort.

**Release:** [v0.2.0-research.1](https://github.com/varunorama/kilnbio-research/releases/tag/v0.2.0-research.1). This is a versioned distribution for approved **noncommercial research**. Read [permitted uses and access](ACCESS.md) and the [license](LICENSE). For commercial or industry-sponsored projects, contact Varun B. Dwaraka before use.

| Starting material | What you can do |
|---|---|
| Bulk RNA, ATAC, ChIP, WGBS/RRBS and 10x gene-expression FASTQs | Validate sample sheets and launch pinned nf-core workflows with explicit references, containers and run provenance. RNA includes an RSEM preset and a separate featureCounts route for aligned reads. |
| Counts, methylated/total cytosine counts, peaks or protein abundance matrices | Fit supported differential models; compare gene calls from DESeq2, edgeR and limma-voom on the same count input; compare methods on log2 protein abundance; inspect feature eligibility, sign and rank agreement. |
| Single-cell count matrices, Oxford Nanopore modification data, biomodal reports or methylation arrays | Run route-specific QC, pseudobulk analysis, modification-count import or compatible biomarkers. Chemistry, model features, required weights and input format are explicit; routes are not interchangeable. |
| Called genomic variants plus molecular assays | Extract panel events from VCF and associate them with specimen- or donor-linked traits; this route starts from calls, not DNA FASTQ. |
| Multiple assays or predictions from the same cohort | Audit specimen joins and sample loss; compare matched assay sets or predictions; evaluate frozen candidate directions against held-out evidence. Agreement alone is not independent replication. |

### Get started

1. Obtain approval and the ZIP plus matching checksum from the maintainer.
2. In Linux or WSL2, verify the file with `sha256sum -c kilnbio-research-0.2.0-research.1.zip.sha256`, then unzip it.
3. From the extracted directory, run:

```bash
bash install.sh "$HOME/.local/share/kilnbio-research" no-runtime
KILN="$HOME/.local/share/kilnbio-research/bin/kilnbio"
"$KILN" --help
"$KILN" assay-audit reach examples/assay_roster.tsv --out /tmp/kilnbio-reach-example.json
```

The last command needs only Python and a new output path. To run real data, see [installation, inputs and dependencies](docs/USER_GUIDE.md), the [command index](CLI.txt), and [worked public-data recipes](recipes/README.md). The recipes acquire published inputs and calculate outputs on your machine.

For the accompanying study, cite the exact tag above. The four primary worked examples are glioma, Quartet, K562 and lung adenocarcinoma. Additional [biological controls](benchmarks/biological/PROTOCOL.md) and the [peak-caller comparison protocol](benchmarks/chromatin/PROTOCOL.md) are also included. These entry points and the supported preprocessing routes do not all have the same extent of biological validation. Later releases receive new tags.

KilnBio supplies the launcher, audits, method comparisons and research recipes. The underlying aligners, statistical packages, nf-core workflows, models and public data have their own authors and terms. Cite those used in your analysis as well as [this release](CITATION.cff); see [third-party notices](THIRD_PARTY_NOTICES.txt). Research use only.
