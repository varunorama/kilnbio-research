# Running a research release

Start with the [automatic installation guide](INSTALLATION.md) to download Python, R packages and the sequencing runtime for your selected profiles. The environment choices below also describe manual setup for an existing installation.

Research use only; not validated for clinical diagnosis or patient management. See [intended use](RESEARCH_USE.md), [configured study orchestration](STUDY_RUNNER.md) and the [cross-omic integration path](INTEGRATING_OMICS.md).

Keep the unpacked release separate from inputs, reference data, scratch space and new results. Record the release tag in your analysis. Use the ZIP and SHA-256 from the same delivery; compare the digest with the one in the maintainer's approval or a separate trusted channel.

### Choose an environment

- The Bash launcher runs on Linux or WSL2. Use Python 3.12 for the worked recipes. The lightweight installer (`bash install.sh "$HOME/.local/share/kilnbio-research" no-runtime`) copies this version into its own directory. It does not alter your shell profile.
- For FASTQ workflows, install a working Docker, Apptainer or Singularity and provide the reference files and disk space appropriate to the assay. On x86-64 Linux, the installer mode `with-runtime` also downloads Java 21 and Nextflow 25.10.4. Container images and workflow code are fetched when needed.
- For direct R methods, use a project-specific R library. The reviewed baseline is R 4.6.2 and Bioconductor 3.23. From the release directory run `Rscript --vanilla environments/setup_r.R "$HOME/.local/share/kilnbio-research/R-library"`; set `KILNBIO_R_LIBRARY` to that path and, if necessary, `KILNBIO_RSCRIPT` to your Rscript executable. The setup script records resolved versions; [tested versions](../environments/R-tested-versions.csv) are supplied for comparison.
- Python dependencies vary by command. Set `KILNBIO_PYTHON` to your virtual environment's Python. [Single-cell requirements](../environments/analysis-requirements.txt) and separate [recipe environments](../recipes/README.md) are included. Clocks and modification callers can require additional upstream tools, weights and permissions; a lightweight install does not install these.

### From a sample sheet to a run

Start from the matching CSV and JSON files under [examples/](../examples). Replace example paths, specify the genome assembly and annotation, and check assay-specific fields (such as ChIP controls, library strandedness or methylation library protocol). A dry run validates the sheet and produces a command without launching Nextflow. Supply actual input files even for this check.

```bash
KILN="$HOME/.local/share/kilnbio-research/bin/kilnbio"
"$KILN" rnaseq my-samples.csv my-params.json runs/rna-plan --profile docker --dry-run
"$KILN" rnaseq my-samples.csv my-params.json runs/rna-analysis --profile docker --cpus 8 --memory-gb 32
```

The run directory records resolved parameters, revision and execution reports. Use separate destinations when changing a scientific comparison. `--resume` is for restarting the same Nextflow run. [versions.json](../versions.json) identifies the pipeline revisions; [CLI.txt](../CLI.txt) lists direct methods. The differences among these input types matter: use unnormalized integer counts for the RNA method comparison and continuous log2 measurements for its protein counterpart. Same-input callers are correlated; method agreement does not create additional biological replicates.

### Check what a second assay costs

In the [synthetic roster](../examples/assay_roster.tsv), three donors have the baseline assay and two have both. The README example should report `donors: 4`, `base_reach: 3` and `combined_reach: 2`. That denominator is the first thing to inspect before comparing a one-assay study to a two-assay subset. For a real study, use `assay-audit pair --help` to check donor, specimen and tissue keys; a consistent label alone does not prove biological identity.

### EpiDiff and third-party models

KilnBio does not distribute the EpiDiff forest or its probe map. Obtain the [pinned author package](https://github.com/aet21/EpiDiff/tree/631dd1ae039992865cb1788e98e9857f120ff1e5) under its terms. Export the sorted probe IDs (`candDiffALL2.v` in the package's `rfTRdiff2` data) as a text file with one ID per line. Obtain the matching [HM450 manifest](https://zhouserver.research.chop.edu/InfiniumAnnotation/HM450/). Then:

```bash
"$KILN" epidiff-coverage make-map --assembly hg38 --ids ids.txt --manifest manifest.tsv.gz --out epi-hg38.tsv
"$KILN" epidiff-coverage panel --assembly hg38 --loci epi-hg38.tsv --bed capture.bed --out coverage.json
```

These commands verify the source files and assembly-specific map. A panel intersection describes potential coverage; it does not validate a differentiation score on sequencing data. For all biomarker routes, confirm the specimen, units, assay chemistry, model terms and feature coverage before interpreting a score.

### Tests and citations

Run `python3 -m unittest discover -s tests` on Linux after unpacking. An R integration test requires the R environment described above. These synthetic tests check software behavior; they are not assay-wide biological validation. Use the [worked recipes](../recipes/README.md) for named scientific examples. Cite the package version, upstream workflows and aligners actually used, statistical methods and the original sources of any data or weights. See [upstream citations](../upstream_citations/) and [third-party notices](../THIRD_PARTY_NOTICES.txt).
