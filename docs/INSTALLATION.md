# Install KilnBio

The managed installer supports Linux x86-64 and Ubuntu in Windows WSL2. It creates an isolated installation below `$HOME/.local/share/kilnbio`. Run setup from the downloaded source directory.

## Prepare your computer

On Ubuntu, install the small bootstrap prerequisites once:

```bash
sudo apt-get update
sudo apt-get install -y python3 curl ca-certificates
```

Python 3.8+ is sufficient to start setup. KilnBio's managed environment uses Python 3.12. On Windows, install WSL2 with `wsl --install -d Ubuntu` in an Administrator PowerShell, restart if Windows requests it, then use the Ubuntu terminal for the Linux commands. [Microsoft WSL installation](https://learn.microsoft.com/en-us/windows/wsl/install)

Allow space for the selected environment and its package download cache. The core environment is small relative to the R/Bioconductor and sequencing profiles. All profiles together can use several GB before input data, container images, references or workflow scratch files. FASTQ storage requirements depend on the assay and study.

## Dependency profiles

```bash
bash setup.sh --plan --profile all
bash setup.sh
bash setup.sh --profile analysis
bash setup.sh --profile analysis --profile sequencing
bash setup.sh --profile all
```

Choose one of the install commands above. `--plan` prints dependencies without installing. Each selected profile includes core.

| Profile | What setup installs |
|---|---|
| `core` (default) | Python 3.12, pip, the KilnBio launcher and a synthetic donor-audit check |
| `bladder` | R 4.6.2, native build tools and checksum-verified limma 3.68.5/statmod 1.5.2 packages |
| `analysis` | Bladder tools, R/Bioconductor packages from `environments/setup_r.R`, and the NumPy/pandas/SciPy requirements for matrix and Olink analyses |
| `singlecell` | Analysis profile plus Scanpy, AnnData and Leiden requirements |
| `sequencing` | Java 21, Nextflow 25.10.4 and a check of the chosen container engine |
| `assistant` | Python dependencies for the optional MCP connector |
| `all` | The profiles above |

R source packages may take substantially longer to build than the core installation. The full Bioconductor profile follows release 3.23 and records the versions resolved by that repository. It uses the existing R setup script. The pinned bladder packages have their own checksums and version checks. Package providers must still supply the recorded versions.

Install to another location with `--prefix /path/to/kilnbio`. Repeat the same command after a transfer failure. Setup retains downloaded packages, verifies installed dependencies and performs a smoke test before switching the command to a completed installation. Dependencies are reused when their specification is unchanged. Different selections create separate environments, so include every profile you want active when rerunning setup.

## Check the installation

```bash
KILN="$HOME/.local/share/kilnbio/bin/kilnbio"
"$KILN" --help
bash setup.sh --check
"$KILN" assay-audit reach examples/assay_roster.tsv --out first-audit.json
```

The example returns four donors, three with the baseline assay and two with both. Choose a new output filename for another run.

Use `export PATH="$HOME/.local/share/kilnbio/bin:$PATH"` to make `kilnbio` available by its short name in your current terminal. You can add that line to your own shell configuration.

The installation receipt is `$HOME/.local/share/kilnbio/installation.json`. It identifies the source content hash, selected profiles, environment bundle and smoke-test result. The bundle records `conda-explicit.txt` and `pip-freeze.txt`. The R analysis library records package versions and session information. Retain these with the study's output receipts.

## Sequencing containers and references

Choose your container engine, then install the sequencing profile:

```bash
bash setup.sh --profile sequencing --container docker
bash setup.sh --check
```

Docker must be installed and running for your account. On a cluster, use your site's Apptainer or Singularity installation and pass `--container apptainer` or `--container singularity` to setup. The matching run option is `--profile docker`, `--profile apptainer` or `--profile singularity`. [Docker installation](https://docs.docker.com/engine/install/) · [Apptainer installation](https://apptainer.org/docs/admin/main/installation.html)

Setup reports a missing or unavailable engine. The installed core tools can run, while `--check` returns a failure for the sequencing profile until its engine is usable. The installer does not need administrator privileges for its managed environments. Installing a host container engine may require your administrator.

The pinned workflows fetch their analysis tool containers at execution. Configure the correct genome assembly, annotation, capture BED and controls in your study. Nextflow resource limits and references are supplied per run. Optional commercial/licensed tools, Nanopore basecaller models, biomarker weights and API credentials have route-specific installation or access steps. See [analysis inputs](USER_GUIDE.md) and [commands](../CLI.txt).

## Run the public bladder recipe with managed R

After `bash setup.sh --profile bladder`, run:

```bash
"$HOME/.local/share/kilnbio/bin/kilnbio" bladder-example --cache ./bladder-inputs --out ./bladder-run
```

This downloads and verifies the published study tables, runs the frozen recipe and records the outputs. A new output directory is required for each calculation. Add `--offline` to reuse the verified input cache.

## Connect an assistant

To register an installed command-line assistant after setting up the assistant profile, run one of:

```bash
"$HOME/.local/share/kilnbio/bin/kilnbio" assistant-install codex --demo
"$HOME/.local/share/kilnbio/bin/kilnbio" assistant-install claude --demo
"$HOME/.local/share/kilnbio/bin/kilnbio" assistant-install antigravity --demo
```

The selected assistant must already be installed and signed in. Registration creates a local synthetic project and connects its MCP server. For a client running on Windows with an analysis host in WSL, or for browser clients, follow the [assistant setup guide](ASSISTANT_SETUP.md).

Managed registration uses the same versioned environment as the command line, including your selected R and sequencing tools. Select every needed profile during setup, for example `--profile analysis --profile assistant`. After updating the installation, repeat registration to connect the client to the new version.

## Other systems and existing environments

Automatic environment provisioning is currently limited to Linux x86-64/WSL2. For macOS Terminal, follow the [Mac installation and verification guide](MACOS.md), which covers Python setup, the decision functions, optional R tools and sequencing limitations. That manual route has not yet been verified on a native Mac. The underlying Bash installation command is:

```bash
bash install.sh "$HOME/.local/share/kilnbio" no-runtime
```

Use `environments/setup_r.R` for the R baseline and the requirements files for the chosen Python route. The bladder walkthrough supports native Windows, macOS and Linux with R 4.6.x. Native Windows users should use WSL2 for the Bash sequencing launcher.

## Troubleshooting

- **Package solve or download fails.** Save the displayed package or URL error, check connectivity and retry the same command. Setup reports failure before replacing the active command.
- **A setup lock remains after interruption.** Check that setup has stopped, then remove only the empty `.setup-lock` directory inside your chosen prefix and rerun.
- **R package compilation fails.** Keep the compiler output. The bladder profile supplies C/C++/Fortran tools, and analysis builds use the same managed environment.
- **Docker is installed but unavailable.** Run `docker info` in the same terminal. Enable the daemon or Docker Desktop's WSL integration and check account permissions.
- **An optional command reports missing tools or weights.** Consult its route guide. The managed profiles list their coverage above.

The package manager is a checksum-verified Micromamba binary. Dependency downloads use conda-forge, Bioconda, CRAN, Bioconductor and the referenced Python package indexes. [Micromamba documentation](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)

## Automatic plotting dependencies

Every profile includes NumPy, pandas and Matplotlib for offline analysis reports.
This adds plotting to the core profile; direct manual installations should also
install `environments/report-requirements.txt`. Analysis-specific R, single-cell
and sequencing dependencies still come from their respective profiles. Reports
run automatically through `kilnbio`, including study/assistant execution.
[View examples](gallery/README.md).
