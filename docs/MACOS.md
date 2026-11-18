# Use KilnBio from macOS Terminal

Start with the Python decision methods and synthetic examples below. They do not need R, Nextflow, Docker or participant data. This is a manual installation route for Intel and Apple silicon Macs; it has not yet been verified by a native macOS installation test. The managed `bash setup.sh` installer currently supports Linux x86-64 and WSL2 only and deliberately rejects macOS.

These instructions apply to the current source repository, including `decision-methods`. The older manuscript release archives do not contain those new functions. Follow the repository's existing [access and license terms](../ACCESS.md).

## Prepare Python and obtain the source

Use Terminal in your normal native architecture. On Apple silicon, keep Python and its packages on ARM64 rather than mixing an Intel Python with ARM64 libraries. The commands use `brew --prefix` instead of assuming either `/opt/homebrew` or `/usr/local`.

Install [Homebrew](https://brew.sh/) if needed, using its official instructions, then run:

```bash
brew install python@3.12 git
git clone https://github.com/varunorama/kilnbio-research.git
cd kilnbio-research

"$(brew --prefix python@3.12)/bin/python3.12" -m venv "$HOME/.venvs/kilnbio"
source "$HOME/.venvs/kilnbio/bin/activate"
python -m pip install --upgrade pip
python -m pip install -r environments/multiomic-requirements.txt
python -m pip check
```

If you already have the repository, enter that checkout instead of cloning over it. Use an existing compatible Python 3.12 installation instead of Homebrew if preferred. Keep package installation inside the virtual environment; the [Homebrew Python guidance](https://docs.brew.sh/Homebrew-and-Python) explains this separation. [Python 3.12 formula](https://formulae.brew.sh/formula/python@3.12)

## Install the command and select its Python

From the repository directory with the environment active:

```bash
bash install.sh "$HOME/.local/share/kilnbio-macos" no-runtime
export KILNBIO_PYTHON="$HOME/.venvs/kilnbio/bin/python"
export PATH="$HOME/.local/share/kilnbio-macos/bin:$PATH"

kilnbio --help
kilnbio decision-methods catalog
```

`no-runtime` copies the toolkit and creates its launcher; it does not install Linux Java binaries or provision R. An existing installed release with different source files is not overwritten. For a later source update, use a new prefix such as `kilnbio-macos-v2` and adjust PATH accordingly.

In a new Terminal session, repeat the `source` and two `export` lines. You can add those lines to your own `~/.zshrc` if you want this environment selected automatically. The installer does not edit shell startup files.

## Verify a small complete example

Run from the source repository so the example paths resolve:

```bash
CHECK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/kilnbio-check.XXXXXX")
kilnbio assay-audit reach examples/assay_roster.tsv --out "$CHECK_DIR/reach.json"
kilnbio decision-methods run examples/decision_rank.json --out "$CHECK_DIR/rank.json"
kilnbio study examples/decision_study.json "$CHECK_DIR/study"
kilnbio study examples/decision_study.json "$CHECK_DIR/study" --resume

python - "$CHECK_DIR" <<'PY'
import json
from pathlib import Path
import sys
p = Path(sys.argv[1])
reach = json.loads((p / "reach.json").read_text())
rank = json.loads((p / "rank.json").read_text())
assert (reach["donors"], reach["base_reach"], reach["combined_reach"]) == (4, 3, 2)
assert rank["result"]["arms"]["combined"]["nominations"] == ["GENE_A", "GENE_C"]
assert (p / "study/rank/result.json").is_file()
print("Synthetic checks passed:", p)
PY

git rev-parse HEAD > "$CHECK_DIR/source-commit.txt"
python -m pip freeze > "$CHECK_DIR/python-packages.txt"
python -c 'import platform; print(platform.platform(), platform.machine())' > "$CHECK_DIR/platform.txt"
```

The inputs are fictional. These checks verify execution and file handoffs, not biological accuracy on your data. Keep the result directory and version records with your study. Every fresh calculation needs a new output path; `--resume` only applies to a completed, unchanged study plan.

For the numerical regression suite, install its additional test dependency and run:

```bash
python -m pip install -r environments/decision-tests-requirements.txt
python -m unittest discover -s tests -p test_decision_methods.py -v
```

The [methodology guide](METHODOLOGY.md) explains which function to choose, its assumptions, and how to interpret its outputs. Python users can also import functions from `downstream.decision_methods` while working from the source checkout.

## Add R only for analyses that require it

The decision-methods example above is Python-only. The bladder reproduction and R differential analyses need their own R dependencies. Install the R version required by the chosen recipe using the matching Intel or Apple silicon package from [CRAN for macOS](https://cran.r-project.org/bin/macosx/).

The frozen bladder recipe checks for R 4.6.x and installs pinned limma 3.68.5 and statmod 1.5.2 into a separate library. Its macOS path compiles source packages, so install Xcode Command Line Tools and the Fortran toolchain specified for your CRAN R build. Follow the [R for macOS toolchain instructions](https://mac.r-project.org/tools/); do not mix incompatible CRAN and Homebrew R toolchains.

After those prerequisites, from the source checkout and active Python environment:

```bash
export KILNBIO_RSCRIPT="$(command -v Rscript)"
export KILNBIO_R_LIBRARY="$HOME/.local/share/kilnbio-macos/bladder-R"
python recipes/bladder/setup.py --rscript "$KILNBIO_RSCRIPT" --library "$KILNBIO_R_LIBRARY"
kilnbio bladder-example --cache ./bladder-inputs --out ./bladder-run
```

The package setup requires a new library directory. The final command downloads the public bladder inputs and runs the recipe; it is separate from the synthetic checks. Other R analyses use `environments/setup_r.R` and the route-specific input guides. This guide does not claim that all R/Bioconductor dependencies have been validated on macOS.

## Raw sequencing is a separate setup

Do not run `setup_linux_runtime.sh` on a Mac. Raw-read workflows additionally require Java 21, Nextflow 25.10.4 (the versions selected by the Linux runtime script), Docker Desktop, references and sufficient storage. See the official [Nextflow installation guide](https://docs.seqera.io/nextflow/install) and [Docker Desktop for Mac](https://docs.docker.com/desktop/setup/install/mac-install/).

Choose Docker Desktop for your Mac's architecture. Some upstream workflow images may require x86-64 execution or emulation on Apple silicon; do not assume every assay container supports ARM64. Establish that compatibility for the chosen workflow before a full run. A Linux x86-64 server is the currently managed route for those workflows, and you can use macOS Terminal to connect to it.

## Common problems

- **`brew` is missing:** finish Homebrew's shell setup instructions and open a new Terminal session.
- **`kilnbio` is missing:** repeat the PATH export, or call `$HOME/.local/share/kilnbio-macos/bin/kilnbio` directly.
- **NumPy or SciPy cannot be imported:** activate the virtual environment, set `KILNBIO_PYTHON` as above, and rerun `python -m pip check`.
- **An architecture error mentions ARM64 or x86_64:** check `uname -m` and `python -c 'import platform; print(platform.machine())'`; rebuild the environment using one consistent architecture.
- **An output already exists:** choose a new result directory. Do not overwrite a previous analysis to bypass its checks.
- **A pinned R package fails to compile:** retain the build output and check the R version and matching compiler toolchain before retrying.

The current verification covers a fresh Linux x86-64/WSL2 core installation, the portable manual-launcher path on Linux, and the numerical method tests. A native Mac test remains outstanding; report the Mac model, macOS version, architecture and failing command when requesting help.
