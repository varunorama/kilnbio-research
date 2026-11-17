# Bladder cancer example

## Current source installation

From a downloaded or cloned current checkout on Linux x86-64/WSL2, install the managed R dependencies and run the example:

```bash
bash setup.sh --profile bladder
"$HOME/.local/share/kilnbio/bin/kilnbio" bladder-example --cache ./bladder-inputs --out ./bladder-run
```

Open `bladder-run/report.json`, `changes.tsv` and `shortlists.tsv` to see the patient counts and how nominations changed. The input download is approximately 220 MB. Use a new output directory for each run and add `--offline` to reuse the verified cache. [Installation details](INSTALLATION.md)

## Archived research.2 instructions


Start here to reproduce the bladder analysis using [the research.2 release](https://github.com/varunorama/kilnbio-research/releases/tag/v0.2.0-research.2). You need only the ZIP, its checksum, this README, public internet access, Python 3.10+ and R 4.6.x (tested: Python 3.12 and R 4.6.2). The public release downloads work without a GitHub account. The analysis starts from deposited TCGA RNA abundance, reverse-phase protein arrays (RPPA) and somatic mutation calls, not FASTQs or mass spectrometry. Allow 1 GB of free disk space. No Python packages, Docker, Nextflow or full KilnBio installation are needed for this example.

1. Start your clock. In the download directory check the named ZIP against its companion checksum: on Linux/WSL run `sha256sum -c kilnbio-research-0.2.0-research.2.zip.sha256`; on Windows compare `(Get-FileHash kilnbio-research-0.2.0-research.2.zip).Hash` with the hexadecimal value in the checksum file. Extract the ZIP and open a terminal in the extracted directory **containing this README**. Select a writable directory for the files generated below. The examples write them beside the extracted release.
2. Make sure `python` (or `python3`) and `Rscript` are on PATH. On Linux and macOS, building the pinned R packages from source also requires the C and Fortran build tools for your R installation (Linux: make, gcc and gfortran with R development headers). Windows uses prebuilt R packages. Run these two commands from the release root, unchanged on Windows, macOS or Linux:

```
python recipes/bladder/setup.py --library ../bladder-R
python recipes/bladder/reproduce.py --r-library ../bladder-R --cache ../bladder-inputs --out ../bladder-run
```

The first command fetches and verifies the two R packages (limma 3.68.5 and statmod 1.5.2) into a new library. The second downloads and verifies seven public study files (~220 MB), audits specimen pairing, fits the three ranking arms and writes results. If R is not on PATH, append `--rscript "path/to/Rscript"` to **both** commands. On a repeat, reuse the library, choose a new output directory and add `--offline` to the second command to verify and reuse the cache. If a transfer is interrupted, rerun the second command with the same cache; if an output directory was already created, give the rerun a new output name. A checksum failure is a stop condition; replace only the identified cache file from the pinned source. For an already installed matching limma/statmod library you may skip setup and pass it to `--r-library`; record that choice.

3. Open `bladder-run/report.json`. Record the number of eligible RNA patients, those also measured by RPPA, and the mutation-positive count in each population. Compare **A→B** (RNA in all eligible patients versus RNA in only the matched patients), then **B→C** (matched RNA versus matched RNA plus protein). The RNA-only and joint lists each have 15 places. Read `changes.tsv` to name which genes enter or leave, and `shortlists.tsv` for their ranks and directions. The sign is relative to the FGFR3 mutation contrast; the joint score is a signed rank, not a p-value. A replacement changes who gets followed up; it does not show that the new target is better. See the [frozen rules](../recipes/bladder/PROTOCOL.md) for missingness, antibody mapping and interpretation.
4. Stop your clock. Use [the handoff form](../recipes/bladder/RETURN_FORM.md) to record your answer, elapsed time (including setup), any help or changes you needed, and any failed step. Identify a gene whose nomination changed and explain which comparison accounts for it. Return the form, `report.json`, `changes.tsv`, `shortlists.tsv`, `receipt.json` and `R-session.txt` to the person who supplied the release through your agreed channel. Do not send the public input matrices or modify the protocol to obtain an answer. A failed run is a useful test result; return the command and error even if no report was produced.

Input sources: [TCGA bladder study](https://pubmed.ncbi.nlm.nih.gov/28988769/), [cBioPortal PanCancer Atlas bladder collection](https://www.cbioportal.org/study/summary?id=blca_tcga_pan_can_atlas_2018), [pinned Datahub files](../recipes/bladder/sources.json). Statistical method: [limma](https://doi.org/10.1093/nar/gkv007).
