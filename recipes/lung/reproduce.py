"""Compare restricting enrollment vs adding protein, then evaluate frozen targets in GEO."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
RELEASE = HERE.parents[1]


def verify(directory, manifest):
    for name, receipt in manifest.items():
        file = directory / name
        with file.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        if digest != receipt["sha256"] or file.stat().st_size != receipt["bytes"]:
            raise ValueError(f"Public source does not match its receipt: {name}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True, type=Path, help="New working and results directory")
    p.add_argument("--cache", type=Path, help="Optional read-only staged discovery/ and external/ inputs; no network")
    p.add_argument("--rscript", default="Rscript")
    p.add_argument("--r-library", type=Path, required=True, help="Installed R library containing limma")
    a = p.parse_args()
    out = a.out.resolve()
    if out.exists():
        p.error("Choose a new output directory")
    if not a.r_library.is_dir():
        p.error("--r-library must be an existing R library directory")
    out.mkdir(parents=True)
    for name in ("fetch_discovery.py","prepare.py","freeze.R","fetch_external.py","external_matrix.py","validate.R","report.py"):
        shutil.copyfile(HERE/name, out/name)
    # Scripts use paths relative to their working copy. Keep the install and any cache unchanged.
    def py(name, *args):
        subprocess.run([sys.executable, str(out/name), *map(str,args)], check=True)
    def r(name):
        env = dict(os.environ, LANGUAGE="C")
        subprocess.run([a.rscript,"--vanilla",str(out/name),str(out),str(a.r_library.resolve())], env=env, check=True)
    def stage(folder, manifest, fetch):
        if a.cache is None:
            py(fetch)
        else:
            for name in manifest:
                src = a.cache.resolve()/folder/name
                target = out/folder/name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(src, target)
        verify(out/folder, manifest)
    discovery = json.loads((HERE/"provenance/discovery_sources.json").read_text())
    external = json.loads((HERE/"provenance/external_sources.json").read_text())["files"]
    stage("discovery", discovery, "fetch_discovery.py")
    py("prepare.py", RELEASE)
    r("freeze.R")
    with (out/"analysis/frozen.tsv").open("rb") as handle:
        frozen = hashlib.file_digest(handle, "sha256").hexdigest()
    (out/"analysis/frozen.sha256").write_text(frozen+"\n", encoding="ascii")
    stage("external", external, "fetch_external.py")
    py("external_matrix.py")
    r("validate.R")
    py("report.py", RELEASE)
    print("Results:", out/"analysis/comparison.json")


if __name__ == "__main__":
    main()
