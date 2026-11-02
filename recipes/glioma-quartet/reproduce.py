"""Reproduce the two lead analyses from checksum-verified public measurements."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CODE = HERE / "analysis"
sys.path.insert(0, str(CODE / "glioma"))
sys.path.insert(0, str(CODE / "public_modeling"))
from evaluate import evaluate
from prepare import erasmus, huashan
from glioma_evidence import FILES, read_table


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def inputs(cache, fetch):
    cache.mkdir(parents=True, exist_ok=True)
    for row in csv.DictReader((HERE / "sources.tsv").open(encoding="utf8"), delimiter="\t"):
        target = cache / row["file"]
        if fetch and not target.exists():
            with tempfile.TemporaryDirectory(dir=cache) as tmp:
                staged = Path(tmp) / "download"
                with urlopen(Request(row["url"], headers={"User-Agent": "Research-reproduction/1.0"}), timeout=180) as src, staged.open("wb") as dst:
                    shutil.copyfileobj(src, dst)
                if digest(staged) != row["sha256"]:
                    raise ValueError(f"Downloaded file differs from the source receipt: {target.name}")
                staged.replace(target)
        if not target.is_file() or digest(target) != row["sha256"]:
            raise ValueError(f"Missing or mismatched public input: {target.name}. See sources.tsv or use --fetch.")


def read(path, **kwargs):
    return pd.read_csv(path, sep="\t", **kwargs)


def glioma(cache, out):
    frozen = out / "nominations_primary.tsv"
    subprocess.run([sys.executable, str(CODE / "glioma/freeze.py"), "--discovery-cache", str(cache), "--out", str(frozen)], check=True)
    if digest(frozen) != "1c9e9095e2569dfb1dcb0bce88f659a7e23c598f2e96d7b7cdf15335cb137b30":
        raise ValueError("Discovery ranking does not reproduce the frozen panel")
    primary = read(frozen)
    _, fc, merged = [read_table(cache, n, True) for n in FILES]
    top = sorted(merged.index, key=lambda g: (merged.at[g, "DPM"], g))[:31]
    joint = pd.DataFrame(dict(panel="joint_overall", rank=range(1, 32), gene=top,
        rna_sign=np.sign(fc.loc[top, "rna"]).astype(int).values,
        protein_sign=np.sign(fc.loc[top, "protein"]).astype(int).values,
        promoter_methylation_sign=np.sign(fc.loc[top, "methylation"]).astype(int).values,
        selection_p=merged.loc[top, "DPM"].values))
    expanded = pd.concat([joint, primary[primary.panel != "joint_only"]], ignore_index=True)
    expanded.to_csv(out / "nominations_exploratory.tsv", sep="\t", index=False)
    needed = set(primary.gene) | set(expanded.gene)
    columns = {"symbol":"gene","hgnc_id":"hgnc_id","entrez_id":"entrez_id","uniprot_ids":"uniprot_ids"}
    hgnc = read(cache / "hgnc_complete_set_2026-10-02.tsv", dtype=str, usecols=list(columns)).rename(columns=columns)
    mapping = hgnc[hgnc.gene.isin(needed)].sort_values("gene")
    if not mapping.gene.is_unique or set(mapping.gene) != needed:
        raise ValueError("Incomplete or ambiguous HGNC identifier mapping")
    all_rows = {k: [] for k in ("confirmation", "differences")}
    for endpoint, prepare in [("RNA", erasmus), ("protein", huashan)]:
        matrix, audit = prepare(cache, mapping)
        matrix_path = out / f"{endpoint}_matrix.tsv"
        matrix.to_csv(matrix_path, sep="\t", index_label="donor")
        audit.to_csv(out / f"{endpoint}_eligibility.tsv", sep="\t", index=False)
        x = read(matrix_path, index_col="donor")
        assert x.index.is_unique and x.columns.is_unique
        for design, panel in (("primary", primary), ("exploratory", expanded)):
            wire = panel.assign(panel=panel.panel.replace({"joint_overall": "joint_only"}))
            detail, counts, diffs = evaluate(x, wire, endpoint, np.random.default_rng(20261101), 2000)
            if design == "exploratory":
                detail["panel"] = detail.panel.replace({"joint_only": "joint_overall"})
                for row in counts:
                    if row["panel"] == "joint_only":
                        row["panel"] = "joint_overall"
            detail.to_csv(out / f"{design}_{endpoint}_genes.tsv", sep="\t", index=False)
            for key, items in [("confirmation", counts), ("differences", diffs)]:
                all_rows[key].extend(dict(design=design, **row) for row in items)
    for key, records in all_rows.items():
        pd.DataFrame(records).to_csv(out / (key + ".tsv"), sep="\t", index=False)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--cache", type=Path, required=True)
    cli.add_argument("--out", type=Path, required=True)
    cli.add_argument("--fetch", action="store_true")
    cli.add_argument("--rscript", default="Rscript")
    args = cli.parse_args()
    cache, out = args.cache.resolve(), args.out.resolve()
    if out.exists():
        cli.error("Please use a new output directory")
    inputs(cache, args.fetch)
    out.mkdir(parents=True)
    print("Verified all public sources. Rebuilding glioma panels, matrices and tests.", flush=True)
    glioma(cache, out)
    print("Fitting the three Quartet comparisons in held-out laboratories.", flush=True)
    env = dict(os.environ, LC_ALL="C", LANGUAGE="C", R_DEFAULT_PACKAGES="datasets,utils,grDevices,graphics,stats,methods")
    subprocess.run([args.rscript, "--vanilla", str(CODE / "quartet/quartet_hidden_replicates.R"), str(cache), str(out / "quartet")], env=env, check=True)
    subprocess.run([sys.executable, str(CODE / "plot_glioma.py"), str(out)], check=True)
    subprocess.run([sys.executable, str(CODE / "plot_quartet.py"), str(out / "quartet")], check=True)
    print("Completed both analyses in", out)


if __name__ == "__main__":
    main()
