"""Audit published IDH glioma association summaries and compare two different claims."""
import argparse
import hashlib
from pathlib import Path
import shutil
import tempfile
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control


BASE = ("https://raw.githubusercontent.com/reimandlab/DPM_publication_code/"
        "fd8f2bfdc9e090654eb44cc942d89ad2136c6ea3/Figure_5/input_data/intermediate_files/")
FILES = {
    "figure5_pvals.tsv": "24c0e55a6154bc6dbaedbae7907d573efd729209530b0ba0b508be6ca4d984a3",
    "figure5_fcs.tsv": "5529918e97c4e85ec8e4cc425f763b6e6752a066b1296857fdcacef5301fcd7e",
    "figure5_merged_pvalues.tsv": "93f8cbcfe6b37c6aac4dccf5b616b68aa5d9022418516fd96b444ea93dddbf5a",
}
LAYERS = ["rna", "protein", "methylation"]


def read_table(cache, filename, offline):
    path = cache / filename
    if not path.exists():
        if offline:
            raise FileNotFoundError(path)
        with tempfile.TemporaryDirectory(dir=cache) as tmp:
            stage = Path(tmp) / filename
            with urlopen(BASE + filename, timeout=60) as resp, stage.open("wb") as out:
                shutil.copyfileobj(resp, out)
            if hashlib.sha256(stage.read_bytes()).hexdigest() != FILES[filename]:
                raise ValueError(f"Unexpected source for {filename}")
            stage.replace(path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != FILES[filename]:
        raise ValueError(f"Checksum failed for {filename}")
    # Source header excludes the gene-ID column, as in R write.table().
    return pd.read_csv(path, sep="\t", index_col=0)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--cache", type=Path, required=True)
    cli.add_argument("--out", type=Path, required=True)
    cli.add_argument("--offline", action="store_true")
    args = cli.parse_args()
    if args.out.exists():
        cli.error("Choose a new output directory")
    args.cache.mkdir(exist_ok=True, parents=True)
    p, effects, merged = [read_table(args.cache, file, args.offline) for file in FILES]
    if list(p) != LAYERS or list(effects) != LAYERS or list(merged) != ["Brown", "DPM"]:
        raise ValueError("Unexpected assay columns")
    if not p.index.is_unique or p.index.hasnans or not p.index.equals(effects.index) or not p.index.equals(merged.index):
        raise ValueError("Unexpected gene identifiers or alignment")
    # Ten published RNA log-fold changes are -Inf. Their signs are defined,
    # but the magnitudes are unbounded and are flagged in the output.
    if effects.isna().any().any() or not all(t.ge(0).all().all() and t.le(1).all().all() for t in (p, merged)):
        raise ValueError("Invalid published summary statistics")
    if p.isna().any().any() or merged.isna().any().any():
        raise ValueError("Unexpected missing probability")
    # The partial-conjunction null allows one true effect. This two-sided
    # test uses no assumed direction or independence among the three assays.
    pc2 = np.minimum(1, 2 * np.sort(p.to_numpy(), axis=1)[:, 1])
    result = pd.concat([p.add_prefix("p_"), effects.add_prefix("source_log2fc_"), merged.add_prefix("p_")], axis=1)
    result["unbounded_source_effect"] = ~np.isfinite(effects).all(axis=1)
    result["p_at_least_two"] = pc2
    for method in ("bh", "by"):
        for layer in (*LAYERS, "DPM", "at_least_two"):
            result[f"{method}_{layer}"] = false_discovery_control(result[f"p_{layer}"], method=method)
    single = result[["bh_" + layer for layer in LAYERS]].min(axis=1) <= 0.05
    dpm = result.bh_DPM <= 0.05
    result["joint_BH_only"] = dpm & ~single
    result["all_three_raw_p_above_005"] = p.min(axis=1) > 0.05
    result["full_direction_agreement"] = (np.sign(effects.rna) == np.sign(effects.protein)) & (
        np.sign(effects.rna) == -np.sign(effects.methylation)) & (effects.rna != 0)
    counts = pd.Series({
        "genes_in_shared_universe": len(result), "any_individual_BH": int(single.sum()),
        "published_DPM_after_BH": int(dpm.sum()), "DPM_BH_no_individual_BH": int(result.joint_BH_only.sum()),
        "DPM_BH_all_three_raw_p_above_005": int((dpm & result.all_three_raw_p_above_005).sum()),
        "at_least_two_BH": int((result.bh_at_least_two <= 0.05).sum()),
        "at_least_two_BY": int((result.by_at_least_two <= 0.05).sum()),
        "DPM_BH_only_also_at_least_two_BY": int((result.joint_BH_only & (result.by_at_least_two <= 0.05)).sum()),
    }, name="genes")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.out.parent) as tmp:
        stage = Path(tmp) / "glioma"
        stage.mkdir()
        counts.to_csv(stage / "counts.tsv", sep="\t", index_label="metric")
        ordered = result.sort_values("p_DPM", kind="stable")
        ordered.to_csv(stage / "genes.tsv", sep="\t", index_label="gene")
        ordered.loc[ordered.joint_BH_only].to_csv(stage / "joint_only.tsv", sep="\t", index_label="gene")
        stage.rename(args.out)
    print(counts.to_string())


if __name__ == "__main__":
    main()
