"""Leave-one-lab-out test of a single concurrent reference on measured Quartet assays."""
import argparse
import hashlib
from pathlib import Path
import shutil
import tempfile
from urllib.request import urlopen

import numpy as np
import pandas as pd


FILES = {
    "metadata": (39433825, "dca5a516abbd81203adff88d0df5057f"),
    "RNA": (39433822, "ac5047aa080964bb4e08e08565475955"),
    "protein": (39433819, "e5737b2bc28c5bb09993617371f7e67b"),
    "metabolite": (39433852, "7ac78e34868f73d8197668f7d4974c45"),
}
CASES = ["D5", "F7", "M8"]


def input_path(cache, name, offline):
    file_id, checksum = FILES[name]
    path = cache / f"quartet_{name}.csv"
    if not path.exists():
        if offline:
            raise FileNotFoundError(path)
        with tempfile.TemporaryDirectory(dir=cache) as tmp:
            staged = Path(tmp) / name
            with urlopen(f"https://ndownloader.figshare.com/files/{file_id}", timeout=90) as source, staged.open("wb") as dest:
                shutil.copyfileobj(source, dest)
            if hashlib.md5(staged.read_bytes()).hexdigest() != checksum:
                raise ValueError(f"Invalid upstream file for {name}")
            staged.replace(path)
    if hashlib.md5(path.read_bytes()).hexdigest() != checksum:
        raise ValueError(f"Unexpected checksum for {name}")
    return path


def assess(assay, path, meta):
    rows = []
    m = meta[(meta.datatype.str.lower() == assay.lower()) & (meta.rep == 1)].copy()
    counts = m.groupby(["batch", "sample"]).size().unstack()
    if set(counts.columns) != {*CASES, "D6"} or counts.shape[0] != 15 or not (counts == 1).all().all():
        raise ValueError("Expected a crossed one-run-per-material design in each batch")
    if not m.groupby("batch").lab.nunique().eq(1).all():
        raise ValueError("Batch maps to multiple labs")
    data = pd.read_csv(path, index_col=0)
    if not data.index.is_unique or not m.library.is_unique or not set(m.library).issubset(data.columns):
        raise ValueError("Nonunique or unknown library/feature")
    samples = data[m.library].to_numpy(dtype=float).T
    if not np.isfinite(samples).all():
        raise ValueError("Expected publisher's finite, preprocessed log2 matrix")
    ref = samples[m["sample"].eq("D6").to_numpy()]
    lookup = {batch: i for i, batch in enumerate(m.loc[m["sample"] == "D6", "batch"])}
    ref_rows = [lookup[b] for b in m.batch]
    labs = m.lab.to_numpy()
    identity = m["sample"].to_numpy()
    for method, matrix in [("as_published", samples), ("one_D6_bridge", samples - ref[ref_rows])]:
        for lab in sorted(set(labs)):
            train = (labs != lab) & np.isin(identity, CASES)
            test = (labs == lab) & np.isin(identity, CASES)
            # Each training lab contributes equally to each donor template.
            templates = []
            for donor in CASES:
                per_lab = [matrix[train & (labs == other) & (identity == donor)].mean(axis=0)
                           for other in sorted(set(labs[train]))]
                templates.append(np.mean(per_lab, axis=0))
            templates = np.array(templates)
            for ix in np.flatnonzero(test):
                distances = np.sqrt(np.mean((templates - matrix[ix]) ** 2, axis=1))
                donor_ix = CASES.index(identity[ix])
                rows.append(dict(assay=assay, method=method, held_out_lab=lab,
                                 batch=m.iloc[ix].batch, material=identity[ix], features=matrix.shape[1],
                                 correct=np.argmin(distances) == donor_ix,
                                 same_material_RMSE=distances[donor_ix]))
    return rows


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--cache", required=True, type=Path)
    cli.add_argument("--out", required=True, type=Path)
    cli.add_argument("--offline", action="store_true")
    a = cli.parse_args()
    if a.out.exists():
        cli.error("Use a new output directory")
    a.cache.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(input_path(a.cache, "metadata", a.offline))
    rows = []
    for assay in ("RNA", "protein", "metabolite"):
        rows.extend(assess(assay, input_path(a.cache, assay, a.offline), meta))
    details = pd.DataFrame(rows)
    by_lab = details.groupby(["assay", "method", "held_out_lab"]).agg(
        accuracy=("correct", "mean"), mean_RMSE=("same_material_RMSE", "mean")).reset_index()
    summary = by_lab.groupby(["assay", "method"]).agg(labs=("held_out_lab", "size"),
        mean_lab_accuracy=("accuracy", "mean"), mean_lab_RMSE=("mean_RMSE", "mean"))
    a.out.mkdir(parents=True)
    details.to_csv(a.out / "held_out_libraries.tsv", sep="\t", index=False)
    by_lab.to_csv(a.out / "held_out_labs.tsv", sep="\t", index=False)
    summary.to_csv(a.out / "summary.tsv", sep="\t")
    print(summary.to_string())


if __name__ == "__main__":
    main()
