"""Donor-held-out RNA and ATAC erythroid stage benchmark on GSE115684 source tables."""
import argparse
import csv
import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


FILES = {
    "RNA": ("GSE115678_RNAseq_rawGeneCounts.txt.gz", "338f6d862926dc1ae9889fc9a9782972e84965c8b1d01a572db524a6d101aa65"),
    "ATAC": ("GSE115672_ery_only.counts.tsv.gz", "22856b9d37f8a5ef715a76e9e2431fe8b3da68c76ee5edea477a1aa1e72a8390"),
    "peaks": ("GSE115672_ery_only.bed.gz", "b46954f44139c62bddf8a14c3e0d467b08549e6897bc5e0304cf5e86df1588e1"),
}


def read_input(root, kind):
    filename, expected = FILES[kind]
    path = root/filename
    with path.open("rb") as f:
        if hashlib.file_digest(f, "sha256").hexdigest() != expected:
            raise ValueError(f"Unexpected source file: {filename}")
    if kind == "peaks":
        return pd.read_csv(path, sep="\t", header=None)
    return pd.read_csv(path, sep="\t")


def counts(root, kind):
    x = read_input(root, kind)
    if kind == "RNA":
        # Keep all rows for library sizes. Pool symbols shared by multiple rows.
        x = x.groupby("genes", dropna=False, sort=False).sum()
    else:
        peaks = read_input(root, "peaks")
        if len(peaks) != len(x) or peaks.iloc[:,:3].duplicated().any():
            raise ValueError("Peak BED and counts disagree")
        x.index = peaks.iloc[:,:3].astype(str).agg(":".join,axis=1)
    labels = []
    for col in x:
        m = re.fullmatch(r"Donor(\d+)_P([1-8])_Rep\d+(_ATAC)?",col)
        if not m:
            raise ValueError(f"Unknown sample label {col}")
        labels.append((int(m[1]), int(m[2])))
    if x.isna().any().any() or (x < 0).any().any() or not all(pd.api.types.is_integer_dtype(t) for t in x.dtypes):
        raise ValueError("Expected nonnegative raw counts")
    samples = x.T
    samples.index = pd.MultiIndex.from_tuples(labels,names=["donor","stage"])
    pooled = samples.groupby(level=["donor","stage"]).sum()
    if (pooled.sum(axis=1) <= 0).any():
        raise ValueError("Empty donor-stage count library")
    return np.log1p(pooled.div(pooled.sum(axis=1),axis=0)*1e6)


def predict(train, y, held):
    # All means, scales, marginal rankings and coefficients use training donors.
    mean, sd = train.mean(axis=0), train.std(axis=0)
    eligible = sd > 1e-8
    z = (train[:,eligible] - mean[eligible])/sd[eligible]
    y0 = y - y.mean()
    order = np.argsort(-np.abs(z.T @ y0), kind="stable")[:100]
    z = z[:,order]
    held_z = ((held[:,eligible] - mean[eligible])/sd[eligible])[:,order]
    # Fixed ridge penalty 10 in the dual representation. No tuning on the held donor.
    return y.mean() + held_z @ z.T @ np.linalg.solve(z @ z.T + 10*np.eye(len(y)), y0)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--query-counts",type=Path,help="Project a feature-by-sample raw count TSV")
    mode.add_argument("--list-features",action="store_true",help="Write the exact reference feature universe")
    p.add_argument("--assay",choices=["RNA","ATAC"])
    args = p.parse_args()
    if bool(args.assay) != bool(args.query_counts or args.list_features):
        p.error("Supply --assay with --query-counts or --list-features")
    if args.assay:
        reference = counts(args.inputs,args.assay)
        if args.list_features:
            reference.columns.to_series(index=None).rename("feature").to_csv(args.out,index=False,mode="x")
            return
        with args.query_counts.open(encoding="utf-8-sig",newline="") as f:
            header = next(csv.reader(f,delimiter="\t"))
        if len(header)<2 or header[0]!="feature" or len(set(header))!=len(header) or any(not s.strip() for s in header):
            raise ValueError("Expected feature followed by unique sample names")
        x = pd.read_csv(args.query_counts,sep="\t",index_col="feature",keep_default_na=False)
        if not x.index.is_unique or set(x.index)!=set(reference.columns):
            raise ValueError("Match the entire source feature universe for comparable library sizes")
        x = x.loc[reference.columns].T
        if not all(pd.api.types.is_integer_dtype(t) for t in x.dtypes) or (x<0).any().any() or (x.sum(axis=1)<=0).any():
            raise ValueError("Expected nonnegative integer counts and a nonempty library for every sample")
        q = np.log1p(x.div(x.sum(axis=1),axis=0)*1e6)
        estimates = predict(reference.to_numpy(),reference.index.get_level_values("stage").to_numpy(),q.to_numpy())
        out = pd.DataFrame({"sample":x.index,"assay":args.assay,"erythroid_reference_stage":estimates,
                            "outside_P1_P8":(estimates<1)|(estimates>8)})
        out.to_csv(args.out,sep="\t",index=False,mode="x")
        return
    predictions = []
    for kind in ("RNA","ATAC"):
        x = counts(args.inputs, kind)
        donors, stages = (x.index.get_level_values(i).to_numpy() for i in ("donor","stage"))
        for donor in np.unique(donors):
            held = donors == donor
            pred = predict(x.to_numpy()[~held],stages[~held],x.to_numpy()[held])
            for stage, value in zip(stages[held],pred):
                predictions.append((donor,stage,kind,value,np.median(stages[~held])))
    df = pd.DataFrame(predictions,columns=["donor","stage","assay","estimate","training_median"])
    wide = df.pivot(index=["donor","stage","training_median"], columns="assay",values="estimate")
    if wide.isna().any().any():
        raise ValueError("Unmatched donor-stage predictions across assays")
    combined = wide.mean(axis=1).rename("estimate").reset_index().assign(assay="equal_mean")
    df = pd.concat([df,combined],ignore_index=True)
    rows = []
    for (assay,donor), g in df.groupby(["assay","donor"]):
        rows.append((assay,donor,len(g),spearmanr(g.stage,g.estimate).statistic,
                     abs(g.stage-g.estimate).mean(),abs(g.stage-g.training_median).mean()))
    stats = pd.DataFrame(rows,columns=["assay","donor","n_stages","spearman","stage_MAE","median_baseline_MAE"])
    args.out.mkdir(parents=True,exist_ok=False)
    df.to_csv(args.out/"predictions.tsv",sep="\t",index=False)
    stats.to_csv(args.out/"by_donor.tsv",sep="\t",index=False)
    print(stats.to_string(index=False))


if __name__ == "__main__":
    main()
