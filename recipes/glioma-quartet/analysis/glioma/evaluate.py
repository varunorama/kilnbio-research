"""Test the precommitted panel hypotheses on patient-level external measurements."""
import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr
from scipy.stats import false_discovery_control, mannwhitneyu


ROOT = Path(__file__).resolve().parent
PANEL_SHA = "1c9e9095e2569dfb1dcb0bce88f659a7e23c598f2e96d7b7cdf15335cb137b30"


def rank_test(mut, wt, wm, ww, sign):
    # Row 0 is the actual cohort. Other rows are bootstrap donor multiplicities.
    levels = np.unique(np.r_[mut[np.isfinite(mut)], wt[np.isfinite(wt)]])
    a = wm @ (mut[:, None] == levels)
    b = ww @ (wt[:, None] == levels)
    n, m = a.sum(axis=1), b.sum(axis=1)
    ties = a + b
    total = n + m
    pairs = n * m
    u = (a * (b.cumsum(axis=1) - b/2)).sum(axis=1)
    if sign < 0:
        u = pairs - u
    p = np.ones(len(wm))
    auc = np.full(len(wm), np.nan)
    eligible = (n >= max(5, len(mut)/2)) & (m >= max(5, len(wt)/2)) & (sign != 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        sd = np.sqrt(pairs/12 * (total + 1 - (ties**3 - ties).sum(axis=1)/(total*(total-1))))
    valid = eligible & (sd > 0)
    p[valid] = ndtr(-(u[valid] - pairs[valid]/2 - .5)/sd[valid])
    auc[eligible] = u[eligible]/pairs[eligible]
    if valid[0]:
        check = mannwhitneyu(mut[np.isfinite(mut)] * sign, wt[np.isfinite(wt)] * sign,
                             alternative="greater", method="asymptotic", use_continuity=True)
        np.testing.assert_allclose(p[0], check.pvalue, rtol=1e-10, atol=1e-14)
    return p, auc, eligible


def evaluate(x, panel, endpoint, rng, reps):
    groups = [x[x.mutant == status] for status in (1,0)]
    weights = []
    for g in groups:
        n = len(g)
        weights.append(np.vstack([np.ones(n),rng.multinomial(n, np.ones(n)/n, size=reps)]))
    col = {"joint_only":endpoint.lower()+"_sign", "RNA_only":"rna_sign", "protein_only":"protein_sign"}
    panel = panel.copy()
    panel["sign"] = [r[col[r["panel"]]] for _,r in panel.iterrows()]
    keys = list(dict.fromkeys(zip(panel.gene,panel.sign)))
    results = []
    for gene, sign in keys:
        v = [g[gene].to_numpy(dtype=float) for g in groups]
        results.append(rank_test(*v, *weights, sign))
    p, auc, eligible = (np.column_stack(xx) for xx in zip(*results))
    adj = false_discovery_control(p, axis=1)
    hits = (adj <= .05) & (auc > .5)
    index = {k:i for i,k in enumerate(keys)}
    ix = [index[k] for k in zip(panel.gene,panel.sign)]
    panel["observed_mut"] = [int(groups[0][g].count()) for g in panel.gene]
    panel["observed_wt"] = [int(groups[1][g].count()) for g in panel.gene]
    for name, vals in [("testable",eligible),("directional_AUC",auc),("directional_p",p),("BH",adj),("confirmed",hits)]:
        panel[name] = vals[0,ix]
    n_hits = {}
    report = []
    for name, rows in panel.groupby("panel",sort=False):
        positions = [index[k] for k in zip(rows.gene,rows.sign)]
        n_hits[name] = hits[:,positions].sum(axis=1)
        report.append(dict(endpoint=endpoint, panel=name, n_mut=len(groups[0]), n_wt=len(groups[1]),
                           family_size=len(keys), nominated=len(rows), testable=int(rows.testable.sum()),
                           direction=int((rows.directional_AUC > .5).sum()),
                           nominal=int(((rows.directional_p<=.05)&(rows.directional_AUC>.5)).sum()),
                           confirmed=int(n_hits[name][0])))
    contrasts = []
    for name in ("RNA_only","protein_only"):
        diff = n_hits["joint_only"] - n_hits[name]
        lo,hi = np.quantile(diff[1:],[.025,.975],method="inverted_cdf")
        contrasts.append(dict(endpoint=endpoint, comparator=name, difference=int(diff[0]),
                              stability_low=int(lo), stability_high=int(hi), resamples=reps))
    return panel, report, contrasts


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--prepared", type=Path, required=True)
    cli.add_argument("--out", type=Path, required=True)
    a=cli.parse_args()
    fp=ROOT / "frozen_panels.tsv"
    if hashlib.sha256(fp.read_bytes().replace(b"\r\n", b"\n")).hexdigest() != PANEL_SHA:
        raise ValueError("Panel checksum mismatch")
    panel=pd.read_csv(fp,sep="\t")
    a.out.mkdir(exist_ok=False,parents=True)
    scores, contrasts = [], []
    for endpoint in ("RNA","protein"):
        x=pd.read_csv(a.prepared/(endpoint+".tsv"),sep="\t",index_col="donor")
        if not x.index.is_unique or not x.columns.is_unique or set(x.mutant.unique()) != {0,1}:
            raise ValueError("Unexpected patient matrix")
        detail, rows, differences = evaluate(x,panel,endpoint,np.random.default_rng(20261101),2000)
        detail.to_csv(a.out/(endpoint+"_genes.tsv"),index=False,sep="\t")
        scores.extend(rows)
        contrasts.extend(differences)
    for filename, records in [("confirmation",scores),("differences",contrasts)]:
        frame=pd.DataFrame(records)
        frame.to_csv(a.out/(filename+".tsv"),sep="\t",index=False)
        print(frame.to_string(index=False))


if __name__=="__main__":
    main()
