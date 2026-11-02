"""Exploratory K562 assay comparison from checksum-verified published measurements."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

BASE = Path(__file__).resolve().parent
INPUTS = BASE / "inputs"
COORD = ["chrom", "chromStart", "chromEnd"]
KEY = COORD + ["measuredGeneEnsemblId"]
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}


def get_pairs(source):
    df = pd.read_csv(source, sep="\t")
    assert pd.api.types.is_bool_dtype(df.Regulated) and pd.api.types.is_bool_dtype(df.ValidConnection)
    eligible = df.loc[df.CellType.eq("K562") & df.ValidConnection & df.chrom.isin(AUTOSOMES)
                      & df.chrom.eq(df.chrTSS) & df.distanceToTSS.gt(0) & df["H3K27ac.RPM"].ge(0)]
    eligible = eligible.dropna(subset=["Regulated", "measuredGeneEnsemblId"])
    eligible = eligible.sort_values(KEY + ["Dataset", "name"])
    valid = eligible.groupby(KEY).Regulated.transform("nunique").eq(1)
    return eligible.loc[valid].drop_duplicates(KEY).copy()


def union_ranges(regions):
    intervals = {}
    for chrom, group in regions.groupby("chrom"):
        merged = []
        for start, end in group[["chromStart", "chromEnd"]].sort_values("chromStart").itertuples(index=False,name=None):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        intervals[chrom] = np.asarray(merged).T
    return intervals


def load_sites(file, ranges):
    pieces = []
    names = ["chrom", "start", "strand", "depth", "percent"]
    for chunk in pd.read_csv(file,sep="\t",header=None,usecols=[0,1,5,9,10],chunksize=400000):
        chunk.columns = names
        for chrom, part in chunk.groupby("chrom", sort=False):
            if chrom not in ranges:
                continue
            assert part.strand.isin(["+", "-", "."]).all()
            positions = part.start.to_numpy() - part.strand.eq("-").to_numpy(dtype=int)
            starts, ends = ranges[chrom]
            ix = np.searchsorted(starts, positions, side="right") - 1
            keep = (ix >= 0) & (positions < ends[ix]) & (part.depth.to_numpy() > 0)
            if np.any(keep):
                x = part.loc[keep,["chrom", "depth", "percent"]].copy()
                assert x.percent.between(0,100).all()
                x["pos"] = positions[keep]
                x["m"] = x.depth * (x.percent / 100)
                pieces.append(x[["chrom","pos","depth","m"]])
    agg = pd.concat(pieces,ignore_index=True).groupby(["chrom","pos"],as_index=False)[["depth","m"]].sum()
    return agg.loc[agg.depth.ge(5)]


def features(regions, input_directory):
    ranges = union_ranges(regions)
    sources = ["ENCFF867JRG","ENCFF721JMB"]
    reps = []
    for acc in sources:
        print("Counting qualifying CpGs",acc,flush=True)
        reps.append(load_sites(input_directory / (acc+".bed.gz"),ranges))
    sites = reps[0].merge(reps[1],on=["chrom","pos"],suffixes=("_1","_2"),validate="one_to_one")
    peaks = pd.read_csv(input_directory/"ENCFF855PCP.bed.gz",sep="\t",header=None,usecols=[0,1,2,6])
    peaks.columns=["chrom","start","end","signal"]
    assert peaks.signal.ge(0).all()
    rows = []
    for chrom,g in regions.groupby("chrom",sort=False):
        p=peaks[peaks.chrom==chrom].sort_values("start")
        starts,ends,sigs=(p[k].to_numpy() for k in ["start","end","signal"])
        prior_max=np.maximum.accumulate(ends)
        c=sites[sites.chrom==chrom].sort_values("pos")
        pos=c.pos.to_numpy()
        arrays=c[["depth_1","depth_2","m_1","m_2"]].to_numpy()
        cum=np.vstack([np.zeros(4),np.cumsum(arrays,axis=0)])
        for st,en in g[["chromStart","chromEnd"]].itertuples(index=False,name=None):
            a=np.searchsorted(prior_max,st,side="right");b=np.searchsorted(starts,en)
            signals=sigs[a:b][ends[a:b]>st]
            j,k=np.searchsorted(pos,[st,en])
            n=k-j
            r1,r2,m1,m2=cum[k]-cum[j]
            rows.append([chrom,st,en,signals.max() if signals.size else 0.,n,
                         (r1+r2)/(2*n) if n else np.nan,
                         (m1/r1 + m2/r2)/2 if n else np.nan])
    return pd.DataFrame(rows,columns=COORD+["atac","cpgs","depth","methylation"])


def predict(train,test,fields):
    vals=np.zeros(len(test))
    for chrom,idx in test.groupby("chrom").indices.items():
        sub=train.loc[train.chrom.ne(chrom)]
        assert sub.Regulated.nunique()==2
        model=make_pipeline(StandardScaler(),LogisticRegression(C=1,max_iter=500))
        model.fit(sub[fields],sub.Regulated)
        vals[idx]=model.predict_proba(test.iloc[idx][fields])[:,1]
    return vals


def run():
    cli=argparse.ArgumentParser()
    cli.add_argument("--benchmark",required=True,type=Path)
    cli.add_argument("--source-dir",type=Path,default=INPUTS)
    cli.add_argument("--verified-features",type=Path)
    args=cli.parse_args()
    # Inputs were pinned in the original retrieval manifest.
    orig=json.loads((args.benchmark/"sources.json").read_text())
    upstream=json.loads((BASE/"sources.json").read_text())
    for row in ([] if args.verified_features else upstream):
        with (args.source_dir/row["file"]).open("rb") as f: assert hashlib.file_digest(f,"md5").hexdigest()==row["md5"]
    frames=[]
    for name in ["enhancer-crispr.tsv.gz","enhancer-heldout.tsv.gz"]:
        src=args.benchmark/name
        with src.open("rb") as f:assert hashlib.file_digest(f,"sha256").hexdigest()==orig[name]["sha256"]
        frames.append(get_pairs(src))
    regions=pd.concat([df[COORD] for df in frames],ignore_index=True).drop_duplicates()
    if args.verified_features:
        obs=pd.read_csv(args.verified_features,sep="\t")
        assert not obs.duplicated(COORD).any()
        assert len(regions.merge(obs,on=COORD,validate="one_to_one"))==len(regions)
    else:
        obs=features(regions,args.source_dir)
    assert len(obs)==len(regions)
    obs.to_csv(BASE/"measured_elements.tsv.gz",sep="\t",index=False)
    train,test=[df.merge(obs,on=COORD,validate="many_to_one").sort_values(KEY).reset_index(drop=True) for df in frames]
    for d in (train,test):
        for a,b in [("d","distanceToTSS"),("h","H3K27ac.RPM"),("a","atac"),("c","cpgs"),("v","depth")]:
            d[a]=np.log1p(d[b])
    summary=[]; cover=[]
    for part,df in [("training",train),("external",test)]:
        for status,mask in [("eligible",np.ones(len(df),dtype=bool)),("common",df.cpgs.ge(3))]:
            sub=df[mask]
            cover.append(dict(split=part,universe=status,pairs=len(sub),positives=int(sub.Regulated.sum()),elements=len(sub[COORD].drop_duplicates()),genes=sub.measuredGeneEnsemblId.nunique()))
    pd.DataFrame(cover).to_csv(BASE/"coverage.tsv",sep="\t",index=False)
    comb=["d","a","h"]; ctrl=comb+["c","v"]
    groups={"distance":["d"],"distance_ATAC":["d","a"],"distance_ChIP":["d","h"],"chromatin":comb,"chromatin_coverage":ctrl,"chromatin_coverage_methylation":ctrl+["methylation"]}
    def collect(t,tag,model,fields,tr):
        vals=predict(tr,t,fields)
        order=np.argsort(-vals,kind="stable")
        hits=t.Regulated.to_numpy()
        summary.append(dict(universe=tag,model=model,n=len(t),positives=int(hits.sum()),AP=average_precision_score(hits,vals),confirmed_top25=int(hits[order[:25]].sum())))
        return vals
    collect(test,"eligible","chromatin",comb,train)
    tr=train[train.cpgs.ge(3)].reset_index(drop=True)
    te=test[test.cpgs.ge(3)].reset_index(drop=True)
    scores={name:collect(te,"common",name,fs,tr) for name,fs in groups.items()}
    pred=te[KEY+["measuredGeneSymbol","Regulated","Dataset","cpgs","depth","methylation","atac","H3K27ac.RPM"]].assign(**scores)
    pred.to_csv(BASE/"predictions.tsv.gz",sep="\t",index=False)
    g1,g2="chromatin_coverage","chromatin_coverage_methylation"
    ranks={g:np.argsort(-scores[g],kind="stable")[:25] for g in [g1,g2]}
    touched=np.union1d(ranks[g1],ranks[g2])
    slots=pred.iloc[touched][KEY+["measuredGeneSymbol","Regulated"]].copy()
    for g in ranks:slots[g+"_rank"]=pd.Series(np.arange(1,26),index=ranks[g]).reindex(touched).to_numpy()
    slots.to_csv(BASE/"followup_slots.tsv",sep="\t",index=False)
    y=te.Regulated.to_numpy(); genes,codes=np.unique(te.measuredGeneEnsemblId,return_inverse=True)
    rng=np.random.default_rng(731)
    values=[]
    for i in range(2000):
        weights=np.bincount(rng.integers(len(genes),size=len(genes)),minlength=len(genes))[codes]
        if not weights[y].sum():continue
        values.append(average_precision_score(y,scores[g2],sample_weight=weights)-average_precision_score(y,scores[g1],sample_weight=weights))
    report=dict(primary_comparison=[g1,g2],AP_gain=average_precision_score(y,scores[g2])-average_precision_score(y,scores[g1]),conditional_interval=np.quantile(values,[.025,.975]).tolist(),valid_draws=len(values),top25_shared=len(set(ranks[g1])&set(ranks[g2])),test_genes=len(genes),seed=731)
    pd.DataFrame(summary).to_csv(BASE/"results.tsv",sep="\t",index=False)
    (BASE/"summary.json").write_text(json.dumps(report,indent=2))
    print(pd.DataFrame(cover).to_string(index=False),pd.DataFrame(summary).to_string(index=False),report,sep="\n",flush=True)


if __name__=="__main__":
    run()
