"""Score locally installed Biolearn DNAm models on complete profiles, with matrix-level QC."""
import argparse
import csv
import gzip
import hashlib
import json
from importlib.metadata import version
from importlib.util import find_spec
from pathlib import Path
import tempfile

import torch  # Load Torch before the other Windows native libraries.
import numpy as np
import pandas as pd
from biolearn.data_library import GeoData
from biolearn.model import model_definitions
from biolearn.model_gallery import ModelGallery
from biolearn.util import get_data_file

from clocks import local_linear_models


def catalog():
    names = sorted(local_linear_models() | {"GrimAgeV1", "GrimAgeV2", "AltumAge", "EpiTOC2", "MiAge",
                                            "GPAge10", "GPAge30", "GPAge71", "GPAgeA", "GPAgeB", "GPAgeC"})
    models, required = {}, {}
    gallery = ModelGallery()
    for name in names:
        model = gallery.get(name, imputation_method="none")
        models[name] = model
        required[name] = set(model.methylation_sites())
        if name == "DunedinPACE":
            required[name].update(pd.read_csv(get_data_file("DunedinPACE_Gold_Means.csv"), index_col=0).index)
    return models, required


def read_matrix(stream, wanted, *, implicit_probe=False):
    """Scan every measurement. Keep only requested rows. The first field of each row is its probe ID."""
    header = next(csv.reader([stream.readline().decode("utf-8-sig")], delimiter="\t"))
    if not header:
        raise ValueError("Empty matrix")
    explicit = {"probe",""} | ({"TargetID","ID_REF","ID","ProbeID"} if implicit_probe else set())
    if implicit_probe and header[0] not in explicit:
        samples = header
    elif header[0] in explicit:
        samples = header[1:]
    else:
        raise ValueError("Expected first field 'probe' in the TSV header")
    if (not samples or any(not s or s != s.strip() for s in samples)
            or len(set(samples)) != len(samples) or "probe" in samples):
        raise ValueError("Invalid or duplicate sample identifiers")
    fields = ["probe", *samples]
    counts = np.zeros(len(samples), dtype=np.int64)
    bad = counts.copy()
    totals = np.zeros(len(samples))
    mins = np.full(len(samples),np.inf)
    maxes = -mins
    seen, kept = set(), []
    # Explicit names and float parsing bound memory even for large cohorts.
    chunks = pd.read_csv(stream, sep="\t", header=None, names=fields, index_col=0,
                         chunksize=max(256,min(4096,400000//len(samples))),
                         dtype={k:float for k in samples}, keep_default_na=True)
    for chunk in chunks:
        ids = chunk.index
        if ids.hasnans or not ids.is_unique or seen.intersection(ids):
            raise ValueError("Missing or duplicate probe identifier")
        seen.update(ids)
        x = chunk.to_numpy(dtype=float,copy=False)
        finite = np.isfinite(x)
        counts += finite.sum(axis=0)
        bad += (np.isinf(x) | (finite & ((x<0)|(x>1)))).sum(axis=0)
        totals += np.where(finite,x,0).sum(axis=0)
        mins = np.minimum(mins,np.min(np.where(finite,x,np.inf),axis=0))
        maxes = np.maximum(maxes,np.max(np.where(finite,x,-np.inf),axis=0))
        selected = chunk.loc[ids.isin(wanted)]
        if not selected.empty:
            kept.append(selected.copy())
    if not seen:
        raise ValueError("Empty matrix")
    qc = pd.DataFrame(dict(probes=len(seen), finite=counts, invalid_beta=bad,
                           missing_fraction=1-counts/len(seen), minimum=np.where(counts,mins,np.nan),
                           maximum=np.where(counts,maxes,np.nan), mean_beta=totals/np.maximum(1,counts),
                           eligible=(bad==0)&(counts>0)&(mins<maxes)), index=samples)
    qc.loc[counts==0,"mean_beta"] = np.nan
    return (pd.concat(kept) if kept else pd.DataFrame(columns=samples,dtype=float)), qc


def predict_complete(name, model, beta, meta):
    """Preserve the sample axis for Biolearn 0.9.1 models that squeeze singleton predictions."""
    ids = beta.columns
    if len(ids) != 1 or not (name == "AltumAge" or name.startswith("GPAge")):
        return model.predict(GeoData(meta, dnam=beta))
    # Each of these fixed models predicts each column independently. A second copy
    # avoids the upstream scalar result. Only the original sample is returned.
    twin = str(ids[0]) + "_shape_check"
    paired_beta = beta.assign(**{twin: beta.iloc[:, 0]})
    paired_meta = pd.concat([meta, meta.rename(index={ids[0]: twin})])
    result = model.predict(GeoData(paired_meta, dnam=paired_beta))
    if list(result.index) != [ids[0], twin] or not np.allclose(result.iloc[0], result.iloc[1]):
        raise ValueError(f"Inconsistent singleton predictions from {name}")
    return result.loc[ids]


def score_panel(beta, qc, metadata, models, required):
    if not beta.columns.equals(qc.index) or not qc.index.is_unique or not beta.index.is_unique:
        raise ValueError("Beta columns and QC sample identifiers must agree in the same order")
    if not metadata.index.is_unique or set(metadata.index)!=set(qc.index):
        raise ValueError("Metadata must contain one record for every input sample")
    meta = metadata.reindex(qc.index).copy()
    meta["age"] = pd.to_numeric(meta.get("age", pd.Series(index=meta.index,dtype=float)),errors="coerce")
    meta["sex"] = meta.get("sex", pd.Series(index=meta.index,dtype=str)).map({"Female":0,"Male":1})
    scores = pd.DataFrame(index=qc.index)
    coverage = []
    for name, model in models.items():
        engine_ready = not name.startswith("GPAge") or find_spec("GPy") is not None
        needed = required[name]
        present = beta.index.intersection(needed)
        complete = np.zeros(len(qc),dtype=bool)
        eligible = complete.copy()
        if len(needed)==len(present):
            values = beta.loc[present]
            complete = np.isfinite(values.to_numpy()).all(axis=0)
            eligible = complete & qc.eligible.to_numpy()
            if name in ("GrimAgeV1","GrimAgeV2"):
                eligible &= (meta.age.ge(0) & np.isfinite(meta.age) & meta.sex.notna()).to_numpy()
            if not engine_ready:
                eligible[:] = False
            if eligible.any():
                ids = qc.index[eligible]
                result = predict_complete(name, model, values.loc[:,ids].copy(), meta.loc[ids].copy())
                if name in ("GrimAgeV1","GrimAgeV2"):
                    # AgeAccelGrim is a cohort-fitted residual. Report only the frozen score.
                    result = result[["DNAmGrimAge"]]
                if not result.index.equals(ids) or not np.isfinite(result.to_numpy()).all():
                    raise ValueError(f"Invalid or misaligned result from {name}")
                for col in result:
                    key = name if result.shape[1]==1 else f"{name}.{col}"
                    scores[key] = result[col].reindex(qc.index)
        spec = model_definitions[name]
        coverage.append(dict(model=name, engine="ready" if engine_ready else "install optional GPy", required=len(needed), present=len(present), complete_profiles=int(complete.sum()),
                             scored_profiles=int(eligible.sum()), intended_tissue=spec.get("tissue"), endpoint=spec.get("output"),
                             source=spec.get("source")))
    return scores, pd.DataFrame(coverage)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("matrix",type=Path)
    p.add_argument("--metadata",type=Path,help="TSV with sample, optional age in years and sex (Female/Male)")
    p.add_argument("--out",required=True,type=Path)
    a = p.parse_args()
    if a.out.exists():
        p.error("Choose a new output directory")
    models, requirements = catalog()
    opener = gzip.open if a.matrix.suffix==".gz" else open
    with opener(a.matrix,"rb") as stream:
        beta, qc = read_matrix(stream,set.union(*requirements.values()))
    meta = (pd.read_csv(a.metadata,sep="\t",dtype={"sample":str}).set_index("sample",verify_integrity=True)
            if a.metadata else pd.DataFrame(index=qc.index))
    scores,coverage = score_panel(beta,qc,meta,models,requirements)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=a.out.parent) as temp:
        staged = Path(temp)/"panel"
        staged.mkdir()
        qc.to_csv(staged/"sample_qc.tsv",sep="\t",index_label="sample")
        scores.to_csv(staged/"scores.tsv",sep="\t",index_label="sample")
        coverage.to_csv(staged/"coverage.tsv",sep="\t",index=False)
        versions = {k:version(k) for k in ("biolearn","pandas","numpy","torch")}
        if find_spec("GPy") is not None:
            versions["GPy"] = version("GPy")
        (staged/"versions.json").write_text(json.dumps(versions, indent=2))
        provenance = {}
        for label, path in (("matrix", a.matrix), ("metadata", a.metadata)):
            if path:
                digest = hashlib.sha256()
                with path.open("rb") as stream:
                    for part in iter(lambda: stream.read(1 << 20), b""):
                        digest.update(part)
                provenance[label] = dict(file=path.name, sha256=digest.hexdigest())
        (staged/"inputs.json").write_text(json.dumps(provenance, indent=2))
        staged.rename(a.out)


if __name__ == "__main__":
    main()
