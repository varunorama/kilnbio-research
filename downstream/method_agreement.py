"""Describe agreement between prespecified analyses of the same contrast and feature family."""
import argparse
import csv
import math
from pathlib import Path

CONTEXT=("study","contrast","feature_family","feature_unit","assembly","channel")
def agreement(values,alpha):
    pos=[];neg=[];tested=[];missing=[]
    for meth,(effect,q) in values.items():
        if effect is None or q is None:
            missing.append(meth);continue
        if not 0<=q<=1:raise ValueError(f"{meth}: adjusted P outside [0,1]")
        tested.append(meth)
        if q<=alpha:
            if effect>0:pos.append(meth)
            elif effect<0:neg.append(meth)
    support=max(len(pos),len(neg));planned=len(values)
    status=("opposite_calls" if pos and neg else "all_up" if len(pos)==planned else
            "all_down" if len(neg)==planned else "partial_up" if pos else
            "partial_down" if neg else "no_calls" if tested else "untestable")
    return [planned,len(tested),len(pos),len(neg),support, support/planned,
            support/len(tested) if tested else "",status,",".join(pos),",".join(neg),",".join(missing)]
def val(s):
    if s in ("NA",""):return None
    v=float(s)
    if not math.isfinite(v):raise ValueError("Expected finite value, blank or NA")
    return v
def read_table(stream):
    reader = csv.DictReader(stream, delimiter="\t")
    fields = reader.fieldnames or []
    if len(set(fields)) != len(fields):
        raise ValueError("Duplicate input header")
    rows = list(reader)
    if not rows or any(None in row or None in row.values() for row in rows):
        raise ValueError("Empty or ragged input")
    return fields, rows

def read(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return read_table(handle)

def run(matrix,manifest,alpha):
    return compare_tables(read(matrix), read(manifest), alpha)

def compare_tables(matrix, manifest, alpha):
    if not 0<alpha<1:raise ValueError("alpha must be between zero and one")
    fields,methods=manifest
    if set(fields)!={"method",*CONTEXT}:raise ValueError("Invalid method manifest columns")
    names=[r["method"] for r in methods]
    if len(names)<2 or len(set(names))!=len(names) or any(not n.strip() or n.strip()!=n or "," in n for n in names):raise ValueError("Specify at least two distinct method names")
    for k in CONTEXT:
        entries={r[k] for r in methods}
        if len(entries)!=1 or not next(iter(entries)).strip():raise ValueError("Incompatible or missing method metadata: "+k)
    fields,rows=matrix
    if not {"feature",*[n+s for n in names for s in ("_effect","_q")]}<=set(fields):raise ValueError("Incomplete method columns")
    ids=[r["feature"] for r in rows]
    if len(set(ids))!=len(ids) or any(not i.strip() or i.strip()!=i for i in ids):raise ValueError("Require one row per feature in the declared family")
    out=[]
    for r in rows:
        vs={n:(val(r[n+"_effect"]),val(r[n+"_q"])) for n in names}
        if any(q is not None and (e is None or not 0<=q<=1) for e,q in vs.values()):raise ValueError("Adjusted P without effect or outside [0,1]")
        out.append([r["feature"],*agreement(vs,alpha),alpha,methods[0]["contrast"],methods[0]["feature_family"]])
    return ["feature","planned","tested","up_calls","down_calls","same_direction_calls","support_planned","support_tested","status","up_methods","down_methods","untested_methods","alpha","contrast","feature_family"],out
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("matrix",type=Path);p.add_argument("manifest",type=Path);p.add_argument("--alpha",type=float,default=.05);p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    try:
        header,rows=run(a.matrix,a.manifest,a.alpha)
        with a.out.open("x",encoding="utf-8",newline="") as f:
            w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(header);w.writerows(rows)
    except (ValueError,OSError) as e:p.exit(2,f"{p.prog}: {e}\n")
if __name__=="__main__":main()
