"""Audit assay availability, specimen pairing and two prespecified prediction policies."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import random
from statistics import fmean

def read_table(path, required):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not set(required) <= set(fields):
            raise ValueError(f"{path}: missing or duplicate columns; require {required}")
        rows = list(reader)
    if not rows or any(None in r or None in r.values() for r in rows):
        raise ValueError(f"{path}: empty or malformed table")
    return rows

def identifier(value):
    if not value or value != value.strip() or any(c in value for c in "\r\n\t"):
        raise ValueError("Identifiers and labels must be nonempty, with no surrounding whitespace")
    return value

def unique(rows, key):
    ids = [identifier(r[key]) for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError(f"Repeated {key}; supply one independent row per {key}")
    return dict(zip(ids, rows))

def flag(s):
    if s not in ("0", "1"):
        raise ValueError("Assay flags must be 0 or 1; unknown status needs separate review")
    return int(s)

def reach(rows):
    unique(rows, "donor")
    counts = {"both":0, "base_only":0, "add_only":0, "neither":0}
    for r in rows:
        counts[{(1,1):"both",(1,0):"base_only",(0,1):"add_only",(0,0):"neither"}[
            flag(r["base_ok"]),flag(r["add_ok"])]] += 1
    return dict(donors=len(rows), patterns=counts, combined_reach=counts["both"],
                base_reach=counts["both"]+counts["base_only"])

def number(s):
    if s=="": return None
    n=float(s)
    if not math.isfinite(n): raise ValueError("Use an empty field for missing values")
    return n

def paired(errors):
    if not errors: return dict(n=0,base_mae=None,policy_mae=None,gain=None,interval=None)
    x,y=zip(*errors)
    ds=[a-b for a,b in errors]
    # Fixed predictions; resample donors jointly, never independently by arm.
    rng=random.Random(17042)
    interval=None
    if len(ds)>1:
        boot=sorted(fmean(rng.choices(ds,k=len(ds))) for _ in range(2000))
        interval=[boot[49],boot[1949]]
    return dict(n=len(ds),base_mae=fmean(x),policy_mae=fmean(y),gain=fmean(ds),interval=interval)

def policies(roster, scores):
    out=reach(roster)
    scores=unique(scores,"donor")
    if set(scores)!={r["donor"] for r in roster}:
        raise ValueError("Prediction roster differs; retain even unscored donors as blank rows")
    all_errors=[]; both=[]; fallback=[]; missing_y=0; eligible_unscored=0
    for r in sorted(roster,key=lambda r:r["donor"]):
        pred=scores[r["donor"]]
        y,b,a=(number(pred[k]) for k in ("observed","baseline","augmented"))
        b_ok=flag(r["base_ok"]); both_ok=b_ok and flag(r["add_ok"])
        if (b is not None)!=bool(b_ok) or (a is not None)!=bool(both_ok):
            raise ValueError("Predictions must match declared scoring eligibility; review missing or unexpected predictions")
        if y is None: missing_y+=1; continue
        if not b_ok: eligible_unscored+=1; continue
        e=(abs(y-b),abs(y-(a if both_ok else b)))
        all_errors.append(e); (both if both_ok else fallback).append(e)
    out.update(missing_outcomes=missing_y, labeled_abstentions=eligible_unscored,
               both_scored=paired(both), fallback_only=paired(fallback), all_scored=paired(all_errors))
    return out

def pairs(roster, observations, base, add, tissue):
    donors=unique(roster,"donor")
    if identifier(base)==identifier(add): raise ValueError("Specify distinct assays")
    identifier(tissue)
    unique(observations,"measurement")
    data={d:[] for d in donors}; owners={}; unknown=0
    for r in observations:
        d=r["donor"]; s=r["specimen"]; t=identifier(r["tissue"]); a=identifier(r["assay"])
        if d not in donors or a not in (base,add): raise ValueError("Unexpected donor or assay")
        if s:
            if owners.setdefault(identifier(s),(d,t))!=(d,t): raise ValueError("Specimen maps to conflicting donors or tissues")
        else: unknown+=1
        data[d].append((s,t,a,r["measurement"]))
    stats=dict(donors=len(donors), measurements=len(observations), unknown_specimen_records=unknown,
               both_any_tissue=0,both_requested_tissue=0,both_same_specimen=0,
               unique_one_to_one_donors=0,ambiguous_donors=0)
    accepted=[]; review=[]
    for d,items in sorted(data.items()):
        wanted=[r for r in items if r[1]==tissue]
        stats["both_any_tissue"]+=len({r[2] for r in items})==2
        stats["both_requested_tissue"]+=len({r[2] for r in wanted})==2
        group={}
        for s,t,a,i in wanted:
            if s: group.setdefault(s,{base:[],add:[]})[a].append(i)
        matches={s:g for s,g in sorted(group.items()) if g[base] and g[add]}
        if not matches: continue
        stats["both_same_specimen"]+=1
        if len(matches)==1 and all(len(x)==1 for x in next(iter(matches.values())).values()):
            s,g=next(iter(matches.items()))
            accepted.append(dict(donor=d,specimen=s,base_measurement=g[base][0],additional_measurement=g[add][0]))
        else:
            review.append(dict(donor=d,candidates={s:{a:sorted(ids) for a,ids in g.items()} for s,g in matches.items()}))
    stats.update(unique_one_to_one_donors=len(accepted),ambiguous_donors=len(review))
    return dict(summary=stats,tissue=tissue,base=base,additional=add,unique_pairs=accepted,review_required=review)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest="action",required=True)
    for name in ("reach","policy","pair"):
        s=sub.add_parser(name)
        s.add_argument("roster",type=Path)
        if name=="policy": s.add_argument("scores",type=Path)
        if name=="pair":
            s.add_argument("observations",type=Path)
            for arg in ("base","add","tissue"): s.add_argument("--"+arg, required=True)
        s.add_argument("--out",type=Path,required=True,help="New JSON file (will not overwrite)")
    a=p.parse_args()
    try:
        cols=["donor"] if a.action=="pair" else ["donor","base_ok","add_ok"]
        roster=read_table(a.roster,cols); inputs=[a.roster]
        if a.action=="pair":
            inputs.append(a.observations)
            result=pairs(roster,read_table(a.observations,["donor","tissue","assay","specimen","measurement"]),a.base,a.add,a.tissue)
        elif a.action=="reach": result=reach(roster)
        else:
            inputs.append(a.scores)
            result=policies(roster,read_table(a.scores,["donor","observed","baseline","augmented"]))
        result["inputs_sha256"]=[hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs]
        payload=json.dumps(result,indent=2,allow_nan=False)+"\n"
        with a.out.open("x",encoding="utf-8") as f: f.write(payload)
    except (OSError,ValueError) as e: p.exit(2,f"{p.prog}: {e}\n")

if __name__=="__main__":main()
