"""Link a pretransplant lot measurement to recipient outcomes, preserving serial graft ancestry."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median
import tempfile


def table(path, fields):
    with path.open(encoding="utf-8-sig", newline="") as f:
        rows = csv.DictReader(f,delimiter="\t")
        if rows.fieldnames != fields.split():
            raise ValueError(f"{path.name}: expected columns {fields}")
        for r in rows:
            if None in r or any(v is None or not v.strip() or v != v.strip() for v in r.values()):
                raise ValueError(f"{path.name}: missing or malformed field")
            yield r


def lineage(lots_path,grafts_path):
    lots, lines = {}, {}
    for r in table(lots_path,"lot donor line"):
        lot,donor,line = r["lot"],r["donor"],r["line"]
        if lot in lots or lines.setdefault(line,donor) != donor:
            raise ValueError("Duplicate lot or line assigned to multiple donors")
        lots[lot] = (donor,line)
    passage, parents = {}, defaultdict(set)
    for r in table(grafts_path,"recipient passage source"):
        rec,level,src = r["recipient"],r["passage"],r["source"]
        if level not in {"primary","secondary"} or passage.setdefault(rec,level)!=level or src in parents[rec]:
            raise ValueError("Duplicate graft edge or inconsistent passage")
        parents[rec].add(src)
    ancestry = {}
    for rec,level in passage.items():
        if level == "primary":
            if not parents[rec] <= lots.keys():
                raise ValueError("Primary graft requires documented lot IDs")
            ancestry[rec] = parents[rec]
    for rec,level in passage.items():
        if level == "secondary":
            if any(passage.get(s)!="primary" for s in parents[rec]):
                raise ValueError("Secondary graft requires primary recipient IDs")
            ancestry[rec] = set().union(*(ancestry[s] for s in parents[rec]))
    if not lots or not passage:
        raise ValueError("Provide enrolled lots and recipients")
    return lots,passage,parents,ancestry


def select(path,lots,passage,pred,outcome):
    # Readouts for one physical specimen share owner, time and compartment.
    specs,seen,units = {},set(),{}
    selected = [defaultdict(lambda:defaultdict(list)),defaultdict(lambda:defaultdict(list))]
    for r in table(path,"specimen owner_type owner day tissue endpoint unit replicate value"):
        spec,kind,owner = r["specimen"],r["owner_type"],r["owner"]
        if kind not in {"lot","recipient"} or owner not in (lots if kind=="lot" else passage):
            raise ValueError("Unknown readout owner")
        day,value = int(r["day"]),float(r["value"])
        if day<0 or not math.isfinite(value):
            raise ValueError("Invalid day or measurement")
        key = (kind,owner,day,r["tissue"])
        if specs.setdefault(spec,key) != key:
            raise ValueError("Specimen mapped to conflicting owners, times or tissues")
        if units.setdefault(r["endpoint"],r["unit"]) != r["unit"]:
            raise ValueError("Mixed units for one endpoint")
        if r["unit"] in {"fraction","percent"} and not 0<=value<=(1 if r["unit"]=="fraction" else 100):
            raise ValueError("Fraction or percent outside its declared range")
        measurement = (spec,r["endpoint"],r["replicate"])
        if measurement in seen:
            raise ValueError("Duplicate specimen/endpoint/technical replicate")
        seen.add(measurement)
        selector=(r["endpoint"],day,r["tissue"])
        if kind=="lot" and selector==pred:
            selected[0][owner][spec].append(value)
        if kind=="recipient" and selector==outcome:
            selected[1][owner][spec].append(value)
    if not all(selected):
        raise ValueError("No measurements for a requested endpoint, time and tissue")
    return [{owner:median([median(rep) for rep in specimens.values()]) for owner,specimens in group.items()}
            for group in selected],units


def write_table(path,fields,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(fields.split())
        w.writerows(rows)


def run(lots_path,grafts_path,readouts_path,pred,outcome,level,out):
    lots,passage,parents,ancestry = lineage(lots_path,grafts_path)
    (predictors,outcomes),units = select(readouts_path,lots,passage,pred,outcome)
    if out.exists():
        raise ValueError("Choose a new output directory")
    eligible = [r for r in sorted(passage) if passage[r]==level]
    if not eligible or not any(r in outcomes for r in eligible):
        raise ValueError("No measured recipient at the selected passage and endpoint")
    rec_rows,lot_rows = [],[]
    enrolled,observed = defaultdict(list),defaultdict(list)
    for rec in eligible:
        sources = sorted(ancestry[rec])
        donors = sorted({lots[lot][0] for lot in sources})
        value = outcomes.get(rec)
        if len(sources)==1:
            enrolled[sources[0]].append(rec)
            if value is not None:
                observed[sources[0]].append(value)
        rec_rows.append([rec,level,json.dumps(sorted(parents[rec])),json.dumps(sources),json.dumps(donors),
                         value,len(sources)==1])
    for lot,(donor,line) in sorted(lots.items()):
        values = observed[lot]
        lot_rows.append([donor,line,lot,predictors.get(lot),units[pred[0]],len(enrolled[lot]),len(values),
                         median(values) if values else None,units[outcome[0]]])
    matched = [r for r in lot_rows if r[3] is not None and r[7] is not None]
    audit = {"predictor":pred,"outcome":outcome,"passage":level,"enrolled_recipients_in_passage":len(eligible),
             "multi_lot_recipients_excluded":sum(len(ancestry[r])!=1 for r in eligible),
             "lots_with_both":len(matched),"distinct_donors_with_both":len({r[0] for r in matched}),
             "technical_replicates":"median within specimen, then median within owner at selected time and tissue",
             "outcome_summary":"median of measured single-lot recipients; unmeasured recipients remain in enrolled counts",
             "input_sha256":{key:hashlib.sha256(path.read_bytes()).hexdigest()
                             for key,path in [("lots",lots_path),("grafts",grafts_path),("readouts",readouts_path)]}}
    out.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out.parent) as tmp:
        stage = Path(tmp)/"result"
        stage.mkdir()
        write_table(stage/"lots.tsv","donor line lot predictor predictor_unit enrolled_recipients measured_recipients median_outcome outcome_unit",lot_rows)
        write_table(stage/"recipients.tsv","recipient passage immediate_sources source_lots human_donors outcome attributable_to_one_lot",rec_rows)
        (stage/"design.json").write_text(json.dumps(audit,indent=2)+"\n",encoding="utf-8")
        stage.rename(out)
    return audit


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for field in ["lots","grafts","readouts","out"]:
        p.add_argument("--"+field,required=True,type=Path)
    for role in ["predictor","outcome"]:
        p.add_argument("--"+role,required=True,metavar=("ENDPOINT","DAY","TISSUE"),nargs=3)
    p.add_argument("--passage",required=True,choices=["primary","secondary"])
    a=p.parse_args()
    pred,outcome=[(v[0],int(v[1]),v[2]) for v in (a.predictor,a.outcome)]
    print(json.dumps(run(a.lots,a.grafts,a.readouts,pred,outcome,a.passage,a.out),indent=2))


if __name__=="__main__":
    main()
