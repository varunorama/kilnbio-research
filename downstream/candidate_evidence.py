"""Report replication by assay and cohort; keep cited perturbations as separate evidence."""
import argparse
import hashlib
import json
from pathlib import Path
from assay_audit import read_table, identifier, flag, number

def keyed(rows, cols):
    d={}
    for r in rows:
        key=tuple(identifier(r[c]) for c in cols)
        if key in d: raise ValueError(f"Duplicate key {key}")
        d[key]=r
    return d

def summarize(nominations, endpoints, observations, references=(), alpha=0.05):
    if not 0<alpha<1: raise ValueError("alpha must be between zero and one")
    nom=keyed(nominations,("gene","assay"))
    design=keyed(endpoints,("assay","cohort"))
    obs=keyed(observations,("gene","assay","cohort"))
    genes=sorted({g for g,a in nom})
    if {a for g,a in nom}!={a for a,c in design}:raise ValueError("Nominated and validation assays differ")
    wanted={(g,a,c) for g,a in nom for a2,c in design if a==a2}
    if set(obs)!=wanted:raise ValueError("Provide every nominated endpoint, including untestable rows")
    reports={g:dict(gene=g,endpoints=[],perturbation_reports=[]) for g in genes}
    for g,a,c in sorted(wanted):
        direction=nom[g,a]["direction"]
        if direction not in ("-1","1"):raise ValueError("Frozen direction must be -1 or 1")
        r=obs[g,a,c]; test=flag(r["testable"])
        e,q=number(r["effect"]),number(r["q"])
        if test and (e is None or q is None or not 0<=q<=1):raise ValueError("Testable rows need finite effects and adjusted P values")
        if not test and (q is not None or e is not None):raise ValueError("Leave effect and q blank for untestable rows")
        status=("unmeasured" if not test else "confirmed" if e*int(direction)>0 and q<=alpha else "not_confirmed")
        reports[g]["endpoints"].append(dict(assay=a,cohort=c,expected_direction=int(direction),effect=e,adjusted_p=q,status=status))
    for k,r in keyed(references,("gene","source","model","intervention","readout")).items():
        if k[0] not in reports or not k[1].startswith(("https://","doi:","PMID:")):raise ValueError("Use nominated genes and traceable citation identifiers")
        identifier(r["finding"])
        reports[k[0]]["perturbation_reports"].append(dict(r))
    for r in reports.values():
        hits=[e for e in r["endpoints"] if e["status"]=="confirmed"]
        r["confirmed_assays"]=sorted({e["assay"] for e in hits})
        r["confirmed_cohort_labels"]=sorted({e["cohort"] for e in hits})
    return dict(nominated_genes=len(genes),alpha=alpha,candidates=list(reports.values()))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ("nominations","endpoints","observations"):p.add_argument(n,type=Path)
    p.add_argument("--perturbations",type=Path)
    p.add_argument("--alpha",type=float,default=.05)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    try:
        cols=[("nominations","gene assay direction"),("endpoints","assay cohort"),("observations","gene assay cohort effect q testable"),("perturbations","gene source model intervention readout finding")]
        tables={key:read_table(getattr(a,key),fields.split()) if getattr(a,key) else [] for key,fields in cols}
        result=summarize(tables["nominations"],tables["endpoints"],tables["observations"],tables["perturbations"],a.alpha)
        result["inputs_sha256"]={key:hashlib.sha256(getattr(a,key).read_bytes()).hexdigest() for key,_ in cols if getattr(a,key)}
        with a.out.open("x",encoding="utf-8") as h:json.dump(result,h,indent=2,allow_nan=False);h.write("\n")
    except (OSError,ValueError) as e:p.exit(2,f"{p.prog}: {e}\n")

if __name__=="__main__":main()
