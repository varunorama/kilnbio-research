"""Make assay-audit inputs from open GDC STAR-counts and 450K beta file metadata."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import urllib.parse
import urllib.request

API="https://api.gdc.cancer.gov/"
BARCODE=re.compile(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]-[0-9]{2}[A-Z]-[A-Z0-9]{4}-[A-Z0-9]{2}")
def get(route, params):
    url=API+route+("?" + urllib.parse.urlencode(params) if params else "")
    with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":"KilnBio-GDC-metadata"}),timeout=90) as resp:
        return json.load(resp)
def query(tag,project):
    cond=[("project.project_id" if tag=="cases" else "cases.project.project_id",project)]
    if tag!="cases":
        cond += [("access","open")]
        cond += [("data_type","Gene Expression Quantification"),("analysis.workflow_type","STAR - Counts")] if tag=="rna" else [("data_type","Methylation Beta Value"),("platform","Illumina Human Methylation 450")]
    fields="case_id,submitter_id,project.project_id" if tag=="cases" else "file_id,file_name,data_type,analysis.workflow_type,platform,access,cases.submitter_id,associated_entities.entity_type,associated_entities.entity_submitter_id,associated_entities.case_id"
    filters={"op":"and","content":[{"op":"in","content":{"field":k,"value":[v]}} for k,v in cond]}
    return {"filters":json.dumps(filters),"fields":fields,"size":500,"sort":("case_id" if tag=="cases" else "file_id")+":asc"}
def records(raw,tag,project,download):
    params=query(tag,project)
    if download:
        pages=[];n=None;offset=0
        while n is None or offset<n:
            data=get("cases" if tag=="cases" else "files",dict(params,**{"from":offset}))["data"]
            if n is not None and n!=data["pagination"]["total"]:raise ValueError("Source changed during pagination")
            n=data["pagination"]["total"];hits=len(data["hits"])
            if not hits and n:raise ValueError("Empty page before end")
            blob=(json.dumps(data,sort_keys=True)+"\n").encode()
            (raw/f"{tag}.{offset}.json").write_bytes(blob)
            pages.append(dict(offset=offset,hits=hits,sha256=hashlib.sha256(blob).hexdigest()));offset+=hits
        (raw/f"{tag}.receipt.json").write_text(json.dumps(dict(params=params,total=n,pages=pages),indent=2),encoding="utf-8")
    receipt=json.loads((raw/f"{tag}.receipt.json").read_text(encoding="utf-8"))
    if receipt["params"]!=params:raise ValueError("Cache does not match requested project/assay query")
    rows=[]
    for page in receipt["pages"]:
        blob=(raw/f"{tag}.{int(page['offset'])}.json").read_bytes()
        if hashlib.sha256(blob).hexdigest()!=page["sha256"]:raise ValueError("Page checksum mismatch")
        data=json.loads(blob)
        if len(rows)!=page["offset"] or len(data["hits"])!=page["hits"]:raise ValueError("Incorrect pagination receipt")
        rows+=data["hits"]
    if len(rows)!=receipt["total"] or len({r["id"] for r in rows})!=len(rows):raise ValueError("Incomplete or duplicate page data")
    return rows

def convert(cases,assays):
    identifiers={c["case_id"]:c["submitter_id"] for c in cases}
    if len(identifiers)!=len(cases) or len(set(identifiers.values()))!=len(cases):raise ValueError("Ambiguous case identifiers")
    obs=[];review=[]
    for tag,files in assays.items():
        for f in files:
            entities=f.get("associated_entities",[]);cs=f.get("cases",[])
            e=entities[0] if len(entities)==1 else {}
            barcode=e.get("entity_submitter_id","")
            if (len(cs)!=1 or e.get("entity_type")!="aliquot" or not BARCODE.fullmatch(barcode)
                    or identifiers.get(e.get("case_id"))!=barcode[:12] or cs[0].get("submitter_id")!=barcode[:12]):
                review.append((tag,f["file_id"],"ambiguous_or_unverified_file_link"));continue
            obs.append((barcode[:12],barcode[:16],barcode[13:15],tag,f["file_id"]))
    return sorted(identifiers.values()),obs,review

def write(path,headers,rows):
    with path.open("x",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(headers);w.writerows(rows)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("project");p.add_argument("out",type=Path)
    p.add_argument("--cache",type=Path,help="Replay a previous raw metadata folder, verifying query and page checksums")
    a=p.parse_args()
    if not re.fullmatch(r"TCGA-[A-Z0-9]+",a.project):p.error("Expected TCGA project code")
    a.out.mkdir(parents=True,exist_ok=False)
    raw=a.cache if a.cache else a.out/"raw"
    if not a.cache:
        raw.mkdir();(raw/"api-status.json").write_text(json.dumps(get("status",{}),indent=2),encoding="utf-8")
    cases,rna,dnam=[records(raw,tag,a.project,not a.cache) for tag in ("cases","rna","dnam")]
    ids,obs,review=convert(cases,{"rna":rna,"dnam":dnam})
    write(a.out/"donors.tsv",["donor"],[(d,) for d in ids])
    write(a.out/"observations.tsv","donor specimen tissue assay measurement".split(),obs)
    write(a.out/"review.tsv",["assay","measurement","reason"],review)
    print(f"{len(ids)} cases; {len(obs)} verified file links; {len(review)} unresolved file links. See review.tsv.")
if __name__=="__main__":main()
