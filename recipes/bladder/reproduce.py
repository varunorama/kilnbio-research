"""Run the blinded bladder handoff from published matrices; no third-party Python packages."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
import time

from fetch import sha, stage

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CLASSES = {"Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site",
           "Translation_Start_Site", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins"}


def read_rows(path, comments=False):
    with path.open(encoding="utf-8-sig", newline="") as f:
        rows = csv.DictReader((s for s in f if not (comments and s.startswith("#"))), delimiter="\t")
        if not rows.fieldnames or len(set(rows.fieldnames)) != len(rows.fieldnames):
            raise ValueError(f"Empty or duplicate headers: {path}")
        yield from rows


def table(path, fields, rows):
    with path.open("x", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(fields)
        w.writerows(rows)


def columns(path):
    with path.open(encoding="utf-8-sig") as f:
        h = f.readline().rstrip("\r\n").split("\t")
    if len(set(h)) != len(h):
        raise ValueError(f"Duplicate sample IDs or columns in {path}")
    return h


def prepare(cache, out):
    rfile = cache / "data_mrna_seq_v2_rsem.txt"
    pfile = cache / "data_rppa.txt"
    primary = {}
    for row in read_rows(cache / "data_clinical_sample.txt", True):
        sid, donor = row["SAMPLE_ID"], row["PATIENT_ID"]
        if re.fullmatch(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01", sid):
            if sid in primary or donor != sid[:12] or row["SAMPLE_TYPE"] != "Primary":
                raise ValueError(f"Ambiguous primary specimen {sid}")
            primary[sid] = donor
    if len(primary) != len(set(primary.values())):
        raise ValueError("More than one primary specimen per donor")
    lines = (cache / "case_lists/cases_sequenced.txt").read_text().splitlines()
    [case_list] = [l.split(":", 1)[1] for l in lines if l.startswith("case_list_ids:")]
    profiled = set(case_list.split()) & primary.keys()
    mutated = {r["Tumor_Sample_Barcode"][:15] for r in read_rows(cache / "data_mutations.txt", True)
               if r["Hugo_Symbol"] == "FGFR3" and r["Variant_Classification"] in CLASSES}
    allrna = sorted(profiled & set(columns(rfile)))
    paired = sorted(set(allrna) & set(columns(pfile)))
    eligible = set(paired)
    if not eligible or set(allrna) == eligible:
        raise ValueError("Cannot estimate both contrasts without matched and RNA-only patients")
    groups = {}
    for key, ids in [("RNA_all", allrna), ("matched", paired), ("RNA_only", sorted(set(allrna)-eligible))]:
        n = sum(s in mutated for s in ids)
        groups[key] = {"patients": len(ids), "positive": n, "no_qualifying_call": len(ids)-n}
        if key != "RNA_only" and min(n, len(ids)-n) < 3:
            raise ValueError(f"Insufficient group size: {key}")
    table(out/"roster.tsv", ["specimen","donor","fgfr3","base_ok","add_ok"],
          ((s,primary[s],int(s in mutated),1,int(s in eligible)) for s in allrna))
    table(out/"observations.tsv", ["donor","specimen","tissue","assay","measurement"],
          ((primary[s],s,"primary",a,a+"_"+s) for s in allrna for a in
           (("RNA","protein") if s in eligible else ("RNA",))))
    audit = [sys.executable, str(ROOT/"downstream/assay_audit.py")]
    subprocess.run(audit+["pair",str(out/"roster.tsv"),str(out/"observations.tsv"),"--base","RNA","--add","protein",
                         "--tissue","primary","--out",str(out/"pairing.json")],check=True)
    subprocess.run(audit+["reach",str(out/"roster.tsv"),"--out",str(out/"reach.json")],check=True)
    verified = json.loads((out/"pairing.json").read_text())
    keys = sorted((r["specimen"],r["base_measurement"],r["additional_measurement"]) for r in verified["unique_pairs"])
    if verified["review_required"] or keys != [(s,"RNA_"+s,"protein_"+s) for s in paired]:
        raise ValueError("Audited pairs differ from exact sample matches; review required")

    pdata = list(read_rows(pfile))
    probes = []
    for r in pdata:
        label = r["Composite.Element.REF"]
        gene, probe = (s.strip() for s in label.split("|",1))
        reason = ("ambiguous_symbol" if not re.fullmatch(r"[A-Za-z0-9-]+",gene) else
                  "contrast_gene" if gene == "FGFR3" else
                  "modified_or_specific_probe" if re.search(r"[_-]p[STY]\d|cleav|acetyl|phosph|methyl|ubiquit|mutant",probe,re.I) else "")
        probes.append([label, gene, reason])
    nprobe = Counter(g for _,g,reason in probes if not reason)
    wanted = {g for g,n in nprobe.items() if n==1}
    rna, hits = {}, Counter()
    for r in read_rows(rfile):
        g = r["Hugo_Symbol"]
        if g in wanted:
            hits[g] += 1
            rna[g] = [r[s] for s in allrna]
    both = sorted(g for g in wanted if hits[g]==1)
    for entry in probes:
        _,g,why = entry
        entry[2] = why or ("multiple_unmodified_probes" if nprobe[g]!=1 else
                            "missing_or_nonunique_RNA" if hits[g]!=1 else "mapped")
    table(out/"antibody_map.tsv", ["source","gene","mapping_status"],probes)
    pm = {g:[r[s] for s in paired] for (_,g,why),r in zip(probes,pdata) if why=="mapped"}
    table(out/"rna_matrix.tsv",["gene"]+allrna,([g]+rna[g] for g in both))
    table(out/"protein_matrix.tsv",["gene"]+paired,([g]+pm[g] for g in both))
    return {"groups":groups,"clinical_primary":len(primary),"mutation_profiled_primary":len(profiled),
            "mapped_before_coverage_filter":len(both)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",required=True,type=Path,help="New results directory")
    p.add_argument("--cache",required=True,type=Path,help="Separate public-input cache")
    p.add_argument("--offline",action="store_true",help="Verify all cache files without network access")
    p.add_argument("--rscript",default="Rscript")
    p.add_argument("--r-library",type=Path,help="Optional extra library containing limma")
    a = p.parse_args()
    out,cache = a.out.resolve(),a.cache.resolve()
    if out.exists() or cache==out or cache in out.parents or out in cache.parents:
        p.error("Choose a new output directory separate from the cache")
    tick = time.monotonic()
    start = datetime.now(timezone.utc).isoformat()
    subprocess.run([a.rscript,"--vanilla",str(HERE/"rank.R"),"--check",str(a.r_library.resolve()) if a.r_library else ""],check=True)
    source = stage(cache,a.offline)
    fetched = time.monotonic()
    out.mkdir(parents=True)
    report = prepare(cache,out)
    prepared = time.monotonic()
    subprocess.run([a.rscript,"--vanilla",str(HERE/"rank.R"),str(out),str(a.r_library.resolve()) if a.r_library else ""],check=True)
    rows = list(read_rows(out/"rankings.tsv"))
    panels = {arm:{r["gene"]:int(r["direction"]) for r in rows if r["arm"]==arm and int(r["rank"])<=15}
              for arm in ("A","B","C")}
    comparisons = {}
    changed = []
    for before,after in ("AB","BC","AC"):
        x,y = panels[before],panels[after]
        if len(x)!=15 or len(y)!=15:
            raise ValueError("Incomplete shortlist")
        shared = x.keys() & y.keys()
        contrasts = {"retained":sorted(shared), "entered":sorted(y.keys()-x.keys()),"exited":sorted(x.keys()-y.keys()),
                     "sign_reversals":sorted(g for g in shared if x[g]*y[g]<0)}
        comparisons[before+"_to_"+after] = contrasts
        changed.extend((before+"_to_"+after,g,x.get(g,""),y.get(g,""),
                        "entered" if g not in x else "exited" if g not in y else
                        "reversed" if x[g]*y[g]<0 else "retained") for g in sorted(x.keys()|y.keys()))
    table(out/"changes.tsv",["comparison","gene","before_direction","after_direction","status"],changed)
    report.update(common_genes=sum(r["arm"]=="A" for r in rows),budget=15,comparisons=comparisons)
    (out/"report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    receipt = {"version":(ROOT/"VERSION").read_text().strip(),"started_utc":start,"python":sys.version,
               "platform":platform.platform(),"protocol_sha256":sha(HERE/"PROTOCOL.md"),"source_manifest":source,
               "offline":a.offline,"seconds":{"dependencies_and_inputs":round(fetched-tick,2),
               "selection_and_pairing":round(prepared-fetched,2),"model_and_report":round(time.monotonic()-prepared,2)},
               "outputs_sha256":{x.name:sha(x) for x in sorted(out.iterdir()) if x.is_file()}}
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print("Finished:",out/"report.json")


if __name__=="__main__":
    main()
