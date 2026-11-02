"""Extract explicit germline calls for an exact panel of normalized diploid autosomal variants."""
import argparse
import csv
import gzip
import math
import re
from pathlib import Path
from assay_audit import read_table, unique

def extract(vcf, panel, samples, assembly, min_dp, min_gq):
    planned=unique(panel,"event"); ids=unique(samples,"sample")
    for r in samples:
        if not r["donor"].strip():raise ValueError("Sample map requires nonempty donor IDs")
    if len({r["donor"] for r in samples})!=len(samples):raise ValueError("Choose one germline sample per donor")
    if min_dp<0 or min_gq<0:raise ValueError("Quality thresholds must be nonnegative")
    keys={}; coord={}
    for event,r in planned.items():
        key=(r["chrom"],int(r["pos"]),r["ref"],r["alt"])
        if (r["assembly"]!=assembly or not re.fullmatch(r"(?:chr)?(?:[1-9]|1[0-9]|2[0-2])",key[0]) or key[1]<1 or
            not all(re.fullmatch("[ACGT]+",a) for a in key[2:]) or key[2]==key[3] or key in keys):raise ValueError("Require unique normalized autosomal biallelic alleles with matching build")
        keys[key]=event;coord.setdefault(key[:2],set()).add(event)
    reader=gzip.open if str(vcf).endswith((".gz",".bgz")) else open
    calls={}; found=set(); header=None; refs=[]
    with reader(vcf,"rt",encoding="utf-8") as f:
        for row in f:
            if row.startswith("##reference="):refs.append(row.strip().split("=",1)[1])
            if row.startswith("##"):continue
            if row.startswith("#CHROM\t"):
                if header is not None:raise ValueError("Repeated VCF header")
                header=row.rstrip("\r\n").split("\t")
                if len(set(header[9:]))!=len(header[9:]) or not ids.keys()<=set(header[9:]):raise ValueError("VCF sample IDs do not match the requested map")
                indices={s:header.index(s,9) for s in ids};continue
            if row.startswith("#"):continue
            if header is None:raise ValueError("Missing VCF sample header")
            c=row.rstrip("\r\n").split("\t")
            pos=(c[0],int(c[1]))
            if pos not in coord:continue
            key=(*pos,c[3],c[4])
            if key not in keys:continue
            event=keys[key]
            if event in found or len(c)!=len(header):raise ValueError("Duplicate targeted VCF allele or malformed row")
            found.add(event); format_keys=c[8].split(":")
            for sample,idx in indices.items():
                fmt=dict(zip(format_keys,c[idx].split(":"))); gt=fmt.get("GT",".")
                reason="";dosage=""
                if c[6]!="PASS":reason="not_PASS"
                elif fmt.get("FT","PASS") not in ("PASS","."):reason="sample_filter"
                elif not re.fullmatch(r"[01][/|][01]",gt):reason="invalid_or_missing_GT"
                else:
                    for threshold,k in [(min_dp,"DP"),(min_gq,"GQ")]:
                        if threshold:
                            try:v=float(fmt.get(k,"."))
                            except ValueError:v=float("nan")
                            if not math.isfinite(v) or v<threshold:reason="low_or_missing_"+k;break
                    if not reason:dosage=int(gt[0])+int(gt[2])
                calls[event,sample]=(dosage,reason)
    if header is None:raise ValueError("Missing VCF header")
    # A reference URL/path is provenance, not a substitute for checking the FASTA.
    canonical={"grch37":"37","hg19":"37","grch38":"38","hg38":"38"}
    declared={v for k,v in canonical.items() if any(k in ref.lower() for ref in refs)}
    want=canonical.get(assembly.lower())
    if want and declared and declared!={want}:raise ValueError("VCF reference header contradicts the requested build")
    result=[]
    for event in planned:
        for sample,r in ids.items():
            d,reason=calls.get((event,sample),("","no_exact_record"))
            result.append([event,"germline",r["donor"],r["specimen"],"unknown" if d=="" else "altered" if d else "reference",d,reason,sample,assembly,";".join(refs)])
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ("vcf","panel","samples"):p.add_argument(n,type=Path)
    p.add_argument("--assembly",required=True);p.add_argument("--min-dp",type=int,default=10);p.add_argument("--min-gq",type=int,default=20);p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    try:
        rows=extract(a.vcf,read_table(a.panel,["event","assembly","chrom","pos","ref","alt"]),read_table(a.samples,["sample","donor","specimen"]),a.assembly,a.min_dp,a.min_gq)
        with a.out.open("x",encoding="utf-8",newline="") as h:
            w=csv.writer(h,delimiter="\t",lineterminator="\n");w.writerow("event origin donor specimen state dosage reason source_sample assembly vcf_reference".split());w.writerows(rows)
    except (OSError,ValueError) as e:p.exit(2,f"{p.prog}: {e}\n")
if __name__=="__main__":main()
