"""Audit coverage of published EpiDiff CpGs. Sequencing input is one KilnBio CpG count import."""
import argparse
from bisect import bisect_right
from collections import defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path


MANIFEST_SHA = {
    "hg19": "88e995655e68b2867105e02924ad9f3ca6d51357e187676053bba935bf1fbc65",
    "hg38": "668a11ea624b3e645aab963dca048bdd4628ae5d297732867c91525c831cf191",
}
IDS_SHA = "daa98b732ee3abc5847e847dc339a41be74ffa07957c287ff2a60ffa83e42711"
MAP_SHA = {"hg19":"2e7385eef980cc2197d42a2153571e9e1cd8cd6cf6b5927feeedf33b4905bc21",
           "hg38":"28c5fe9b79d71324c9d86508fd2b382f1173ea0d878bc9c4d67a49bf310f67cd"}


def checked_map(path, assembly):
    with path.open("rb") as f:
        if hashlib.file_digest(f,"sha256").hexdigest() != MAP_SHA[assembly]:
            raise ValueError("Locus map differs from the pinned assembly; rebuild with make-map")
    return path
PRIMARY_CHROMS = {str(i) for i in range(1,23)} | {"X","Y"}


def make_map(assembly, manifest, ids, out):
    probes = sorted(ids.read_text().splitlines())
    if hashlib.sha256(("\n".join(probes)+"\n").encode()).hexdigest() != IDS_SHA:
        raise ValueError("Expected the 1,298 features of the pinned EpiDiff forest")
    with manifest.open("rb") as f:
        if hashlib.file_digest(f, "sha256").hexdigest() != MANIFEST_SHA[assembly]:
            raise ValueError("Manifest checksum differs from the reviewed assembly")
    wanted = set(probes)
    rows = []
    with gzip.open(manifest, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["probeID"] in wanted:
                pos = int(r["CpG_beg"]) + 1  # BED dyad to reference-forward C in chr,pos,N,X.
                if int(r["CpG_end"]) != pos + 1 or r["MASK_general"] not in {"TRUE","FALSE"}:
                    raise ValueError("Unexpected manifest interval or quality flag")
                rows.append((r["probeID"], r["CpG_chrm"], pos, int(r["MASK_general"] == "TRUE")))
    if sorted(r[0] for r in rows) != probes or len({r[1:3] for r in rows}) != len(probes):
        raise ValueError("Incomplete or ambiguous model-locus map")
    with out.open("x", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["probe","chr","pos","array_mask_general"])
        w.writerows(sorted(rows))


def chromosome(name):
    bare = name[3:] if name.startswith("chr") else name
    return "chr"+bare if bare in PRIMARY_CHROMS else name


def coverage(counts, loci, depth):
    if depth < 1:
        raise ValueError("Minimum depth must be positive")
    targets = {(r["chr"],int(r["pos"])): int(r["array_mask_general"]) for r in loci}
    observed = {}
    for row in counts:
        pos, n, x = (int(row[k]) for k in ("pos","N","X"))
        if pos < 1 or not 0 <= x <= n:
            raise ValueError("Invalid coordinate or modified/total counts")
        key = chromosome(row["chr"]), pos
        if key in targets:
            if key in observed:
                raise ValueError("Repeated target. Supply strand-collapsed counts")
            observed[key] = n
    passed = [k for k,n in observed.items() if n >= depth]
    return {"model_cpgs":len(targets), "min_depth":depth, "at_depth":len(passed),
            "below_depth_with_reads":sum(0 < n < depth for n in observed.values()),
            "no_reads_or_absent":len(targets)-sum(n > 0 for n in observed.values()),
            "array_masked_total":sum(targets.values()), "array_masked_at_depth":sum(targets[k] for k in passed)}


def rrbs_example(path, loci):
    # GSM1548420 reports two-base BED dyads and separate strand-specific X/N.
    with path.open("rb") as f:
        if hashlib.file_digest(f, "sha256").hexdigest() != "bdddc55b81431baa4f3f06bd818311831c159c84cb89e0a2027d3b9aceec5923":
            raise ValueError("Expected published GSM1548420 hg19 RRBS BED")
    targets = {(r["chr"],int(r["pos"])) for r in loci}
    counts = {}
    with gzip.open(path,"rt") as f:
        for line in f:
            chrom, start, end, ratio, _, strand = line.split()
            key = chrom, int(start)+1
            if key not in targets:
                continue
            if int(end) != key[1]+1 or strand not in {"+","-"}:
                raise ValueError("Unexpected GEO dyad or strand")
            x,n = map(int,ratio.strip("'").split("/"))
            strands = counts.setdefault(key,{})
            if strand in strands or not 0 <= x <= n:
                raise ValueError("Repeated strand or invalid source counts")
            strands[strand] = (x,n)
    for (chrom,pos), strands in counts.items():
        yield dict(chr=chrom,pos=pos,X=sum(x for x,n in strands.values()),N=sum(n for x,n in strands.values()))


def panel_coverage(path, loci):
    intervals = defaultdict(list)
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith(("#","track ","browser ")):
                continue
            chrom, start, end = line.split()[:3]
            start,end = int(start),int(end)
            if not 0 <= start < end:
                raise ValueError("Target BED requires nonempty zero-based half-open intervals")
            intervals[chromosome(chrom)].append((start,end))
    if not intervals:
        raise ValueError("Empty target design")
    # Duplicate and overlapping baits count once. A single bait must contain the whole CpG.
    lookup = {}
    for chrom, bounds in intervals.items():
        starts, ends = [], []
        for start,end in sorted(bounds):
            starts.append(start)
            ends.append(max(end,ends[-1]) if ends else end)
        lookup[chrom] = starts,ends
    covered = 0
    for r in loci:
        starts, ends = lookup.get(r["chr"],([],[]))
        beg = int(r["pos"])-1
        i = bisect_right(starts,beg)-1
        covered += i >= 0 and ends[i] >= beg+2
    return {"model_cpgs":len(loci),"dyads_in_design":covered,"dyads_outside_design":len(loci)-covered,
            "endpoint":"design_overlap_only"}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    m = sub.add_parser("make-map")
    m.add_argument("--assembly", choices=MANIFEST_SHA, required=True)
    for arg in ("manifest","ids","out"):
        m.add_argument("--"+arg, type=Path, required=True)
    a = sub.add_parser("audit")
    a.add_argument("--import-dir", type=Path, required=True)
    a.add_argument("--min-depth", type=int, default=10)
    a.add_argument("--out", type=Path, required=True)
    e = sub.add_parser("rrbs-example")
    e.add_argument("--geo-bed",type=Path,required=True)
    e.add_argument("--min-depth",type=int,default=10)
    e.add_argument("--out",type=Path,required=True)
    b = sub.add_parser("panel")
    b.add_argument("--assembly",choices=MANIFEST_SHA,required=True)
    b.add_argument("--bed",type=Path,required=True)
    b.add_argument("--out",type=Path,required=True)
    for subcommand in (a,b,e):
        subcommand.add_argument("--loci",type=Path,required=True,help="Local map for the matching assembly, generated with make-map")
    args = p.parse_args()
    if args.command == "make-map":
        make_map(args.assembly,args.manifest,args.ids,args.out)
        return
    if args.command == "panel":
        with checked_map(args.loci,args.assembly).open() as f:
            result = panel_coverage(args.bed,list(csv.DictReader(f,delimiter="\t")))
        result["assembly"] = args.assembly
        with args.out.open("x",encoding="utf-8") as f:
            f.write(json.dumps(result,indent=2)+"\n")
        return
    if args.command == "rrbs-example":
        assembly, chemistry = "hg19", "bisulfite"
        with checked_map(args.loci,"hg19").open() as f:
            loci = list(csv.DictReader(f,delimiter="\t"))
        result = coverage(rrbs_example(args.geo_bed,loci),loci,args.min_depth)
    else:
        status = json.loads((args.import_dir/"status.json").read_text(encoding="utf-8"))
        cfg = status["config"]
        if (status.get("complete") is not True or status.get("coordinates") != "1-based-cytosine"
                or status.get("context") != "CpG" or cfg.get("modification") != "modC"
                or cfg.get("measurement_reviewed") is not True):
            raise ValueError("Requires completed and reviewed combined-modC CpG count import")
        assembly = {"hg19":"hg19","GRCh37":"hg19","hg38":"hg38","GRCh38":"hg38"}[cfg["assembly"]]
        chemistry = cfg["chemistry"]
        with checked_map(args.loci,assembly).open() as ref, (args.import_dir/"counts.tsv").open() as inp:
            counts = csv.DictReader(inp,delimiter="\t")
            if counts.fieldnames != ["chr","pos","N","X"]:
                raise ValueError("Expected KilnBio chr,pos,N,X count columns")
            result = coverage(counts,csv.DictReader(ref,delimiter="\t"),args.min_depth)
    result.update(assembly=assembly, chemistry=chemistry, endpoint="coverage_only")
    with args.out.open("x",encoding="utf-8") as out:
        out.write(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
