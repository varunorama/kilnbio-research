"""Audit an explicit target design, measured features and optional DNA alignment depth."""
import argparse
from bisect import bisect_right
from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from analysis import digest


PRESETS = {
    "twist-human-methylome-hg38": "82f54b6e17d92c98e07e9204cddeffad9bd4b7882e925e123954013b33b65398",
    "twist-exome-2-hg38": "8fa8f1a9fd5d2dcbefa2b48713e8b242065eb60bef0f9b8482322246ed4d5a77",
}


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8-sig") if path.suffix == ".gz" else path.open(encoding="utf-8-sig")


def reference_lengths(path):
    lengths = {}
    with open_text(path) as handle:
        for line in handle:
            chrom, length, *_ = line.rstrip("\r\n").split("\t")
            if chrom in lengths or int(length) < 1:
                raise ValueError(f"Invalid reference length for {chrom}")
            lengths[chrom] = int(length)
    if not lengths:
        raise ValueError("Empty reference size table")
    return lengths


def check_interval(chrom, start, end, lengths):
    if chrom not in lengths or not 0 <= start < end <= lengths[chrom]:
        raise ValueError(f"Unknown contig or coordinates outside reference: {chrom} {start} {end}")


def bed_rows(path, lengths):
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.rstrip("\r\n").split("\t")
            chrom, left, right = fields[:3]
            start, end = int(left), int(right)
            check_interval(chrom, start, end, lengths)
            yield chrom, start, end, fields[3] if len(fields) > 3 else "."


def target_union(path, lengths):
    intervals = defaultdict(list)
    number = 0
    for chrom, start, end, _ in bed_rows(path, lengths):
        intervals[chrom].append((start, end))
        number += 1
    if not number:
        raise ValueError("Empty target BED")
    result = {}
    for chrom, rows in intervals.items():
        merged = []
        for start, end in sorted(rows):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        result[chrom] = merged
    return result, number


def within(chrom, start, end, regions):
    rows = regions.get(chrom, ())
    i = bisect_right(rows, (start, float("inf"))) - 1
    return i >= 0 and end <= rows[i][1]


def measurement_rows(path, kind, lengths, assembly):
    if kind == "bed":
        yield from bed_rows(path, lengths)
        return
    with open_text(path) as handle:
        if kind == "cpg":
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["chr", "pos", "N", "X"]:
                raise ValueError("Expected a strand-collapsed KilnBio chr,pos,N,X count table")
            last = {}
            for row in reader:
                start = int(row["pos"]) - 1
                if start <= last.get(row["chr"], -1):
                    raise ValueError("CpG positions must be unique and increasing on each contig")
                last[row["chr"]] = start
                n, x = int(row["N"]), int(row["X"])
                if n <= 0 or not 0 <= x <= n:
                    raise ValueError("Expected a covered CpG with valid integer counts")
                check_interval(row["chr"], start, start + 2, lengths)
                yield row["chr"], start, start + 2, f"{x}/{n}"
        else:
            header = False
            builds = {"hg19": "37", "grch37": "37", "hg38": "38", "grch38": "38"}
            for line in handle:
                if line.startswith("##reference="):
                    declared = {v for k, v in builds.items() if k in line.lower()}
                    expected = builds.get(assembly.lower())
                    if expected and declared and declared != {expected}:
                        raise ValueError("VCF reference header contradicts requested assembly")
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM\t"):
                    header = True
                    continue
                if not header:
                    raise ValueError("Missing VCF column header")
                chrom, pos, _, ref, alt, _, filt, *_ = line.rstrip().split("\t")
                if filt != "PASS" or any(c not in "ACGT" for c in ref + alt) or not ref or not alt:
                    continue
                start = int(pos) - 1
                check_interval(chrom, start, start + len(ref), lengths)
                yield chrom, start, start + len(ref), f"{ref}>{alt}"
            if not header:
                raise ValueError("Missing VCF column header")


def depth_summary(args, stage, regions, lengths, bases):
    header = subprocess.check_output(["samtools", "view", "-H", str(args.bam)], text=True)
    sam_lengths = {}
    for line in header.splitlines():
        if line.startswith("@SQ\t"):
            tags = dict(field.split(":", 1) for field in line.split("\t")[1:])
            sam_lengths[tags["SN"]] = int(tags["LN"])
    if any(sam_lengths.get(c) != lengths[c] for c in regions):
        raise ValueError("BAM contigs or lengths disagree with the target reference")
    cmd = ["samtools", "depth", "-s", "-q", str(args.min_bq), "-Q", str(args.min_mq),
           "-b", str(stage / "targets.bed"), str(args.bam)]
    covered = passing = total_depth = 0
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as proc:
        for line in proc.stdout:
            chrom, pos, depth = line.rstrip().split("\t")
            n = int(depth)
            if not within(chrom, int(pos) - 1, int(pos), regions):
                raise ValueError("Samtools emitted an unexpected position outside targets")
            total_depth += n
            covered += int(n > 0)
            passing += int(n >= args.min_depth)
        if proc.wait():
            raise RuntimeError("Samtools depth failed")
    return {"target_bases": bases, "covered_bases": covered, "zero_depth_bases": bases - covered,
            "bases_at_min_depth": passing, "fraction_at_min_depth": passing / bases,
            "mean_depth": total_depth / bases, "min_depth": args.min_depth,
            "min_mapping_quality": args.min_mq, "min_base_quality": args.min_bq,
            "samtools": subprocess.check_output(["samtools", "--version-only"], text=True).strip(),
            "overlapping_mates": "samtools depth -s; count once"}


def run(args):
    out = args.out.resolve()
    if out.exists():
        raise ValueError("Choose a new output directory")
    if args.min_depth < 1 or min(args.min_mq, args.min_bq) < 0:
        raise ValueError("Invalid depth or quality thresholds")
    if (args.measurements or args.bam) and args.data_assembly != args.assembly:
        raise ValueError("Measurement and panel assembly labels must match")
    signature = digest(args.targets)
    if args.preset != "custom" and (args.assembly not in ("hg38", "GRCh38") or signature != PRESETS[args.preset]):
        raise ValueError("Fixed preset requires the original vendor hg38 BED; use custom for a revised design")
    lengths = reference_lengths(args.sizes)
    regions, rows = target_union(args.targets, lengths)
    bases = sum(end-start for spans in regions.values() for start,end in spans)
    report = {"panel": args.name, "preset": args.preset, "assembly": args.assembly, "target_rows": rows,
              "target_union_bases": bases, "target_bed_sha256": signature, "sizes_sha256": digest(args.sizes)}
    if args.baits:
        baits, bait_rows = target_union(args.baits, lengths)
        report["baits"] = {"rows": bait_rows, "union_bases": sum(e-s for v in baits.values() for s,e in v),
                           "sha256": digest(args.baits)}
    out.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".capture-", dir=out.parent))
    try:
        with (stage/"targets.bed").open("w", newline="") as h:
            for chrom in lengths:
                for start,end in regions.get(chrom, ()):
                    h.write(f"{chrom}\t{start}\t{end}\n")
        if args.measurements:
            tested = inside = 0
            with (stage/"features.tsv").open("w", newline="", encoding="utf-8") as h:
                writer = csv.writer(h, delimiter="\t", lineterminator="\n")
                writer.writerow(["chrom","start","end","feature","fully_on_target"])
                for chrom,start,end,name in measurement_rows(args.measurements,args.format,lengths,args.assembly):
                    hit = within(chrom,start,end,regions)
                    tested += 1
                    inside += hit
                    writer.writerow([chrom,start,end,name,int(hit)])
            report["features"] = {"format":args.format,"tested":tested,"fully_on_target":inside,
                                  "outside_or_partial":tested-inside,"input_sha256":digest(args.measurements)}
        if args.bam:
            report["alignment_depth"] = depth_summary(args,stage,regions,lengths,bases)
            report["alignment_depth"]["bam_sha256"] = digest(args.bam)
        (stage/"summary.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
        stage.rename(out)
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--targets", type=Path, required=True, help="BED3+ covered targets (0-based, half-open)")
    p.add_argument("--baits", type=Path, help="Optional separate capture probe intervals")
    p.add_argument("--name", required=True, help="Panel and design version or custom identifier")
    p.add_argument("--preset", choices=["custom", *PRESETS], default="custom")
    p.add_argument("--assembly", required=True)
    p.add_argument("--sizes", type=Path, required=True, help="Matching FASTA index or chromosome sizes")
    p.add_argument("--data-assembly", help="Required for any measurement input")
    p.add_argument("--measurements", type=Path, help="BED, called VCF, or KilnBio CpG counts")
    p.add_argument("--format", choices=["bed","vcf","cpg"], default="bed")
    p.add_argument("--bam", type=Path, help="Coordinate-sorted alignment for DNA coverage")
    p.add_argument("--min-depth", type=int, default=10)
    p.add_argument("--min-mq", type=int, default=20)
    p.add_argument("--min-bq", type=int, default=20)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    try:
        run(args)
    except (OSError,ValueError,RuntimeError,subprocess.CalledProcessError) as e:
        p.exit(2, f"ERROR: {e}\n")


if __name__ == "__main__":
    main()
