#!/usr/bin/env bash
# Independently reproduce native ATAC bin counts with BEDTools.
set -euo pipefail
study=$(realpath "${1:?benchmark directory}")
code=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
image=quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0
docker run --rm -i --network none -v "$study:/b" -v "$code/experiments.tsv:/panels.tsv:ro" "$image" bash -s <<'BASH'
set -euo pipefail
audit=/b/audit
mkdir -p "$audit"
bedtools makewindows -g /b/ref/panel.sizes -w 100 |
    awk '$3-$2==100' |
    bedtools intersect -v -a stdin -b /b/ref/mask.bed |
    awk 'BEGIN {OFS="\t"} {print $1,$2+50,$2+51}' > "$audit/eligible.bed"
printf 'experiment\trep\tcaller\tselected_bins\tdnase_reference_bins\tdnase_overlap_bins\n'
while read -r exp dnase1 dnase2; do
    reference="$audit/$exp.reference.bed"
    bedtools intersect -u -a "$audit/eligible.bed" -b "/b/input/$dnase1.bed.gz" |
        bedtools intersect -u -a stdin -b "/b/input/$dnase2.bed.gz" > "$reference"
    positives=$(wc -l < "$reference")
    for r in 1 2; do
        for caller in MACS2 MACS3 HMMRATAC; do
            name=calls_peaks.narrowPeak
            if [[ $caller == HMMRATAC ]]; then name=calls_accessible_regions.narrowPeak; fi
            peaks="/b/$exp/r$r-$caller/$name"
            selected=$(bedtools intersect -u -a "$audit/eligible.bed" -b "$peaks" | wc -l)
            hits=$(bedtools intersect -u -a "$reference" -b "$peaks" | wc -l)
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$exp" "$r" "$caller" "$selected" "$positives" "$hits"
        done
    done
done < <(awk -F'\t' '$1=="ATAC" {print $3,$9,$10}' /panels.tsv)
BASH
