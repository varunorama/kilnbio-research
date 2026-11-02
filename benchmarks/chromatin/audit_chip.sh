#!/usr/bin/env bash
# Independently count held-out single-end tags in native broad calls.
set -euo pipefail
study=$(realpath "${1:?benchmark directory}")
code=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
docker run --rm -i --network none -v "$study:/b" -v "$code/experiments.tsv:/panels.tsv:ro" \
    quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0 bash -s <<'BASH'
set -euo pipefail
audit=/b/audit/chip
mkdir -p "$audit"
bedtools makewindows -g /b/ref/panel.sizes -w 1000 | awk '$3-$2==1000' |
    bedtools intersect -v -a stdin -b /b/ref/mask.bed > "$audit/eligible.bed"
bedtools merge -i "$audit/eligible.bed" > "$audit/eligible_ranges.bed"
awk 'BEGIN {OFS="\t"} {print $1,$2+500,$2+501}' "$audit/eligible.bed" > "$audit/mid.bed"
printf 'experiment\trep\tcaller\tselected_bins\theldout_treatment_hits\theldout_treatment_total\theldout_input_hits\theldout_input_total\n'
while read -r exp; do
    for lib in control.test rep1 rep2; do
        awk 'BEGIN {OFS="\t"} {base=($6=="+" ? $2+100 : $3-101); if (base>=0) print $1,base,base+1}' "/b/$exp/$lib.bed" |
            bedtools intersect -u -a stdin -b "$audit/eligible_ranges.bed" > "$audit/$lib.bed"
    done
    ninput=$(wc -l < "$audit/control.test.bed")
    for rep in 1 2; do
        tags="$audit/rep$((3-rep)).bed"
        ntag=$(wc -l < "$tags")
        for caller in MACS3 Omnipeak; do
            name=calls_peaks.broadPeak
            if [[ $caller == Omnipeak ]]; then name=calls.broadPeak; fi
            bedtools intersect -u -a "$audit/mid.bed" -b "/b/$exp/r$rep-$caller/$name" |
                awk 'BEGIN {OFS="\t"} {print $1,$2-500,$2+500}' > "$audit/selected.bed"
            nbin=$(wc -l < "$audit/selected.bed")
            bedtools merge -i "$audit/selected.bed" > "$audit/selected_ranges.bed"
            chip=$(bedtools intersect -u -a "$tags" -b "$audit/selected_ranges.bed" | wc -l)
            input=$(bedtools intersect -u -a "$audit/control.test.bed" -b "$audit/selected_ranges.bed" | wc -l)
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$exp" "$rep" "$caller" "$nbin" "$chip" "$ntag" "$input" "$ninput"
        done
    done
done < <(awk -F'\t' 'NR>1 && $1!="ATAC" {print $3}' /panels.tsv)
BASH
