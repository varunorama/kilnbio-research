#!/usr/bin/env bash
set -euo pipefail
base=$(realpath "${1:?benchmark directory}")
id=${2:?ENCODE accession}
kind=${3:?ATAC or ChIP}
[[ $id =~ ^ENCFF[A-Z0-9]+$ && -s "$base/input/$id.bam.md5" ]]
mkdir -p "$base/tags"
target="$base/tags/$id.bed"
exec 9>"$target.lock"
flock 9
[[ ! -f $target ]] || exit 0
sam=(docker run --rm --network none -v "$base:/b:ro" community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd samtools view -@ 2 -q 30 -L /b/ref/panel.bed)
bed=(docker run --rm -i --network none -v "$base/ref:/ref:ro" quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0)
case "$kind" in
    ATAC)
        "${sam[@]}" -f 2 -F 3852 "/b/input/$id.bam" |
            awk 'BEGIN {OFS="\t"} $7=="=" && $9>0 && $9<=2000 {print $3,$4-1,$4-1+$9}' |
            "${bed[@]}" bedtools intersect -v -a stdin -b /ref/mask.bed > "$target.partial"
        ;;
    ChIP)
        "${sam[@]}" -u -F 3844 "/b/input/$id.bam" |
            "${bed[@]}" bash -o pipefail -c 'bedtools bamtobed -i stdin | bedtools intersect -v -a stdin -b /ref/mask.bed' |
            awk 'BEGIN {OFS="\t"} {print $1,$2,$3,".",$5,$6}' > "$target.partial"
        ;;
    *) exit 2 ;;
esac
[[ -s "$target.partial" ]]
mv "$target.partial" "$target"
