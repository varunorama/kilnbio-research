#!/usr/bin/env bash
set -euo pipefail
study=$(realpath "${1:?benchmark directory}")
while read -r tag sha; do
    docker pull "$tag"
    repo=${tag%:*}
    docker inspect --format '{{range .RepoDigests}}{{println .}}{{end}}' "$tag" | grep -Fx "$repo@sha256:$sha"
done <<'IMAGES'
quay.io/biocontainers/macs2:2.2.7.1--py38h4a8c8d9_3 b8506329f67a88c4c9c0e433028c25ba5c29b199a8936f629a02c40a4ecf1942
quay.io/biocontainers/macs3:3.1.0--py312h9639c01_0 f0f328c01b1d08f042519de66542446501119cdf976815d087d5e6e192592397
quay.io/biocontainers/macs3:3.0.1--py311h0152c62_3 38967cf3d3217209bcb8baa9f277c53910e1d48d1afd825008ea27848f83c8d7
community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd a55ddea590e567a91df592300a960aa534cfc1bd16e7623e3938ec21f4f3df15
quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0 b0018bd0a10853e19ee92f6d46d8d12f1c41e516845105e1f02c91b4d7b961b1
IMAGES
cd "$study"
jar=omnipeak-1.6.6883.jar
if [[ ! -f $jar ]]; then
    curl -fL --retry 3 "https://download.jetbrains.com/biolabs/omnipeak/$jar" -o "$jar.partial"
    mv "$jar.partial" "$jar"
fi
echo "562e94981fa6064bc3633651349527bdf75ce19292e39fd38dd594c482a35c11  $jar" | sha256sum -c -
