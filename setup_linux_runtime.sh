#!/usr/bin/env bash
# User-local runtime; does not change shell startup files or system Java.
set -Eeuo pipefail
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'This runtime installer currently supports Linux x86_64/WSL2 only.' >&2; exit 2; }
for tool in curl python3 tar sha256sum; do command -v "$tool" >/dev/null || { echo "Missing prerequisite: $tool" >&2; exit 2; }; done
DEST=${1:-"$HOME/kilnbio-runtime"}
mkdir -p "$DEST/downloads" "$DEST/bin" "$DEST/nextflow-home"
DEST=$(cd "$DEST" && pwd)
curl --fail --location --retry 2 --connect-timeout 20 --max-time 90 \
  'https://api.adoptium.net/v3/assets/latest/21/hotspot?architecture=x64&image_type=jre&os=linux&vendor=eclipse' \
  -o "$DEST/downloads/java-release.json"
python3 - "$DEST/downloads/java-release.json" > "$DEST/downloads/java.tsv" <<'PY'
import json,sys
release=json.load(open(sys.argv[1]))[0]
p=release['binary']['package']
print(p['link'],p['checksum'],sep='\t')
PY
IFS=$'\t' read -r URL SHA < "$DEST/downloads/java.tsv"
curl --fail --location --retry 2 --connect-timeout 20 --max-time 300 "$URL" -o "$DEST/downloads/java.tar.gz"
printf '%s  %s\n' "$SHA" "$DEST/downloads/java.tar.gz" | sha256sum -c -
mkdir -p "$DEST/java"
tar -xzf "$DEST/downloads/java.tar.gz" --strip-components=1 -C "$DEST/java"
curl --fail --location --retry 2 --connect-timeout 20 --max-time 90 \
  'https://api.github.com/repos/nextflow-io/nextflow/releases/tags/v25.10.4' -o "$DEST/downloads/nextflow-release.json"
python3 - "$DEST/downloads/nextflow-release.json" > "$DEST/downloads/nextflow.tsv" <<'PY'
import json,sys
r=json.load(open(sys.argv[1]))
assert not r['prerelease'] and not r['draft']
asset=next(a for a in r['assets'] if a['name'].endswith('-dist'))
print(asset['browser_download_url'],r['tag_name'],sep='\t')
PY
IFS=$'\t' read -r URL VERSION < "$DEST/downloads/nextflow.tsv"
curl --fail --location --retry 2 --connect-timeout 20 --max-time 300 "$URL" -o "$DEST/bin/nextflow"
printf '%s  %s\n' '60db869a7618783e3ee60183a2b0c6c50ffc67221eabec60800b904991786aaa' "$DEST/bin/nextflow" | sha256sum -c -
chmod +x "$DEST/bin/nextflow"
sha256sum "$DEST/bin/nextflow" > "$DEST/downloads/nextflow.sha256"
printf 'export JAVA_HOME=%q\nexport NXF_HOME=%q\nexport PATH=%q:"$PATH"\n' \
 "$DEST/java" "$DEST/nextflow-home" "$DEST/bin:$DEST/java/bin" > "$DEST/env.sh"
source "$DEST/env.sh"
java -version
nextflow -version
echo "Runtime installed. Source $DEST/env.sh in the Linux shell. Nextflow release: $VERSION"
