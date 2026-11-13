#!/usr/bin/env bash
# Requires the scientific Python and R dependencies; no network or private inputs.
set -Eeuo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
[[ $# == 1 ]] || { echo 'Usage: bash ci/run_platform_checks.sh NEW_OUTPUT' >&2; exit 2; }
"${KILNBIO_PYTHON:-python3}" "$ROOT/tests/integration_platform_extensions.py"
"${KILNBIO_RSCRIPT:-Rscript}" --vanilla "$ROOT/tests/integration_multiomic_layers.R" "$ROOT" "$1"
"${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/multiomic_meta.py" "$1/association/manifest.tsv" --out "$1/joint"
