#!/usr/bin/env bash
# Bootstrap the managed installation with the host's standard-library Python.
set -Eeuo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
command -v python3 >/dev/null || { echo 'Install Python 3.8+ first. On Ubuntu: sudo apt-get install python3 curl ca-certificates' >&2; exit 2; }
exec python3 "$ROOT/environments/install.py" "$@"
