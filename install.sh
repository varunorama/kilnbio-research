#!/usr/bin/env bash
# Install a versioned user-local copy; never changes global PATH/startup files.
set -Eeuo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PREFIX=${1:-"$HOME/.local/share/kilnbio"}
MODE=${2:-with-runtime}
[[ "$MODE" == with-runtime || "$MODE" == no-runtime ]] || { echo 'Second argument: with-runtime or no-runtime' >&2; exit 2; }
command -v python3 >/dev/null || { echo 'Python 3 is required' >&2; exit 2; }
VERSION=$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).read_text().strip())' "$ROOT/VERSION")
DEST=$(python3 - "$ROOT" "$PREFIX" "$VERSION" <<'PY'
import hashlib,pathlib,shutil,sys
root=pathlib.Path(sys.argv[1]).resolve();prefix=pathlib.Path(sys.argv[2]).expanduser().resolve();version=sys.argv[3]
dest=prefix/'releases'/version
if dest==root or root in dest.parents: raise SystemExit('Install outside the source directory')
allowed=['VERSION', 'LICENSE', 'CITATION.cff', 'CLI.txt', 'README.md', 'README.txt', 'THIRD_PARTY_NOTICES.txt', 'ACCESS.md', 'versions.json', 'run.sh', 'prepare.py', 'install.sh', 'setup_linux_runtime.sh', 'docs', 'upstream_licenses', 'upstream_citations', 'schemas', 'examples', 'environments', 'workflows', 'downstream', 'longread', 'modifications', 'singlecell', 'aging', 'tests', 'benchmarks', 'check_methylation_calls.py', 'check_sequencing_outputs.py', 'recipes', 'ci', 'evidence', 'assistants', '.claude/skills/kilnbio']
for name in allowed:
    src=root/name
    if not src.exists(): continue
    entries=src.rglob('*') if src.is_dir() else [src]
    for f in entries:
        if not f.is_file() or '__pycache__' in f.parts or f.suffix=='.pyc': continue
        out=dest/f.relative_to(root)
        if out.exists() and out.read_bytes()!=f.read_bytes():
            raise SystemExit(f'Existing release differs at {out}; choose a new prefix/version')
        out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,out)
print(dest)
PY
)
PREFIX=$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).expanduser().resolve())' "$PREFIX")
if [[ "$MODE" == with-runtime ]]; then bash "$DEST/setup_linux_runtime.sh" "$PREFIX/runtime"; fi
mkdir -p "$PREFIX/bin"
printf '#!/usr/bin/env bash\nset -Eeuo pipefail\nif [[ -f %q ]]; then source %q; fi\nexec bash %q "$@"\n' \
 "$PREFIX/runtime/env.sh" "$PREFIX/runtime/env.sh" "$DEST/run.sh" > "$PREFIX/bin/kilnbio"
chmod +x "$PREFIX/bin/kilnbio"
echo "Installed research release $VERSION. Command: $PREFIX/bin/kilnbio"
echo 'Configure your container runtime and run tiny validation tests on this machine before using real data.'
