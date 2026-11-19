#!/usr/bin/env bash
# KilnBio bulk short-read sequencing launcher. No project or client defaults.
set -Eeuo pipefail
trap 'printf "ERROR: launcher failed at line %s\n" "$LINENO" >&2' ERR
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Reports are part of a successful command, including commands started by study/MCP.
# Do not wrap study itself: each child command must generate its own report.
if [[ ${KILNBIO_REPORT_WRAPPED:-} != "${1:-}" ]]; then
  case "${1:-}" in
    variant-annotations|alphafold|vep-annotate|alphagenome|multiomic-layers|genomic-associate|differential-compare|compare-assays|featurecounts|spikein-factors|de|dm|scrna|pseudobulk|import-mods|ont|wgbs|rrbs|atacseq|chipseq|rnaseq|scrnaseq|exome|ch-somatic|olink-import|ch-screen|ch-import|capture-panel|genomic-panel|method-agreement|assay-audit|candidate-evidence|compare-predictions|clock-panel|clocks|assay-clocks|epidiff-coverage|erythroid-stage|graft-readouts|bladder-example|decision-methods)
      exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/auto_reports.py" "$@" ;;
  esac
fi
# Direct scientific modules share the installed entry point. They do not call nf-core.
case "${1:-}" in
  decision-methods) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/decision_methods.py" "$@" ;;
  assistant-serve) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/assistants/kilnbio_mcp.py" "$@" ;;
  assistant-install) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/assistants/install_client.py" "$@" ;;
  bladder-example)
    shift
    R_OPTIONS=()
    if [[ -n ${KILNBIO_R_LIBRARY:-} ]]; then R_OPTIONS+=(--r-library "$KILNBIO_R_LIBRARY"); fi
    exec "${KILNBIO_PYTHON:-python3}" "$ROOT/recipes/bladder/reproduce.py" --rscript "${KILNBIO_RSCRIPT:-Rscript}" "${R_OPTIONS[@]}" "$@"
    ;;
  study) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/study.py" "$@" ;;
  variant-annotations) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/variant_annotations.py" "$@" ;;
  alphafold) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/alphafold_lookup.py" "$@" ;;
  vep-annotate) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/vep_annotate.py" "$@" ;;
  alphagenome) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/alphagenome_query.py" "$@" ;;
  olink-import) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/olink_import.py" "$@" ;;
  multiomic-meta) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/multiomic_meta.py" "$@" ;;
  circos) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/circos_report.py" "$@" ;;
  pathway-compare) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/pathway_compare.py" "$@" ;;
  multiomic-layers) shift; exec "${KILNBIO_RSCRIPT:-Rscript}" --vanilla "$ROOT/downstream/multiomic_layers.R" "$@" ;;
  ch-screen) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/clonal_hematopoiesis.py" screen "$@" ;;
  ch-import) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/clonal_hematopoiesis.py" import "$@" ;;
  capture-panel) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/capture_panel.py" "$@" ;;
  genomic-panel) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/genomic_panel.py" "$@" ;;
  genomic-associate) shift; exec "${KILNBIO_RSCRIPT:-Rscript}" --vanilla "$ROOT/downstream/genomic_associate.R" "$@" ;;
  method-agreement) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/method_agreement.py" "$@" ;;
  differential-compare) shift; exec "${KILNBIO_RSCRIPT:-Rscript}" --vanilla "$ROOT/downstream/differential_compare.R" "$@" ;;
  assay-audit) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/assay_audit.py" "$@" ;;
  candidate-evidence) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/candidate_evidence.py" "$@" ;;
  compare-assays) shift; exec "${KILNBIO_RSCRIPT:-Rscript}" --vanilla "$ROOT/downstream/compare_assays.R" "$@" ;;
  compare-predictions) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/benchmarks/public_modeling/paired_predictions.py" "$@" ;;
  clock-panel) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/aging/clock_panel.py" "$@" ;;
  featurecounts) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/featurecounts.py" "$@" ;;
  spikein-factors) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/spikein.py" "$@" ;;
  de|dm) exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/analysis.py" "$@" ;;
  enrich) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/enrichment.py" "$@" ;;
  scrna|pseudobulk) exec "${KILNBIO_PYTHON:-python3}" "$ROOT/singlecell/analyze.py" "$@" ;;
  import-mods) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/modifications/import_counts.py" "$@" ;;
  clocks) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/aging/clocks.py" "$@" ;;
  assay-clocks) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/aging/assay_clocks.py" "$@" ;;
  epidiff-coverage) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/aging/epidiff_coverage.py" "$@" ;;
  erythroid-stage) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/aging/erythroid_stages.py" "$@" ;;
  graft-readouts) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/downstream/graft_readouts.py" "$@" ;;
  ont) shift; exec "${KILNBIO_PYTHON:-python3}" "$ROOT/longread/ont.py" "$@" ;;
  help|--help|-h) cat "$ROOT/CLI.txt"; exit 0 ;;
esac
if [[ $# -lt 4 ]]; then
  echo "Usage: kilnbio {wgbs|rrbs|atacseq|chipseq|rnaseq|scrnaseq|exome|ch-somatic} samples.csv params.json RUN_DIR [--dry-run] [--resume] [--profile docker|apptainer|singularity] [--cpus N] [--memory-gb N] (source: bash run.sh ...)" >&2
  exit 2
fi
ASSAY=$1; SAMPLES=$2; PARAMS=$3; RUN=$4; shift 4
DRY=0; RESUME=0; PROFILE=docker; CPUS=4; MEMORY=24
if [[ -r /proc/meminfo ]]; then
  MEMORY=$(python3 -c 'import re; n=int(re.search(r"MemTotal:\s+(\d+)",open("/proc/meminfo").read()).group(1)); print(min(24,max(2,int(n/1048576*0.75))))')
fi
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY=1; shift ;;
    --resume) RESUME=1; shift ;;
    --profile|--cpus|--memory-gb)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      case "$1" in --profile) PROFILE=$2;; --cpus) CPUS=$2;; --memory-gb) MEMORY=$2;; esac
      shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
case "$PROFILE" in docker|apptainer|singularity) ;; *) echo 'Unsupported container profile' >&2; exit 2;; esac
command -v python3 >/dev/null || { echo 'Python 3.10+ is required' >&2; exit 2; }
RUN=$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).resolve())' "$RUN")
PREP_EXTRA=()
if (( RESUME )); then PREP_EXTRA+=(--resume); fi
python3 "$ROOT/prepare.py" "$ASSAY" "$SAMPLES" "$PARAMS" "$RUN" --cpus "$CPUS" --memory-gb "$MEMORY" --profile "$PROFILE" "${PREP_EXTRA[@]}"
IFS=$'\t' read -r PIPELINE REVISION < "$RUN/provenance/workflow.tsv"
ATTEMPT="$RUN/provenance/executions/$(date -u +%Y%m%dT%H%M%S)_$$"
mkdir -p "$ATTEMPT"
CMD=(nextflow run "nf-core/$PIPELINE" -r "$REVISION" -profile "$PROFILE"
  -params-file "$RUN/provenance/params.resolved.json" -c "$RUN/provenance/resources.config"
  -work-dir "$RUN/work" -with-report "$ATTEMPT/execution_report.html"
  -with-trace "$ATTEMPT/trace.tsv" -with-timeline "$ATTEMPT/timeline.html")
if (( RESUME )); then CMD+=(-resume); fi
printf '%q ' "${CMD[@]}" | tee "$RUN/provenance/command.sh"; printf '\n'
cp -- "$RUN/provenance/command.sh" "$ATTEMPT/command.sh"
if (( DRY )); then echo 'Validated plan only; no Nextflow tasks launched.'; exit 0; fi
command -v nextflow >/dev/null || { echo 'Install Nextflow and Java; see README.txt' >&2; exit 2; }
command -v "$PROFILE" >/dev/null || { echo "Container runtime missing: $PROFILE" >&2; exit 2; }
if [[ $PROFILE == docker ]]; then docker info >/dev/null; fi
nextflow -version > "$ATTEMPT/nextflow_version.txt"
python3 - "$ATTEMPT/nextflow_version.txt" <<'PY'
import re,sys
s=open(sys.argv[1]).read(); m=re.search(r'version\s+(\d+)\.(\d+)\.(\d+)',s)
if not m or tuple(map(int,m.groups())) != (25,10,4):
    raise SystemExit('This toolkit requires Nextflow 25.10.4 (RNA release minimum); 26.10.4 failed pinned-workflow config resolution in local testing')
PY
# Stable launch directory preserves the Nextflow cache required by -resume.
cd -- "$RUN"
"${CMD[@]}" 2>&1 | tee -a "$RUN/provenance/launcher.log" "$ATTEMPT/launcher.log"
