#!/usr/bin/env bash
# ============================================================
# run_pipeline.sh – Execução do pipeline
# Uso: bash run_pipeline.sh [local|slurm|test] [--resume] [args extras]
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${1:-}" != -* ]]; then
    PROFILE="${1:-local}"
    shift 2>/dev/null || true
else
    PROFILE="local"
fi

# ── Detecta conda/mamba ──────────────────────────────────────
CONDA_CMD=""
for cmd in mamba conda; do
    command -v "$cmd" &>/dev/null && { CONDA_CMD="$cmd"; break; }
done
[ -z "$CONDA_CMD" ] && { echo "ERRO: mamba/conda não encontrado."; exit 1; }

echo "═══════════════════════════════════════════════════════"
echo "  RNASeq Insight – Arabidopsis thaliana"
echo "  Profile: $PROFILE"
echo "  $(date)"
echo "═══════════════════════════════════════════════════════"

LOG_DIR="$SCRIPT_DIR/logs"
mkdir -p "$LOG_DIR"
LOGFILE="$LOG_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"

$CONDA_CMD run -n rnaseq-tools \
  nextflow run "$SCRIPT_DIR/main.nf" \
    -profile         "$PROFILE" \
    -params-file     "$SCRIPT_DIR/params.yaml" \
    -with-report     "$SCRIPT_DIR/results/nextflow_report.html" \
    -with-trace      "$SCRIPT_DIR/results/nextflow_trace.txt" \
    -with-timeline   "$SCRIPT_DIR/results/nextflow_timeline.html" \
    -with-dag        "$SCRIPT_DIR/results/nextflow_dag.html" \
    "$@" 2>&1 | tee "$LOGFILE"

echo ""
echo "═══════════════════════════════════════════════════════"
echo "  Pipeline concluído!"
echo "  Log: $LOGFILE"
echo "  Resultados: $SCRIPT_DIR/results/"
echo ""
echo "  Dashboard interativo:"
echo "  RESULTS_DIR=results Rscript -e \"shiny::runApp('dashboard/app.R', port=3838)\""
echo "═══════════════════════════════════════════════════════"
