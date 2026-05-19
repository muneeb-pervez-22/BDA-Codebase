#!/usr/bin/env bash
# =============================================================================
# run_full_dataset.sh
# =============================================================================
# Runs step3 on every accession in GSE43857_RAW/ (927 files), skipping any
# that already have output in methylation_full_run/.
# Then runs step4 to aggregate all results into results_full/.
#
# Usage:
#   chmod +x run_full_dataset.sh
#   ./run_full_dataset.sh
#
# Optional: pass -j N to control parallel jobs (default: number of CPU cores)
#   ./run_full_dataset.sh -j 4
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="${SCRIPT_DIR}/run_full_dataset.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== Run started: $(date) ==="
ALLC_DIR="${SCRIPT_DIR}/GSE43857_RAW"
REGIONS="${SCRIPT_DIR}/lncRNAs_only.bed"
OUT_DIR="${SCRIPT_DIR}/methylation_full_run"   # separate dir — avoids mixing with old allc_files results
RESULTS_DIR="${SCRIPT_DIR}/results_full"
JOBS="${1:-$(nproc)}"

if [[ "${1:-}" == "-j" ]]; then
    JOBS="${2:?-j requires a number}"
fi

echo "============================================="
echo "  Full-dataset methylation pipeline"
echo "============================================="
echo "  Input dir  : $ALLC_DIR"
echo "  Regions    : $REGIONS"
echo "  Out (step3): $OUT_DIR"
echo "  Out (step4): $RESULTS_DIR"
echo "  Parallel   : $JOBS jobs"
echo ""

mkdir -p "$OUT_DIR" "$RESULTS_DIR"

# ── Build list of files that still need processing ───────────────────────────
PENDING=()
for f in "$ALLC_DIR"/*.tsv; do
    [[ -f "$f" ]] || continue
    accession="$(basename "$f" .tsv)"
    expected="${OUT_DIR}/${accession}_lncRNA_methylation.tsv"
    if [[ ! -f "$expected" ]]; then
        PENDING+=("$f")
    fi
done

TOTAL=$(ls "$ALLC_DIR"/*.tsv 2>/dev/null | wc -l)
ALREADY_DONE=$(( TOTAL - ${#PENDING[@]} ))

echo "[step3] $TOTAL total accessions in $ALLC_DIR"
echo "        $ALREADY_DONE already processed, ${#PENDING[@]} remaining"
echo ""

if [[ ${#PENDING[@]} -eq 0 ]]; then
    echo "        Nothing to do — all accessions already processed."
else
    echo "[step3] Processing ${#PENDING[@]} accessions with $JOBS parallel jobs ..."
    printf '%s\n' "${PENDING[@]}" | parallel -j "$JOBS" \
        python3 "${SCRIPT_DIR}/step3_methylation_context.py" \
            --allc {} \
            --regions "$REGIONS" \
            --out_dir "$OUT_DIR"
    echo ""
    echo "[step3] Done."
fi

# ── Count all outputs before running step4 ───────────────────────────────────
N_OUT=$(ls "$OUT_DIR"/*_lncRNA_methylation.tsv 2>/dev/null | wc -l)
echo ""
echo "[step4] Aggregating $N_OUT accession files → $RESULTS_DIR ..."
python3 "${SCRIPT_DIR}/step4_classify_and_summarise.py" \
    --tsv_dir "$OUT_DIR" \
    --out_dir "$RESULTS_DIR"

echo ""
echo "All done. Results in $RESULTS_DIR/"
