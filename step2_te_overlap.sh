#!/usr/bin/env bash
# =============================================================================
# step2_te_overlap.sh
# =============================================================================
# Finds which lncRNAs overlap with annotated Transposable Elements (TEs),
# and which do not. This is the key first split before methylation analysis:
#
#   lncRNAs
#   ├── TE-overlapping  →  likely to show teM-like (non-CG) methylation
#   └── TE-free         →  intergenic/antisense; our primary region of interest
#                          for novel gbM-like or UM classification
#
# Requires: bedtools >= 2.29
#
# Usage:
#   chmod +x step2_te_overlap.sh
#   ./step2_te_overlap.sh \
#       lncRNAs_only.bed \
#       Panda_AT-TEs_annotation_v1_0.bed \
#       ./overlap_results/
# =============================================================================

set -euo pipefail

# ── Arguments ──────────────────────────────────────────────────────────────
LNCRNA_BED="${1:-lncRNAs_only.bed}"
TE_BED="${2:-Panda_AT-TEs_annotation_v1_0.bed}"
OUTDIR="${3:-overlap_results}"

# ── Validate inputs ─────────────────────────────────────────────────────────
if [[ ! -f "$LNCRNA_BED" ]]; then
    echo "ERROR: lncRNA BED not found: $LNCRNA_BED" >&2; exit 1
fi
if [[ ! -f "$TE_BED" ]]; then
    echo "ERROR: TE BED not found: $TE_BED" >&2; exit 1
fi
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not in PATH" >&2; exit 1; }

mkdir -p "$OUTDIR"

echo "============================================="
echo "  lncRNA × TE Overlap Analysis"
echo "============================================="
echo "  lncRNA BED : $LNCRNA_BED"
echo "  TE BED     : $TE_BED"
echo "  Output dir : $OUTDIR"
echo ""

# ── Strip the header from the TE file (lines starting with '#') ─────────────
# bedtools requires no header. The Panda file uses '#' for the header.
TE_CLEAN="$OUTDIR/TE_no_header.bed"
grep -v "^#" "$TE_BED" > "$TE_CLEAN"
echo "[1/5] Cleaned TE file → $TE_CLEAN"

# ── Sort both BED files ──────────────────────────────────────────────────────
# bedtools intersect -sorted requires lexicographic sort on chr, then numeric on start/end
LNCRNA_SORTED="$OUTDIR/lncRNAs_sorted.bed"
TE_SORTED="$OUTDIR/TE_sorted.bed"

sort -k1,1 -k2,2n "$LNCRNA_BED" > "$LNCRNA_SORTED"
sort -k1,1 -k2,2n "$TE_CLEAN"   > "$TE_SORTED"
echo "[2/5] Sorted input files"

# ── (A) lncRNAs that DO overlap a TE ────────────────────────────────────────
# -u: report each lncRNA only once (even if it overlaps multiple TEs)
# -wa: write original lncRNA entry
# We also capture the TE metadata in a separate file for downstream analysis
TE_OVERLAP="$OUTDIR/lncRNAs_TE_overlapping.bed"
bedtools intersect \
    -a "$LNCRNA_SORTED" \
    -b "$TE_SORTED" \
    -wa -u \
    -sorted \
> "$TE_OVERLAP"
echo "[3/5] TE-overlapping lncRNAs → $TE_OVERLAP"

# Additionally: save the full overlap with TE metadata (which TE family etc.)
TE_OVERLAP_WITH_META="$OUTDIR/lncRNAs_TE_overlap_with_metadata.bed"
bedtools intersect \
    -a "$LNCRNA_SORTED" \
    -b "$TE_SORTED" \
    -wa -wb \
    -sorted \
> "$TE_OVERLAP_WITH_META"
echo "      (with TE metadata) → $TE_OVERLAP_WITH_META"

# ── (B) lncRNAs that do NOT overlap any TE ──────────────────────────────────
# -v: invert the match (anti-join)
TE_FREE="$OUTDIR/lncRNAs_TE_free.bed"
bedtools intersect \
    -a "$LNCRNA_SORTED" \
    -b "$TE_SORTED" \
    -v \
    -sorted \
> "$TE_FREE"
echo "[4/5] TE-free lncRNAs       → $TE_FREE"

# ── Summary statistics ───────────────────────────────────────────────────────
N_TOTAL=$(wc -l < "$LNCRNA_SORTED")
N_OVERLAP=$(wc -l < "$TE_OVERLAP")
N_FREE=$(wc -l < "$TE_FREE")

echo ""
echo "[5/5] Summary"
echo "─────────────────────────────────────────────"
printf "  Total lncRNAs:           %6d\n" "$N_TOTAL"
printf "  TE-overlapping lncRNAs:  %6d  (%.1f%%)\n" \
    "$N_OVERLAP" "$(echo "scale=1; $N_OVERLAP * 100 / $N_TOTAL" | bc)"
printf "  TE-free lncRNAs:         %6d  (%.1f%%)\n" \
    "$N_FREE"    "$(echo "scale=1; $N_FREE * 100 / $N_TOTAL" | bc)"
echo "─────────────────────────────────────────────"
echo ""
echo "Next step: run step3_methylation_context.py on each accession's"
echo "bisulfite file, using lncRNAs_TE_free.bed and/or lncRNAs_TE_overlapping.bed"
