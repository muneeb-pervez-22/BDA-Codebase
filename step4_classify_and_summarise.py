#!/usr/bin/env python3
"""
step4_classify_and_summarise.py
────────────────────────────────
Aggregates the per-accession methylation TSVs (from step3) across all
accessions, then classifies each lncRNA into one of three epigenetic states
per accession, mirroring the logic used by Shahzad et al. for coding genes:

    UM   — Unmethylated: no significant CG, CHG, or CHH methylation
    gbM  — Gene-Body-Methylation-like: CG methylation present, CHG/CHH absent
    teM  — TE-like methylation: non-CG (CHG or CHH) methylation present

It also computes:
  • Per-lncRNA population conservation (% accessions in each state)
  • Context-level summary statistics across all accessions and lncRNAs
  • A comparison table analogous to Fig. 1 of the reference paper

Usage:
    python step4_classify_and_summarise.py \\
        --tsv_dir  methylation_per_accession/ \\
        --out_dir  results/

Outputs:
    results/lncRNA_per_accession_states.tsv  ─ gene_id × accession state matrix
    results/lncRNA_conservation.tsv          ─ per-lncRNA UM/gbM/teM frequencies
    results/context_summary.tsv             ─ population-level context fractions
"""

import argparse
import os
from pathlib import Path

import pandas as pd
import numpy as np


# ─── Classification thresholds (match Shahzad et al. / Choi et al. 2020) ────
# A region is considered methylated in a context if the methylation fraction
# exceeds the threshold AND there are at least MIN_SITES covered sites.

CG_THRESHOLD  = 0.10   # ≥10% CG methylation → CG-methylated
NON_CG_THRESHOLD = 0.05  # ≥5% CHG or CHH → non-CG methylated
MIN_SITES     = 3       # minimum covered cytosines in that context


def classify_state(row: pd.Series) -> str:
    """
    Classify one lncRNA in one accession into UM / gbM / teM.

    Priority order (matches the paper's segmentation logic):
      1. If CHG or CHH fraction exceeds threshold → teM
         (non-CG methylation is the hallmark of TE silencing)
      2. Elif CG fraction exceeds threshold → gbM-like
      3. Else → UM
    """
    chg_ok  = (row["n_CHG"] >= MIN_SITES) and (row["frac_CHG"] >= NON_CG_THRESHOLD)
    chh_ok  = (row["n_CHH"] >= MIN_SITES) and (row["frac_CHH"] >= NON_CG_THRESHOLD)
    cg_ok   = (row["n_CG"]  >= MIN_SITES) and (row["frac_CG"]  >= CG_THRESHOLD)

    if chg_ok or chh_ok:
        return "teM"
    if cg_ok:
        return "gbM"
    return "UM"


def load_all_accessions(tsv_dir: str) -> dict[str, pd.DataFrame]:
    """
    Load every *_lncRNA_methylation.tsv produced by step3.
    Returns a dict: accession_name → DataFrame.
    """
    files = sorted(Path(tsv_dir).glob("*_lncRNA_methylation.tsv"))
    if not files:
        raise FileNotFoundError(f"No TSV files found in {tsv_dir}. "
                                "Run step3 first.")
    print(f"Found {len(files)} accession files.")
    data = {}
    for f in files:
        accession = f.stem.replace("_lncRNA_methylation", "")
        df = pd.read_csv(f, sep="\t")
        data[accession] = df
    return data


def build_state_matrix(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Build a gene_id × accession matrix of epigenetic states (UM/gbM/teM/NA).
    NA = not enough covered sites in this accession to make a call.
    """
    all_genes = sorted(
        set(gene for df in data.values() for gene in df["gene_id"])
    )
    print(f"Total lncRNA loci across all accessions: {len(all_genes):,}")

    state_rows = {}
    for accession, df in data.items():
        # classify each row
        df = df.copy()
        df["state"] = df.apply(classify_state, axis=1)
        state_series = df.set_index("gene_id")["state"]
        state_rows[accession] = state_series

    matrix = pd.DataFrame(state_rows, index=all_genes)
    matrix.index.name = "gene_id"
    # Any gene not observed in a given accession stays NaN → "NA"
    matrix = matrix.fillna("NA")
    return matrix


def compute_conservation(matrix: pd.DataFrame) -> pd.DataFrame:
    """
    For each lncRNA, compute:
      - n_called: accessions with a non-NA call
      - pct_UM, pct_gbM, pct_teM: % of called accessions in each state
    """
    records = []
    for gene_id, row in matrix.iterrows():
        called = row[row != "NA"]
        n = len(called)
        if n == 0:
            records.append({"gene_id": gene_id, "n_called": 0,
                             "pct_UM": np.nan, "pct_gbM": np.nan,
                             "pct_teM": np.nan})
            continue
        vc = called.value_counts()
        records.append({
            "gene_id":  gene_id,
            "n_called": n,
            "pct_UM":   vc.get("UM",  0) / n * 100,
            "pct_gbM":  vc.get("gbM", 0) / n * 100,
            "pct_teM":  vc.get("teM", 0) / n * 100,
        })
    return pd.DataFrame(records).sort_values("pct_gbM", ascending=False)


def compute_context_summary(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Compute population-wide mean methylation fraction per context (CG/CHG/CHH)
    across all lncRNAs and all accessions.
    This is the lncRNA equivalent of the teM/gbM/SNP variance partitioning
    shown in Fig. 2 of the reference paper — here we start simpler:
    just mean methylation fraction.
    """
    all_fracs = {"CG": [], "CHG": [], "CHH": []}
    for df in data.values():
        for ctx in ["CG", "CHG", "CHH"]:
            col = f"frac_{ctx}"
            vals = df[col].dropna().values
            all_fracs[ctx].extend(vals.tolist())

    rows = []
    for ctx, vals in all_fracs.items():
        arr = np.array(vals)
        rows.append({
            "context":     ctx,
            "n_obs":       len(arr),
            "mean_frac":   np.nanmean(arr),
            "median_frac": np.nanmedian(arr),
            "std_frac":    np.nanstd(arr),
            "pct_>0":      (arr > 0).mean() * 100,
            "pct_>0.1":    (arr > 0.1).mean() * 100,
        })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--tsv_dir",  required=True,
                        help="Directory with per-accession TSVs from step3")
    parser.add_argument("--out_dir",  default="results",
                        help="Output directory")
    args = parser.parse_args()

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    print("[1/4] Loading per-accession methylation files ...")
    data = load_all_accessions(args.tsv_dir)

    print("[2/4] Building lncRNA × accession state matrix ...")
    matrix = build_state_matrix(data)
    matrix_path = os.path.join(args.out_dir, "lncRNA_per_accession_states.tsv")
    matrix.to_csv(matrix_path, sep="\t")
    print(f"      → {matrix_path}")

    print("[3/4] Computing population conservation per lncRNA ...")
    conservation = compute_conservation(matrix)
    cons_path = os.path.join(args.out_dir, "lncRNA_conservation.tsv")
    conservation.to_csv(cons_path, sep="\t", index=False)
    print(f"      → {cons_path}")

    # Quick summary to screen
    print("\n  ── Population-level state summary ──────────────────")
    print(f"  lncRNAs predominantly UM   (>90% UM):  "
          f"{(conservation['pct_UM'] >= 90).sum():>5,}")
    print(f"  lncRNAs with gbM-like     (>50% gbM):  "
          f"{(conservation['pct_gbM'] >= 50).sum():>5,}")
    print(f"  lncRNAs with teM-like     (>50% teM):  "
          f"{(conservation['pct_teM'] >= 50).sum():>5,}")
    print("  ─────────────────────────────────────────────────────\n")

    print("[4/4] Computing population-wide context fraction summary ...")
    ctx_summary = compute_context_summary(data)
    ctx_path = os.path.join(args.out_dir, "context_summary.tsv")
    ctx_summary.to_csv(ctx_path, sep="\t", index=False)
    print(f"      → {ctx_path}\n")

    print("  ── Methylation context summary (mean fraction) ──────")
    print(ctx_summary.to_string(index=False, float_format="{:.4f}".format))
    print("  ─────────────────────────────────────────────────────")
    print("\nDone. Key question to answer from results/:")
    print("  Do lncRNA regions show CG-only (gbM-like) or non-CG (teM-like)")
    print("  methylation, and how does their frequency compare to coding genes")
    print("  (55% gbM / 12% teM / 33% UM in the reference paper)?")


if __name__ == "__main__":
    main()
