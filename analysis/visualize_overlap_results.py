#!/usr/bin/env python3
"""
visualize_overlap_results.py
────────────────────────────
Analytical and plotting script for lncRNA × transposable-element (TE) overlap
results, optionally joined with the full methylation conservation output.

This script is meant to be run after:
  1. step1_prepare_beds.py
  2. step2_te_overlap.sh
  3. optionally run_full_dataset.sh / step4_classify_and_summarise.py

It answers questions such as:
  • How many lncRNAs overlap annotated TEs?
  • How long are TE-overlapping vs TE-free lncRNAs?
  • Are TE-overlapping lncRNAs enriched for teM-like methylation states?
  • Are TE-free lncRNAs mostly UM or do they contain gbM-like candidates?
  • Which TE metadata values are most common among overlapping lncRNAs?

Example:
    python analysis/visualize_overlap_results.py \
        --lncrna-bed lncRNAs_only.bed \
        --te-overlap-bed overlap_results/lncRNAs_TE_overlapping.bed \
        --te-free-bed overlap_results/lncRNAs_TE_free.bed \
        --overlap-meta-bed overlap_results/lncRNAs_TE_overlap_with_metadata.bed \
        --conservation-tsv results_full/lncRNA_conservation.tsv \
        --out-dir analysis_overlap_results

Outputs:
    analysis_overlap_results/
      overlap_counts.tsv
      lncRNA_overlap_status.tsv
      methylation_by_te_status.tsv                      [if conservation supplied]
      dominant_methylation_state_by_te_status.tsv       [if conservation supplied]
      top_te_metadata_values.tsv                        [if overlap metadata supplied]
      *.png figures
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


BED6_COLS = ["chr", "start", "end", "gene_id", "score", "strand"]


def read_bed6(path: str | Path) -> pd.DataFrame:
    """Read a headerless BED6 file into a DataFrame."""
    df = pd.read_csv(path, sep="\t", header=None, names=BED6_COLS, dtype={"chr": str})
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["length"] = df["end"] - df["start"]
    return df


def make_overlap_status(
    lncrna_bed: str | Path,
    te_overlap_bed: str | Path,
    te_free_bed: Optional[str | Path] = None,
) -> pd.DataFrame:
    """Create one row per lncRNA with TE-overlap status."""
    all_lnc = read_bed6(lncrna_bed)
    overlap_ids = set(read_bed6(te_overlap_bed)["gene_id"])

    all_lnc["te_status"] = np.where(
        all_lnc["gene_id"].isin(overlap_ids),
        "TE-overlapping",
        "TE-free",
    )

    # Optional sanity check against the anti-join file from step2_te_overlap.sh.
    if te_free_bed is not None and Path(te_free_bed).exists():
        free_ids = set(read_bed6(te_free_bed)["gene_id"])
        inferred_free = set(all_lnc.loc[all_lnc["te_status"] == "TE-free", "gene_id"])
        mismatch = inferred_free.symmetric_difference(free_ids)
        if mismatch:
            print(
                f"WARNING: TE-free BED disagrees with inferred TE-free set for "
                f"{len(mismatch)} lncRNA IDs. Proceeding with status inferred "
                f"from --te-overlap-bed."
            )

    return all_lnc


def classify_dominant_methylation_state(row: pd.Series) -> str:
    """
    Collapse pct_UM / pct_gbM / pct_teM into a human-readable dominant state.

    These rules mirror the current project narrative:
      • >90% UM = predominantly unmethylated
      • >50% gbM = predominantly gbM-like
      • >50% teM = predominantly teM-like
      • otherwise = mixed/intermediate
    """
    if row["pct_UM"] >= 90:
        return "predominantly UM (>90%)"
    if row["pct_gbM"] >= 50:
        return "predominantly gbM-like (>50%)"
    if row["pct_teM"] >= 50:
        return "predominantly teM-like (>50%)"
    return "mixed/intermediate"


def load_conservation(path: str | Path, min_called: Optional[int] = None) -> pd.DataFrame:
    """Load results_full/lncRNA_conservation.tsv and optionally filter low-call loci."""
    cons = pd.read_csv(path, sep="\t")
    required = {"gene_id", "n_called", "pct_UM", "pct_gbM", "pct_teM"}
    missing = required - set(cons.columns)
    if missing:
        raise ValueError(f"Conservation file is missing required columns: {sorted(missing)}")

    if min_called is not None:
        before = len(cons)
        cons = cons[cons["n_called"] >= min_called].copy()
        print(f"Filtered conservation table by n_called >= {min_called}: {len(cons)}/{before} loci retained")

    cons["dominant_methylation_state"] = cons.apply(classify_dominant_methylation_state, axis=1)
    return cons


def save_overlap_counts(status_df: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    counts = (
        status_df["te_status"]
        .value_counts()
        .rename_axis("te_status")
        .reset_index(name="n_lncRNAs")
    )
    counts["percent"] = counts["n_lncRNAs"] / counts["n_lncRNAs"].sum() * 100
    counts.to_csv(out_dir / "overlap_counts.tsv", sep="\t", index=False)
    return counts


def plot_overlap_counts(counts: pd.DataFrame, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(counts["te_status"], counts["n_lncRNAs"])
    ax.set_ylabel("Number of lncRNAs")
    ax.set_title("lncRNA overlap with annotated TEs")
    for i, row in counts.iterrows():
        ax.text(i, row["n_lncRNAs"], f"{row['n_lncRNAs']:,}\n({row['percent']:.1f}%)", ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(out_dir / "overlap_counts_bar.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.pie(counts["n_lncRNAs"], labels=counts["te_status"], autopct="%1.1f%%")
    ax.set_title("TE-overlapping vs TE-free lncRNAs")
    fig.tight_layout()
    fig.savefig(out_dir / "overlap_counts_pie.png", dpi=200)
    plt.close(fig)


def plot_length_by_status(status_df: pd.DataFrame, out_dir: Path) -> None:
    groups = [grp["length"].values for _, grp in status_df.groupby("te_status")]
    labels = [name for name, _ in status_df.groupby("te_status")]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.boxplot(groups, labels=labels, showfliers=False)
    ax.set_ylabel("lncRNA length (bp)")
    ax.set_title("lncRNA length by TE-overlap status")
    fig.tight_layout()
    fig.savefig(out_dir / "lncRNA_length_by_te_status_boxplot.png", dpi=200)
    plt.close(fig)


def join_conservation(status_df: pd.DataFrame, cons: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    joined = status_df.merge(cons, on="gene_id", how="left")
    joined.to_csv(out_dir / "methylation_by_te_status.tsv", sep="\t", index=False)

    state_counts = (
        joined.dropna(subset=["dominant_methylation_state"])
        .groupby(["te_status", "dominant_methylation_state"])
        .size()
        .reset_index(name="n_lncRNAs")
    )
    total_by_status = state_counts.groupby("te_status")["n_lncRNAs"].transform("sum")
    state_counts["percent_within_te_status"] = state_counts["n_lncRNAs"] / total_by_status * 100
    state_counts.to_csv(out_dir / "dominant_methylation_state_by_te_status.tsv", sep="\t", index=False)
    return joined


def plot_state_by_status(joined: pd.DataFrame, out_dir: Path) -> None:
    plot_df = joined.dropna(subset=["dominant_methylation_state"]).copy()
    if plot_df.empty:
        print("No methylation conservation rows available for plotting.")
        return

    counts = pd.crosstab(plot_df["te_status"], plot_df["dominant_methylation_state"])
    fractions = counts.div(counts.sum(axis=1), axis=0) * 100

    fig, ax = plt.subplots(figsize=(9, 5))
    bottom = np.zeros(len(fractions))
    x = np.arange(len(fractions.index))
    for col in fractions.columns:
        vals = fractions[col].values
        ax.bar(x, vals, bottom=bottom, label=col)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(fractions.index, rotation=0)
    ax.set_ylabel("Percent of lncRNAs within TE-status group")
    ax.set_title("Dominant methylation state by TE-overlap status")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    fig.savefig(out_dir / "dominant_methylation_state_by_te_status_stacked.png", dpi=200)
    plt.close(fig)

    for col in ["pct_UM", "pct_gbM", "pct_teM"]:
        groups = [grp[col].dropna().values for _, grp in plot_df.groupby("te_status")]
        labels = [name for name, _ in plot_df.groupby("te_status")]
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.boxplot(groups, labels=labels, showfliers=False)
        ax.set_ylabel(f"{col} across called accessions")
        ax.set_title(f"{col} by TE-overlap status")
        fig.tight_layout()
        fig.savefig(out_dir / f"{col}_by_te_status_boxplot.png", dpi=200)
        plt.close(fig)


def infer_metadata_column(meta: pd.DataFrame, start_col: int = 6) -> int:
    """
    Heuristically choose a TE metadata column to plot.

    The overlap metadata file is produced by `bedtools intersect -wa -wb`, so the
    first six columns are lncRNA BED6 and the remaining columns are TE metadata.
    Because the Panda TE BED can be used with or without a header, this function
    picks a non-numeric TE-side column with a moderate number of unique values.
    """
    best_col = start_col
    best_score = -np.inf
    for col in meta.columns[start_col:]:
        s = meta[col].dropna().astype(str)
        if s.empty:
            continue
        numeric_fraction = pd.to_numeric(s, errors="coerce").notna().mean()
        n_unique = s.nunique()
        if numeric_fraction > 0.8:
            continue
        if n_unique < 2:
            continue
        # Prefer categorical-ish columns, not unique IDs.
        unique_ratio = n_unique / len(s)
        score = -abs(unique_ratio - 0.05) - (0.01 * max(n_unique - 200, 0))
        if score > best_score:
            best_score = score
            best_col = col
    return int(best_col)


def analyze_overlap_metadata(
    overlap_meta_bed: str | Path,
    out_dir: Path,
    te_metadata_col: Optional[int] = None,
    top_n: int = 20,
) -> None:
    """Summarize and plot a selected TE metadata column from the full overlap file."""
    meta = pd.read_csv(overlap_meta_bed, sep="\t", header=None, dtype=str)
    if meta.empty:
        print("Overlap metadata file is empty; skipping TE metadata plot.")
        return

    if te_metadata_col is None:
        te_metadata_col = infer_metadata_column(meta, start_col=6)
        print(f"Inferred TE metadata column for plotting: column {te_metadata_col} (0-indexed in full overlap file)")

    if te_metadata_col not in meta.columns:
        raise ValueError(f"Requested --te-metadata-col {te_metadata_col}, but file has columns 0..{meta.columns.max()}")

    counts = (
        meta[te_metadata_col]
        .fillna("NA")
        .astype(str)
        .value_counts()
        .head(top_n)
        .rename_axis("metadata_value")
        .reset_index(name="n_overlaps")
    )
    counts.insert(0, "metadata_column_0_indexed", te_metadata_col)
    counts.to_csv(out_dir / "top_te_metadata_values.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(9, max(4, 0.3 * len(counts))))
    y = np.arange(len(counts))
    ax.barh(y, counts["n_overlaps"])
    ax.set_yticks(y)
    ax.set_yticklabels(counts["metadata_value"])
    ax.invert_yaxis()
    ax.set_xlabel("Number of lncRNA–TE overlap records")
    ax.set_title(f"Top TE metadata values from column {te_metadata_col}")
    fig.tight_layout()
    fig.savefig(out_dir / "top_te_metadata_values.png", dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--lncrna-bed", default="lncRNAs_only.bed", help="All lncRNAs BED6 from step1")
    parser.add_argument("--te-overlap-bed", default="overlap_results/lncRNAs_TE_overlapping.bed", help="TE-overlapping lncRNAs BED6 from step2")
    parser.add_argument("--te-free-bed", default="overlap_results/lncRNAs_TE_free.bed", help="TE-free lncRNAs BED6 from step2")
    parser.add_argument("--overlap-meta-bed", default="overlap_results/lncRNAs_TE_overlap_with_metadata.bed", help="Full lncRNA × TE overlap file from step2")
    parser.add_argument("--conservation-tsv", default="results_full/lncRNA_conservation.tsv", help="Optional methylation conservation table from step4")
    parser.add_argument("--min-called", type=int, default=None, help="Optional minimum n_called filter for conservation table, e.g. 649 for 70% of 927")
    parser.add_argument("--te-metadata-col", type=int, default=None, help="0-indexed column in full overlap metadata file to summarize; inferred if omitted")
    parser.add_argument("--top-n", type=int, default=20, help="Number of TE metadata values to plot")
    parser.add_argument("--out-dir", default="analysis_overlap_results", help="Directory for tables and figures")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    status_df = make_overlap_status(args.lncrna_bed, args.te_overlap_bed, args.te_free_bed)
    status_df.to_csv(out_dir / "lncRNA_overlap_status.tsv", sep="\t", index=False)

    counts = save_overlap_counts(status_df, out_dir)
    plot_overlap_counts(counts, out_dir)
    plot_length_by_status(status_df, out_dir)

    if args.conservation_tsv and Path(args.conservation_tsv).exists():
        cons = load_conservation(args.conservation_tsv, min_called=args.min_called)
        joined = join_conservation(status_df, cons, out_dir)
        plot_state_by_status(joined, out_dir)
    else:
        print("No conservation TSV found/provided; skipping methylation-by-TE-status plots.")

    if args.overlap_meta_bed and Path(args.overlap_meta_bed).exists():
        analyze_overlap_metadata(args.overlap_meta_bed, out_dir, args.te_metadata_col, args.top_n)
    else:
        print("No overlap metadata BED found/provided; skipping TE metadata plots.")

    print(f"Done. Wrote overlap analysis tables and figures to: {out_dir}")


if __name__ == "__main__":
    main()
