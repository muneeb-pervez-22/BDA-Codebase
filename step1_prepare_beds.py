#!/usr/bin/env python3
"""
step1_prepare_beds.py
─────────────────────
Separates lncRNAs (CUFF_NC) from protein-coding genes (CUFF_PC) in the
combined annotation BED file, and emits a clean, sorted 6-column BED file
ready for bedtools operations.

Usage:
    python step1_prepare_beds.py \
        --input  lncRNA_annotations.bed \
        --out_nc  lncRNAs_only.bed \
        --out_pc  protein_coding_only.bed

Output columns (BED6):
    chr  start  end  gene_id  score  strand
"""

import argparse
import pandas as pd
import sys
from pathlib import Path


BED6_COLS = ["chr", "start", "end", "gene_id", "score", "strand"]


def load_annotation_bed(path: str) -> pd.DataFrame:
    """
    Load the combined annotation BED.
    Handles the header line gracefully regardless of whether it starts with '#'.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=0,          # first non-comment line is the header
        dtype=str,
    )
    # Rename columns to standard BED6 names (order is fixed in this file)
    df.columns = BED6_COLS
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    return df


def filter_by_prefix(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Return rows whose gene_id starts with the given prefix, sorted by position."""
    mask = df["gene_id"].str.startswith(prefix)
    subset = df[mask].copy()
    subset = subset.sort_values(["chr", "start", "end"]).reset_index(drop=True)
    return subset


def write_bed(df: pd.DataFrame, path: str) -> None:
    """Write a BED6 file without a header (bedtools expects no header)."""
    df.to_csv(path, sep="\t", header=False, index=False)
    print(f"  → Wrote {len(df):,} records to {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",  required=True, help="Combined annotation BED (with header)")
    parser.add_argument("--out_nc", default="lncRNAs_only.bed",
                        help="Output BED for lncRNAs only (CUFF_NC)")
    parser.add_argument("--out_pc", default="protein_coding_only.bed",
                        help="Output BED for protein-coding genes (CUFF_PC)")
    args = parser.parse_args()

    print(f"[1/3] Loading {args.input} ...")
    df = load_annotation_bed(args.input)
    print(f"      Total entries: {len(df):,}")

    print(f"[2/3] Filtering by gene type ...")
    lnc = filter_by_prefix(df, "CUFF_NC")
    pc  = filter_by_prefix(df, "CUFF_PC")
    print(f"      lncRNAs (CUFF_NC): {len(lnc):,}")
    print(f"      Protein-coding (CUFF_PC): {len(pc):,}")

    # Sanity check: flag entries that fall into neither category
    unclassified = df[~df["gene_id"].str.startswith(("CUFF_NC", "CUFF_PC"))]
    if len(unclassified):
        print(f"  WARNING: {len(unclassified)} unclassified entries (neither NC nor PC):",
              file=sys.stderr)
        print(unclassified["gene_id"].unique()[:10], file=sys.stderr)

    print(f"[3/3] Writing output BED files ...")
    write_bed(lnc, args.out_nc)
    write_bed(pc,  args.out_pc)

    print("\nDone. Next step: run step2_te_overlap.sh")


if __name__ == "__main__":
    main()
