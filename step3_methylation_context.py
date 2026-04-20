#!/usr/bin/env python3
"""
step3_methylation_context.py
─────────────────────────────
For ONE bisulfite sequencing accession (one .allc or tab-separated
cytosine-report file from GSE43857), compute per-lncRNA methylation
fractions in each sequence context (CG, CHG, CHH).

This script is designed to be run in a simple loop (or via GNU parallel /
a cluster job array) over all 927 accession files:

    # Serial loop example:
    for f in allc_files/*.allc; do
        python step3_methylation_context.py \\
            --allc    "$f" \\
            --regions lncRNAs_only.bed \\        # or TE_free / TE_overlapping
            --out_dir methylation_per_accession/
    done

    # GNU parallel (much faster):
    ls allc_files/*.allc | parallel -j 8 \\
        python step3_methylation_context.py \\
            --allc {} \\
            --regions lncRNAs_only.bed \\
            --out_dir methylation_per_accession/

──────────────────────────────────────────────────────────────────────────────
Input allc format (GSE43857 / 1001 methylomes standard):
    chr  pos  strand  mc_class  mc_count  total  methylated
    Chr1  1    +       CGT       5         10     1
    ...
    mc_class is the trinucleotide context; we collapse to:
        CG   → mc_class starts with "CG"
        CHG  → mc_class starts with "C" and index-2 is G (e.g. CAG, CTG, CCG)
        CHH  → mc_class starts with "C" and index-2 is not G (e.g. CAA, CTT)

Output (per accession):
    A tab-separated file with one row per lncRNA:
    gene_id  chr  start  end  strand
    n_CG_sites  CG_mc  CG_total  CG_fraction
    n_CHG_sites CHG_mc CHG_total CHG_fraction
    n_CHH_sites CHH_mc CHH_total CHH_fraction
──────────────────────────────────────────────────────────────────────────────
"""

import argparse
import sys
import os
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np


# ── Context classification ───────────────────────────────────────────────────

def classify_context(mc_class: str) -> str:
    """
    Map a trinucleotide mc_class to CG / CHG / CHH.
    CG:  mc_class[0:2] == 'CG'
    CHG: mc_class[0] == 'C', mc_class[2] == 'G' (and not CG)
    CHH: mc_class[0] == 'C', mc_class[2] != 'G'
    Returns '' for anything that doesn't start with C (shouldn't happen).
    """
    mc = mc_class.upper()
    if not mc.startswith("C") or len(mc) < 3:
        return ""
    if mc[1] == "G":
        return "CG"
    if mc[2] == "G":
        return "CHG"
    return "CHH"


# ── Load BED regions ─────────────────────────────────────────────────────────

def load_regions(bed_path: str) -> pd.DataFrame:
    """Load a headerless BED6 file (output of step1/step2)."""
    df = pd.read_csv(
        bed_path, sep="\t", header=None,
        names=["chr", "start", "end", "gene_id", "score", "strand"],
        dtype={"chr": str, "start": int, "end": int,
               "gene_id": str, "score": str, "strand": str},
    )
    return df.sort_values(["chr", "start"]).reset_index(drop=True)


# ── Load allc file ────────────────────────────────────────────────────────────

def load_allc(allc_path: str) -> pd.DataFrame:
    """
    Load a bisulfite allc file.
    Handles both plain-text and gzip-compressed (.allc.gz) files.
    """
    allc_col_names = ["chr", "pos", "strand", "mc_class",
                      "mc_count", "total", "methylated"]
    kwargs = dict(
        sep="\t",
        header=0,  # Change this from None to 0 to skip the header row
        names=allc_col_names,
        dtype={"chr": str, "pos": int, "strand": str,
               "mc_class": str, "mc_count": int,
               "total": int, "methylated": int},
    )
    if allc_path.endswith(".gz"):
        return pd.read_csv(allc_path, compression="gzip", **kwargs)
    return pd.read_csv(allc_path, **kwargs)


# ── Per-region methylation aggregation ──────────────────────────────────────

def aggregate_methylation(regions: pd.DataFrame,
                          allc: pd.DataFrame) -> pd.DataFrame:
    """
    For each lncRNA region, sum up mc_count and total across all cytosines
    that fall within [start, end) for each context (CG, CHG, CHH).

    Uses a chromsome-by-chromosome approach with binary search to stay
    memory-efficient on the large allc files.

    Returns a DataFrame with one row per region.
    """
    # Pre-classify contexts in allc
    allc = allc.copy()
    allc["context"] = allc["mc_class"].map(classify_context)
    allc = allc[allc["context"] != ""]   # drop any malformed rows

    # Group allc by chromosome for fast lookup
    allc_by_chr = {chrom: grp for chrom, grp in allc.groupby("chr")}

    records = []
    contexts = ["CG", "CHG", "CHH"]

    for _, region in regions.iterrows():
        # Strip "Chr" prefix if it exists to match the TSV files
        chrom = region["chr"].replace("Chr", "") 
        start = region["start"]
        end   = region["end"]

        row = {
            "gene_id": region["gene_id"],
            "chr":     chrom,
            "start":   start,
            "end":     end,
            "strand":  region["strand"],
        }

        if chrom not in allc_by_chr:
            # No cytosine data on this chromosome for this accession
            for ctx in contexts:
                row[f"n_{ctx}"]       = 0
                row[f"mc_{ctx}"]      = 0
                row[f"total_{ctx}"]   = 0
                row[f"frac_{ctx}"]    = np.nan
        else:
            chrom_allc = allc_by_chr[chrom]
            # Boolean slice — efficient because allc is sorted
            in_region = chrom_allc[
                (chrom_allc["pos"] >= start) &
                (chrom_allc["pos"] <  end)
            ]

            for ctx in contexts:
                ctx_data = in_region[in_region["context"] == ctx]
                mc    = ctx_data["mc_count"].sum()
                total = ctx_data["total"].sum()
                row[f"n_{ctx}"]     = len(ctx_data)
                row[f"mc_{ctx}"]    = int(mc)
                row[f"total_{ctx}"] = int(total)
                row[f"frac_{ctx}"]  = (mc / total) if total > 0 else np.nan

        records.append(row)

    return pd.DataFrame(records)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--allc",    required=True,
                        help="Path to the bisulfite allc file for one accession "
                             "(plain or .gz)")
    parser.add_argument("--regions", required=True,
                        help="BED6 file of lncRNA regions (no header, "
                             "from step1/step2)")
    parser.add_argument("--out_dir", default="methylation_per_accession",
                        help="Directory to write per-accession output TSV")
    parser.add_argument("--min_sites", type=int, default=3,
                        help="Minimum number of covered cytosine sites required "
                             "to keep a region (default: 3, matching the paper)")
    args = parser.parse_args()

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    # Derive accession name from the allc filename (e.g. Col-0.allc → Col-0)
    accession = Path(args.allc).stem.replace(".allc", "")
    out_path  = os.path.join(args.out_dir, f"{accession}_lncRNA_methylation.tsv")

    print(f"[{accession}] Loading regions from {args.regions} ...")
    regions = load_regions(args.regions)
    print(f"[{accession}]   {len(regions):,} lncRNA regions")

    print(f"[{accession}] Loading allc from {args.allc} ...")
    allc = load_allc(args.allc)
    print(f"[{accession}]   {len(allc):,} cytosine records")

    print(f"[{accession}] Aggregating methylation per region ...")
    result = aggregate_methylation(regions, allc)

    # Filter: drop regions with fewer than min_sites covered in ANY context
    # (mirrors the coverage filter in Shahzad et al.)
    before = len(result)
    result = result[
        (result["n_CG"]  >= args.min_sites) |
        (result["n_CHG"] >= args.min_sites) |
        (result["n_CHH"] >= args.min_sites)
    ].reset_index(drop=True)
    print(f"[{accession}]   Retained {len(result):,}/{before:,} regions "
          f"(≥{args.min_sites} covered cytosine sites in at least one context)")

    result.to_csv(out_path, sep="\t", index=False)
    print(f"[{accession}] Written → {out_path}")


if __name__ == "__main__":
    main()
