# lncRNA Methylation Context Pipeline
## Arabidopsis thaliana — Extending Shahzad et al. 2025 to non-coding regions

---

### Background

Shahzad et al. (*Nature Plants*, 2025) showed that in *Arabidopsis* coding genes:
- **gbM** (CG-only methylation) occurs in ~55% of genes, promotes transcription
- **teM** (non-CG methylation — CHG + CHH) occurs in ~12% of genes, silences them
- **UM** (unmethylated) occurs in ~33% of genes

These three states are **largely independent** and reside in different gene types.

This pipeline asks: **Do lncRNA loci follow similar rules?**  
Specifically, the professor's first question is:
> *In lncRNA regions, which methylation contexts (CG, CHG, CHH) are present,
> at what frequencies, and do they resemble gbM-like or teM-like patterns?*

---

### Pipeline Overview

```
lncRNA_annotations.bed   Panda_AT-TEs_annotation_v1_0.bed   GSE43857 allc files
        │                          │                                 │
        ▼                          │                                 │
[ step1_prepare_beds.py ]          │                                 │
  Separate CUFF_NC / CUFF_PC       │                                 │
        │                          │                                 │
        ▼                          ▼                                 │
[ step2_te_overlap.sh ]                                              │
  bedtools intersect: lncRNAs ∩ TEs                                  │
  → lncRNAs_TE_overlapping.bed                                       │
  → lncRNAs_TE_free.bed                                              │
        │                                                            │
        └──────────────────────────────────────────────┐            │
                                                        ▼            ▼
                                              [ step3_methylation_context.py ]
                                                Run once per accession (×927)
                                                Extracts CG/CHG/CHH fractions
                                                per lncRNA per accession
                                                        │
                                                        ▼
                                              [ step4_classify_and_summarise.py ]
                                                Classify each lncRNA per accession:
                                                  UM / gbM-like / teM-like
                                                Compute population conservation
                                                Produce context summary tables
```

---

### Setup

```bash
# Create conda environment
conda create -n lncrna_meth python=3.11 pandas numpy bedtools -c conda-forge -c bioconda
conda activate lncrna_meth
```

---

### Step-by-step Usage

#### Step 1 — Separate lncRNAs from protein-coding genes

```bash
python step1_prepare_beds.py \
    --input  lncRNA_annotations.bed \
    --out_nc lncRNAs_only.bed \
    --out_pc protein_coding_only.bed
```

Output:
- `lncRNAs_only.bed` — 11,071 CUFF_NC entries, sorted, no header (BED6)
- `protein_coding_only.bed` — 23,676 CUFF_PC entries

---

#### Step 2 — Find lncRNAs overlapping TEs

```bash
chmod +x step2_te_overlap.sh
./step2_te_overlap.sh \
    lncRNAs_only.bed \
    Panda_AT-TEs_annotation_v1_0.bed \
    overlap_results/
```

Output (in `overlap_results/`):
- `lncRNAs_TE_overlapping.bed` — lncRNAs with ≥1 bp overlap with a TE
- `lncRNAs_TE_free.bed` — lncRNAs with no TE overlap (intergenic lncRNAs)
- `lncRNAs_TE_overlap_with_metadata.bed` — overlapping lncRNAs + TE family info

**Why this split matters:**
lncRNAs overlapping TEs are expected to carry teM (non-CG methylation) because
TE-silencing machinery deposits CHG/CHH methylation throughout the TE body.
TE-free lncRNAs are the interesting case: do they acquire CG-only (gbM-like)
methylation? Or are they predominantly unmethylated?

---

#### Step 3 — Extract methylation context per lncRNA, per accession

Download GSE43857 allc files from GEO (one file per accession, ~927 files).
Each file is named like `allc_ColO.tsv.gz` and has 7 columns:

```
chr   pos   strand  mc_class  mc_count  total  methylated
Chr1  1     +       CGT       5         10     1
```

Run step 3 on each file. For a cluster (SLURM/PBS), submit as a job array.
For a local machine with GNU parallel:

```bash
# Make output directory
mkdir -p methylation_per_accession/

# Run across all accessions in parallel (8 cores)
ls path/to/allc_files/*.allc.gz | parallel -j 8 \
    python step3_methylation_context.py \
        --allc {} \
        --regions lncRNAs_only.bed \
        --out_dir methylation_per_accession/

# To run only on TE-free lncRNAs:
ls path/to/allc_files/*.allc.gz | parallel -j 8 \
    python step3_methylation_context.py \
        --allc {} \
        --regions overlap_results/lncRNAs_TE_free.bed \
        --out_dir methylation_per_accession_TE_free/
```

Each accession produces one file:
`methylation_per_accession/{accession}_lncRNA_methylation.tsv`

Columns: `gene_id, chr, start, end, strand,`  
`n_CG, mc_CG, total_CG, frac_CG,`  
`n_CHG, mc_CHG, total_CHG, frac_CHG,`  
`n_CHH, mc_CHH, total_CHH, frac_CHH`

---

#### Step 4 — Classify and summarise across all accessions

```bash
python step4_classify_and_summarise.py \
    --tsv_dir methylation_per_accession/ \
    --out_dir results/
```

**Classification logic (mirrors Shahzad et al.):**

| State   | Rule                                           |
|---------|------------------------------------------------|
| `teM`   | CHG fraction ≥ 5% OR CHH fraction ≥ 5%        |
| `gbM`   | CG fraction ≥ 10% AND no significant non-CG   |
| `UM`    | No context exceeds its threshold               |

Outputs:
- `results/lncRNA_per_accession_states.tsv` — (gene_id × accession) state matrix
- `results/lncRNA_conservation.tsv` — per-lncRNA % UM/gbM/teM across accessions
- `results/context_summary.tsv` — population-wide mean CG/CHG/CHH fractions

---

### Key Comparison to the Reference Paper

| Metric                    | Coding genes (Shahzad 2025) | lncRNAs (this study) |
|---------------------------|----------------------------:|---------------------:|
| % gbM (>90% conservation) | 41%                         | **TBD**              |
| % UM                      | 33%                         | **TBD**              |
| % teM                     | 12%                         | **TBD**              |
| CG fraction (gbM genes)   | ~0.25–0.50                  | **TBD**              |
| CHH fraction (teM genes)  | high (>0.10)                | **TBD**              |

The critical question: are lncRNAs **more UM** than coding genes?
Do they show **independent gbM vs teM** (like coding genes), or is their
methylation predominantly driven by overlapping TE contamination?

---

### Notes on Scalability

- **step3** is the bottleneck: it processes 927 accession × 11,071 lncRNAs.
  Each run takes ~2–5 min depending on file size. Total: ~30 hours serial,
  ~4 hours with 8 cores. **Always use parallel execution.**
- **step4** loads all TSVs into memory simultaneously; with 927 files × 11K
  rows each, peak RAM usage is ~4–8 GB. Run on a machine with ≥16 GB RAM.
- The scripts are modular: you can substitute a different regions BED
  (TE-free, TE-overlapping, or protein-coding) into step3 without changes.

---

### File Formats Reference

**BED6** (standard, no header):
```
Chr1  3630  5899  CUFF_NC.1  0  +
```

**allc** (GSE43857):
```
Chr1  1  +  CGT  5  10  1
```
Column 4 (mc_class) encodes the trinucleotide context:
- Starts with `CG` → CG context
- Starts with `C`, ends with `G` (e.g. `CAG`, `CTG`) → CHG context
- Starts with `C`, ends with non-G (e.g. `CAA`, `CTT`) → CHH context
