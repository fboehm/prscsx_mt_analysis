"""
scripts/format_pan_ukbb.py
--------------------------
Convert one Pan-UKBB flat-file (.tsv.bgz) into per-population summary
statistics files in PRScsx_mt format:

    SNP   A1   A2   BETA   SE

Called by Snakemake rule `format_sumstats`.

Inputs  (via snakemake object)
------
  snakemake.input.bgz     – Pan-UKBB compressed flat file for one phenotype
  snakemake.input.snpinfo – PRScsx HM3 SNP info file (snpinfo_mult_1kg_hm3)

Outputs (via snakemake object)
-------
  One file per population: {out_dir}/{trait}_{POP}.txt

Params
------
  snakemake.params.pops    – list of uppercase population codes, e.g. ["EUR","EAS",...]
  snakemake.params.out_dir – output directory

Wildcard
--------
  snakemake.wildcards.trait – trait name used in output filenames
"""

import os
import sys

import pandas as pd

# ---------------------------------------------------------------------------
# Inputs / params from Snakemake
# ---------------------------------------------------------------------------
bgz_file = snakemake.input.bgz
snpinfo  = snakemake.input.snpinfo
pops     = list(snakemake.params.pops)
out_dir  = snakemake.params.out_dir
trait    = snakemake.wildcards.trait

# Pan-UKBB uses lowercase population codes in column names
POP_MAP = {
    "AFR": "afr", "AMR": "amr", "CSA": "csa",
    "EAS": "eas", "EUR": "eur", "MID": "mid", "SAS": "sas",
}

ATGC = {"A", "T", "G", "C"}

# ---------------------------------------------------------------------------
# 1. Load HM3 SNP set for filtering
# ---------------------------------------------------------------------------
print(f"[{trait}] Loading HM3 SNP set from {snpinfo} ...")
hm3_snps = set()
with open(snpinfo) as fh:
    next(fh)  # skip header
    for line in fh:
        hm3_snps.add(line.split()[1])   # column 1 = rsID
print(f"[{trait}]   {len(hm3_snps):,} HM3 SNPs")

# ---------------------------------------------------------------------------
# 2. Identify which columns to read
# ---------------------------------------------------------------------------
fixed_cols = ["rsid", "ref", "alt"]
pop_cols   = []
for pop in pops:
    pk = POP_MAP[pop]
    pop_cols += [f"beta_{pk}", f"se_{pk}", f"low_confidence_{pk}"]

usecols = fixed_cols + pop_cols

# ---------------------------------------------------------------------------
# 3. Read the compressed file
# bgzip is gzip-compatible; pandas can read it directly.
# Use chunked reading to keep peak RAM manageable for large files.
# ---------------------------------------------------------------------------
print(f"[{trait}] Reading {bgz_file} (this may take a few minutes) ...")

chunk_size = 500_000
chunks = []

reader = pd.read_csv(
    bgz_file,
    sep="\t",
    compression="gzip",
    usecols=usecols,
    na_values=["NA", "nan", ""],
    low_memory=False,
    chunksize=chunk_size,
)

for i, chunk in enumerate(reader):
    # Filter to valid rsIDs
    chunk = chunk[chunk["rsid"].str.startswith("rs", na=False)]
    # Filter to single-nucleotide ATGC variants (no indels)
    chunk = chunk[chunk["ref"].isin(ATGC) & chunk["alt"].isin(ATGC)]
    # Filter to HM3 SNPs early to reduce memory
    chunk = chunk[chunk["rsid"].isin(hm3_snps)]
    if not chunk.empty:
        chunks.append(chunk)
    if (i + 1) % 10 == 0:
        n_kept = sum(len(c) for c in chunks)
        print(f"[{trait}]   chunk {i+1}: {n_kept:,} HM3 SNPs kept so far")

df = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=usecols)
print(f"[{trait}] {len(df):,} HM3 SNPs after initial filters")

if df.empty:
    print(f"[{trait}] WARNING: no SNPs remain after filtering — check phenocode and column names")
    # Still create empty output files so Snakemake doesn't fail
    os.makedirs(out_dir, exist_ok=True)
    for pop in pops:
        out_file = os.path.join(out_dir, f"{trait}_{pop}.txt")
        with open(out_file, "w") as fh:
            fh.write("SNP\tA1\tA2\tBETA\tSE\n")
    sys.exit(0)

# ---------------------------------------------------------------------------
# 4. Deduplicate rsIDs (keep first occurrence)
# ---------------------------------------------------------------------------
df = df.drop_duplicates(subset="rsid", keep="first")
print(f"[{trait}] {len(df):,} SNPs after deduplication")

# ---------------------------------------------------------------------------
# 5. Write one file per population
# ---------------------------------------------------------------------------
os.makedirs(out_dir, exist_ok=True)

for pop in pops:
    pk       = POP_MAP[pop]
    out_file = os.path.join(out_dir, f"{trait}_{pop}.txt")

    # Apply per-population filters
    b_col  = f"beta_{pk}"
    se_col = f"se_{pk}"
    lc_col = f"low_confidence_{pk}"

    # Treat missing low_confidence as True (exclude)
    lc = df[lc_col].fillna(True)

    mask = (
        df[b_col].notna()  &
        df[se_col].notna() &
        (df[se_col] > 0)   &
        (~lc.astype(bool))
    )
    pop_df = df.loc[mask, ["rsid", "alt", "ref", b_col, se_col]].copy()
    pop_df.columns = ["SNP", "A1", "A2", "BETA", "SE"]

    pop_df.to_csv(out_file, sep="\t", index=False)

    # Approximate median N from SE: N ~ 1/median(SE)^2 (rough, continuous traits)
    if len(pop_df) > 0:
        med_n_approx = int(round(1.0 / pop_df["SE"].median() ** 2))
    else:
        med_n_approx = 0

    print(
        f"[{trait}] {pop}: {len(pop_df):,} SNPs written to {out_file} "
        f"(median-SE-implied N ≈ {med_n_approx:,})"
    )
