# Analysis Plan: Multi-Trait Multi-Ancestry PRS Using PRS-CSx-MT

**Reference:** He et al. (2026), *Nature Genetics*.
"Multi-trait and multi-ancestry genetic analysis of comorbid lung diseases and traits improves genetic discovery and polygenic risk prediction."

---

## Overview

This plan adapts the PRSxtra framework from He et al. using the local **PRS-CSx-MT** implementation (`PRScsx_mt/PRScsx_mt.py`).

**Key departure from He et al.:** He et al. used a two-step approach — MTAG (multi-trait) followed by PRS-CSx (multi-ancestry) — because PRS-CSx can only handle one trait at a time.  PRS-CSx-MT models cross-trait covariance and cross-ancestry LD *jointly in a single MCMC*, making the MTAG pre-processing step unnecessary and, if applied, would double-count the cross-trait information.  The pipeline here therefore drops Phase 1 and feeds raw GWAS summary statistics directly into PRS-CSx-MT.

| Phase | He et al. method | This plan |
|-------|-----------------|-----------|
| 1 | MTAG (multi-trait boost per ancestry) | **Skipped** — PRS-CSx-MT handles this internally |
| 2 | PRS-CSx per trait (multi-ancestry shrinkage) | **PRS-CSx-MT** (multi-trait + multi-ancestry jointly) |
| 3 | Ridge regression to combine 39 candidate scores | Ridge regression to combine 32 candidate scores |

**Target traits (following He et al.):**

| Trait category | Phenotypes |
|---------------|-----------|
| Disease | Asthma, COPD, Lung cancer |
| Spirometry | FEV1, FVC, FEV1/FVC |
| Smoking | Smoking status (ever/never), Cigarettes per day |

**Ancestry groups:** AFR, AMR, EAS, EUR (matching the paper's reference panels).

**Validation cohort:** All of Us Research Program (independent of GWAS discovery cohorts).

---

## Phase 1 — Data Acquisition and Harmonization

### Goal
Obtain and harmonize ancestry-specific GWAS summary statistics for all 8 traits, ready for direct input to PRS-CSx-MT.

### Inputs

| Trait | Source |
|-------|--------|
| Asthma, COPD | Global Biobank Meta-analysis Initiative (GBMI) |
| Lung cancer | LC-GWMA (lung cancer multipopulation genome-wide meta-analysis) |
| FEV1, FVC, FEV1/FVC | GWAS Catalog accessions GCST90705067–GCST90705072 (He et al. EAS); Pan-UKBB or LF-GWMA (EUR/AFR/AMR) |
| Smoking status, Cig/day | GSCAN (GWAS & Sequencing Consortium of Alcohol and Nicotine use) |

### Steps

1. **Download GWAS summary statistics** for each trait × ancestry combination (32 files: 8 traits × 4 ancestries). Prioritize ancestry-specific summary statistics where available; use multi-ancestry meta-analysis results for ancestries with no trait-specific GWAS.

2. **Harmonize to PRS-CSx-MT input format:**
   - Required columns: `SNP A1 A2 BETA SE N`
   - For binary traits (asthma, COPD, lung cancer): convert log-OR to `BETA`, with corresponding `SE`
   - Filter to HapMap3 SNPs only (using `snpinfo_mult_1kg_hm3`); ~1.1 M SNPs
   - Align effect alleles to the HapMap3 reference strand
   - Use GRCh37 coordinates throughout
   - Script: `scripts/format_gwas_sumstats.py` (adapt from `scripts/format_pan_ukbb.py`)

3. **Quality-check each file:**
   - Confirm genomic inflation factor (λGC) is reasonable
   - Verify effect direction of well-known sentinel variants (e.g., SERPINA1 for COPD, FTO/MC4R region for BMI-related traits)
   - Note: exclude traits with mean χ² ≤ 1.02 for a given ancestry, as these provide insufficient signal (following He et al.'s MTAG inclusion criterion, applied here as a QC filter)

### Outputs
`data/sumstats/formatted/{TRAIT}_{POP}.txt` — 32 harmonized summary statistic files.

---

## Phase 2 — Multi-Trait Multi-Ancestry Shrinkage with PRS-CSx-MT

### Goal
Jointly estimate posterior SNP effect sizes across all 8 traits and 4 ancestry groups in a single MCMC run, producing 32 candidate polygenic scores.

### Rationale for skipping MTAG
MTAG boosts per-trait summary statistics by borrowing information from correlated traits, and He et al. used it specifically because PRS-CSx cannot model multiple traits simultaneously.  PRS-CSx-MT incorporates the same cross-trait information sharing internally, via a multi-trait continuous shrinkage prior, making MTAG a redundant pre-processing step.  Running MTAG before PRS-CSx-MT would inflate the apparent cross-trait correlations by encoding them twice — once in the boosted summary statistics and again in the MCMC joint model.

### Inputs
- Harmonized GWAS summary statistics from Phase 1 (32 files)
- 1KG LD reference panels (`data/ld_ref/ldblk_1kg_{pop}/`)
- HapMap3 SNP info file (`data/ld_ref/snpinfo_mult_1kg_hm3`)
- Target `.bim` file (All of Us HM3 SNP subset)

### Configuration (`config_respiratory.yaml`)

```yaml
ref_dir:  "data/ld_ref"
bim_prefix: "data/aou/aou_hm3"   # replace with real AoU BIM prefix

traits:      ["asthma", "copd", "lung_cancer", "fev1", "fvc", "fev1_fvc",
              "smoking", "cpd"]
populations: ["EUR", "EAS", "AFR", "AMR"]

# GWAS sample sizes per trait × population (fill from actual GWAS manifests)
n_gwas:
  asthma:
    EUR: 1768106
    EAS: 100000    # update from GBMI manifest
    AFR: 100000
    AMR: 100000
  copd:
    EUR: 1370418
    EAS: 100000
    AFR: 100000
    AMR: 100000
  lung_cancer:
    EUR: 70156
    EAS: 10000     # update from LC-GWMA
    AFR: 5000
    AMR: 5000
  fev1:
    EUR: 562003
    EAS: 129685    # KCPS-II + TWB
    AFR: 10000
    AMR: 10000
  fvc:
    EUR: 562008
    EAS: 129685
    AFR: 10000
    AMR: 10000
  fev1_fvc:
    EUR: 561866
    EAS: 129685
    AFR: 10000
    AMR: 10000
  smoking:
    EUR: 3371039
    EAS: 100000
    AFR: 100000
    AMR: 100000
  cpd:
    EUR: 782050
    EAS: 50000
    AFR: 50000
    AMR: 50000

phi:      null    # let algorithm estimate; set 1e-2 for a faster initial run
n_iter:   1000
n_burnin: 500
thin:     5
n_jobs:   8       # chromosome-level parallelism

# Sample overlap correction
# All GBMI traits share participants across ancestries; set rho_pheno once
# phenotypic correlations between traits are estimated in training data.
rho_pheno: null   # e.g. 0.3 for lung function traits; null = no correction
n_overlap: null   # set if exact overlapping sample counts are known

sumstats_fmt_dir: "data/sumstats/formatted"
out_dir:          "results/respiratory/prscsx_mt"
```

### Steps

1. **Run PRS-CSx-MT** via an updated Snakefile (`Snakefile_respiratory`):
   ```bash
   snakemake -s Snakefile_respiratory --cores 8
   ```
   The core rule calls:
   ```bash
   python PRScsx_mt/PRScsx_mt.py \
     --ref_dir=data/ld_ref \
     --bim_prefix=data/aou/aou_hm3 \
     --sst_file='formatted/asthma_EUR.txt,asthma_EAS.txt,...;copd_EUR.txt,...;...' \
     --n_gwas='1768106,100000,...;1370418,...;...' \
     --pop=EUR,EAS,AFR,AMR \
     --out_dir=results/respiratory/prscsx_mt \
     --out_name=respiratory \
     --n_iter=1000 --n_burnin=500 --n_jobs=8
   ```

2. **Verify outputs.** PRS-CSx-MT writes one posterior weight file per trait × ancestry × chromosome.
   Expected: 8 traits × 4 ancestries × 22 chromosomes = 704 files.
   After concatenating chromosomes per trait × ancestry: **32 candidate scores**.

3. **Concatenate chromosome-level weight files** per trait × ancestry into genome-wide scoring files:
   ```bash
   # Example for asthma, EUR
   cat results/respiratory/prscsx_mt/respiratory_EUR_asthma_pst_eff_a1_b0.5_*_chr*.txt \
     > results/respiratory/scores/respiratory_EUR_asthma_genome.txt
   ```

### Outputs
`results/respiratory/scores/respiratory_{POP}_{TRAIT}_genome.txt` — 32 genome-wide posterior weight files.

---

## Phase 3 — Regularization via Ridge Regression

### Goal
Linearly combine the 32 candidate PRSs into a single final score that maximally predicts each disease outcome.

### Study population: All of Us Research Program

**Inclusion criteria** (matching He et al.):
- Whole-genome sequencing data, v7 release
- Self-reported sex and date of birth
- For COPD and lung cancer: self-reported smoking status required
- Exclude individuals with homozygous SERPINA1 variants (rs6647, rs709932, rs28929474) from COPD analysis

**Phenotype ascertainment** (ICD-9/ICD-10 codes, following He et al. Supplementary Tables 36–40):
- Asthma: ICD codes + self-reported family history
- COPD: ICD codes + spirometry data where available
- Lung cancer: ICD codes

**Train/validation split:** randomly assign 70% of QC-passed participants to training, 30% to held-out validation (stratified by ancestry group).

### Steps

1. **Score all All of Us participants** on each of the 32 candidate PRSs using PLINK2:
   ```bash
   plink2 \
     --bfile data/aou/aou_hm3 \
     --score results/respiratory/scores/respiratory_EUR_asthma_genome.txt 1 2 3 \
     --score-col-nums 3 \
     --out results/scores/cand_EUR_asthma
   ```
   Repeat for each trait × ancestry combination (32 calls, or parallelize).

2. **Standardize** each candidate PRS to mean 0, SD 1 using training-set statistics.

3. **Fit ridge regression** (R `glmnet`) separately for each disease outcome (asthma, COPD, lung cancer):
   ```r
   library(glmnet)

   # x: n_train × 32 matrix of standardized candidate PRSs
   # y: binary disease status (0/1)
   fit <- cv.glmnet(x, y, alpha = 0,        # alpha=0 → ridge
                    family = "binomial",
                    nfolds = 10,
                    type.measure = "auc")

   best_lambda <- fit$lambda.min
   coefs <- coef(fit, s = best_lambda)      # 32 weights + intercept
   ```

4. **Compute PRSxtra** in the validation set as the linear combination of candidate scores weighted by ridge coefficients.

5. **Evaluate performance** in the held-out 30%:
   - Primary metric: AUC (DeLong's two-sided test for comparisons)
   - Stratify by ancestry: AFR, AMR, EAS, EUR (and MID, SAS if sample sizes permit)
   - Benchmark comparisons:
     - Single-trait, ancestry-matched PRS (using PRS-CS on each trait × ancestry separately)
     - PRSxa: PRS-CSx per trait, ridge regression across ancestries only
     - PRSmix+ (multi-trait, single-ancestry library)
   - Report odds ratios per SD of PRSxtra
   - Evaluate in clinical subgroups: smokers vs. never-smokers (COPD/lung cancer); with vs. without family history (asthma)

6. **Assess disease exacerbation prediction** (secondary outcome):
   - COPD exacerbations: ICD-coded hospitalization and ER visits
   - Asthma exacerbations

### Outputs
- Ridge regression coefficient files: `results/ridge/PRSxtra_{DISEASE}_coefs.txt`
- Validation AUC summary table (by disease × ancestry)
- OR plots by PRS decile
- Score distribution and prevalence-by-decile figures

---

## Summary of Key Departures from He et al.

| Decision | He et al. | This plan | Rationale |
|----------|-----------|-----------|-----------|
| Multi-trait step | MTAG (separate, per ancestry) | **None** | PRS-CSx-MT models cross-trait info internally |
| Multi-ancestry step | PRS-CSx (separate, per trait) | **PRS-CSx-MT** (jointly) | Unified model avoids double-counting |
| Candidate scores | 39 (incl. meta-ancestry scores) | 32 (4 ancestries × 8 traits) | Simpler; meta-ancestry step optional |
| LD reference | 1KG | 1KG (same) | — |
| Target cohort | All of Us v7 | All of Us v7 (same) | — |
| Regularization | Ridge (glmnet) | Ridge (glmnet, same) | — |

---

## File Structure

```
test_claude/
├── config_respiratory.yaml          # new config for this analysis
├── Snakefile_respiratory            # new Snakefile (adapts Snakefile_real)
├── data/
│   ├── ld_ref/                      # 1KG LD panels (already downloaded)
│   ├── aou/
│   │   └── aou_hm3.bim              # export from All of Us
│   └── sumstats/
│       ├── raw/                     # downloaded GWAS files
│       └── formatted/               # harmonized, HM3-filtered files
├── results/
│   ├── respiratory/
│   │   ├── prscsx_mt/               # Phase 2 chromosome-level weight files
│   │   └── scores/                  # genome-wide weight files (32 files)
│   └── ridge/                       # Phase 3 coefficients + PRSxtra scores
├── scripts/
│   ├── format_gwas_sumstats.py      # harmonize raw GWAS → PRS-CSx-MT format
│   ├── compute_scores.sh            # PLINK2 scoring (32 candidate PRSs)
│   └── ridge_regression.R           # Phase 3 regularization
└── PRScsx_mt/                       # existing software
```

---

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| PRS-CSx-MT | local (`PRScsx_mt/`) | Phase 2 joint shrinkage |
| PLINK2 | ≥2.0 | Compute individual scores |
| R `glmnet` | ≥4.0 | Phase 3 ridge regression |
| Python 3 | ≥3.8 | Formatting scripts |
| Snakemake | ≥7.0 | Workflow management |

---

## Milestones

1. **Data acquisition** — Download and harmonize GWAS summary statistics for all 8 traits × 4 ancestries; export AoU HM3 `.bim` file.
2. **Phase 2** — Run PRS-CSx-MT; verify and concatenate weight files (32 genome-wide scores).
3. **Phase 3 training** — Score All of Us training set on 32 candidate PRSs; fit ridge regression per disease.
4. **Phase 3 validation** — Evaluate PRSxtra vs. benchmarks in held-out set, stratified by ancestry.
5. **Write-up** — AUC tables, OR plots, ancestry-stratified comparisons.
