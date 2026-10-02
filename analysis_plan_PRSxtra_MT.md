# Analysis Plan: PRS-CSx-MT vs PRSxtra in All of Us

**Reference:** He et al. (2026), *Nature Genetics*.
"Multi-trait and multi-ancestry genetic analysis of comorbid lung diseases and traits improves genetic discovery and polygenic risk prediction."

**Pipeline:** `Snakefile_prsxtra` (Stage A, summary statistics) and `aou/` (Stage B, inside the All of Us Researcher Workbench). The runbook is `aou/README.md`.

---

## Question

Does PRS-CSx-MT, with or without ridge stacking, predict better in All of Us than **PRSxtra**, He et al.'s MTAG → PRS-CSx → ridge pipeline? This is the real-data counterpart of the simulation in `sim_sweep/compare_prsxtra.py` (results: `sim_sweep/analysis_prsxtra/report.md`). There, PRS-CSx-MT + ridge beat PRSxtra in every grid cell, by +0.009 (EUR) and +0.014 (EAS) in corr(PRS, y).

## Methods compared

| Method | Multi-trait step | Multi-ancestry step | Tuning in All of Us |
|---|---|---|---|
| **PRSxtra** (He et al.) | MTAG, separately within each ancestry | PRS-CSx per trait | ridge over all trait × ancestry scores |
| **PRS-CSx-MT + ridge** | PRS-CSx-MT (joint) | PRS-CSx-MT (joint) | ridge over all trait × ancestry scores |
| **PRS-CSx-MT** | PRS-CSx-MT (joint) | PRS-CSx-MT (joint) | none: the score matching the target ancestry and trait |
| PRSxa (benchmark) | none | PRS-CSx | ridge over the target trait's ancestry scores |
| PRS-CSx, MTAG → PRS-CSx | none / MTAG | PRS-CSx | none |

PRS-CSx-MT gets raw GWAS, never MTAG output. It models cross-trait sharing itself, so MTAG first would count that information twice. PRSxtra and PRS-CSx-MT + ridge have the same number of candidate scores (traits × ancestries), the same tuning participants and the same CV folds. Their difference isolates MTAG + PRS-CSx versus the joint model.

Not reproduced from He et al.: their meta-ancestry (`--meta`) candidate scores (39 vs 32 candidates). They could be added as extra columns in both families.

## Trait sets

| Set | Candidate traits | Outcomes in All of Us | GWAS sources |
|---|---|---|---|
| Lipids | LDL, HDL, TG | LDL, HDL, TG (EHR labs) | Pan-UKBB (EUR, EAS, AFR, AMR) |
| Respiratory (He et al.) | asthma, COPD, lung cancer, FEV1, FVC, FEV1/FVC, smoking, CPD | asthma, COPD, lung cancer | GBMI; LC-GWMA; GCST90705067–72 (He et al. EAS) + Pan-UKBB/LF-GWMA; GSCAN |

Ancestries: AFR, AMR, EAS, EUR (1KG LD panels). PRS-CSx-MT fits the full trait × ancestry grid, so every trait needs a GWAS in every ancestry used. Drop an ancestry rather than leave a gap.

---

## Stage A: summary statistics (outside All of Us)

`snakemake -s Snakefile_prsxtra --configfile config_prsxtra_{lipids,respiratory}.yaml`

1. **Harmonize** each GWAS to `SNP A1 A2 BETA SE` on HM3 SNPs (`scripts/format_pan_ukbb.py` and `scripts/format_gwas_sumstats.py`).
   - Matching is by rsID or by GRCh37/38 position, using a GRCh38 liftover of HM3 (`scripts/liftover_hm3.py`).
   - Odds ratios are converted to log-OR. For binary traits, N is the effective N = 4/(1/cases + 1/controls).
   - QC table: `qc/sumstats_qc.tsv` (mean χ², λGC).
2. **PRS-CSx-MT:** all traits × ancestries jointly, one job per chromosome.
3. **PRS-CSx:** per trait, on the raw GWAS.
4. **MTAG** within each ancestry (`scripts/run_mtag_pop.py`), using the real MTAG software with that ancestry's LD scores.
   - Traits with mean χ² ≤ 1.02 are passed through unchanged, as in He et al.
   - The MTAG output goes to PRS-CSx as z-scores with MTAG's GWAS-equivalent N.
5. **PRS-CSx** per trait on the MTAG output.
6. **Assemble** `candidates/candidates.tsv.gz` (three families × traits × ancestries) and `manifest.json`.

All fits share the LD reference, SNP set (`bim_prefix`), phi (auto), MCMC length and per-chromosome seeds. For the final run, set `bim_prefix` to `aou_hm3_matched` from Stage B step 1, so fits use exactly the SNPs that can be scored.

## Stage B: All of Us Researcher Workbench

0. **Extract** the HM3 positions from the GRCh38 ACAF genotype files (`aou/00_extract_hm3.sh`).
1. **Match** HM3 rsIDs to All of Us variants by chr:pos:alleles (`aou/01_match_variants.py`).
2. **Score** all candidates in one plink2 pass (`aou/02_score_candidates.py`).
3. **Extract phenotypes, covariates and ancestry** (`aou/03_extract_phenotypes.py`).
   - Lipids: per-person median of LOINC lab values, LDL divided by 0.7 for statin users, TG log-transformed.
   - Diseases: ICD-9/10 codes. **Replace the default code lists with He et al. Supplementary Tables 36–40**, and add their smoking-status and SERPINA1 criteria.
   - Ancestry groups and PCs come from the All of Us `ancestry_preds.tsv`.
4. **Stack and evaluate** (`aou/04_stack_evaluate.R`):
   - Within each ancestry: 70/30 tuning/validation split, stratified by case status.
   - Ridge: `cv.glmnet(alpha = 0)`, 10-fold, `lambda.min`; covariates (age, sex, PC1–10) unpenalized.
   - Continuous outcomes: incremental R² over covariates, with paired-bootstrap contrasts.
   - Binary outcomes: AUC, with paired DeLong contrasts, and OR per SD.
   - Optional repeated splits as a robustness check (not independent replicates).
   - Counts below 20 are masked before export.

## Primary and secondary contrasts

Primary, per ancestry × outcome:
- PRS-CSx-MT + ridge − PRSxtra
- PRS-CSx-MT − PRSxtra

Secondary:
- PRSxtra − PRSxa: what MTAG adds within the PRSxtra pipeline
- PRS-CSx-MT − PRS-CSx: the multi-trait gain without tuning
- PRS-CSx-MT + ridge − PRS-CSx-MT: the stacking gain

The simulation predicts small but consistent gains for PRS-CSx-MT. They should be largest in under-represented ancestries (EAS in the simulation) and where traits share causal variants.

## Milestones

1. Download the respiratory GWAS; fill column names and case/control counts in `config_prsxtra_respiratory.yaml`. Set up MTAG (Python 2.7) and per-ancestry LD scores.
2. Run Stage A for both trait sets with the HM3 proxy BIM.
3. Workbench steps 0–1; export `aou_hm3_matched.bim`; rerun Stage A on it.
4. Workbench steps 2–4 for both trait sets.
5. Write-up: tables by ancestry × outcome; comparison with the simulation.
