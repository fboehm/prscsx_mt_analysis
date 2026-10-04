# PRS-CSx-MT vs PRSxtra in All of Us

Real-data version of the simulation comparison in `sim_sweep/compare_prsxtra.py`.
It compares PRS-CSx-MT, untuned and with ridge stacking, against PRSxtra
(He et al. 2026, *Nat Genet*: MTAG → PRS-CSx → ridge). Two trait sets:

| config | traits (candidate scores) | outcomes evaluated in All of Us |
|---|---|---|
| `config_prsxtra_lipids.yaml` | LDL, HDL, TG (Pan-UKBB) | LDL, HDL, TG (EHR labs) |
| `config_prsxtra_respiratory.yaml` | asthma, COPD, lung cancer, FEV1, FVC, FEV1/FVC, smoking, CPD | asthma, COPD, lung cancer |

Each family has one score per trait × ancestry plus a cross-ancestry META
score per trait (PRS-CSx `--meta`), matching He et al.'s candidate set.

## Methods

Every method gets the same LD reference, SNP set, phi setting, MCMC length,
tuning participants and CV folds.

| method | candidate scores | tuning |
|---|---|---|
| `prscsx_mt` | PRS-CSx-MT, ancestry- and trait-matched | none |
| `prscsx_mt_ridge` | all PRS-CSx-MT scores (traits × ancestries) | ridge |
| `prsxtra` | all MTAG → PRS-CSx scores (traits × ancestries) | ridge |
| `prsxa` | PRS-CSx scores of the target trait (all ancestries) | ridge |
| `prscsx`, `mtag_prscsx` | ancestry- and trait-matched single scores | none |

The main contrasts are `prscsx_mt_ridge − prsxtra` (same tuning data, same
number of candidates) and `prscsx_mt − prsxtra`. The second one gives PRSxtra
the tuning sample.

## Stage A: summary statistics → candidate weights (outside All of Us)

```bash
# lipids: Pan-UKBB is downloaded automatically
snakemake -s Snakefile_prsxtra --configfile config_prsxtra_lipids.yaml --cores 8

# respiratory: download the public GWAS first (GBMI, He et al. EAS lung function,
# Shrine et al. 2023, GSCAN, MVP + Biobank Japan lung cancer; ~12 GB)
python scripts/download_respiratory_gwas.py          # --list shows the plan
snakemake -s Snakefile_prsxtra --configfile config_prsxtra_respiratory.yaml --cores 8

# on Slurm
snakemake -s Snakefile_prsxtra --configfile <config> --profile sim_sweep/profiles/slurm
```

Prerequisites:

- **Python 3** with numpy, scipy, h5py, pandas and pyliftover. Sibling repos
  `../PRScsx` and `../PRScsx_mt`.
- **MTAG** (`git clone https://github.com/JonJala/mtag tools/mtag`). It needs
  **Python 2.7**: `conda env create -f tools/mtag/environment.yml`, then point
  `mtag.python` at that env's interpreter.
- **LD scores for each ancestry for MTAG** (`mtag.ld_ref_panel`). MTAG ships EUR
  only (`eur_w_ld_chr/`). For EAS, AFR and AMR, supply chromosome-split ldsc
  LD scores, e.g. computed from 1KG with ldsc, or the Pan-UKBB per-ancestry LD
  scores split by chromosome.
- **Independence:** confirm that no GWAS includes All of Us participants.

Outputs in `out_dir`:
- `candidates/candidates.tsv.gz`: rsID-keyed weights. Columns are
  `<family>__<trait>__<POP>`, where family is `prscsx_mt`, `mtag_prscsx` or `prscsx`.
- `candidates/manifest.json`: column list, fit settings, MTAG summary, QC and the
  `aou:` evaluation settings.
- `qc/sumstats_qc.tsv`: N, mean χ² and λGC per GWAS.
- `mtag/summary_<POP>.tsv`: which traits went through MTAG, and their N_eff.

## Stage B: inside the All of Us Researcher Workbench

Upload `candidates.tsv.gz`, `manifest.json`, `data/ld_ref/hm3_hg38.tsv` and
this `aou/` folder to the workspace bucket. Then, in a notebook terminal (R with
data.table, glmnet, pROC and jsonlite; Python with pandas; plink2):

```bash
# 0. HM3 subset of the GRCh38 genotypes. Copy your CDR's ACAF plink_bed
#    files locally first (gsutil -u $GOOGLE_PROJECT cp ...).
bash aou/00_extract_hm3.sh 'acaf/chr{CHR}' hm3_hg38.tsv hm3_aou

# 1. match HM3 rsIDs to All of Us variants (chr:pos:alleles, GRCh38)
python aou/01_match_variants.py --pvar hm3_aou.pvar --hm3_hg38 hm3_hg38.tsv --out_dir match

# 2. score everyone on every candidate
python aou/02_score_candidates.py --pfile hm3_aou --candidates candidates.tsv.gz \
    --variant_map match/variant_map.tsv --out_prefix scores/lipids

# 3. phenotypes + covariates + genetic ancestry/PCs (BigQuery)
python aou/03_extract_phenotypes.py --trait_set lipids \
    --ancestry_preds <CDR>/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv \
    --out pheno_lipids.tsv

# 4. ridge stacking + evaluation, per ancestry x outcome
Rscript aou/04_stack_evaluate.R --scores scores/lipids.scores.tsv.gz \
    --pheno pheno_lipids.tsv --manifest manifest.json --out_dir eval_lipids
```

Repeat steps 2–4 with the respiratory files.

### Recommended: refit on the All of Us SNP set

`01_match_variants.py` writes `match/aou_hm3_matched.bim`. It is a list of
variants with no participant data, so it can leave the Workbench. Download it,
set `bim_prefix` in both Stage A configs to its prefix, and rerun Stage A. Every
method is then fit on exactly the SNPs that can be scored, instead of the
all-HM3 proxy.

### Evaluation details (`04_stack_evaluate.R`)

- **Split:** within each ancestry, a 70/30 tuning/validation split, stratified
  by case status for binary outcomes. Set `aou.n_splits` > 1 to repeat it.
- **Ridge:** `cv.glmnet(alpha = 0)`, 10-fold CV, `lambda.min`. Covariates are
  unpenalized and the score is the PRS part of the fitted predictor.
- **Continuous outcomes:** incremental R² over covariates, with paired-bootstrap
  CIs and p-values for the contrasts. **Binary outcomes:** AUC of the score with
  DeLong CIs and paired DeLong tests, plus OR per SD adjusted for covariates.
- **Ancestries without a GWAS** (e.g. MID, SAS) get the ridge methods only;
  `--pop_map MID=EUR` would give them untuned scores as well.
- **Disclosure:** all outputs are aggregate. Cells with any count below 20 are
  flagged in `cells.tsv` and their counts masked, per All of Us dissemination
  policy.

### Before the final analysis

- **Phenotype codes:** disease outcomes use He et al.'s concept IDs
  (Supplementary Tables 36–40, in `aou/he2026_codes.tsv`) with child codes
  included. Not yet reproduced: their SERPINA1 exclusions for COPD and the
  smoking-status subgroups.
- **Statin flag:** the flag is "ever exposed", not "on a statin at the time of
  the measurement". The LDL adjustment divides by 0.7 for flagged participants.
- **Age:** for lipids, age is taken at the latest measurement. For diseases it
  is taken at `--ref_date`; set that to your CDR cutoff date.

## Tested

A synthetic end-to-end run covered every step except BigQuery (step 3) and
real MTAG. It used chromosome 22, 600 SNPs, EUR and EAS, three traits in three
raw-file layouts, and a stand-in MTAG using the same estimator. The run went
from raw GWAS files through Stage A, then steps 0, 1, 2 and 4 on a VCF with
GRCh38-style IDs and mixed REF/ALT orientation. Liftover was checked against
the real UCSC chain.
