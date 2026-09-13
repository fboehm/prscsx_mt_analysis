# prscsx_mt_analysis

Analysis and simulation code for the **PRS-CSx-MT methods paper** — simulation
sweeps, the four-method benchmark, and the All of Us real-data analysis plan.
This repo holds the code that *uses* PRS-CSx-MT; the method itself lives in a
separate repo.

## Dependencies on sibling repos

This repo does **not** vendor the method or baseline code. It imports them from
two sibling repositories:

| Repo | Provides |
|------|----------|
| [`PRScsx_mt`](../PRScsx_mt) | The method: `PRScsx_mt.py`, `mcmc_gtb_mt.py`, `parse_genet_mt.py`, `gigrnd.py` |
| [`PRScsx`](../PRScsx) | The original PRS-CSx baseline: `mcmc_gtb.py`, `parse_genet.py` |

By default they are located as **siblings** of this repo (i.e. all three cloned
into one common directory). Override with environment variables if they live
elsewhere:

```bash
export PRSCSX_MT_DIR=/path/to/PRScsx_mt
export PRSCSX_DIR=/path/to/PRScsx
```

### Cloning to a cluster

```bash
cd /shared/project           # one common parent directory
git clone <url>/PRScsx_mt.git
git clone <url>/PRScsx.git
git clone <url>/prscsx_mt_analysis.git
# now PRScsx_mt/, PRScsx/, and prscsx_mt_analysis/ are siblings — imports resolve
```

Python needs **numpy, scipy, h5py**. Locally that is `/usr/bin/python3.12` (the
default conda `python3` lacks scipy). On the cluster, set `PYTHON` in
`sim_sweep/submit_jobs.sh` to an interpreter that has them.

## Running the simulation sweep (Slurm)

```bash
cd sim_sweep
python generate_job_list.py                                   # scenario x seed job list
sbatch --array=1-$(wc -l < job_list.txt) submit_jobs.sh       # run the array
python collect_results.py                                     # -> all_results.csv, deltas.csv
```

The `rg_frac_grid` scenario group is the crux experiment: it crosses genetic
correlation (`rg`) with the shared-causal fraction (`frac_shared_causal`) under
both an auto-φ and a matched fixed-φ arm, so the multi-trait benefit can be
decomposed into φ-estimation pooling vs. genuine cross-trait borrowing.

### Running the same sweep with Snakemake (Slurm)

`sim_sweep/Snakefile` wraps the harness above in a Snakemake DAG: one job per
`(scenario, seed)` replicate (calling `run_one_replicate.py`), then a `collect`
step (`collect_results.py`). It tracks and retries individual replicates instead
of resubmitting the whole array, and requests per-replicate Slurm resources that
scale with the scenario's `n_trait × n_pop`.

```bash
cd sim_sweep
snakemake -s Snakefile -n                       # dry run: preview the DAG
snakemake -s Snakefile --cores 8                # run locally

# On Slurm (Snakemake 7.x): edit profiles/slurm/config.yaml first
# (set slurm_account / slurm_partition), then:
mkdir -p logs/slurm
snakemake -s Snakefile --profile profiles/slurm
```

Sweep parameters (replicate count, MCMC length, interpreter, sibling-repo paths,
resource scaling) live in `sim_sweep/config.yaml`. The Snakefile header documents
the Snakemake ≥ 8 `--executor slurm` invocation (same resource names, no edits).

## Layout

- `sim_sweep/` — Slurm simulation harness (the main pipeline)
- `simulate_mt.py`, `simulate_indiv.py` — simulators
- `compare_methods.py`, `run_e2e_test.py`, `run_underrep_rescue.py`, `test_*.py` — standalone comparisons/tests
- `scripts/`, `Snakefile*`, `config*.yaml` — older Snakemake pipeline (produced the pilot `results/`)
- `analysis_plan_PRSxtra_MT.md` — All of Us real-data analysis plan
- `results/` — pilot simulation outputs
