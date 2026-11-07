#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# SLURM array job: PRS-CSx-MT simulation sweep
#
# Usage (after running generate_job_list.py):
#   sbatch --array=1-$(wc -l < job_list.txt) submit_jobs.sh
#
# Or submit the whole sweep at once:
#   python generate_job_list.py && sbatch --array=1-$(wc -l < job_list.txt) submit_jobs.sh
#
# PRS-CSx-MT vs PRSxtra comparison (compare_prsxtra.py) on the rg x frac grid:
#   python generate_job_list.py --group rg_frac_grid --out job_list_prsxtra.txt
#   sbatch --job-name=prsxtra --export=ALL,MODE=prsxtra \
#          --array=1-$(wc -l < job_list_prsxtra.txt) submit_jobs.sh
# (add PHI=1e-2 to --export for the fixed-phi arm; results go to results_prsxtra/)
# Summarize afterwards: python compare_prsxtra.py --collect results_prsxtra
# ─────────────────────────────────────────────────────────────────────────────

#SBATCH --job-name=prscsx_mt_sweep
#SBATCH --time=03:00:00   # 2.5x the untuned cost: a 4-value a_grid runs 4 MCMC
                          # fits in the auto-phi arm and 1 in the fixed arm. The
                          # Snakemake path asks 138 min for the same work.
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
# Each task runs 2 phi arms x (per-trait PRS-CSx + PRS-CSx-MT). The two requests
# above are flat for every scenario, because this path has no per-scenario
# scaling -- unlike the Snakemake profile, which at 5000/2000 asks 55 min / 3.0 GB
# for 45 of the 49 scenarios and 115 min / 4.5 GB for the one with
# n_trait x n_pop = 10. The flat 2 h / 5G covers that worst case with margin:
# the 2026-11-06 pilot measured mean 13.2 min, max 29.3 min and peak max_rss
# 0.8 GB over all 49 scenarios. The old 03:00:00 / 8G was sized for 1000/500.
#SBATCH --output=logs/%x_%a.out    # %x = job name, so sweep/prsxtra logs don't collide
#SBATCH --error=logs/%x_%a.err

# ── cluster-specific setup ────────────────────────────────────────────────────
# Uncomment and edit as needed for your cluster:
# module load python/3.9
# source /path/to/venv/bin/activate
#
# The runner needs numpy + scipy + h5py. Set PYTHON to an interpreter that has
# them if `python3` on the compute nodes does not (locally that is python3.12,
# NOT the default conda python3). Example:
#   export PYTHON=/usr/bin/python3.12
PYTHON="${PYTHON:-/mmfs1/home/jacks.local/frederick.boehm/.conda/envs/prscsx/bin/python}"

# ── resolve paths ─────────────────────────────────────────────────────────────
# Slurm runs a *copy* of this script from its spool dir, so ${BASH_SOURCE[0]}
# doesn't point here — use the submit directory (where `sbatch` was run) instead,
# falling back to the script location when run outside Slurm.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# ── grid settings ───────────────────────────────────────────────────────
# Read from config.yaml, which the Snakemake path already uses, so the two launch
# paths cannot drift apart again (they did: this script hardcoded 1000/500 after
# the 2026-11-06 chain-length check moved the sweep to 5000/2000). The literals
# are only a fallback for an unreadable or reformatted config. Override either
# per submission:  N_ITER=1000 N_BURNIN=500 sbatch --array=... submit_jobs.sh
# Reads scalars of any type: strips a trailing comment and surrounding quotes,
# so a_grid: "1,0.25,0.05,0.01" comes back as 1,0.25,0.05,0.01.
_cfg() {
    sed -n "s/^$1:[[:space:]]*\(.*\)/\1/p" "${SCRIPT_DIR}/config.yaml" | head -1 \
      | sed -e 's/[[:space:]]*#.*$//' -e 's/^"\(.*\)"$/\1/' -e "s/^'\(.*\)'$/\1/" \
            -e 's/[[:space:]]*$//'
}
N_SNP="${N_SNP:-$(_cfg n_snp)}";          N_SNP="${N_SNP:-2000}"
N_ITER="${N_ITER:-$(_cfg n_iter)}";       N_ITER="${N_ITER:-5000}"
N_BURNIN="${N_BURNIN:-$(_cfg n_burnin)}"; N_BURNIN="${N_BURNIN:-2000}"
# These four were never passed at all -- the script relied on run_one_replicate's
# own defaults happening to equal the config. That held until n_test moved to
# 10000 on 2026-11-07, so read them explicitly now.
SIM_MODE="${SIM_MODE:-$(_cfg sim_mode)}";   SIM_MODE="${SIM_MODE:-realistic}"
N_REF="${N_REF:-$(_cfg n_ref)}";            N_REF="${N_REF:-500}"
N_TEST="${N_TEST:-$(_cfg n_test)}";         N_TEST="${N_TEST:-10000}"
LD_HETERO="${LD_HETERO:-$(_cfg ld_hetero)}"; LD_HETERO="${LD_HETERO:-0.15}"
A_GRID="${A_GRID:-$(_cfg a_grid)}";         A_GRID="${A_GRID:-1}"
TUNE_FRAC="${TUNE_FRAC:-$(_cfg tune_frac)}"; TUNE_FRAC="${TUNE_FRAC:-0.5}"
# The PRSxtra track keeps its own cohort size, pinned at what it already ran at.
PRSXTRA_N_TEST="${PRSXTRA_N_TEST:-$(_cfg prsxtra_n_test)}"
PRSXTRA_N_TEST="${PRSXTRA_N_TEST:-5000}"
PRSXTRA_TUNE_FRAC="${PRSXTRA_TUNE_FRAC:-$(_cfg prsxtra_tune_frac)}"
PRSXTRA_TUNE_FRAC="${PRSXTRA_TUNE_FRAC:-0.5}"
echo "Grid: n_snp=${N_SNP} n_iter=${N_ITER} n_burnin=${N_BURNIN}"
echo "Sim:  sim_mode=${SIM_MODE} n_ref=${N_REF} n_test=${N_TEST} ld_hetero=${LD_HETERO}"
echo "a:    a_grid=${A_GRID} tune_frac=${TUNE_FRAC}"

# MODE selects the runner: 'sweep' (run_one_replicate.py, default) or
# 'prsxtra' (compare_prsxtra.py). Each mode has its own job list and results dir.
MODE="${MODE:-sweep}"
case "${MODE}" in
    sweep)
        JOB_LIST="${SCRIPT_DIR}/job_list.txt"
        RESULTS_DIR="${SCRIPT_DIR}/results" ;;
    prsxtra)
        JOB_LIST="${SCRIPT_DIR}/job_list_prsxtra.txt"
        RESULTS_DIR="${SCRIPT_DIR}/results_prsxtra" ;;
    *)
        echo "ERROR: unknown MODE='${MODE}' (expected sweep or prsxtra)" >&2
        exit 1 ;;
esac

mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${RESULTS_DIR}"

# ── read this task's (scenario, seed) ─────────────────────────────────────────
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${JOB_LIST}")
SCENARIO_NAME=$(echo "${LINE}" | awk '{print $1}')
SEED=$(echo "${LINE}" | awk '{print $2}')

echo "Job ${SLURM_ARRAY_TASK_ID} [${MODE}]: scenario=${SCENARIO_NAME}  seed=${SEED}"
echo "Start: $(date)"

# ── run ───────────────────────────────────────────────────────────────────────
if [ "${MODE}" = "prsxtra" ]; then
    "${PYTHON}" "${SCRIPT_DIR}/compare_prsxtra.py" \
        --scenario_name  "${SCENARIO_NAME}"               \
        --seed           "${SEED}"                         \
        --scenarios_file "${SCRIPT_DIR}/scenarios.json"   \
        --out_dir        "${RESULTS_DIR}"                  \
        --phi            "${PHI:-auto}"                    \
        --n_snp          "${N_SNP}"                        \
        --n_iter         "${N_ITER}"                       \
        --n_burnin       "${N_BURNIN}"                     \
        --n_ref          "${N_REF}"                        \
        --n_test         "${PRSXTRA_N_TEST}"               \
        --ld_hetero      "${LD_HETERO}"                    \
        --tune_frac      "${PRSXTRA_TUNE_FRAC}"
else
    "${PYTHON}" "${SCRIPT_DIR}/run_one_replicate.py" \
        --scenario_name  "${SCENARIO_NAME}"               \
        --seed           "${SEED}"                         \
        --scenarios_file "${SCRIPT_DIR}/scenarios.json"   \
        --out_dir        "${RESULTS_DIR}"                  \
        --n_snp          "${N_SNP}"                        \
        --n_iter         "${N_ITER}"                       \
        --n_burnin       "${N_BURNIN}"                     \
        --sim_mode       "${SIM_MODE}"                     \
        --n_ref          "${N_REF}"                        \
        --n_test         "${N_TEST}"                       \
        --ld_hetero      "${LD_HETERO}"                    \
        --a_grid         "${A_GRID}"                       \
        --tune_frac      "${TUNE_FRAC}"
fi

EXIT_CODE=$?
echo "End: $(date)  exit=${EXIT_CODE}"
exit ${EXIT_CODE}
