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
#SBATCH --time=03:00:00           # each task now runs 2 phi arms x (per-trait PRS-CSx + PRS-CSx-MT)
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
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
        --n_snp          2000                              \
        --n_iter         1000                              \
        --n_burnin       500
else
    "${PYTHON}" "${SCRIPT_DIR}/run_one_replicate.py" \
        --scenario_name  "${SCENARIO_NAME}"               \
        --seed           "${SEED}"                         \
        --scenarios_file "${SCRIPT_DIR}/scenarios.json"   \
        --out_dir        "${RESULTS_DIR}"                  \
        --n_snp          2000                              \
        --n_iter         1000                              \
        --n_burnin       500
fi

EXIT_CODE=$?
echo "End: $(date)  exit=${EXIT_CODE}"
exit ${EXIT_CODE}
