#!/bin/bash
# Chain-length check: 4 scenarios x 5 seeds x {PRS-CSx, PRS-CSx-MT} x 4 chain lengths = 160 tasks.
#
# PRSCSX_MT_DIR must be a PRScsx_mt checkout on branch learned-lambda-psi, e.g.
#   git -C ../../../PRScsx_mt fetch && git worktree add ../../../PRScsx_mt-learned-lambda learned-lambda-psi
# Then, from this directory (sim_sweep/chain_length/):
#   PRSCSX_MT_DIR=$(realpath ../../../PRScsx_mt-learned-lambda) \
#     sbatch --export=ALL --array=1-$(wc -l < jobs.txt) submit_chain_length.sh
# Re-submit only missing tasks:  sbatch --export=ALL --array=$(python summarize_chain_length.py --missing) ...
# Summarize:                     python summarize_chain_length.py
#
# Rough cost (local timings): 20000 iterations ~50-60 min per task, 1000 ~3 min; ~60 core-hours total.

#SBATCH --job-name=chainlen
#SBATCH --account=337t79
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.err

PYTHON="${PYTHON:-/mmfs1/home/jacks.local/frederick.boehm/.conda/envs/prscsx/bin/python}"
DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
OUT="${OUT:-${DIR}/out}"
: "${PRSCSX_MT_DIR:?set PRSCSX_MT_DIR to a PRScsx_mt checkout on branch learned-lambda-psi}"
export PRSCSX_MT_DIR OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "${DIR}/logs" "${OUT}"

read -r SC SEED METHOD NITER NBURN <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DIR}/jobs.txt")"
echo "task ${SLURM_ARRAY_TASK_ID}: ${SC} seed=${SEED} ${METHOD} ${NITER}/${NBURN}  MT=$(git -C "${PRSCSX_MT_DIR}" rev-parse --short HEAD 2>/dev/null)"
echo "start $(date)"
"${PYTHON}" "${DIR}/chain_length.py" "${SC}" "${SEED}" "${METHOD}" "${NITER}" "${NBURN}" "${OUT}"
rc=$?
echo "end $(date) exit=${rc}"
exit ${rc}
