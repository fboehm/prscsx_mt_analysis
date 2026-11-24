#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# SLURM array: the EUR overshoot check, one task per line of jobs.txt (20).
#
# From this directory, on the login node (not plain sbatch; see launch.sh):
#   ./launch.sh -p quickq
# Afterwards:
#   python3 summarize.py
#
# Sizing as the scale check (same n_snp / chain / a_grid 1): ~53 min locally,
# ~0.9 GB peak, cluster ~1.16x faster.
# ─────────────────────────────────────────────────────────────────────────────
#SBATCH --job-name=mvcs_overshoot
#SBATCH --time=02:00:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.err

set -euo pipefail
PYTHON="${PYTHON:-/mmfs1/home/jacks.local/frederick.boehm/.conda/envs/prscsx/bin/python}"
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$HERE"
export SWEEP_DIR="$(cd ../.. && pwd)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Compute nodes have no git: launch.sh verifies a clean checkout on the login
# node and passes the SHA in; re-read HEAD from .git's files here to catch a
# checkout that moved after submission.
: "${SCALE_MT_SHA:?submit with ./launch.sh, not sbatch directly}"
MT_DIR="${PRSCSX_MT_DIR:-$(cd "$SWEEP_DIR/../.." && pwd)/PRScsx_mt}"
head_ref=$(cat "$MT_DIR/.git/HEAD")
if [[ "$head_ref" == ref:* ]]; then
    ref=${head_ref#ref: }
    have=$(cat "$MT_DIR/.git/$ref" 2>/dev/null || awk -v r="$ref" '$2==r{print $1}' "$MT_DIR/.git/packed-refs")
else
    have=$head_ref
fi
if [ "$have" != "$SCALE_MT_SHA" ]; then
    echo "PRScsx_mt HEAD is now ${have:0:7}; launched against ${SCALE_MT_SHA:0:7}" >&2
    exit 2
fi

read -r SC SEED < <(sed -n "${SLURM_ARRAY_TASK_ID}p" jobs.txt)
mkdir -p runs
"$PYTHON" run_overshoot.py "runs/trace_${SC}_seed${SEED}.json" \
    --scenario_name "$SC" --seed "$SEED" --scenarios_file scenarios_overshoot.json \
    --out_dir runs --n_snp 2000 --n_iter 5000 --n_burnin 2000 \
    --n_ref 500 --n_test 10000 --a_grid 1
