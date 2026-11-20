#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# SLURM array: the mvcs scale check, one task per line of jobs.txt (42).
#
# From this directory:
#   sbatch --array=1-$(wc -l < jobs.txt) submit_scale.sh
# Afterwards:
#   python summarize.py
#
# Sizing: locally a replicate took ~53 min single-threaded (PRS-CSx ~45 of it)
# and peaked at ~0.9 GB RSS; the cluster has run ~1.16x faster than local.
# ─────────────────────────────────────────────────────────────────────────────
#SBATCH --job-name=mvcs_scale
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

# The check is only meaningful on the fixed mvcs sampler. run_one_replicate.py
# finds PRScsx_mt at $PRSCSX_MT_DIR or as a sibling of the analysis repo.
MT_DIR="${PRSCSX_MT_DIR:-$(cd "$SWEEP_DIR/../.." && pwd)/PRScsx_mt}"
want=3fb1b99
have=$(git -C "$MT_DIR" rev-parse --short=7 HEAD)
if [ "$have" != "$want" ] || [ -n "$(git -C "$MT_DIR" status --porcelain --untracked-files=no)" ]; then
    echo "PRScsx_mt at $MT_DIR is $have (dirty: $(git -C "$MT_DIR" status --porcelain --untracked-files=no | wc -l)); need clean $want" >&2
    exit 2
fi

read -r SC SEED < <(sed -n "${SLURM_ARRAY_TASK_ID}p" jobs.txt)
mkdir -p runs
"$PYTHON" run_scale.py "runs/clamp_${SC}_seed${SEED}.json" \
    --scenario_name "$SC" --seed "$SEED" --scenarios_file scenarios_scale.json \
    --out_dir runs --n_snp 2000 --n_iter 5000 --n_burnin 2000 \
    --n_ref 500 --n_test 10000 --a_grid 1
