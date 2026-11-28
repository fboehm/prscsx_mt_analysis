#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# SLURM array: stage 2 (real chr22 LD), one task per line of
# jobs.txt (120).
#
# From this directory, on the login node (not plain sbatch; see launch.sh):
#   ./launch.sh -p quickq
# Afterwards:
#   python3 summarize.py > summary.txt
#
# chr22 from the panel (n_snp comes from it; --n_snp is ignored). check_panel.py projects the runtime.
# ─────────────────────────────────────────────────────────────────────────────
#SBATCH --job-name=mvcs_rld
#SBATCH --time=16:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%a.out
#SBATCH --error=logs/%x_%a.err

set -euo pipefail
PYTHON="${PYTHON:-/mmfs1/home/jacks.local/frederick.boehm/.conda/envs/prscsx/bin/python}"
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$HERE"
export SWEEP_DIR="$(cd ../.. && pwd)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Compute nodes have no git: launch.sh verifies clean checkouts on the login
# node and passes the SHAs in; re-read HEAD from .git's files here to catch a
# checkout that moved after submission.
: "${SCALE_MT_SHA:?submit with ./launch.sh, not sbatch directly}"
: "${SCALE_CS_SHA:?submit with ./launch.sh, not sbatch directly}"
: "${LD_PANEL_DIR:?submit with ./launch.sh, not sbatch directly}"
head_of() {
    local head_ref ref
    head_ref=$(cat "$1/.git/HEAD")
    if [[ "$head_ref" == ref:* ]]; then
        ref=${head_ref#ref: }
        cat "$1/.git/$ref" 2>/dev/null || awk -v r="$ref" '$2==r{print $1}' "$1/.git/packed-refs"
    else
        echo "$head_ref"
    fi
}
for pair in "$PRSCSX_MT_DIR:$SCALE_MT_SHA" "$PRSCSX_DIR:$SCALE_CS_SHA"; do
    dir=${pair%:*}; want=${pair##*:}
    have=$(head_of "$dir")
    if [ "$have" != "$want" ]; then
        echo "$dir HEAD is now ${have:0:7}; launched against ${want:0:7}" >&2
        exit 2
    fi
done

read -r SC SEED < <(sed -n "${SLURM_ARRAY_TASK_ID}p" jobs.txt)
mkdir -p runs
"$PYTHON" ../psi_cap/run_cap.py "runs/trace_${SC}_seed${SEED}.json" \
    --scenario_name "$SC" --seed "$SEED" --scenarios_file scenarios_real_ld.json \
    --out_dir runs --n_snp 10000 --n_iter 5000 --n_burnin 2000 \
    --n_ref 500 --n_test 10000 --a_grid 1
