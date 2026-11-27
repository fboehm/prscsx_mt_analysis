#!/bin/bash
# Run on the LOGIN node (compute nodes have no git). Checks that PRScsx_mt is a
# clean checkout of the fixed sampler and that this analysis checkout has no
# uncommitted edits to the simulator, runner or this folder, then
# submits the array. Extra arguments go to sbatch, e.g.:  ./launch.sh -p quickq
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
WANT=a60d13084db3875abeb79286cc8d37d73bf3ecb9
SWEEP_DIR="$(cd ../.. && pwd)"
REPO_DIR="$(cd "$SWEEP_DIR/.." && pwd)"
MT_DIR="${PRSCSX_MT_DIR:-$(cd "$REPO_DIR/.." && pwd)/PRScsx_mt}"

have=$(git -C "$MT_DIR" rev-parse HEAD)
dirty=$(git -C "$MT_DIR" status --porcelain --untracked-files=no)
if [ "$have" != "$WANT" ] || [ -n "$dirty" ]; then
    echo "PRScsx_mt at $MT_DIR is ${have:0:7}${dirty:+ (dirty)}; need clean ${WANT:0:7} (branch exp/mvcs-psi-cap)" >&2
    exit 2
fi
grep -q 'true_ld_ref' "$REPO_DIR/simulate_realistic.py" \
    || { echo "simulate_realistic.py has no true_ld_ref; pull this branch first" >&2; exit 2; }
if [ -n "$(git -C "$REPO_DIR" status --porcelain --untracked-files=no -- simulate_realistic.py sim_sweep/run_one_replicate.py sim_sweep/mvcs_diagnosis/psi_cap)" ]; then
    echo "uncommitted edits to the simulator, runner or this folder; commit first" >&2
    exit 2
fi
[ -e runs ] && [ -n "$(ls -A runs 2>/dev/null)" ] \
    && { echo "runs/ is not empty; move it aside before relaunching" >&2; exit 2; }
mkdir -p logs runs
echo "launched $(date '+%F %T') analysis $(git -C "$REPO_DIR" rev-parse --short HEAD) PRScsx_mt ${WANT:0:7}" > runs/LAUNCH.txt
sbatch --export=ALL,PRSCSX_MT_DIR="$MT_DIR",SCALE_MT_SHA="$WANT" \
       --array=1-$(wc -l < jobs.txt) "$@" submit_cap.sh
