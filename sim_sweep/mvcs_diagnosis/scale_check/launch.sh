#!/bin/bash
# Run on the LOGIN node (compute nodes have no git). Checks that PRScsx_mt is a
# clean checkout of the fixed sampler, then submits the array. Extra arguments
# go to sbatch, e.g.:  ./launch.sh -p quickq
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
WANT=3fb1b999ebbaaace78b1c4d5e687a21dbe64ce6e
SWEEP_DIR="$(cd ../.. && pwd)"
MT_DIR="${PRSCSX_MT_DIR:-$(cd "$SWEEP_DIR/../.." && pwd)/PRScsx_mt}"

have=$(git -C "$MT_DIR" rev-parse HEAD)
dirty=$(git -C "$MT_DIR" status --porcelain --untracked-files=no)
if [ "$have" != "$WANT" ] || [ -n "$dirty" ]; then
    echo "PRScsx_mt at $MT_DIR is ${have:0:7}${dirty:+ (dirty)}; need clean ${WANT:0:7}" >&2
    exit 2
fi
mkdir -p logs runs
sbatch --export=ALL,PRSCSX_MT_DIR="$MT_DIR",SCALE_MT_SHA="$WANT" \
       --array=1-$(wc -l < jobs.txt) "$@" submit_scale.sh
