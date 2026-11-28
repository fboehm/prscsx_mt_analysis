#!/bin/bash
# Run on the LOGIN node (compute nodes have no git). Checks that PRScsx_mt and
# PRScsx are clean checkouts of their psi_cap branches and that this analysis
# checkout has no uncommitted edits to the simulator, runner, tracer or this
# folder, then submits the array. Extra arguments go to sbatch, e.g.:
#   ./launch.sh -p quickq
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
WANT_MT=a60d13084db3875abeb79286cc8d37d73bf3ecb9   # PRScsx_mt exp/mvcs-psi-cap
WANT_CS=29711ce9cd2bab8d3976b25adcf57251d9d68ef1   # PRScsx    exp/uncapped-psi
SWEEP_DIR="$(cd ../.. && pwd)"
REPO_DIR="$(cd "$SWEEP_DIR/.." && pwd)"
MT_DIR="${PRSCSX_MT_DIR:-$(cd "$REPO_DIR/.." && pwd)/PRScsx_mt}"
CS_DIR="${PRSCSX_DIR:-$(cd "$REPO_DIR/.." && pwd)/PRScsx}"

check() {   # dir want branch
    local have dirty
    have=$(git -C "$1" rev-parse HEAD)
    dirty=$(git -C "$1" status --porcelain --untracked-files=no)
    if [ "$have" != "$2" ] || [ -n "$dirty" ]; then
        echo "$1 is ${have:0:7}${dirty:+ (dirty)}; need clean ${2:0:7} (branch $3)" >&2
        exit 2
    fi
}
check "$MT_DIR" "$WANT_MT" exp/mvcs-psi-cap
check "$CS_DIR" "$WANT_CS" exp/uncapped-psi
grep -q "def load_ld_panel" "$REPO_DIR/simulate_realistic.py" \
    || { echo "simulate_realistic.py has no load_ld_panel; pull this branch first" >&2; exit 2; }
grep -q 'true_ld_ref' "$REPO_DIR/simulate_realistic.py" \
    || { echo "simulate_realistic.py has no true_ld_ref; pull this branch first" >&2; exit 2; }
if [ -n "$(git -C "$REPO_DIR" status --porcelain --untracked-files=no -- simulate_realistic.py sim_sweep/run_one_replicate.py sim_sweep/mvcs_diagnosis/psi_cap/run_cap.py sim_sweep/mvcs_diagnosis/real_ld)" ]; then
    echo "uncommitted edits to the simulator, runner, tracer or this folder; commit first" >&2
    exit 2
fi
LD_PANEL_DIR="${LD_PANEL_DIR:-$REPO_DIR/data/ld_ref}"
for f in snpinfo_mult_1kg_hm3 ldblk_1kg_eur/ldblk_1kg_chr22.hdf5 ldblk_1kg_eas/ldblk_1kg_chr22.hdf5; do
    [ -f "$LD_PANEL_DIR/$f" ] || { echo "LD panel file missing: $LD_PANEL_DIR/$f (set LD_PANEL_DIR)" >&2; exit 2; }
done
[ -e runs ] && [ -n "$(ls -A runs 2>/dev/null)" ] \
    && { echo "runs/ is not empty; move it aside before relaunching" >&2; exit 2; }
mkdir -p logs runs
echo "launched $(date '+%F %T') panel $LD_PANEL_DIR analysis $(git -C "$REPO_DIR" rev-parse --short HEAD) PRScsx_mt ${WANT_MT:0:7} PRScsx ${WANT_CS:0:7}" > runs/LAUNCH.txt
sbatch --export=ALL,PRSCSX_MT_DIR="$MT_DIR",PRSCSX_DIR="$CS_DIR",LD_PANEL_DIR="$LD_PANEL_DIR",SCALE_MT_SHA="$WANT_MT",SCALE_CS_SHA="$WANT_CS" \
       --array=1-$(wc -l < jobs.txt) "$@" submit_real_ld.sh
