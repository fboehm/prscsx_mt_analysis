#!/bin/bash
# Submit one or more scenario groups as a sweep, on the LOGIN node.
#
#   ./launch_sweep.sh <results_dir> <group> [<group> ...] -- [sbatch args]
#   e.g. ./launch_sweep.sh results_h0.1 rg_frac_grid_h0.1 cross_trait_mvcs_h0.1 -- -p quickq
#
# Compute nodes have no git, so run_one_replicate.py used to stamp every result
# 'unknown'. This script checks here, where git exists, that PRScsx_mt and
# PRScsx are clean checkouts of the expected commits, and exports the SHAs for
# the runner to stamp. It also writes the job list next to the results, so
# relaunching a different group cannot clobber it.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
SWEEP_DIR=$(pwd)
PARENT=$(cd ../.. && pwd)

# Expected samplers. PRScsx_mt 3fb1b99 = fixed mvcs (its `mult` path is
# byte-identical to 31d3dc2, which ran sweep v4). PRScsx 4db1f54 = psi guard.
WANT_MT=3fb1b999ebbaaace78b1c4d5e687a21dbe64ce6e
WANT_CSX=4db1f541c0063be067ebf8573434a74d512cc6ff
MT_DIR="${PRSCSX_MT_DIR:-$PARENT/PRScsx_mt}"
CSX_DIR="${PRSCSX_DIR:-$PARENT/PRScsx}"

[ $# -ge 2 ] || { sed -n 2,5p "$0"; exit 1; }

check() {  # dir want label
    local have dirty
    have=$(git -C "$1" rev-parse HEAD)
    dirty=$(git -C "$1" status --porcelain --untracked-files=no)
    if [ "$have" != "$2" ] || [ -n "$dirty" ]; then
        echo "$3 at $1 is ${have:0:7}${dirty:+ (dirty)}; need clean ${2:0:7}" >&2
        exit 2
    fi
}
check "$MT_DIR" "$WANT_MT" PRScsx_mt
check "$CSX_DIR" "$WANT_CSX" PRScsx
[ -z "$(git status --porcelain --untracked-files=no -- . )" ] \
    || { echo "prscsx_mt_analysis/sim_sweep has uncommitted edits; commit first" >&2; exit 2; }

RESULTS_DIR="$SWEEP_DIR/$1"; shift
SC_GROUPS=()
while [ $# -gt 0 ] && [ "$1" != "--" ]; do SC_GROUPS+=("$1"); shift; done
[ "${1:-}" = "--" ] && shift

[ ${#SC_GROUPS[@]} -gt 0 ] || { echo "no scenario groups given" >&2; exit 1; }
JOB_LIST="$RESULTS_DIR/job_list.txt"
[ -e "$JOB_LIST" ] && { echo "$JOB_LIST exists; refusing to overwrite a launched sweep" >&2; exit 2; }
mkdir -p "$RESULTS_DIR" logs
TMP="$SWEEP_DIR/.joblist_$$"; trap 'rm -f "$TMP" "$TMP".g' EXIT
: > "$TMP"
for g in "${SC_GROUPS[@]}"; do   # an unknown group exits nonzero, leaving no job_list.txt
    python3 generate_job_list.py --group "$g" --out "$(basename "$TMP").g" > /dev/null
    cat "$TMP.g" >> "$TMP"
done
mv "$TMP" "$JOB_LIST"
N=$(wc -l < "$JOB_LIST")
{
    echo "launched $(date '+%F %T') from $(hostname)"
    echo "analysis $(git rev-parse HEAD)"
    echo "PRScsx_mt $WANT_MT"
    echo "PRScsx $WANT_CSX"
    echo "groups ${SC_GROUPS[*]}"
    echo "tasks $N"
} > "$RESULTS_DIR/LAUNCH.txt"
echo "$N tasks from: ${SC_GROUPS[*]}  ->  $RESULTS_DIR"

sbatch --job-name="sweep_$(basename "$RESULTS_DIR")" \
       --export=ALL,JOB_LIST="$JOB_LIST",RESULTS_DIR="$RESULTS_DIR",SWEEP_PRSCSX_MT_SHA="$WANT_MT",SWEEP_PRSCSX_SHA="$WANT_CSX",PRSCSX_MT_DIR="$MT_DIR",PRSCSX_DIR="$CSX_DIR" \
       --array=1-"$N" "$@" submit_jobs.sh
