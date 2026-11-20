#!/bin/bash
# Runs every line of jobs.txt ("<scenario> <seed>"), 5 at a time (7 GB VM; ~0.7 GB each), one thread each.
cd "$(dirname "$0")"
export SWEEP_DIR=$(cd ../.. && pwd) OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
one() {
  sc=$1; sd=$2
  /usr/bin/python3.12 run_scale.py runs/clamp_${sc}_seed${sd}.json --scenario_name $sc --seed $sd \
    --scenarios_file scenarios_scale.json --out_dir runs --n_snp 2000 --n_iter 5000 --n_burnin 2000 \
    --n_ref 500 --n_test 10000 --a_grid 1 > logs/${sc}_seed${sd}.log 2>&1; rc=$?
  echo "$(date '+%F %T') exit=$rc $sc $sd" >> logs/done.txt
}
export -f one
xargs -P 5 -L 1 bash -c 'one $0 $1' < jobs.txt
