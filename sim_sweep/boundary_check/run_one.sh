#!/bin/sh
# One boundary-check replicate. Args: <scenario> <seed>. BC_DIR from the environment.
sc="$1"; s="$2"
cd "/mnt/c/Users/Fred/OneDrive - South Dakota State University - SDSU/test_claude/prscsx_mt_analysis/sim_sweep" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/usr/bin/python3.12 run_one_replicate.py \
    --scenario_name "$sc" --seed "$s" \
    --n_snp 2000 --n_iter 5000 --n_burnin 2000 \
    --n_ref 500 --n_test 10000 \
    --a_grid 0.05,0.01,0.005,0.002 \
    --out_dir "$BC_DIR/out" > "$BC_DIR/log_${sc}_${s}.txt" 2>&1
code=$?
echo "$sc $s exit=$code" >> "$BC_DIR/status.txt"
rm -rf "$BC_DIR/out/$sc/seed_$s/data"      # 168 MB each, regenerable from the seed
exit $code
