#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 OMP_THREAD_LIMIT=1
cd /home/jinshengxi/bladder_metabolism
pids=()
for worker in 0 1 2; do
 taskset -c "$worker" Rscript results/llps_C5_cellchat_v2/scripts/02_run_specimens.R "$worker" >> "results/llps_C5_cellchat_v2/00_audit/worker_${worker}.log" 2>&1 &
 pids+=("$!")
done
for pid in "${pids[@]}"; do wait "$pid"; done
