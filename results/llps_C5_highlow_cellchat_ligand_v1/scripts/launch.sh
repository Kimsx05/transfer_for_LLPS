#!/bin/bash
set -euo pipefail
O=/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 OMP_THREAD_LIMIT=1
for w in 0 1 2; do
 taskset -c "$w" Rscript "$O/scripts/02_cellchat.R" "$w" > "$O/02_CellChat_highlow/worker_${w}.log" 2>&1 &
done
wait
