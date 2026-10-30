#!/usr/bin/env bash
set -euo pipefail
cd /home/jinshengxi/bladder_metabolism
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
v3_cpus=$(python -c 'import os; print(",".join(map(str, sorted(os.sched_getaffinity(0))[:8])))')
taskset -pc "$v3_cpus" $$ >/dev/null
exec Rscript "$@"
