#!/bin/bash
export PYTHONPATH=/tmp/c5_esda_deps PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMBA_CACHE_DIR=/tmp/c5_esda_numba
exec taskset -c "$(/usr/bin/python -c 'import os;print(",".join(map(str,sorted(os.sched_getaffinity(0))[:4])))')" /usr/bin/python "$@"
