#!/bin/bash
export PYTHONPATH=/tmp/c5_spatial_pydeps PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
exec taskset -c "$(/usr/bin/python -c 'import os;print(",".join(map(str,sorted(os.sched_getaffinity(0))[:4])))')" /usr/bin/python "$@"
