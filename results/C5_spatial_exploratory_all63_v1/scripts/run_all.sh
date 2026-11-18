#!/bin/bash
set -euo pipefail
cd /home/jinshengxi/bladder_metabolism
export PYTHONPATH=/tmp/c5_spatial_pydeps PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
TASK_CPUS=$(/usr/bin/python -c 'import os;print(",".join(map(str,sorted(os.sched_getaffinity(0))[:4])))')
O=results/C5_spatial_exploratory_all63_v1
runpy(){ taskset -c "$TASK_CPUS" /usr/bin/python "$O/scripts/$1" "${@:2}"; }
runpy score_all.py pilot
runpy batches.py
runpy analyze.py
runpy epithelial_sensitivity.py
runpy gene_attribution.py
taskset -c "$TASK_CPUS" Rscript "$O/scripts/plots.R"
runpy validate_figures.py
runpy finalize.py
