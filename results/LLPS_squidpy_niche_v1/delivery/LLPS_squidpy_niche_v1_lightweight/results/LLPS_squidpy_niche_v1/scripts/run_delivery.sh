#!/bin/bash
set -euo pipefail
O=/home/jinshengxi/bladder_metabolism/results/LLPS_squidpy_niche_v1
for script in 14_stage_audit 13_model_sensitivity; do
 "$O/scripts/run_python.sh" "$O/scripts/$script.py" > "$O/logs/$script.log" 2>&1
done
ulimit -v 97656250
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
/usr/bin/time -v taskset -c 0-3 Rscript "$O/scripts/08_skill_plots.R" > "$O/logs/08_skill_plots.log" 2>&1
for script in 09_report 11_package; do
 "$O/scripts/run_python.sh" "$O/scripts/$script.py" > "$O/logs/$script.log" 2>&1
done
