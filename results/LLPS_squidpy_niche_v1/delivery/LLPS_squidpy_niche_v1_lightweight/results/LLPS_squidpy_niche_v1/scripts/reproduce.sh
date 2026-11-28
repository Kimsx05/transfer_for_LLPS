#!/bin/bash
set -euo pipefail
O="$(cd "$(dirname "$0")/.." && pwd)"
cd "$O/../.."
for script in 00_audit 01_composition 02_cluster 03_native 04_inference 05_spatial_stats 06_clinical 10_supplement 14_stage_audit 13_model_sensitivity 15_adjusted_model_sensitivity 07_figures 16_control_maps 17_final_audit; do
 "$O/scripts/run_python.sh" "$O/scripts/$script.py" > "$O/logs/$script.log" 2>&1
done
ulimit -v 97656250
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
/usr/bin/time -v taskset -c 0-3 Rscript "$O/scripts/08_skill_plots.R" > "$O/logs/08_skill_plots.log" 2>&1
for script in 09_report 11_package; do
 "$O/scripts/run_python.sh" "$O/scripts/$script.py" > "$O/logs/$script.log" 2>&1
done
