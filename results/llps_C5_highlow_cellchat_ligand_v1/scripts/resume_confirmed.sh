#!/bin/bash
set -euo pipefail
O=/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 OMP_THREAD_LIMIT=1
export PYTHONPATH=/tmp/c5_esda_deps:/tmp/c5_spatial_pydeps PYTHONDONTWRITEBYTECODE=1
while [ ! -f "$O/00_input_audit/stromal_specimen_counts.csv" ]; do sleep 10; done
taskset -c 3 python "$O/scripts/01b_complete_zero_samples.py"
bash "$O/scripts/launch.sh" &
cc_pid=$!
taskset -c 3 python "$O/scripts/reuse_identical.py" > "$O/02_CellChat_highlow/reuse.log" 2>&1 &
reuse_pid=$!
wait "$cc_pid"
wait "$reuse_pid"
for stage in 04_summarize 12_communication_depth 14_native_QC 16_cellchat_patient_influence; do taskset -c 3 Rscript "$O/scripts/$stage.R" > "$O/00_input_audit/$stage.log" 2>&1; done
taskset -c 3 python "$O/scripts/15_result_index.py"
taskset -c 3 python "$O/scripts/06_integrate.py"
taskset -c 3 Rscript "$O/scripts/07_plots.R" > "$O/figures/run.log" 2>&1
for stage in 11_independent_ligand_check 13_independent_within_check 09_report 10_validate 09_report; do taskset -c 3 python "$O/scripts/$stage.py" > "$O/00_input_audit/$stage.log" 2>&1; done
