#!/bin/bash
set -euo pipefail
O=/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 OMP_THREAD_LIMIT=1
export PYTHONPATH=/tmp/c5_esda_deps:/tmp/c5_spatial_pydeps PYTHONDONTWRITEBYTECODE=1
for w in 0 1 2; do taskset -c 3 Rscript -e "source('$O/scripts/02_cellchat.R')" "$w" > "$O/00_input_audit/checkpoint_recheck_$w.log" 2>&1; done
for stage in 14_native_QC 03_ligand 08_extra_sensitivity 04_summarize 12_communication_depth 16_cellchat_patient_influence; do
 echo "START $stage"
 taskset -c 3 Rscript -e "source('$O/scripts/$stage.R')" > "$O/00_input_audit/${stage}_final.log" 2>&1
 echo "DONE $stage"
done
for stage in 15_result_index 06_integrate; do taskset -c 3 python "$O/scripts/$stage.py"; done
taskset -c 3 Rscript -e "source('$O/scripts/07_plots.R')" > "$O/figures/run_final.log" 2>&1
for stage in 11_independent_ligand_check 13_independent_within_check 09_report 10_validate 17_questions; do taskset -c 3 python "$O/scripts/$stage.py" > "$O/00_input_audit/${stage}_final.log" 2>&1; done
echo 'POSTPROCESS COMPLETE'
