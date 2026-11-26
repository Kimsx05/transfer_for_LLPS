#!/bin/bash
set -euo pipefail
O=/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1
export PYTHONPATH=/tmp/c5_esda_deps:/tmp/c5_spatial_pydeps PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 OMP_THREAD_LIMIT=1
while ! grep -q 'LIGAND ASSOCIATIONS COMPLETE' "$O/05_LLPS_ligand_correlation/run_final.log"; do
 if grep -q 'Execution halted' "$O/05_LLPS_ligand_correlation/run_final.log"; then echo 'LIGAND FAILED'; exit 1; fi
 sleep 15
done
taskset -c 3 Rscript "$O/scripts/05_continuous.R" > "$O/07_sensitivity_analysis/continuous.log" 2>&1
taskset -c 3 Rscript "$O/scripts/08_extra_sensitivity.R" > "$O/07_sensitivity_analysis/extra.log" 2>&1
while true; do
 n=$(find "$O/02_CellChat_highlow" -name COMPLETE.csv | wc -l)
 e=$(find "$O/02_CellChat_highlow" -name ERROR.csv | wc -l)
 if [ "$e" -gt 0 ]; then echo 'CELLCHAT FAILED'; exit 1; fi
 if [ "$n" -eq 55 ]; then break; fi
 sleep 15
done
taskset -c 3 Rscript "$O/scripts/04_summarize.R" > "$O/03_patient_level_comparison/run.log" 2>&1
taskset -c 3 python "$O/scripts/06_integrate.py" > "$O/06_NicheNet_spatial_integration/run.log" 2>&1
taskset -c 3 Rscript "$O/scripts/07_plots.R" > "$O/figures/run.log" 2>&1
taskset -c 3 python "$O/scripts/09_report.py" > "$O/tables/report.log" 2>&1
taskset -c 3 python "$O/scripts/10_validate.py" > "$O/00_input_audit/final_QC.log" 2>&1
echo 'REMAINING PIPELINE COMPLETE'
