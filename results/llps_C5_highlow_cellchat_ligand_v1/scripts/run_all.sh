#!/bin/bash
set -euo pipefail
O=/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 OMP_THREAD_LIMIT=1
export PYTHONPATH=/tmp/c5_esda_deps:/tmp/c5_spatial_pydeps PYTHONDONTWRITEBYTECODE=1
if [ ! -f "$O/00_input_audit/C5_expression_small.rds" ]; then taskset -c 3 Rscript -e "source('$O/scripts/01_input.R')"; fi
if [ ! -f "$O/00_input_audit/stromal_identity_final.csv.gz" ]; then taskset -c 3 Rscript -e "source('$O/scripts/01c_stromal_identity.R')"; fi
taskset -c 3 python "$O/scripts/01b_complete_zero_samples.py"
taskset -c 3 python "$O/scripts/reuse_identical.py" --audit-only
# Native inference: three single-threaded R processes plus one lightweight monitor.
# No simultaneous helper analyses. COMPLETE checkpoints are skipped.
taskset -c 3 python "$O/scripts/supervise_native.py"
# All following stages are sequential; R source() parses the entire script before execution.
for stage in 14_native_QC 03_ligand 05_continuous 08_extra_sensitivity 04_summarize 12_communication_depth 16_cellchat_patient_influence 18_database_audit; do
 taskset -c 3 Rscript -e "source('$O/scripts/$stage.R')"
done
for stage in 15_result_index 06_integrate; do taskset -c 3 python "$O/scripts/$stage.py"; done
taskset -c 3 Rscript -e "source('$O/scripts/07_plots.R')"
for stage in 11_independent_ligand_check 13_independent_within_check 09_report 10_validate 17_questions; do taskset -c 3 python "$O/scripts/$stage.py"; done

# Refresh file manifest and build a light local delivery archive (no automatic upload).
taskset -c 3 python "$O/scripts/19_delivery.py"
