#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
TASK_RUN=mp10_cellchat_high_definition_sensitivity_v1
exec 9>"$TASK_RUN/objects/run.lock"
flock -n 9 || { echo 'Another execution holds this bundle lock' >&2; exit 1; }
for script in 01_prepare 02_infer 10_expression_validation 03_summarize 06_evidence 05_validate 07_boundary_checks 04_replot 11_finalize; do
  Rscript "$TASK_RUN/scripts/$script.R" > "$TASK_RUN/logs/$script.log" 2>&1
done
/usr/bin/python3 "$TASK_RUN/scripts/08_render_audit.py"
/usr/bin/python3 "$TASK_RUN/scripts/09_report.py"
/usr/bin/python3 "$TASK_RUN/scripts/08_render_audit.py"
