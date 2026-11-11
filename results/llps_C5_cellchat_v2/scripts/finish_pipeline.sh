#!/bin/bash
set -euo pipefail
cd /home/jinshengxi/bladder_metabolism
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 OMP_THREAD_LIMIT=1
# The entire helper tree stays on CPU 3; inference uses CPU 0,1,2.
if [[ "${1:-}" != "--pinned" ]]; then exec taskset -c 3 bash "$0" --pinned; fi
out=results/llps_C5_cellchat_v2
while true; do
 Rscript "$out/scripts/02_run_specimens.R" --export-only >> "$out/00_audit/export_repair_continuous.log" 2>&1
 count=$(python - <<'PY'
from pathlib import Path
p=Path('results/llps_C5_cellchat_v2')
print(sum(len(list(d.glob('*/SUCCESS.csv'))) for d in p.glob('0[123]_threshold*')))
PY
 )
 if [[ "$count" == "177" ]]; then break; fi
 sleep 30
done
Rscript "$out/scripts/03_summarize.R" > "$out/00_audit/summaries.log" 2>&1
Rscript "$out/scripts/04_compare_plot.R" > "$out/00_audit/comparison_figures.log" 2>&1
Rscript "$out/scripts/05_validation.R" > "$out/00_audit/validation.log" 2>&1
python "$out/scripts/06_report.py" > "$out/00_audit/report.log" 2>&1
