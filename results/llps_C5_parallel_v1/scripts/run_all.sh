#!/usr/bin/env bash
set -euo pipefail
cd /home/jinshengxi/bladder_metabolism
s=results/llps_C5_parallel_v1/scripts
bash "$s/run_limited.sh" "$s/01_composition.R"
bash "$s/run_limited.sh" "$s/02_stage.R"
python "$s/03_bulk_inputs.py"
bash "$s/run_limited.sh" "$s/03_bulk_survival.R"
bash "$s/run_limited.sh" "$s/03b_reference_comparison.R"
bash "$s/run_limited.sh" "$s/04_markers.R"
bash "$s/run_limited.sh" "$s/04c_volcano.R"
bash "$s/run_limited.sh" "$s/05_ORA.R"
ORA_DIRECTION=down bash "$s/run_limited.sh" "$s/05_ORA.R"
bash "$s/run_limited.sh" "$s/05b_ORA_figures.R"
bash "$s/run_limited.sh" "$s/06_GSEA.R"
bash "$s/run_limited.sh" "$s/07_validate.R"
python "$s/08_figure_QC.py"
python "$s/09_report.py"
