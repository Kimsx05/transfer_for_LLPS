#!/usr/bin/env bash
set -euo pipefail
cd /home/jinshengxi/bladder_metabolism
D=results/C5_Endothelial_GSEA_nominalNicheNet_v1/scripts
bash "$D/run_R.sh" "$D/01_GSEA_sets_LOPO.R" > "$D/01_run.log" 2>&1
bash "$D/run_R.sh" "$D/03_identity_audit.R" > "$D/03_identity.log" 2>&1
bash "$D/run_R.sh" "$D/03b_receptor_audit.R" > "$D/03b_receptor.log" 2>&1
bash "$D/run_R.sh" "$D/02_NicheNet.R" > "$D/02_run.log" 2>&1
bash "$D/run_R.sh" "$D/02b_target_influence.R" > "$D/02b_target.log" 2>&1
bash "$D/run_python.sh" "$D/06_freeze_new_axes.py"
bash "$D/run_python.sh" "$D/06_light_spatial.py" > "$D/06_light_spatial.log" 2>&1
bash "$D/run_python.sh" "$D/04_summarize.py" > "$D/04_summary.log" 2>&1
bash "$D/run_R.sh" "$D/05_plots.R" > "$D/05_plots.log" 2>&1
bash "$D/run_python.sh" "$D/07_report.py" > "$D/07_report.log" 2>&1
bash "$D/run_python.sh" "$D/08_validate.py" > "$D/08_validate.log" 2>&1
