#!/bin/bash
set -euo pipefail
cd /home/jinshengxi/bladder_metabolism
O=results/C5_fivegene_LLPS77_spatial_v1
bash "$O/scripts/run_python.sh" "$O/scripts/prepare.py"
bash "$O/scripts/run_python.sh" -c "import pandas as pd,subprocess; o='$O'; ids=','.join(pd.read_csv(o+'/00_input_audit/pilot_prespecified.tsv',sep='\t').section_id);subprocess.run(['bash',o+'/scripts/run_python.sh',o+'/scripts/score_new.py',ids],check=True)"
bash "$O/scripts/run_python.sh" "$O/scripts/reproduce_pilot.py"
bash "$O/scripts/run_python.sh" "$O/scripts/launch_all.py"
bash "$O/scripts/run_python.sh" "$O/scripts/analyze_new.py"
bash "$O/scripts/run_python.sh" "$O/scripts/threshold_rings.py"
bash "$O/scripts/run_python.sh" "$O/scripts/AUCell_rings.py"
bash "$O/scripts/run_python.sh" "$O/scripts/zero_diagnostics.py"
bash "$O/scripts/run_python.sh" "$O/scripts/summarize.py"
bash "$O/scripts/run_python.sh" "$O/scripts/supplement.py"
bash "$O/scripts/run_python.sh" -c "import subprocess;subprocess.run(['Rscript','$O/scripts/plots_new.R'],check=True)"
bash "$O/scripts/run_python.sh" "$O/scripts/report.py"
bash "$O/scripts/run_python.sh" "$O/scripts/interpret_results.py"
bash "$O/scripts/run_python.sh" "$O/scripts/finish_additional.py"
bash "$O/scripts/run_python.sh" -c "import subprocess;subprocess.run(['Rscript','$O/scripts/additional_plots.R'],check=True)"
bash "$O/scripts/run_python.sh" "$O/scripts/validate_final.py"
