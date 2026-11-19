#!/bin/bash
set -euo pipefail
cd /home/jinshengxi/bladder_metabolism
o=results/C5_spatial_colocalization_esda_v1
bash "$o/scripts/run_python.sh" "$o/scripts/analyze.py" audit
bash "$o/scripts/run_python.sh" "$o/scripts/audit_inputs.py"
bash "$o/scripts/run_python.sh" "$o/scripts/analyze.py" pilot
bash "$o/scripts/run_R.sh" "$o/scripts/lee.R" pilot
bash "$o/scripts/run_python.sh" "$o/scripts/analyze.py" all
bash "$o/scripts/run_R.sh" "$o/scripts/lee.R" all
bash "$o/scripts/run_python.sh" "$o/scripts/translation.py"
bash "$o/scripts/run_python.sh" "$o/scripts/summarize.py"
bash "$o/scripts/run_R.sh" "$o/scripts/plots.R" maps
bash "$o/scripts/run_R.sh" "$o/scripts/plots.R" summary
bash "$o/scripts/run_R.sh" "$o/scripts/shared_legends.R"
bash "$o/scripts/run_python.sh" "$o/scripts/report.py"
bash "$o/scripts/run_python.sh" "$o/scripts/finalize.py"
