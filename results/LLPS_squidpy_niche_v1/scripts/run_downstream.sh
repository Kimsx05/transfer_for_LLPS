#!/bin/bash
set -euo pipefail
O=/home/jinshengxi/bladder_metabolism/results/LLPS_squidpy_niche_v1
for script in 03_native 04_inference 05_spatial_stats 06_clinical; do
 "$O/scripts/run_python.sh" "$O/scripts/$script.py" > "$O/logs/$script.log" 2>&1
done
