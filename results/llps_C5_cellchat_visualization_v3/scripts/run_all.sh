#!/bin/bash
set -euo pipefail
O=/home/jinshengxi/bladder_metabolism/results/llps_C5_cellchat_visualization_v3
for s in 01_prepare 02_expression 03_network_figures 04_pathway_LR 05_expression_figures 06_nichenet_figures 07_integration; do
 bash "$O/scripts/run_limited.sh" Rscript -e "source('$O/scripts/$s.R')"
done
bash "$O/scripts/run_limited.sh" python "$O/scripts/08_legacy_audit.py"
bash "$O/scripts/run_limited.sh" Rscript -e "source('$O/scripts/09_spatial_context.R')"
for s in 10_report 11_validate 12_delivery; do bash "$O/scripts/run_limited.sh" python "$O/scripts/$s.py"; done
