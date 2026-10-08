#!/usr/bin/env bash
set -euo pipefail
CC_RESULT_ROOT='/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7'
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# Installation is a separate user-approved step, never performed by this runner.
Rscript -e "source('$CC_RESULT_ROOT/scripts/common.R'); stopifnot(requireNamespace('CellChat',quietly=TRUE),packageVersion('CellChat')=='2.1.2')"
Rscript "$CC_RESULT_ROOT/scripts/test_consensus_semantics.R" > "$CC_RESULT_ROOT/00_audit/test_consensus_semantics.log" 2>&1
# Run one complete biological sample first; it establishes the common DB/API audit.
if [ ! -f "$CC_RESULT_ROOT/02_cellchat_per_sample/S001/status.tsv" ]; then
  Rscript "$CC_RESULT_ROOT/scripts/03_cellchat_per_sample.R" 1 > "$CC_RESULT_ROOT/02_cellchat_per_sample/S001/runner.log" 2>&1
fi
Rscript -e "source('$CC_RESULT_ROOT/scripts/common.R'); m<-fread(file.path(ROOT,'01_metadata/per_sample_input_manifest.tsv')); todo<-which(!file.exists(file.path(dirname(m\$input_path),'status.tsv'))); writeLines(as.character(todo),file.path(ROOT,'00_audit/pending_sample_indices.txt'))"
if ! cat "$CC_RESULT_ROOT/00_audit/pending_sample_indices.txt" | xargs -r -n 1 -P 8 bash -c 'Rscript "$1/scripts/03_cellchat_per_sample.R" "$2" > "$1/02_cellchat_per_sample/runner_${2}.log" 2>&1' _ "$CC_RESULT_ROOT"; then
  printf '%s\n' 'A sample worker ended nonzero; retain sample outcomes and recover completed native objects before cohort aggregation.' >> "$CC_RESULT_ROOT/00_audit/worker_command_status.log"
fi
Rscript "$CC_RESULT_ROOT/scripts/04_consensus.R" > "$CC_RESULT_ROOT/03_cellchat_consensus/04_consensus.log" 2>&1
Rscript "$CC_RESULT_ROOT/scripts/10_sensitivity_defined_totals.R" > "$CC_RESULT_ROOT/03_cellchat_consensus/10_sensitivity.log" 2>&1
Rscript "$CC_RESULT_ROOT/scripts/05_figures.R" > "$CC_RESULT_ROOT/09_figures/05_figures.log" 2>&1
Rscript "$CC_RESULT_ROOT/scripts/06_report.R" > "$CC_RESULT_ROOT/00_audit/06_report.log" 2>&1
