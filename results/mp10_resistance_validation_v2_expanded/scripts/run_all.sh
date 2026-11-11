#!/usr/bin/env bash
set -euo pipefail
export R_PROFILE=/dev/null R_PROFILE_USER=/dev/null OMP_THREAD_LIMIT=1 RCPP_PARALLEL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
B=$(cd "$(dirname "$0")/.." && pwd)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
# Frozen input verification; original audit/freeze scripts are not rerun over an existing freeze.
/usr/bin/python3 "$B/scripts/check_freeze.py"
Rscript "$B/scripts/15_bulk.R" > "$B/logs/15_bulk_replay.log" 2>&1
/usr/bin/python3 "$B/scripts/18_model_drug.py" > "$B/logs/18_model_drug_replay.log" 2>&1
/usr/bin/python3 "$B/scripts/24_endpoint_fix_and_pairs.py" > "$B/logs/24_pairs_replay.log" 2>&1
Rscript "$B/scripts/19_clinical.R" > "$B/logs/19_clinical_replay.log" 2>&1
/usr/bin/python3 "$B/scripts/16_epi_sparse.py" > "$B/logs/16_epi_replay.log" 2>&1
/usr/bin/python3 "$B/scripts/21_epi_worker.py" > "$B/logs/21_cycle_replay.log" 2>&1
Rscript "$B/scripts/20_epi_pseudobulk.R" > "$B/logs/20_epiPB_replay.log" 2>&1
/usr/bin/python3 "$B/scripts/35_stage_audit.py"
/usr/bin/python3 "$B/scripts/23_epi_statistics.py" > "$B/logs/23_epi_replay.log" 2>&1
/usr/bin/python3 "$B/scripts/27_measured_sources.py"
Rscript "$B/scripts/30_replay.R"
/usr/bin/python3 "$B/scripts/33_completion_tables.py"
/usr/bin/python3 "$B/scripts/32_provenance_finish.py"
/usr/bin/python3 "$B/scripts/26_summary.py"
/usr/bin/python3 "$B/scripts/31_audit_final.py"
/usr/bin/python3 "$B/scripts/26_summary.py"
if [[ -f "$B/00_plan_audit/official_full_QC.json" ]]; then /usr/bin/python3 "$B/scripts/39_full_completion_check.py"; fi
Rscript "$B/scripts/28_figures.R" > "$B/logs/28_figures_replay.log" 2>&1
/usr/bin/python3 "$B/scripts/37_repeat_identity_qc.py"
/usr/bin/python3 "$B/scripts/29_report.py"

/usr/bin/python3 "$B/scripts/34_artifact_checks.py"
