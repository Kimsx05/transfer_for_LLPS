#!/usr/bin/env bash
set -euo pipefail
B=$(cd "$(dirname "$0")/.." && pwd)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 OMP_THREAD_LIMIT=1
export R_PROFILE=/dev/null R_PROFILE_USER=/dev/null
hour=$(TZ=Asia/Shanghai date +%H)
if [[ "${MP10_RUN_NOW:-0}" != "1" ]] && ((10#$hour >= 7 && 10#$hour < 23)); then echo 'Outside23:00-07:00CST; launch through night_gate.py'; exit 3; fi
/usr/bin/python3 "$B/scripts/check_freeze.py"
/usr/bin/python3 "$B/scripts/21_epi_worker.py" --official > "$B/logs/21_official_full.log" 2>&1
/usr/bin/python3 "$B/scripts/38_official_method_qc.py" > "$B/logs/38_official_method_qc.log" 2>&1
/usr/bin/python3 "$B/scripts/35_stage_audit.py"
/usr/bin/python3 "$B/scripts/23_epi_statistics.py" > "$B/logs/23_epi_night.log" 2>&1
/usr/bin/python3 "$B/scripts/33_completion_tables.py"
/usr/bin/python3 "$B/scripts/26_summary.py" > "$B/logs/26_night.log" 2>&1
/usr/bin/python3 "$B/scripts/31_audit_final.py" > "$B/logs/31_night_validation.log" 2>&1
if [[ -f "$B/00_plan_audit/official_full_QC.json" ]]; then /usr/bin/python3 "$B/scripts/39_full_completion_check.py"; fi
Rscript "$B/scripts/28_figures.R" > "$B/logs/28_night_figures.log" 2>&1
/usr/bin/python3 "$B/scripts/37_repeat_identity_qc.py"
/usr/bin/python3 "$B/scripts/29_report.py" > "$B/logs/29_night_report.log" 2>&1

/usr/bin/python3 "$B/scripts/34_artifact_checks.py" > "$B/logs/34_night_artifacts.log" 2>&1
