#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONPATH="$ROOT/.pydeps" PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4 NUMBA_NUM_THREADS=4
export MPLCONFIGDIR="$ROOT/.mplconfig" NUMBA_CACHE_DIR="$ROOT/.numba_cache"
exec /usr/bin/python "$ROOT/scripts/limited_exec.py" "$@"
