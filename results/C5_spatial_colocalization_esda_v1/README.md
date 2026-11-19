# C5 spatial colocalization (esda) v1

Full exploratory analysis of 63 sections / 265,002 spots. Frozen Five-gene-associated program; formal C5-specific localization remains NO-GO.

See [comprehensive report](COMPREHENSIVE_REPORT.md) and [all-section gallery](GALLERY.html).

Reproduce from the project root with `bash results/C5_spatial_colocalization_esda_v1/scripts/run_pipeline.sh`. Uses isolated `/tmp/c5_esda_deps` (version pins in requirements.txt); R uses this result's `.Rlib`. Original parquet scores, raw data and historical CellChat files are read-only. Every task uses at most four affinity cores, BLAS/OpenMP/Numba/R single-thread; no nested parallelism. seed=123.

A section checkpoint in 00_input_audit is written only after all Python outputs are complete. Resume skips completed sections. Changing data, parameters, or statistical definitions requires a new version directory. Permutation simulation arrays and duplicate expression objects are not saved.

Statistical scope: global primary 252 tests; local BH per section over four targets and all valid spots. Package native min-tail P preserved; conservative doubled P for main BH. Isolates/missing observations are explicitly excluded. Exact rings 1–3 are not cumulative. Known-patient inference uses 29 audited patients only; unknown identities stay section-level descriptive. Local HH is a within-section relative program/neighbor abundance pattern, not a validated C5 niche or communication.

If the temporary Python dependency folder has been removed, reinstall with `python -m pip install --target /tmp/c5_esda_deps -r results/C5_spatial_colocalization_esda_v1/scripts/requirements.txt`. The isolated R library and original install logs remain in this output directory. Existing input hashes are verified on resume; changed inputs require a new result version.

Final checks: 1,665 passed; 134 PNG and 403 PDF figures (including standalone panels/legends), 63 local-statistic Parquets. PNG 300 dpi, sampled PDF vector/font checks, exact score/spot-id preservation and per-section local BH rechecked.
