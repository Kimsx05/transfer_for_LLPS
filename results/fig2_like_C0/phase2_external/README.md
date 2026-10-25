# C0 external validation

Current state: **STOP_ABUNDANCE_CLINICAL_ANALYSIS**. Read `../summary/phase2_summary.md` for results and limitations. Stage/KM/Cox were deliberately not run after failed validation. This is not a negative survival-association result.

Run paths are relative to `/home/jinshengxi/bladder_metabolism`. Existing valid results need not be recomputed. Reproduction order if intentionally rebuilding:

1. `scripts/01_download_public.py`, then `scripts/01b_download_public_additional.py`: one cached copy of GEO/Xena public inputs and method sources. Download status, URLs and hashes are retained.
2. `scripts/02_aggregate_reference.R`: frozen full-atlas raw counts → sparse sample×type counts; no reclustering or score recomputation. Main atlas uses qs2 format; frozen stromal object uses legacy qs format.
3. `scripts/03_cohort_audit.py`: raw metadata, endpoint/stage definitions and actual platform gene coverage → tables/; expression matrices remain compressed source files.
4. `scripts/04_validate_MuSiC.R`: unchanged official MuSiC numeric kernels; patient/dataset isolation, fixed validation plan, all predictions and accuracy metrics. `Rlib/nnls` is the isolated dependency; `resources/MuSiC_*.R` and license are preserved.
5. `scripts/05_validation_plots.R`: four independent plots (PDF + 350 dpi PNG).
6. `scripts/06_final_checks.R`: numerical invariants, patient independence of DEG sets, original proportion definition and stop-gate compliance.

Phase 1B scripts and report are in `../phase1B_review/`. Original Phase 1 files were preserved byte-for-byte. The definitive absolute-path artifact list is `../summary/phase1B_phase2_file_index.tsv`; figure links are in `../summary/phase1B_phase2_figures.md`.

No clinical pipeline should use `MuSiC_full_reference_statistics.rds` as a validated abundance reference. It is retained to inspect the failed reference and reproduce validation only.
