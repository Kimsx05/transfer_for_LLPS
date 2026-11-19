# C5–Endothelial NicheNet spatial v1

Status: analysis complete; RESPONSE_MODE=CURATED_EXPLORATORY. Formal C5-specific spatial localization NO-GO retained.

Start with COMPREHENSIVE_REPORT.md, FIGURE_GALLERY.html, tables/final_axis_summary.tsv and 06_CellChat_integration/complete_expression_activity_LR_target_evidence.tsv.

21 tumor patients meet frozen sampling criteria, but no qualifying positive endothelial response genes. Three standard programs evaluated; activation term too small. 246 broad ligands; five main axes and three ECM controls; all63 sections and265002 spots. Only2 sections achieve full matching balance (one known patient); no spatial-translation test survives global FDR.

Reproduction: run from /home/jinshengxi/bladder_metabolism using scripts/run_R.sh or scripts/run_python.sh. See numerical scripts and frozen protocol; 03d_GO extends the two Hallmark programs once, before spatial selection. 04_spatial.py resumes completed sections using DONE files. 09_balance_boundary_fix.py documents the identical-constant-covariate SMD correction (current04 already contains the corrected definition). 10/11 are standard ECM-membership program sensitivities. 07_report.py generates report; 12_finalize.py appends interpretations and builds gallery; do not repeatedly append without regenerating07 first.

Source checksums, resource versions and validation are in00_input_audit. Raw objects and historical results are read-only. No CellChat, deconvolution, Five scoring or global spatial-statistic rerun. Seed123; serial BLAS/R/Python with affinity capped at4 CPUs. Local existing human v2 NicheNet resources only.

All figures use project plotting skill styles, vector PDF and300dpi PNG. Full/panel/legend PDFs are supplied when there is a legend. Tables provide plotting values. Raw spatial expressions are log1p(CP10k); maps explicitly show within-section percentiles. No causal or cell-intrinsic interpretation.
