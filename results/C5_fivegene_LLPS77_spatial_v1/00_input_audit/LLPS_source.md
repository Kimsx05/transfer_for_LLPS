# LLPS-77 provenance and directionality

Original file: `/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/LLPS_77genes.rdata`. Direct R inspection: one object `LLPS_77genes`, character vector length77. Exactly77 unique gene symbols. No direction, effect size, or UP/DOWN fields are attached. Absence of direction metadata is not evidence that all genes promote LLPS.

Single-cell scripts `scripts/config.R` and `scripts/01_run_AUCell_chunk.R` in that source directory use this same file/object, frozen as `senior_LLPS77`, as one unsigned AUCell gene set. The exported list exactly matches `results/tables/LLPS77_matched_genes.txt` (all77). Prior input audit used exact case-sensitive matching without symbol replacement. Literature selection rationale is not independently recoverable from the source vector; no invented provenance or experimentally validated LLPS-activity interpretation.

Prior spatial coverage report reused: `results/C5_spatial_analysis_v1/00_input_audit/LLPS77_coverage/section_coverage.tsv`. New var-only audit confirms covered symbols per slice. Features absent from an H5AD are omitted from the set denominator; measured but zero features remain in the full expression universe. Common36 sensitivity retains only46.8% of the original list: reduced-program sensitivity, not an interchangeable corrected LLPS77 score. Five/A26 have no overlap. Core15 overlaps COL1A1, hence a separate 76-gene nominal LLPS variant.

All old scores and all four abundance sets are reused through the existing per-section Parquets. Original subtype identifiability PASS_WITH_WARNING is retained. Fibroblast and Mural remain separate. No new cell2location model or inferred patient/tissue label.
