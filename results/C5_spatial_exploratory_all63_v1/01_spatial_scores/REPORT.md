# All-section C5-associated scores

63 Parquet files preserve 265,002 original spot IDs, metadata, all four abundance tables, raw UCell/AUCell scores and within-section percentiles. Signatures are frozen in `../00_inputs_and_QC/frozen_signatures.json`.

UCell maxRank=1500; AUCell ceil(5% of each section's full feature universe), keepZeroesAsNA=TRUE, normAUC default TRUE; seed=123. Absent features are excluded from the signature denominator and explicitly recorded, not silently imputed as zero. A26 has 26/26 feature coverage in all sections. Single-marker MMP11 scoring is diagnostic only.

Original 7 pilot five/three-gene scores reproduce exactly. See the comprehensive report and per-section detection summary for actual detection and method limitations. These are mixed-spot transcriptional scores, not C5 cell abundance.
