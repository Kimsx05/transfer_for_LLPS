# C5-associated spatial exploration (all 63 sections)

Exploratory application authorized after the v2 formal localization NO-GO. Scores are mixed-spot transcriptional signatures, never C5 cell abundance. Historical results and raw inputs are read-only. No AUROC gate is applied.

## Frozen analysis choices
- seed 123; CPU affinity limited to four available CPUs; BLAS/OpenMP/NumExpr and R/BiocParallel serial. Each worker reads one section at a time; at most three internally serial section workers.
- Frozen A26 and core15 from v2; original five and marker3 from v1; MMP11 is a single-marker diagnostic. Lists in `00_inputs_and_QC`. C2C4 program uses the previously defined single-cell reference contrast formula, fixed before spatial scoring; it cannot distinguish C2 from C4.
- Full section gene universe, UCell maxRank=1500, AUCell ceil(0.05*number of genes), keepZeroesAsNA=TRUE, normAUC default TRUE. Only present signature genes are scored; missing genes and changed denominators are explicitly recorded. Never mix absent features into zeros without disclosure. Raw UCell/AUCell and within-section midrank percentiles retained.
- Pilot: the original 7 sections selected by depth representativeness, solely technical acceptance (IDs, scores, prior score concordance), followed by all remaining sections.
- Array adjacency: six Visium hex offsets (0,+/-2),(+/-1,+/-1), only observed spots. No connections across missing spots. No physical distance conversion.
- Moran's I: binary undirected graph, 199 within-section label permutations, one-sided positive autocorrelation, BH across 63 sections separately for each score. Permutations describe spatial pattern, not patient-level biological inference. Raw A26 and rank-residual A26 (Fibroblast, nCount, nFeature, Epithelial adjustment) assessed.
- High-score sets: strictly greater than section 90th percentile; ties can reduce size. Sensitivity 80th/95th percentiles. Components and actual selected fraction reported. Operational clustering flag: positive Moran I with BH<0.05 and at least 10 selected spots. This is not proof of a discrete niche.
- Fibroblast and depth: Spearman and rank partial correlations; matched high vs non-high controls using rank Fibroblast and rank UMI, caliper 0.1 on each, greedy one-to-one without replacement. Balance and match rate reported; no spot-level biological P values.
- Epithelial regions: major mean fraction >=0.3, sensitivity >=0.2/0.4. Fractions refer to cell2location major means, not malignant cell labels. Graph distance rings 0,1,2,3,4+; disconnected components without epithelial seeds remain missing. q05 support is summarized, not interpreted as a calibrated posterior probability.
- Endothelial/Myeloid: graph rings around high-score and matched-control seed sets, descriptive mean log1p abundance contrasts. Rings may overlap across conditions; this is exploratory association, not independent cell interaction evidence.
- Independent unit: aggregate effects within known patient and dataset; unknown patient sections remain descriptive. Dataset-specific results and leave-one-dataset-out known-patient summaries; no pooled spot P values. Tissue labels only from existing audited workbook mapping.
- Coordinate maps: original low-resolution pixel coordinates, upper origin, no rotation/mirroring/rescaling. HE has provenance/geometric checks only, not independent histological validation. Primary maps without HE; auxiliary HE view explicitly labelled provisional.
- Raw signature maps fixed [0,1] across sections, percentiles separately labelled. Abundance/depth/expression use shared global 99th-percentile upper bounds, explicitly clipped. Spatial colors are not evidence of calibrated cross-dataset score comparability.

## Reproduce / resume
Run `bash scripts/run_python.sh scripts/score_all.py pilot`, then `... score_all.py remaining` from the project root using full relative script paths. Completed section Parquets are reused. Subsequent scripts produce spatial statistics and figures without retaining duplicate expression matrices. Logs and source hashes are retained.

The complete resumable entry point is `bash results/C5_spatial_exploratory_all63_v1/scripts/run_all.sh` from the project root. Scoring batches use up to three independent, internally serial section workers within the same four-CPU affinity; later stages run sequentially. SciPy Matrix Market writer parallelism is explicitly set to one. No nested algorithmic parallelism is used. Completed Parquets publish atomically. The combined C2/C4 program is a single-cell-derived contrast, not a validated competitor cell label.

For contribution diagnostics, A26-high UCell rank weights max(1500 - average expression rank, 0) are summarized by gene; these are score weights, not RNA fractions or causal effects. No signature is modified in response to them.

## Additional diagnostic sensitivity
The initial Fibroblast/UMI matching left epithelial composition imbalance. A documented post-diagnostic sensitivity adds epithelial abundance as a third rank-matching variable (same 0.1 per-variable caliper), and computes distance partial correlations with/without epithelial adjustment. This is not a preplanned confirmatory test or signature tuning. All matched sections and an explicitly defined balanced subset (match rate >=0.5; all three absolute SMD <=0.1) are reported; neither subset is selected for favorable effects.

Primary coordinate maps crop the display window to observed spots with 3% padding, shared across panels of that section. HE auxiliary maps retain the full image window. Original coordinate orientation is preserved; cropping does not assert a calibrated common physical scale across sections. Dataset colors use the plotting skill's fixed seven-color palette mapped once to all seven source labels.
