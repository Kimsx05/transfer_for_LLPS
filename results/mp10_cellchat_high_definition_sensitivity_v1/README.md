# MP10 High-definition sensitivity (fine CellChat)

## Run and redraw

Run from the project directory:

```bash
bash mp10_cellchat_high_definition_sensitivity_v1/run_all.sh
```

The script uses three CellChat worker processes with one BLAS/OMP/MKL thread each. Analysis and rendering are sequential outside inference. During this delivered run, one additional single-threaded preparation/plot/QC process was used alongside the three inference workers (at most four active computational threads).

To redraw solely from result tables (no CellChat, Seurat or expression access):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
Rscript mp10_cellchat_high_definition_sensitivity_v1/scripts/04_replot.R
```

To reaggregate compact unfiltered probability/P caches without inference:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript mp10_cellchat_high_definition_sensitivity_v1/scripts/03_summarize.R
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript mp10_cellchat_high_definition_sensitivity_v1/scripts/06_evidence.R
```

R libraries are the existing project R 4.5/4.4 libraries. No software or database was upgraded. The delivered run completed 21 networks, 61 figures (244 files), 2,476 analysis checks and 61 figure/table consistency checks. All PDF files passed independent text rendering and vector checks; PNG files have 300 dpi metadata. Contact sheets and representative full-size figures were visually reviewed.

Plot helpers are copied from the installed `jins-routine-plotting` skill; actual Arial, Cairo PDF and ragg PNG are required. The new analysis scripts and small C++ accelerator are included. Original expression and baseline outputs remain in place; no large expression or full CellChat object is duplicated.

## Frozen population and sampling

The baseline is `c2lfinal_20261108T123140Z_0b5fd64f`. Its exact final-reference membership supplies the 30 final cell labels and major groups. Original eligibility is `epi_state_eligible=TRUE`, Tumor, yielding 40 specimens and 116,448 Epi. `MP10_cell`, the Top20 genes and the original score definition are unchanged. Spatial expression, spot scores, inferred spatial abundances and `sample_MP10_group` do not define these states.

R `quantile(type=7)` is computed once from every eligible Epi before sampling. L50 is `<=Q50`; H50/H10/H05 are respectively `>Q50/>Q90/>Q95`. Ties are never split. C10 requires each specimen's H05 and L50 counts to exceed 10; C30 requires both to exceed 30. C10 contains 38 specimens (37 recorded patient identifiers); C30 contains 30 specimens (30 recorded identifiers). These identifiers are derived/provisional, not verified independent biological patients. BC6 and BC8 are excluded from C10 by H05 counts of 8 and 7.

A single R Mersenne-Twister random priority (`seed=123`, Inversion, Rejection) is assigned to original Epi sorted by `sample_uid, barcode`, independently of their score. Within each pool the smallest priorities are selected. Per specimen `k_s=min(100,n_H05,n_L50)`. All three definitions have the same k_s and exactly identical Low barcodes/order. Their High memberships can overlap naturally; H05 members are not forced into H50/H10. C10 contains 2,651 cells per Epi state and C30 2,479 per state. Different k_s values mean specimens are not exactly equally weighted in inference.

TME cells are the original baseline ALL fine capped sample members, including pre-inference members belonging to small nodes. They are shared once, in baseline order. Dataset inputs are ordered subsets of the unified C10 ALL members. C30 is an ordered specimen subset of C10; nothing is resampled. `tables/membership` preserves candidate pools, priorities, selection order, full input order and modeled order. All excluded/unused Epi stay outside the network; no EpiOther node is created.

## Expression, domain and inference

All calculations use nonnegative single-cell normalized RNA:data. `raw.use=TRUE` means CellChat uses its unprojected signaling expression, not raw UMI counts. The old cache lacked IFNLR1 and IL2RA; both were extracted from the original same-normalization RNA:data and saved in a tiny two-gene supplement. No library size or normalization was recomputed. Original full features, required/auxiliary coverage and actual all-cell expression-zero status are recorded separately.

The original CellChat 1.6.1 human DB hash matches baseline. Of 1,939 definitions, 1,933 have every mandatory ligand/receptor/complex subunit available. Six definitions are not computable due to missing original RNA features: HLA-C_KIR2DL2, HLA-C_KIR2DS1, HLA-F_KIR3DS1, RAET1F_KLRK1, HLA-DRB3_CD4, HLA-DRB4_CD4. They are outside the inference domain, never zero-filled. Native 1.6.1 cofactor handling intersects with available genes; its handling and per-LR auxiliary coverage are saved. Features present with zero expression remain in the fixed domain.

No condition-specific overexpression filtering is run. The complete frozen LR data.frame is passed via `LR.use` and assigned to `object@LR$LRsig`; array dimension 3 is asserted identical to its row order. All networks use triMean, raw.use=TRUE, population.size=FALSE, nboot=100, seed.use=123. The numerically verified original type-7 triMean C++ acceleration is retained. All other inference follows the installed 1.6.1 function. CellChat's data_max is recorded per network. Absolute W values across independently scaled networks are not calibrated measurements.

The legacy 1,875 LR identifiers are read from the actual ALL fine domain table; all are in the new domain. Secondary summaries restrict the very same inference to this intersection. They do not rebuild networks or reconstruct the domain from supported rows. File, database, full-feature, LR, actual expression-content and ordered membership SHA256 hashes are saved. Cache reuse requires matching actual members, labels, parameters, accelerator and expression/domain hashes, not a matching filename.

## Two different count checks

Modeling nodes have total sampled count >10 before inference, using identical nodes across the three High definitions. Each network is inferred once. Postfilters min.cells=10 and 30 remove nodes with total counts <=threshold from results; they do not remove those cells and refit. Compact caches preserve all original probabilities and internal P values on modeled nodes.

C30 is different: it changes specimen eligibility, then refits three ALL networks using the C10-frozen cell subsets. C30 dataset coverage is delivered without additional dataset inference. Min10/min30 and C10/C30 comparisons explicitly recompute backgrounds on their common evaluable nodes.

All supported results use `probability>0 & internal_P<0.05`. Installed 1.6.1 pathway aggregation includes P=0.05 by default, whereas overall aggregation excludes it. Our derived pathway and overall tables both use the strict mask. Original arrays are never modified. Internal permutation P is not a specimen-level High/Low differential test; cells, nested cohorts, schemes and postfilters are not independent biological replicates.

## Effects, backgrounds and semantics

For both Epi→TME (output) and TME→Epi (input), all 29 TME targets are retained. W sums supported LR probabilities, R sums all raw probabilities, and C counts supported LR. Raw differences and ratios are reported alongside background-adjusted effects. For target j, the background contains all other nodes in its explicitly named evaluable TME universe, including computed zero-valued edges:

`r_B=sum(H_background)/sum(L_background)`; `RR=(H/L)/r_B`; `E=H-r_B*L`.

There are no pseudocounts. Low=0 makes the corresponding ratio undefined, but its computed H/L values stay intact and positive High-only supported links are flagged. No background, background Low=0, background High=0, insufficient nodes and unmodeled nodes have distinct statuses. Only finite RR>1 is called relative enhancement; RR=1 is background-like. Ranking uses W RR, with R as support-threshold sensitivity.

`background_id` hashes the ordered universe before removing the target; each row also stores actual target-excluded background members and both sums. `background_comparison_map.tsv` maps native, min10/min30, C10/C30, and ALL/dataset comparable contexts to these identities. Repeated contexts with identical universes share one effects record. Keys are cohort, scope, high, min_cells, lr_domain, background_id, direction, target and metric (plus pathway/LR when applicable). Context tables add a comparison_context column.

`E_path=H_path-r_B*L_path` uses the overall target-specific background ratio, so pathway contributions sum to E. Focus candidates are the predeclared six plus positive ALL top3 unions. Up to five positive H05 C10 ALL W pathways are selected once per target/direction and displayed for all three definitions. Local positive pathways do not override a nonenhanced overall RR. LR tables explicitly identify ligand sender and receptor receiver; interaction names alone do not establish upstream biological causation.

## Files and zero recovery

- `tables/specimen_quantiles_coverage.tsv`: all original 40 specimens, descriptive distribution statistics, actual pool fractions, k_s and exclusion reasons.
- `tables/epi_full_scores_pools.tsv.gz`: original Epi score/pool/random-priority trace.
- `tables/network_plan.tsv`, `network_status.tsv`, `node_counts.tsv`, `membership/`: input order/hash, runtime completion and modeled coverage.
- `config/domain_provenance.json`, `tables/LR_coverage_audit.tsv`, `fixed_LR_domain.tsv`, `gene_universe_coverage.tsv.gz`, `DB_*`: expression and domain provenance.
- `objects/networks/*.qs`: only unfiltered numeric prob/P, counts, LR metadata and hashes, no expression-bearing CellChat objects.
- `tables/epi_tme_edges.tsv.gz`, `native_effects_rankings.tsv.gz`, `background_adjusted_effects.tsv.gz`, `comparable_background_effects.tsv.gz`: complete W/R/C, statuses, backgrounds, effects and rankings.
- `tables/epi_tme_pathways.tsv.gz`, `pathway_background_contributions.tsv.gz`: all targets and all computed pathways, including computed zeros and NE.
- `tables/lr/*__raw_sparse.tsv.gz`: all Epi↔TME LR entries where probability!=0 or P!=1. On modeled node pairs and only within fixed_LR_domain, omitted entries are exactly probability=0 and P=1. Outside that domain, unmodeled nodes and postfiltered nodes are NE, never zero. Threshold/domain/background joins are defined by the edge and domain tables.
- `scripts/recover_lr.R`: callable table-only restoration with these rules, including postfilter and legacy-domain restrictions.
- `tables/focus_LR_background_contributions.tsv.gz`, `focus_top10_positive_LR_contributions.tsv.gz`: explicit raw/supported and overall-background LR contributions.
- `tables/ALL_high_definition_changes.tsv`, `dataset_direction_support.tsv`, `C10_C30_same_background_comparison.tsv`, `LR_domain_sensitivity.tsv`: descriptive robustness checks.
- `figures/`, `plot_data/`, `config/plot_config.json`, `config/figure_manifest.tsv`: vector PDF, 300 dpi PNG, separate panel and legend PDFs where applicable, each figure's direct data and fixed color scales.
- `qc/`: input, native-boundary, support, summation, sparse-restoration, numerical and rendered-output checks.

Distribution ECDFs preserve zero mass and use identical score axes. Dataset ECDFs average specimen ECDFs with equal specimen weights. Heatmaps share log2(RR) limits [-1.5,1.5], which cover all finite effects in this run without clipping. Pathway facet labels explicitly mark overall RR<=1 and retain three decimal places near the no-enrichment boundary. A cell labeled `0` is a true RR zero; `NE` is not evaluated; `L0` means target Low=0, `B0` undefined background ratio, and `NB` no valid background. Exact numeric effects remain in direct plot data. No arbitrary pseudocount makes these cells finite.

## Interpretation limits

The three High definitions are prespecified sensitivity analyses and are all retained. H05 is a high-score tail, not a validated active-cell class; the remaining 95% are not presumed noise. Inference uses pooled cells and fixed triMean; it can reveal distribution-dependent communication estimates, not establish a biological activity boundary. Recorded patient relationships, sample composition, small-node coverage, LR-domain changes, internal support and per-network scaling are considered separately in the accompanying results report. No final candidate or activity threshold is frozen automatically.

Direct plot TSV files can contain quoted multiline display labels. Read them with an explicit tab delimiter, e.g. `data.table::fread(path, sep="\t", header=TRUE)` or `pandas.read_csv(path, sep="\t")`.
