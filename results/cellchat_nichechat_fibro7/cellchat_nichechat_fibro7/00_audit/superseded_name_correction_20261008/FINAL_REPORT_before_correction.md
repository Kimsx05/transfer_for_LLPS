# Fibro7 cohort communication report

CELLCHAT_STATUS = PARTIAL
NICHECHAT_STATUS = AUDIT_ONLY
OVERALL_STATUS = PARTIAL

## Exact frozen inputs

- Full Seurat atlas: `/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/results/objects/major_annotated_with_LLPS77_AUC.qs` (QS2 format despite .qs suffix).
- Final frozen fibro assignment: `/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order/global_res012_object.qs` (standard QS format).
- New copy with cellchat_group: `/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/01_metadata/major_frozen_with_cellchat_group.qs2`.
- Frozen inputs and historical outputs are preserved. Existing cell identities, C0–C6, LLPS, major annotations, reductions, clustering and QC were not changed. Input SHA256 and metadata/reduction preservation checks passed.

| metric | value |
| --- | --- |
| Seurat_version | 5.5.1 |
| SeuratObject_version | 5.4.0 |
| CellChat_version | NOT_INSTALLED |
| total_cells | 353623 |
| biological_samples | 59 |
| patient_number | UNKNOWN: no patient field |
| biological_sample_field | orig.ident |
| patient_uid_field | ABSENT |
| dataset_field | dataset |
| major_field | major_celltype |
| fibro_subtype_field | subcluster in frozen fibro object |
| RNA_assay_structure | Assay5; counts/data dgCMatrix 71683x353623; scale.data matrix 3000x353623 |
| normalized_layer | RNA/data |
| normalization | LogNormalize |
| scale_factor | 10000 |
| primary_included_cells | 344377 |
| excluded_unresolved_cycling | 7324 |
| excluded_previous_fibro_decisions | 1922 |

## Final groups and eligibility

Fibro_C0, Fibro_C1, Fibro_C2, Fibro_C3, Fibro_C4, Fibro_C5, Fibro_C6, Epithelial, Endothelial, Myeloid, T_NK, B, Plasma, Mast
- Final 14 groups only; generic Fibro/Pericyte are absent. Cycling lineage unresolved: 7,324 excluded. Frozen earlier stromal exclusions: 1,922 inherited exactly from historical barcode lists. Final fibro cells: 20,146.
- 59 biological samples in 5 datasets. orig.ident is one-to-one with dataset; no sample identifiers collide across datasets. Biological sample replication is retained in every input.
- No patient_uid or patient identifier is present. Patients are UNKNOWN; long-table patient_uid values are NA. Sample-level replication cannot establish patient independence.

| final_cellchat_group | total_cells | n_samples_present | n_samples_eligible | max_sample_cells |
| --- | --- | --- | --- | --- |
| Fibro_C0 | 10717 | 43 | 20 | 2774 |
| Fibro_C1 | 5280 | 39 | 22 | 1715 |
| Fibro_C2 | 1305 | 47 | 10 | 443 |
| Fibro_C3 | 913 | 39 | 7 | 380 |
| Fibro_C4 | 811 | 39 | 7 | 159 |
| Fibro_C5 | 716 | 37 | 10 | 109 |
| Fibro_C6 | 404 | 25 | 3 | 168 |
| Epithelial | 130804 | 57 | 49 | 9601 |
| Endothelial | 8163 | 54 | 33 | 1497 |
| Myeloid | 34536 | 59 | 55 | 2398 |
| T_NK | 121222 | 55 | 50 | 6149 |
| B | 18622 | 51 | 39 | 4610 |
| Plasma | 8374 | 57 | 30 | 1884 |
| Mast | 2510 | 51 | 20 | 670 |

## CellChat method and execution

- Parameters: min.cells=30; type=truncatedMean; trim=0.1; population.size=FALSE; nboot=100; seed.use=123. RNA/data LogNormalize, scale.factor=10000. No LLPS/reduction input.
- Database: CellChatDB.human; Secreted Signaling + ECM-Receptor + Cell-Cell Contact. All samples use the same DB, standard pipeline and parameters.
- Approved-installation target: official CellChat v2.1.2, commit a0d3b2d231d46c8787177fffeac908270c253747. Actual installed/inference version is recorded in the preflight/API/session audit after installation.
- Attempted samples=0; SUCCESS=0; FAILED=0; not attempted=59.
- Inference has NOT run: CellChat and four direct dependencies are missing. User installation approval is pending. Audit/input preparation and consensus semantic checks are completed; prepared inference/figure scripts have only been syntax checked.

## Explicit sample-edge consensus

- source and target >=30 cells -> structurally eligible. Undetected communication in a successful eligible sample -> weight=0. Ineligible -> NA. Inference failure -> NA plus FAILED, never weight=0.
- For complete evaluable edges: n_eligible=sum(eligible); n_detected=sum(weight>0 over eligible); prevalence=n_detected/n_eligible; consensus_strength=median(all eligible weights including real zeros). QC includes mean, IQR, median nonzero and maximum.
- Primary display rule is exactly n_eligible>=3 and prevalence>=0.30. Non-retained display entries are zero; all audit rows remain. Edges passing this rule can still have median=0 and thus no visible line.
- If an eligible sample failed, complete-cohort median/prevalence are unknown. Such edges retain n_missing_eligible, observed median, and prevalence bounds; incomplete edges are not displayed. This avoids both false zeros and silently redefining the eligible denominator.
- Structural eligible sample-edge observations=2523.
- Retained edges, prevalence distribution and strongest signals are NOT AVAILABLE before inference.

## Strongest overall relationships

NOT AVAILABLE: inference not completed.

## Strongest fibro outgoing relationships

NOT AVAILABLE: inference not completed.

## Strongest fibro incoming relationships

NOT AVAILABLE: inference not completed.

## NicheChat implementation and cross-method comparison

- NicheChat implementation NOT IDENTIFIED; exact method, parameters, input compatibility, LR output and per-sample API are unknown. See 05_nichechat_method_audit/NicheChat_method_audit.md.
- NicheChat inference samples=0; coverage=0; inference failures=0 (not attempted). No method is substituted under the NicheChat name.
- Concordant signals: NOT EVALUABLE. Discordant signals: NOT EVALUABLE. No BOTH_SUPPORTED/CELLCHAT_ONLY/NICHECHAT_ONLY classification is inferred from an unrun method.
- Figures C and D cannot be generated without a validated NicheChat implementation. If CellChat inference completes first, its feature table can be supplied with NicheChat columns=NA and explicit NOT_EVALUABLE status.

## LLPS descriptive integration and sensitivity

- Frozen fibro subtype LLPS77_AUC summary is joined only AFTER communication inference. No score recomputation, subtype reassignment or cell-level LLPS/communication significance test.
- Primary strength remains raw per-sample net$weight. Sensitivity divides sample weights by that sample total eligible network weight; then repeats zero-inclusive aggregation and checks ranking, fibro targets, topology and relative C0–C6 patterns.
- Sensitivity results and LLPS integration are not yet computed.

## Limitations and coverage flags

- Fibro_C6 is eligible in only 3 of 59 samples; C3/C4 in 7 each. All frozen subtypes remain unchanged. Coverage flags are descriptive, not automatic deletion/merging rules.
- The historical marker audit flags C1/C2 identity uncertainty and C6 cycling state. Communication inference cannot repair identity uncertainty, confirm spatial contact, or establish causality.
- Patient mapping is absent; repeated specimens from a patient cannot be audited. Dataset effects and differing depth may affect probability scales. Normalized-weight sensitivity addresses scale, without changing the primary definition.
- Unrun or failed analyses are not treated as negative biological evidence. No Normal/NMIBC/MIBC analysis was performed.

## Final artifact inventory

- Every existing final table/figure/audit and its absolute path is listed in 00_audit/final_artifact_manifest.tsv. Scripts, per-sample expression inputs and downloaded package source are listed separately in 00_audit/all_output_file_manifest.tsv.
- Primary Figure A: NOT GENERATED.
- Supplemental circle: NOT GENERATED.
- Figure B: NOT GENERATED.
- Figures C/D: NOT GENERATED (NicheChat audit only).

CELLCHAT_STATUS = PARTIAL
NICHECHAT_STATUS = AUDIT_ONLY
OVERALL_STATUS = PARTIAL
