# Fibro7 cohort communication report

CELLCHAT_STATUS = PASS
NICHECHAT_STATUS = NOT_APPLICABLE
NICHENET_STATUS = READY
OVERALL_STATUS = PASS

Status scope: this round requires completed sample-replicated CellChat plus READY NicheNet pre-flight. NicheNet ligand activity is intentionally deferred and is not claimed as PASS.

## Exact frozen inputs

- Full Seurat atlas: `/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/results/objects/major_annotated_with_LLPS77_AUC.qs` (QS2 format despite .qs suffix).
- Final frozen fibro assignment: `/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order/global_res012_object.qs` (standard QS format).
- New copy with cellchat_group: `/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/01_metadata/major_frozen_with_cellchat_group.qs2`.
- Frozen inputs and historical outputs are preserved. Existing cell identities, C0–C6, LLPS, major annotations, reductions, clustering and QC were not changed. Input SHA256 and metadata/reduction preservation checks passed.

| metric | value |
| --- | --- |
| Seurat_version | 5.5.1 |
| SeuratObject_version | 5.4.0 |
| CellChat_version | 2.1.2 |
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
- 59 biological samples in 5 datasets. Each orig.ident maps to exactly one dataset; no sample identifiers collide across datasets. Biological sample replication is retained in every input.
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
- Export-only sample repairs from complete inference objects: 28. Original failed-attempt status/log records are retained.
- Database: CellChatDB.human; Secreted Signaling + ECM-Receptor + Cell-Cell Contact. All samples use the same DB, standard pipeline and parameters.
- Installed official source: official CellChat v2.1.2, commit a0d3b2d231d46c8787177fffeac908270c253747. Actual installed/inference version is recorded in the preflight/API/session audit after installation.
- Attempted samples=59; SUCCESS=59; FAILED=0; not attempted=0.
- Per-sample complete objects, runtimes, errors and warnings are listed in 02_cellchat_per_sample/sample_run_status.tsv. Early data.table export/runner errors were retained in attempt archives and repaired from complete native objects; no inferential parameters were changed and no failed inference was recoded as zero. The first sample was rerun once before object preservation was added; subsequent export repairs reused the complete native objects.

## Explicit sample-edge consensus

- source and target >=30 cells -> structurally eligible. Undetected communication in a successful eligible sample -> weight=0. Ineligible -> NA. Inference failure -> NA plus FAILED, never weight=0.
- For complete evaluable edges: n_eligible=sum(eligible); n_detected=sum(weight>0 over eligible); prevalence=n_detected/n_eligible; consensus_strength=median(all eligible weights including real zeros). QC includes mean, IQR, median nonzero and maximum.
- Primary display rule is exactly n_eligible>=3 and prevalence>=0.30. Non-retained display entries are zero; all audit rows remain. Edges passing this rule can still have median=0 and thus no visible line.
- If an eligible sample failed, complete-cohort median/prevalence are unknown. Such edges retain n_missing_eligible, observed median, and prevalence bounds; incomplete edges are not displayed. This avoids both false zeros and silently redefining the eligible denominator.
- Structural eligible sample-edge observations=2523.
- Rule-retained edges=182; positive edges drawn=182; incomplete edges=0.

## Strongest overall relationships

| source | target | n_eligible | n_detected | prevalence | consensus_strength |
| --- | --- | --- | --- | --- | --- |
| Fibro_C5 | Fibro_C6 | 3 | 3 | 1 | 6.08067577935694 |
| Fibro_C6 | Fibro_C5 | 3 | 3 | 1 | 5.86299115831768 |
| Fibro_C5 | Fibro_C1 | 9 | 9 | 1 | 5.76827723055659 |
| Fibro_C5 | Fibro_C5 | 10 | 10 | 1 | 5.64355737707704 |
| Fibro_C5 | Fibro_C4 | 4 | 4 | 1 | 5.49134815395411 |
| Fibro_C1 | Fibro_C6 | 3 | 3 | 1 | 5.45771835944738 |
| Fibro_C6 | Fibro_C1 | 3 | 3 | 1 | 5.36351657944872 |
| Fibro_C6 | Fibro_C6 | 3 | 3 | 1 | 5.31888635508987 |
| Fibro_C1 | Fibro_C5 | 9 | 9 | 1 | 5.05196285890301 |
| Fibro_C4 | Fibro_C5 | 4 | 4 | 1 | 4.89111894428053 |
| Fibro_C5 | Fibro_C3 | 5 | 5 | 1 | 4.63761758811123 |
| Fibro_C0 | Fibro_C4 | 5 | 5 | 1 | 4.56670794801729 |

## Strongest fibro outgoing relationships

| source | target | n_eligible | n_detected | prevalence | consensus_strength |
| --- | --- | --- | --- | --- | --- |
| Fibro_C5 | Fibro_C6 | 3 | 3 | 1 | 6.08067577935694 |
| Fibro_C6 | Fibro_C5 | 3 | 3 | 1 | 5.86299115831768 |
| Fibro_C5 | Fibro_C1 | 9 | 9 | 1 | 5.76827723055659 |
| Fibro_C5 | Fibro_C5 | 10 | 10 | 1 | 5.64355737707704 |
| Fibro_C5 | Fibro_C4 | 4 | 4 | 1 | 5.49134815395411 |
| Fibro_C1 | Fibro_C6 | 3 | 3 | 1 | 5.45771835944738 |
| Fibro_C6 | Fibro_C1 | 3 | 3 | 1 | 5.36351657944872 |
| Fibro_C6 | Fibro_C6 | 3 | 3 | 1 | 5.31888635508987 |
| Fibro_C1 | Fibro_C5 | 9 | 9 | 1 | 5.05196285890301 |
| Fibro_C4 | Fibro_C5 | 4 | 4 | 1 | 4.89111894428053 |
| Fibro_C5 | Fibro_C3 | 5 | 5 | 1 | 4.63761758811123 |
| Fibro_C0 | Fibro_C4 | 5 | 5 | 1 | 4.56670794801729 |

## Strongest fibro incoming relationships

| source | target | n_eligible | n_detected | prevalence | consensus_strength |
| --- | --- | --- | --- | --- | --- |
| Fibro_C5 | Fibro_C6 | 3 | 3 | 1 | 6.08067577935694 |
| Fibro_C6 | Fibro_C5 | 3 | 3 | 1 | 5.86299115831768 |
| Fibro_C5 | Fibro_C1 | 9 | 9 | 1 | 5.76827723055659 |
| Fibro_C5 | Fibro_C5 | 10 | 10 | 1 | 5.64355737707704 |
| Fibro_C5 | Fibro_C4 | 4 | 4 | 1 | 5.49134815395411 |
| Fibro_C1 | Fibro_C6 | 3 | 3 | 1 | 5.45771835944738 |
| Fibro_C6 | Fibro_C1 | 3 | 3 | 1 | 5.36351657944872 |
| Fibro_C6 | Fibro_C6 | 3 | 3 | 1 | 5.31888635508987 |
| Fibro_C1 | Fibro_C5 | 9 | 9 | 1 | 5.05196285890301 |
| Fibro_C4 | Fibro_C5 | 4 | 4 | 1 | 4.89111894428053 |
| Fibro_C5 | Fibro_C3 | 5 | 5 | 1 | 4.63761758811123 |
| Fibro_C0 | Fibro_C4 | 5 | 5 | 1 | 4.56670794801729 |

## NicheNet pre-flight and planned mechanism prioritization

- The user corrected NicheChat to NicheNet/nichenetr on 2026-10-08. NICHECHAT_STATUS=NOT_APPLICABLE; previous name-audit requirements are revoked. Historical audit files are archived, not used as computation prerequisites.
| metric | value |
| --- | --- |
| package_version | 2.2.1.1 |
| source_commit | 66f90d5eeafef280b2b2f339b3fd70ffec1781dd |
| resource_record | https://zenodo.org/records/7074291 |
| ligand_target_dimensions | 33354x1226 |
| LR_unique_pairs | 4986 |
| geneset_oi | UNDEFINED |
| ligand_activity_run | FALSE |

| resource_component | n_resource_symbols | n_exact_RNA_matches | resource_match_fraction |
| --- | --- | --- | --- |
| ligand_target_targets | 33354 | 23248 | 0.69700785512982 |
| LR_ligands | 1226 | 1181 | 0.963295269168026 |
| LR_receptors | 1067 | 1043 | 0.977507029053421 |

- Official saeyslab/nichenetr is used for ligand prioritization against a receiver transcriptional program; no second 14x14 network is constructed.
- Current NicheNet state: READY. See 05_nichenet_preflight/NicheNet_preflight_audit.md for installed version, commit, resource checksum and symbol compatibility.
- Sender: frozen Fibro_C0–C6. Candidate receivers: Epithelial, Myeloid, Endothelial, T_NK and other major populations with stable CellChat support.
- Per-sample sender ligand and receiver receptor expression are audited from the unchanged RNA/data layer. Expression is detection>0; sample/group eligibility stays >=30; no expression cutoff is silently promoted to biological proof.
- geneset_oi is NOT DEFINED. Ligand activity, ligand-target inference and cross-method high-confidence mechanism calls are NOT RUN. A documented, sample-aware receiver program is required first; no stage comparison or arbitrary top-DEG set is used.
- CellChat-prioritized fibro→receiver pairs are saved as follow-up candidates, with their exact LR complex definitions and expression gates. This is preparation, not evidence of NicheNet ligand activity.
- Expected future candidate evidence: stable CellChat support + sender ligand expression + receiver receptor expression + high NicheNet activity + interpretable downstream receiver targets. These evidence types remain separate.

## LLPS descriptive integration and sensitivity

- Frozen fibro subtype LLPS77_AUC summary is joined only AFTER communication inference. No score recomputation, subtype reassignment or cell-level LLPS/communication significance test.
- Primary strength remains raw per-sample net$weight. Sensitivity divides weights by each sample eligible total, omitting the single zero-total sample from the normalized denominator only; its raw zeros remain in the primary network.
- Sensitivity materially qualifies subtype ranking: raw total outgoing/incoming strengths place C5 first; normalized outgoing places C0 first and C5 third, with C1/ C5 epithelial–endothelial target ordering changes. Do not claim unconditional C5 dominance. Epithelial and Endothelial remain key fibro receiver relationships; retain C0 and C5 in later receiver-program design.
- Sensitivity report: `/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/03_cellchat_consensus/CellChat_consensus_sensitivity.md`.

## Limitations and coverage flags

- Fibro_C6 is eligible in only 3 of 59 samples; C3/C4 in 7 each. All frozen subtypes remain unchanged. Coverage flags are descriptive, not automatic deletion/merging rules.
- The historical marker audit flags C1/C2 identity uncertainty and C6 cycling state. Communication inference cannot repair identity uncertainty, confirm spatial contact, or establish causality.
- Patient mapping is absent; repeated specimens from a patient cannot be audited. Dataset effects and differing depth may affect probability scales. Normalized-weight sensitivity addresses scale, without changing the primary definition.
- Unrun or failed analyses are not treated as negative biological evidence. No Normal/NMIBC/MIBC analysis was performed.

## Final artifact inventory

- Final validation: 00_audit/FINAL_VALIDATION.txt.
- Every existing final table/figure/audit and its absolute path is listed in 00_audit/final_artifact_manifest.tsv. Scripts, per-sample expression inputs and downloaded package source are listed separately in 00_audit/all_output_file_manifest.tsv.
- Primary Figure A: 04_cellchat_figures/Figure_A_cohort_consensus_circle.{pdf,png,svg}
- Supplemental circle: 04_cellchat_figures/Figure_A_fibro_centered_circle.{pdf,png,svg}
- Figure B: 09_figures/Figure_B_fibro_major_consensus_heatmap.{pdf,png,svg}
- NicheNet ligand-activity/mechanism figures: DEFERRED (receiver geneset_oi required).

CELLCHAT_STATUS = PASS
NICHECHAT_STATUS = NOT_APPLICABLE
NICHENET_STATUS = READY
OVERALL_STATUS = PASS
