# CELL2LOCATION FINAL VALIDATION GATE REPORT

## Decision

**PASS_WITH_WARNING**

## Frozen inputs

- Final sections: 63 / 63 complete and valid
- Final reference: 300,033 cells, 17,987 genes, 30 factors
- Reference model: `/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/reference_training_v2/11_final_spatial_reference/30type_regression_model/model/model.pt`
- Estimated signatures: `/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/reference_training_v2/11_final_spatial_reference/30type_regression_model/learned_signatures_30type.csv.gz`
- RCTD unrestricted weights: 63 / 63 sections found

## Reference QC

- Median signature-vs-empirical Spearman: 0.841
- Minimum signature-vs-empirical Spearman: 0.741
- Sparse raw-count empirical factor means and frozen estimated signatures were used. No fit/train call was made.

## Identifiability

- High-similarity review pairs: 48 (screening only: Spearman >=0.80 OR cosine >=0.95 OR Top50 specificity-marker Jaccard >=0.30)
- LOW_IDENTIFIABILITY factors: 3
- INTERPRET_WITH_CAUTION factors: 24
- RARE_CONTEXT_SPECIFIC factors: 0

Top review pairs:

| factor_A | factor_B | spearman | cosine | jaccard50 |
| --- | --- | --- | --- | --- |
| C1QC TAM | SPP1 TAM | 0.9486 | 0.9474 | 0.0000 |
| CD4 Helper/Conventional | CD4 Naive/Memory-like | 0.9436 | 0.9891 | 0.0204 |
| Arterial EC | Capillary EC | 0.9364 | 0.9684 | 0.0000 |
| CD4 Helper/Conventional | Tfh/CXCL13-like | 0.9363 | 0.9867 | 0.0000 |
| CD8 GZMK/Effector-memory-like | NK | 0.9354 | 0.9878 | 0.0000 |
| CD4 Helper/Conventional | CD8 GZMK/Effector-memory-like | 0.9319 | 0.9889 | 0.0204 |
| CD4 Naive/Memory-like | CD8 GZMK/Effector-memory-like | 0.9221 | 0.9924 | 0.0417 |
| SMC | rCAF | 0.9134 | 0.9547 | 0.0000 |
| CD4 Helper/Conventional | Treg | 0.9097 | 0.9877 | 0.0000 |
| CD8 Exhausted | CD8 GZMK/Effector-memory-like | 0.9081 | 0.9880 | 0.0000 |
| CD8 Exhausted | Tfh/CXCL13-like | 0.9067 | 0.9835 | 0.0000 |
| Tfh/CXCL13-like | Treg | 0.9020 | 0.9756 | 0.0000 |
| CD4 Helper/Conventional | CD8 Exhausted | 0.9020 | 0.9871 | 0.0000 |
| CD4 Naive/Memory-like | Tfh/CXCL13-like | 0.9000 | 0.9729 | 0.0000 |
| cDC1 | cDC2 | 0.8968 | 0.9506 | 0.0000 |
| Pericyte | rCAF | 0.8877 | 0.9300 | 0.0000 |
| CD4 Helper/Conventional | NK | 0.8802 | 0.9808 | 0.0000 |
| CD8 Exhausted | NK | 0.8801 | 0.9856 | 0.0000 |
| CD8 GZMK/Effector-memory-like | Tfh/CXCL13-like | 0.8783 | 0.9746 | 0.0000 |
| SPP1 TAM | cDC2 | 0.8777 | 0.7281 | 0.0000 |

## Spatial posterior audit

Pairwise correlations were calculated only where both factors had support (each >=10 spots at abundance >=0.01 and total abundance >1). High-spot overlap used supported section-specific top 20% and top 10%. Fixed-ratio behavior used joint-support log abundance-ratio SD/IQR.

## Major-level c2l vs RCTD unrestricted

- Exact barcode matching was used.
- c2l absolute abundance was converted to per-spot proportion.
- Mural remains explicit in the mapping audit and was combined with Fibroblast only to align to the RCTD 8-major taxonomy.
- Majors with negative median Spearman: 0

| major | n_sections | median_rho | IQR_rho | fraction_rho_ge_0_3 | fraction_rho_ge_0_5 |
| --- | --- | --- | --- | --- | --- |
| B | 58 | 0.1368 | 0.2249 | 0.1379 | 0.0517 |
| Endothelial | 63 | 0.4363 | 0.2298 | 0.8095 | 0.3651 |
| Epithelial | 63 | 0.8201 | 0.2313 | 0.9841 | 0.9048 |
| Fibroblast | 63 | 0.8399 | 0.1676 | 1.0000 | 0.9524 |
| Mast | 60 | 0.1291 | 0.1220 | 0.0833 | 0.0167 |
| Myeloid | 63 | 0.4792 | 0.3344 | 0.7778 | 0.4762 |
| Plasma | 62 | 0.6097 | 0.5372 | 0.7258 | 0.6129 |
| T_NK | 59 | 0.2487 | 0.2200 | 0.4237 | 0.0678 |

## Downstream interpretation

Use `09_final_status/DOWNSTREAM_USAGE_TABLE.csv`. Recommendations are advisory; no subtype was merged or removed. High reference similarity alone was not treated as failure when discriminating markers, non-collinear posterior patterns, or context-specific recurrence supported separation.

## MP10 blinding and hard stop

MP10 was not used during cell2location validation or factor identifiability assessment.

No niche clustering, MP10 enrichment, CellChat, reference training, spatial training, annotation change, or downstream biological analysis was performed.
