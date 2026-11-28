# LLPS spatial niche clustering v1



**Analysis scope: LLPS only. C5-associated transcriptional signal is not validated C5 CAF abundance; LLPS-associated activity is not physical phase separation.**



## 1. Included material

63 sections; 265,002 spots; 7 source groups (some source identities remain UNVERIFIED). 29 confirmed patients across 30 sections. No inferred patient identities or stages. All audited spots retained.

## 2. R1/R2/R3

| radius | mean_size | mean_median | min_size | max_size |
|---|---|---|---|---|
| 1 | 6.612 | 6.984 | 1 | 7 |
| 2 | 17.31 | 18.68 | 1 | 19 |
| 3 | 32.7 | 35.71 | 1 | 37 |

| radius | mean_coverage | complete_fraction | under_half_fraction |
|---|---|---|---|
| 1 | 0.9446 | 0.814 | 0.0177 |
| 2 | 0.911 | 0.6722 | 0.0378 |
| 3 | 0.8838 | 0.5583 | 0.05764 |

Exact array offsets ensure no cross-section edges or interpolation. Squidpy grid discrepancies: 0 extra and 0 missing directed edges corrected to exact array adjacency. R2/R3 use graph-reachable cumulative neighborhoods, counting centers once. Larger radius averages more tissue and increases smoothing.

## 3. Continuous abundance versus dominant-type native method

| radius | n_niches_full | n_niches_pilot | silhouette | calinski_harabasz | davies_bouldin | scope |
|---|---|---|---|---|---|---|
| 1 | 353 | 329 | 0.3932 | 1532 | 0.9851 | full pooled native labels evaluated on balanced pilot |
| 2 | 532 | 493 | 0.5414 | 352.9 | 1.003 | full pooled native labels evaluated on balanced pilot |
| 3 | 531 | 464 | 0.3148 | 196.3 | 1.024 | full pooled native labels evaluated on balanced pilot |

| radius | n_unique_profiles | fraction_spots_in_profiles_with_at_least_30_copies | n_clusters | median_cluster_size | n_clusters_below_100_spots | fraction_spots_in_clusters_below_100 |
|---|---|---|---|---|---|---|
| 1 | 3661 | 0.9383 | 353 | 122 | 156 | 0.03271 |
| 2 | 34048 | 0.7086 | 532 | 81 | 305 | 0.06143 |
| 3 | 74576 | 0.5338 | 531 | 66 | 344 | 0.06177 |

The native dominant-type control produces hundreds of clusters in these fixed settings. This is severe overpartitioning; high silhouette in a discrete composition space is not evidence of superior ecological validity. Repeated quantized profiles and KNN ties are diagnosed explicitly. It is not recommended as a biological niche taxonomy.

| model_a | model_b | ARI | NMI |
|---|---|---|---|
| major_best_niche | native_R1 | 0.1448 | 0.3559 |
| major_best_niche | native_R2 | 0.1301 | 0.3315 |
| major_best_niche | native_R3 | 0.1467 | 0.322 |
| subtype_best_niche | native_R1 | 0.0815 | 0.3288 |
| subtype_best_niche | native_R2 | 0.07341 | 0.3143 |
| subtype_best_niche | native_R3 | 0.09725 | 0.307 |
| recommended_niche | native_R1 | 0.1448 | 0.3559 |
| recommended_niche | native_R2 | 0.1301 | 0.3315 |
| recommended_niche | native_R3 | 0.1467 | 0.322 |

Native Squidpy1.6.5 runs on pooled all-spot dominant composition with scale=False; main method is continuous section-balanced training plus unified mapping. Thus comparison includes training-balance and representation differences. Native cumulative binary radius graphs supplied at distance1 avoid multiple path weighting. The installed1.6.5 missing-neighbor padding bug is corrected with sparse profile counts, validated on a toy and all input profiles. Mixing loss is quantified in dominant_information_loss.tsv.

## 4–5. Niche counts and recommendation

| role | model | n_niches | silhouette | min_fraction | selection_score |
|---|---|---|---|---|---|
| major_best | major_R1_raw_k30_res0.2 | 9 | 0.1989 | 0.01968 | 0.1989 |
| subtype_best | subtype_R1_raw_k30_res0.2 | 16 | -0.0287 | 0.01175 | -0.1487 |
| recommended | major_R1_raw_k30_res0.2 | 9 | 0.1989 | 0.01968 | 0.1989 |
| R1_representative | major_R1_raw_k30_res0.2 | 9 | 0.1989 | 0.01968 | 0.1989 |
| R2_representative | major_R2_raw_k30_res0.2 | 11 | 0.1838 | 0.01365 | 0.1838 |
| R3_representative | major_R3_raw_k30_res0.2 | 13 | 0.1173 | 0.006587 | -0.01272 |

The best model means best among evaluated candidates, not an exhaustive optimum. Recommendation fixed using composition-only silhouette with prespecified tiny/overfragmented-cluster penalties. No signature, clinical variable, communications or spatial coordinate entered selection. See model_selection_report.md for exact staged search and mapping. Candidate quality indexes, full pooled validation and low-confidence mappings are retained. Seed stability isolates Leiden randomness on a fixed pilot and fixed KNN, not repeated biological sampling.

| model | section_balanced_mean_confidence | section_balanced_low_confidence_fraction |
|---|---|---|
| major_R1_raw_k30_res0.2 | 0.9629 | 0.02619 |
| major_R2_raw_k30_res0.2 | 0.9628 | 0.02688 |
| major_R3_raw_k30_res0.2 | 0.963 | 0.02697 |
| subtype_R1_raw_k30_res0.2 | 0.9555 | 0.03374 |

| model | min_ARI | mean_ARI | min_NMI |
|---|---|---|---|
| major_R1_raw_k30_res0.2 | 0.5365 | 0.6896 | 0.6975 |
| major_R2_raw_k30_res0.2 | 0.5411 | 0.6434 | 0.7395 |
| major_R3_raw_k30_res0.2 | 0.5307 | 0.6125 | 0.7564 |
| subtype_R1_raw_k30_res0.2 | 0.5044 | 0.6127 | 0.7328 |

| model_a | model_b | scope | ARI | NMI |
|---|---|---|---|---|
| major_R1_raw_k30_res0.2 | major_best__pooled | all_spots_balanced_vs_pooled | 0.376 | 0.6334 |
| subtype_R1_raw_k30_res0.2 | subtype_best__pooled | all_spots_balanced_vs_pooled | 0.2222 | 0.532 |

| model | diagnostic | seed | n_niches | ARI_vs_two_iterations_seed123 | ARI_vs_converged_seed123 | NMI_vs_converged_seed123 |
|---|---|---|---|---|---|---|
| major_R1_raw_k30_res0.2 | Leiden_until_convergence | 123 | 9 | 0.9508 | 1 | 1 |
| major_R1_raw_k30_res0.2 | Leiden_until_convergence | 17 | 10 | 0.574 | 0.6016 | 0.7666 |
| major_R1_raw_k30_res0.2 | Leiden_until_convergence | 321 | 11 | 0.8239 | 0.8643 | 0.8742 |
| major_R1_raw_k30_res0.2 | Leiden_until_convergence | 2026 | 11 | 0.6502 | 0.6837 | 0.8012 |
| major_R1_raw_k30_res0.2 | Kmeans_seed_sensitivity | 123 | 5 | NA | 1 | 1 |
| major_R1_raw_k30_res0.2 | Kmeans_seed_sensitivity | 17 | 5 | NA | 0.9918 | 0.9861 |
| major_R1_raw_k30_res0.2 | Kmeans_seed_sensitivity | 321 | 5 | NA | 0.999 | 0.9974 |
| major_R1_raw_k30_res0.2 | Kmeans_seed_sensitivity | 2026 | 5 | NA | 0.997 | 0.9941 |
| subtype_R1_raw_k30_res0.2 | Leiden_until_convergence | 123 | 17 | 0.8983 | 1 | 1 |
| subtype_R1_raw_k30_res0.2 | Leiden_until_convergence | 17 | 15 | 0.5384 | 0.531 | 0.7641 |
| subtype_R1_raw_k30_res0.2 | Leiden_until_convergence | 321 | 15 | 0.6926 | 0.7297 | 0.8553 |
| subtype_R1_raw_k30_res0.2 | Leiden_until_convergence | 2026 | 17 | 0.6866 | 0.6994 | 0.8467 |

## 6. Biological characterization

| niche | name | top_celltype | top_fraction |
|---|---|---|---|
| 0 | Niche 0 | T_NK | 0.2564 |
| 1 | Niche 1 | Mural | 0.2583 |
| 2 | Epithelial-rich niche | Epithelial | 0.6527 |
| 3 | Epithelial-rich niche | Epithelial | 0.4443 |
| 4 | Niche 4 | Epithelial | 0.3058 |
| 5 | Myeloid-rich niche | Myeloid | 0.6049 |
| 6 | Epithelial-rich niche | Epithelial | 0.4132 |
| 7 | Mural-rich niche | Mural | 0.6607 |
| 8 | Niche 8 | Epithelial | 0.246 |

Names are descriptive of section-balanced central composition. Separate central/neighborhood and pooled/section/patient-balanced tables prevent presenting only the clustering input. No label is evidence of malignant status or a functional mechanism.

## 7. Dataset effect

| model | niche | n_spots | fraction | n_sections | n_datasets | max_dataset_share | dominant_dataset | max_section_share |
|---|---|---|---|---|---|---|---|---|
| recommended_niche | 0 | 66666 | 0.2516 | 59 | 7 | 0.432 | GSE238145 | 0.09395 |
| recommended_niche | 1 | 21802 | 0.08227 | 58 | 7 | 0.3438 | GSE238145 | 0.1589 |
| recommended_niche | 2 | 49987 | 0.1886 | 50 | 6 | 0.3242 | CanDis__UNVERIFIED | 0.0841 |
| recommended_niche | 3 | 21420 | 0.08083 | 22 | 6 | 0.3966 | NPJ__UNVERIFIED | 0.1877 |
| recommended_niche | 4 | 72097 | 0.2721 | 61 | 7 | 0.3281 | CanDis__UNVERIFIED | 0.149 |
| recommended_niche | 5 | 5722 | 0.02159 | 22 | 5 | 0.5952 | GSE238145 | 0.2466 |
| recommended_niche | 6 | 15335 | 0.05787 | 39 | 7 | 0.7414 | GSE238145 | 0.5001 |
| recommended_niche | 7 | 4255 | 0.01606 | 24 | 6 | 0.5835 | GSE285715 | 0.3239 |
| recommended_niche | 8 | 7718 | 0.02912 | 19 | 5 | 0.449 | CanDis__UNVERIFIED | 0.4426 |

Dataset and section dominance are descriptive warnings, not proof of a removable technical batch. No Harmony or forced mixing used. Known original subtype identifiability limits remain relevant to subtype niches.

## 8. Fragmentation and smoothing

| model | mean_singleton_fraction | mean_largest_fraction | mean_same_label_edges |
|---|---|---|---|
| R1_representative_niche | 0.1448 | 0.6136 | 0.8976 |
| R2_representative_niche | 0.08892 | 0.6648 | 0.9068 |
| R3_representative_niche | 0.08407 | 0.7209 | 0.9315 |
| kmeans_niche | 0.03181 | 0.6454 | 0.8989 |
| major_best_niche | 0.1448 | 0.6136 | 0.8976 |
| native_R1 | 0.6886 | 0.5415 | 0.5798 |
| native_R2 | 0.6848 | 0.5816 | 0.5975 |
| native_R3 | 0.6842 | 0.6155 | 0.6674 |
| recommended_niche | 0.1448 | 0.6136 | 0.8976 |
| subtype_best_niche | 0.1588 | 0.6417 | 0.8844 |

Spatial continuity is auxiliary: averaging itself induces it. Original missing tissue and disconnected tissue islands also contribute to fragmentation. Connected-component size distributions and all-section maps must be visually reviewed for boundary merging; automated metrics cannot exclude over-smoothing.

## 9–10. C5 and LLPS enrichment

Verified-patient paired contrasts (within-section z, min20 per niche and rest; then same-patient mean):

| niche | feature | n_units | effect | ci_low | ci_high | FDR | positive_fraction |
|---|---|---|---|---|---|---|---|
| 0 | Endothelial_abundance__z | 24 | -0.3541 | -0.5782 | -0.1075 | 0.04225 | 0.25 |
| 0 | LLPS77_UCell__z | 24 | 0.5332 | 0.2183 | 0.8444 | 0.002279 | 0.8333 |
| 0 | five_UCell__z | 24 | 0.3445 | 0.01856 | 0.666 | 0.1968 | 0.625 |
| 1 | Endothelial_abundance__z | 19 | -0.01348 | -0.2985 | 0.2703 | 1 | 0.5263 |
| 1 | LLPS77_UCell__z | 19 | 2.043 | 1.71 | 2.447 | 3.052e-05 | 1 |
| 1 | five_UCell__z | 19 | 0.2558 | -0.04434 | 0.5575 | 0.3913 | 0.5789 |
| 2 | Endothelial_abundance__z | 14 | 0.4731 | 0.1239 | 0.783 | 0.04225 | 0.7143 |
| 2 | LLPS77_UCell__z | 14 | -0.8441 | -1.092 | -0.5259 | 0.002279 | 0.07143 |
| 2 | five_UCell__z | 14 | -0.4016 | -0.6474 | -0.1615 | 0.04297 | 0.1429 |
| 3 | Endothelial_abundance__z | 4 | -0.07781 | -0.4995 | 0.3721 | 1 | 0.5 |
| 3 | LLPS77_UCell__z | 4 | -0.9968 | -1.504 | -0.433 | 0.1429 | 0 |
| 3 | five_UCell__z | 4 | -0.8977 | -1.202 | -0.4463 | 0.25 | 0 |
| 4 | Endothelial_abundance__z | 26 | -0.1764 | -0.3803 | 0.02383 | 0.1864 | 0.3846 |
| 4 | LLPS77_UCell__z | 26 | -0.1235 | -0.3406 | 0.1007 | 0.1899 | 0.3077 |
| 4 | five_UCell__z | 26 | 0.1414 | -0.06073 | 0.3386 | 0.318 | 0.6154 |
| 5 | Endothelial_abundance__z | 15 | -0.782 | -0.8853 | -0.6914 | 0.0004883 | 0 |
| 5 | LLPS77_UCell__z | 15 | -0.4855 | -0.732 | -0.2367 | 0.006714 | 0.2667 |
| 5 | five_UCell__z | 15 | -0.662 | -0.8708 | -0.453 | 0.0009766 | 0.06667 |
| 6 | Endothelial_abundance__z | 8 | 1.435 | 1.007 | 1.841 | 0.03125 | 1 |
| 6 | LLPS77_UCell__z | 8 | -0.4676 | -0.8817 | 0.003692 | 0.125 | 0.125 |
| 6 | five_UCell__z | 8 | 0.1087 | -0.3019 | 0.5258 | 0.8438 | 0.5 |
| 7 | Endothelial_abundance__z | 4 | -0.08555 | -0.8446 | 1.296 | 1 | 0.25 |
| 7 | LLPS77_UCell__z | 4 | 2.221 | 1.906 | 2.516 | 0.1429 | 1 |
| 7 | five_UCell__z | 4 | -0.705 | -1.486 | 0.07573 | 0.4286 | 0.25 |
| 8 | Endothelial_abundance__z | 1 | 1.059 | NA | NA | NA | 1 |
| 8 | LLPS77_UCell__z | 1 | -0.2521 | NA | NA | NA | 0 |
| 8 | five_UCell__z | 1 | 0.7177 | NA | NA | NA | 1 |

All-section sensitivity (sections are not guaranteed independent patients):

| niche | feature | n_units | effect | ci_low | ci_high | FDR | positive_fraction |
|---|---|---|---|---|---|---|---|
| 0 | Endothelial_abundance__z | 53 | -0.1484 | -0.4033 | 0.1245 | 0.08926 | 0.3019 |
| 0 | LLPS77_UCell__z | 53 | 0.8137 | 0.553 | 1.093 | 3.441e-06 | 0.8302 |
| 0 | five_UCell__z | 53 | 0.5785 | 0.3292 | 0.8669 | 0.0005493 | 0.717 |
| 1 | Endothelial_abundance__z | 44 | 0.07399 | -0.1868 | 0.363 | 0.9493 | 0.5227 |
| 1 | LLPS77_UCell__z | 44 | 1.824 | 1.554 | 2.09 | 5.627e-11 | 0.9773 |
| 1 | five_UCell__z | 44 | 0.4061 | 0.1809 | 0.6491 | 0.004254 | 0.6818 |
| 2 | Endothelial_abundance__z | 36 | 0.03928 | -0.2537 | 0.3161 | 0.7704 | 0.5 |
| 2 | LLPS77_UCell__z | 36 | -1.075 | -1.351 | -0.8259 | 1.153e-08 | 0.05556 |
| 2 | five_UCell__z | 36 | -0.6354 | -0.8243 | -0.4537 | 6.403e-06 | 0.1667 |
| 3 | Endothelial_abundance__z | 13 | 0.1703 | -0.3049 | 0.6617 | 0.7704 | 0.5385 |
| 3 | LLPS77_UCell__z | 13 | -1.215 | -1.529 | -0.8874 | 0.0005493 | 0 |
| 3 | five_UCell__z | 13 | -0.9593 | -1.195 | -0.7308 | 0.0005493 | 0 |
| 4 | Endothelial_abundance__z | 57 | -0.05374 | -0.2485 | 0.1505 | 0.5625 | 0.386 |
| 4 | LLPS77_UCell__z | 57 | -0.04775 | -0.2872 | 0.1995 | 0.8831 | 0.4912 |
| 4 | five_UCell__z | 57 | 0.2376 | 0.00878 | 0.4418 | 0.003643 | 0.7193 |
| 5 | Endothelial_abundance__z | 18 | -0.7576 | -0.849 | -0.6671 | 6.866e-05 | 0 |
| 5 | LLPS77_UCell__z | 18 | -0.4392 | -0.6412 | -0.2302 | 0.002369 | 0.2222 |
| 5 | five_UCell__z | 18 | -0.6282 | -0.8207 | -0.432 | 0.0002403 | 0.1111 |
| 6 | Endothelial_abundance__z | 16 | 1.241 | 0.7976 | 1.62 | 0.0009613 | 0.875 |
| 6 | LLPS77_UCell__z | 16 | 0.1573 | -0.2428 | 0.5478 | 0.5942 | 0.5 |
| 6 | five_UCell__z | 16 | 0.5004 | 0.1991 | 0.782 | 0.01683 | 0.75 |
| 7 | Endothelial_abundance__z | 12 | -0.1321 | -0.6944 | 0.4846 | 0.7704 | 0.4167 |
| 7 | LLPS77_UCell__z | 12 | 2.15 | 1.883 | 2.4 | 0.0008789 | 1 |
| 7 | five_UCell__z | 12 | -0.1144 | -0.5264 | 0.264 | 1 | 0.4167 |
| 8 | Endothelial_abundance__z | 5 | 0.4534 | -0.3121 | 1.137 | 0.5625 | 0.6 |
| 8 | LLPS77_UCell__z | 5 | -0.7684 | -1.437 | -0.1082 | 0.1607 | 0.2 |
| 8 | five_UCell__z | 5 | 0.02164 | -0.6113 | 0.61 | 1 | 0.4 |

Raw-score and percentile results, min10/20/30 and fibroblast P50/P75/P90 subsets remain in full tables; no significant-only filtering.

## 11. Fibroblast adjustment

| niche | feature | n_units | effect | ci_low | ci_high | FDR | positive_fraction |
|---|---|---|---|---|---|---|---|
| 0 | A26_UCell | 20 | 0.2235 | 0.01697 | 0.4358 | 0.1292 | 0.65 |
| 0 | LLPS77_UCell | 20 | -0.01507 | -0.2105 | 0.2057 | 0.5217 | 0.35 |
| 0 | five_UCell | 20 | 0.112 | -0.093 | 0.3215 | 0.5958 | 0.55 |
| 1 | A26_UCell | 15 | -0.2444 | -0.4284 | -0.05008 | 0.1292 | 0.2667 |
| 1 | LLPS77_UCell | 15 | 1.184 | 0.961 | 1.399 | 0.0004272 | 1 |
| 1 | five_UCell | 15 | -0.3541 | -0.4999 | -0.2114 | 0.004272 | 0.1333 |
| 2 | A26_UCell | 13 | -0.2349 | -0.4063 | -0.0827 | 0.1196 | 0.3077 |
| 2 | LLPS77_UCell | 13 | -0.1985 | -0.265 | -0.1143 | 0.008545 | 0.1538 |
| 2 | five_UCell | 13 | -0.1708 | -0.2995 | -0.03424 | 0.0752 | 0.2308 |
| 3 | A26_UCell | 3 | -0.3088 | -0.3608 | -0.2777 | 0.25 | 0 |
| 3 | LLPS77_UCell | 3 | -0.4041 | -0.9033 | -0.144 | 0.2917 | 0 |
| 3 | five_UCell | 3 | -0.3062 | -0.3436 | -0.2562 | 0.2917 | 0 |
| 4 | A26_UCell | 24 | 0.08271 | -0.06847 | 0.2284 | 0.25 | 0.6667 |
| 4 | LLPS77_UCell | 24 | -0.2146 | -0.3593 | -0.06773 | 0.02452 | 0.3333 |
| 4 | five_UCell | 24 | 0.1367 | -0.005577 | 0.2714 | 0.08551 | 0.75 |
| 5 | A26_UCell | 8 | 0.08763 | 0.01608 | 0.1707 | 0.25 | 0.75 |
| 5 | LLPS77_UCell | 8 | 0.2155 | 0.05805 | 0.3845 | 0.0957 | 0.75 |
| 5 | five_UCell | 8 | 0.08777 | -0.007171 | 0.1665 | 0.1367 | 0.75 |
| 6 | A26_UCell | 4 | -0.4547 | -0.7571 | -0.1499 | 0.25 | 0.25 |
| 6 | LLPS77_UCell | 4 | -0.4674 | -1.254 | -0.05491 | 0.175 | 0 |
| 6 | five_UCell | 4 | -0.1518 | -0.3075 | -0.04285 | 0.175 | 0 |
| 7 | A26_UCell | 2 | -0.4722 | NA | NA | NA | 0 |
| 7 | LLPS77_UCell | 2 | 1.883 | NA | NA | NA | 1 |
| 7 | five_UCell | 2 | -0.4974 | NA | NA | NA | 0 |
| 8 | A26_UCell | 0 | NA | NA | NA | NA | NA |
| 8 | LLPS77_UCell | 0 | NA | NA | NA | NA | NA |
| 8 | five_UCell | 0 | NA | NA | NA | NA | NA |

Section-specific shared-support regression controls log absolute and relative Fibroblast abundance with quadratic terms. Low support or ill-conditioning returns NA; residual plots alone are not adjusted niche-effect tests. Exact eligibility and all failures are in niche_fibroblast_adjusted_results.tsv. Adjustment is associational, with remaining tissue/technical confounding.

| signature | n_overlap | overlap_genes | sensitivity |
|---|---|---|---|
| five | 0 | NA | not applicable |
| A26 | 0 | NA | not applicable |

## 12. Joint C5–LLPS–Endothelial niche

| model | C5_Endothelial_ids | C5_LLPS_Endothelial_ids |
|---|---|---|
| R2_representative_niche | NA | NA |
| R3_representative_niche | 4 | NA |
| kmeans_niche | NA | NA |
| native_R1 | NA | NA |
| native_R2 | NA | NA |
| native_R3 | NA | NA |
| recommended_niche | NA | NA |
| subtype_best_niche | NA | NA |

C5 unadjusted enrichment is partition-dependent: R2/R3, subtype and Kmeans controls can contain positive niches despite the negative primary screen. Across all five continuous models, no positive verified-patient C5 adjusted contrast passes the stated BH/CI screen, whereas LLPS adjusted positives recur. R3 has an unadjusted C5/Endothelial pair-only pattern; none of the eight fixed models has a triple-positive pattern. This distinction prevents calling every C5-associated pattern absent based only on the main taxonomy.

| model | feature | positive_ids | negative_ids | adjusted_positive_ids | n_estimable_niches | n_niches |
|---|---|---|---|---|---|---|
| R2_representative_niche | Endothelial_abundance__z | NA | 2;5 | not evaluated | 8 | 11 |
| R2_representative_niche | LLPS77_UCell__z | 0;2 | 1;4;6;7 | 0 | 8 | 11 |
| R2_representative_niche | five_UCell__z | 2 | 4;5;6 | NA | 8 | 11 |
| R3_representative_niche | Endothelial_abundance__z | 4 | 5;6;7 | not evaluated | 8 | 13 |
| R3_representative_niche | LLPS77_UCell__z | 0;5 | 3;7 | 0 | 8 | 13 |
| R3_representative_niche | five_UCell__z | 4 | NA | NA | 8 | 13 |
| kmeans_niche | Endothelial_abundance__z | 1;3 | 0;4 | not evaluated | 5 | 5 |
| kmeans_niche | LLPS77_UCell__z | 0;2 | 1;3;4 | 2;4 | 5 | 5 |
| kmeans_niche | five_UCell__z | 0 | 1;4 | NA | 5 | 5 |
| native_R1 | Endothelial_abundance__z | 0 | 205;44;82 | not evaluated | 96 | 230 |
| native_R1 | LLPS77_UCell__z | 2;43;48 | 0 | not evaluated | 96 | 230 |
| native_R1 | five_UCell__z | 48;80 | NA | not evaluated | 96 | 230 |
| native_R2 | Endothelial_abundance__z | 0 | 293 | not evaluated | 99 | 308 |
| native_R2 | LLPS77_UCell__z | 1;217;255;7 | 0 | not evaluated | 99 | 308 |
| native_R2 | five_UCell__z | 217;7 | 0 | not evaluated | 99 | 308 |
| native_R3 | Endothelial_abundance__z | NA | 18;309 | not evaluated | 71 | 263 |
| native_R3 | LLPS77_UCell__z | 1;60 | 0 | not evaluated | 71 | 263 |
| native_R3 | five_UCell__z | 1;24 | 0 | not evaluated | 71 | 263 |
| recommended_niche | Endothelial_abundance__z | 2;6 | 0;5 | not evaluated | 8 | 9 |
| recommended_niche | LLPS77_UCell__z | 0;1 | 2;5 | 1 | 8 | 9 |
| recommended_niche | five_UCell__z | NA | 2;5 | NA | 8 | 9 |
| subtype_best_niche | Endothelial_abundance__z | NA | 4 | not evaluated | 9 | 16 |
| subtype_best_niche | LLPS77_UCell__z | 0 | 3;6 | 0 | 9 | 16 |
| subtype_best_niche | five_UCell__z | 2 | 3;6 | NA | 9 | 16 |

Under the explicit descriptive screen of positive verified-patient effects, nominal95%CI above0 and BH FDR<0.05 for all three within-section standardized features, jointly supported niche IDs are []. C5–Endothelial pair-only supported IDs: []. Fibroblast and Endothelial composition helped define niches, so their enrichment is partly definitional rather than independent evidence. This screen is not an independent confirmatory joint test and does not establish a fibrovascular mechanism. Absence of qualifying niches is a negative result for this screen, not proof of absence.

## 13. Across-section and source reproducibility

Per-section Pearson and Spearman correlations are saved for center and neighborhood scales, each niche and whole section, before/after Fibroblast adjustment. Fisher-z mean correlations and bootstrap section/patient intervals use biological-unit aggregation. Dataset-specific direction summaries retain heterogeneity; section and patient inference do not cover identical source populations.

| feature_a | feature_b | adjustment | n_units | r | ci_low | ci_high | positive_fraction |
|---|---|---|---|---|---|---|---|
| LLPS77_UCell | Endothelial_abundance | fibroblast_adjusted | 29 | -0.2517 | -0.3184 | -0.1865 | 0.03448 |
| LLPS77_UCell | Endothelial_abundance | raw | 29 | 0.00638 | -0.09074 | 0.1093 | 0.5517 |
| five_UCell | Endothelial_abundance | fibroblast_adjusted | 29 | -0.2176 | -0.2788 | -0.1547 | 0.06897 |
| five_UCell | Endothelial_abundance | raw | 29 | 0.1732 | 0.09394 | 0.2614 | 0.7586 |
| five_UCell | LLPS77_UCell | fibroblast_adjusted | 29 | 0.2218 | 0.1434 | 0.2981 | 0.8621 |
| five_UCell | LLPS77_UCell | raw | 29 | 0.3448 | 0.2373 | 0.4489 | 0.931 |

Additional source sensitivity for adjusted LLPS in main Niche1:

| dataset | niche | feature | mean_effect | n_sections | positive_fraction |
|---|---|---|---|---|---|
| CanDis__UNVERIFIED | 1 | LLPS77_UCell | 1.043 | 5 | 1 |
| GSE171351 | 1 | LLPS77_UCell | 1.172 | 1 | 1 |
| GSE238145 | 1 | LLPS77_UCell | 1.235 | 7 | 1 |
| GSE285715 | 1 | LLPS77_UCell | 0.121 | 3 | 0.3333 |
| GSE319536 | 1 | LLPS77_UCell | 1.136 | 7 | 1 |
| Inhousedata__UNVERIFIED | 1 | LLPS77_UCell | 0.9083 | 5 | 1 |
| NPJ__UNVERIFIED | 1 | LLPS77_UCell | 0.9376 | 7 | 1 |

GSE285715 has weaker/mixed adjusted directions and unresolved patient identities. Reused LLPS77 scores also have source-specific gene coverage; within-section normalization does not recover absent genes.

| dataset | min_genes | max_genes | n_sections |
|---|---|---|---|
| CanDis__UNVERIFIED | 75 | 75 | 13 |
| GSE171351 | 77 | 77 | 4 |
| GSE238145 | 45 | 76 | 10 |
| GSE285715 | 56 | 65 | 6 |
| GSE319536 | 75 | 76 | 16 |
| Inhousedata__UNVERIFIED | 76 | 76 | 7 |
| NPJ__UNVERIFIED | 73 | 76 | 7 |

## 14. Niche adjacency

R1,1000 within-section label permutations,seed123. Observed, expected,SD,Z,O/E saved for every section; absent pairs and zero-variance Z are NA. Section/patient/source summaries prioritize O/E and direction. Spatially defined niches inherently autocorrelate, so these maps describe organization and cannot independently demonstrate cell communication.

## 15. Clinical associations

| dataset | tissue | source_sample_type | n_sections | n_patients |
|---|---|---|---|---|
| CanDis__UNVERIFIED | UNRESOLVED | NA | 13 | 0 |
| GSE171351 | Tumor | MIBC | 4 | 4 |
| GSE238145 | Tumor | MIBC | 10 | 9 |
| GSE285715 | Normal | Normal | 3 | 0 |
| GSE285715 | Tumor | MIBC | 3 | 0 |
| GSE319536 | UNRESOLVED | NA | 16 | 16 |
| Inhousedata__UNVERIFIED | UNRESOLVED | NA | 7 | 0 |
| NPJ__UNVERIFIED | UNRESOLVED | NA | 7 | 0 |

| dataset | T_stage | n_sections | n_verified_patients |
|---|---|---|---|
| GSE171351 | T2 | 3 | 3 |
| GSE171351 | T3 | 1 | 1 |
| GSE285715 | T3 | 2 | 0 |
| GSE285715 | T4 | 1 | 0 |

Original workbook re-audit recovered explicit T stages for7 exact-accession sections (GSE171351:T2 three,T3 one; GSE285715:T3 two,T4 one with unknown patients). Within-dataset independent groups are too small for reliable stage inference; descriptive tables retained.

Verified labels only; dataset-stratified description. Formal tests require >=5 known independent patients in both groups. Compositional closure prevents treating niche fractions as independent outcomes. No new stage inference or unsupported compositional modeling.

Eligible nominal tests: 0; rows without sufficient verified groups remain descriptive.

## 16. Figure candidates

Recommended: section-balanced central composition; major/subtype UMAP with source QC; score-blind representative spatial maps; patient/section contrasts with confidence intervals; O/E summary heatmap. Representative sections selected by niche diversity per source then coverage, never by C5/LLPS. Additional Kmeans/R2/R3 and native-control spatial QC maps use these same representatives; native previews color only the five largest groups and label the rest as pooled other clusters. All63 spatial outputs retained. UMAP plotting samples50000 spots only for readability; full embeddings and labels saved.

## 17. Negative, unstable or unidentifiable findings

K-means K5 has much higher seed stability than the finer Leiden taxonomy. It should be reviewed as a coarser alternative, not treated as confirmation of identical niche boundaries. Native dominant profiles overpartition into hundreds of groups and are not a credible standalone niche taxonomy here. Subtype silhouette is negative for the retained candidate; fine subtype niche identities are especially uncertain.

Report all negative contrasts, zero/low signature values, low-support adjusted comparisons, mapping uncertainty, cross-scale disagreement, source-dominated niches and unresolved clinical groups. Known C5 specificity limitations prohibit C5 cell-localization claims. Unverified cohorts and only29 verified patients limit generalization. No iid spot-level significance is used.

## 18. Next steps

Within the frozen main taxonomy, prioritize Niche1 for independent LLPS-associated transcriptional validation (mixed central composition with relatively prominent Mural cells); its C5-associated signal is not a supported positive adjusted effect. Niche6 is an Endothelial-enriched comparison region, not an established C5/LLPS fibrovascular niche; source imbalance must be checked. Niche7 has a strong Mural component but too few verified patients for strong claims. Do not equate reference Matrix_CAF or other CAF subtypes with validated C5 cells.

After manual spatial/histology review, prioritize reproducible composition niches with adequate cross-source and patient support. Consider targeted expression or spatial validation separating fibroblast abundance from transcriptional state. Niche DEG/pseudobulk, GO/KEGG/GSEA, CellChat and NicheNet were intentionally not run in this round.

## Evidence decisions

| domain | rating | reason |
|---|---|---|
| Niche typing reliability | CONDITIONAL GO | Seed stability and balanced/pooled agreement; representation/scale uncertainty remains explicit. |
| C5-associated enrichment | CONDITIONAL GO | Unadjusted positives depend on partition (4 continuous models); positive adjusted patient contrasts: 0. No supported fibroblast-independent enrichment when this count is0. |
| LLPS-associated enrichment | GO | Positive fibroblast-adjusted transcriptional contrasts recur across 5 continuous models; this is not proof of physical phase separation or independent replication across methods. |
| C5–Endothelial common enrichment | CONDITIONAL GO | Unadjusted pair-only pattern appears in R3_representative_niche; no robust fibroblast-independent joint evidence. |
| Triple C5–LLPS–Endothelial enrichment | NO-GO | No tested fixed model passes the triple-positive verified-patient screen when all IDs are empty; not a proof of biological absence. |

Overall: **CONDITIONAL GO for exploratory spatial organization**, contingent on histological review, subtype identifiability and mapping/balance stability. A negative C5/Endothelial hypothesis does not invalidate spatial niche analysis.

## Reproducibility and limits

4-core CPU affinity; numerical thread caps4; address-space limit100GB. Seed123 except seed-stability checks. Exact input IDs, composition and graph QC, scripts, frozen software and per-stage resource logs retained. Original inputs read-only; no expression matrix or full H&E copied. Environment and figures are excluded from biological inference. See logs for dependency compatibility fixes and any actual execution failures; outputs must pass final_validation.json before delivery.