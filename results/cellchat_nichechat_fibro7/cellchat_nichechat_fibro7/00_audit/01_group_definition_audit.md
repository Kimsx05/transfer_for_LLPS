# CellChat group-definition audit

- Authoritative full atlas: /home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/results/objects/major_annotated_with_LLPS77_AUC.qs
- Authoritative final fibro7 assignment: /home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order/global_res012_object.qs
- No existing single full-atlas object contains the final fibro7 subtype field. Exact barcode join adds cellchat_group only to a new analysis copy.
- Actual fields: sample=orig.ident; dataset=dataset; major=major_celltype; final fibro=subcluster (frozen fibro object).
- Cells=353623; biological samples=59; primary included=344377.
- patient_uid and other patient identifiers are ABSENT. Patient count is unknown; long tables retain patient_uid=NA, without inferring patients from samples.
- Frozen C0–C6 cells=20146; assignment equals clean_snn_res.0.12, Idents and exported frozen cell labels. Major annotation, samples, dataset and LLPS77_AUC match the parent atlas exactly.
- Previously excluded stromal cells=1922 = 1759 (old C4/C5/C6) + 163 (old clean-v2 C5). Exact union equals all atlas fibro/pericyte barcodes absent from final fibro7. These historical exclusions are inherited, with no new QC, no reassignment, and no generic Fibro/Pericyte group.
- Major Cycling cells excluded=7324. Frozen annotation does not resolve lineage. Fibro C6 is retained as requested; it is distinct from excluded major Cycling.
- RNA input is the full atlas RNA/data LogNormalize layer; commands verify scale.factor=10000. No score/reduction/integrated input.
- CellChat 2.1.2 is now installed and verified (2026-10-08); exact parameters, API and session records are in this audit directory.
- All 14 groups remain in every count/matrix framework. Counts below 30 are structurally ineligible; absence is never recoded as biological zero.

- Group order: Fibro_C0, Fibro_C1, Fibro_C2, Fibro_C3, Fibro_C4, Fibro_C5, Fibro_C6, Epithelial, Endothelial, Myeloid, T_NK, B, Plasma, Mast

## Group coverage

    final_cellchat_group total_cells n_samples_present n_samples_eligible
                  <char>       <int>             <int>              <int>
 1:             Fibro_C0       10717                43                 20
 2:             Fibro_C1        5280                39                 22
 3:             Fibro_C2        1305                47                 10
 4:             Fibro_C3         913                39                  7
 5:             Fibro_C4         811                39                  7
 6:             Fibro_C5         716                37                 10
 7:             Fibro_C6         404                25                  3
 8:           Epithelial      130804                57                 49
 9:          Endothelial        8163                54                 33
10:              Myeloid       34536                59                 55
11:                 T_NK      121222                55                 50
12:                    B       18622                51                 39
13:               Plasma        8374                57                 30
14:                 Mast        2510                51                 20
    max_sample_cells
               <int>
 1:             2774
 2:             1715
 3:              443
 4:              380
 5:              159
 6:              109
 7:              168
 8:             9601
 9:             1497
10:             2398
11:             6149
12:             4610
13:             1884
14:              670
