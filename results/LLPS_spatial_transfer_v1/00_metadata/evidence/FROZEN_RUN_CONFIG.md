# Frozen cell2location v2 30-type full-run config

REFERENCE:
- reference path: /home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/reference_training_v2/11_final_spatial_reference/30type_regression_model/model/model.pt
- signature path: /home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/reference_training_v2/11_final_spatial_reference/30type_regression_model/learned_signatures_30type.csv.gz
- cells: 300033
- genes: 17987
- batches: 58
- labels (30):
- Arterial EC
- C1QC TAM
- CD4 Helper/Conventional
- CD4 Naive/Memory-like
- CD8 Exhausted
- CD8 GZMK/Effector-memory-like
- Capillary EC
- Classical Mono
- Epithelial
- GC B
- LAMP3 DC
- Lymphatic EC
- Mast
- Matrix_CAF
- NK
- Naive/Memory B
- Neutrophil
- Pericyte
- Plasma
- SMC
- SPP1 TAM
- Tfh/CXCL13-like
- Treg
- Venous EC
- cDC1
- cDC2
- iCAF
- pDC
- rCAF
- vCAF

SPATIAL:
- eligible sections = 63
- first 10 queue sections:
- CanDis_2026_Sample1
- GSE171351_A1
- GSE238145_B10
- GSE285715_N1
- GSE319536_BC01
- Inhousedata_T4384654
- NPJ_sample1
- CanDis_2026_Sample10
- CanDis_2026_Sample11
- CanDis_2026_Sample12

MODEL:
- N_cells_per_location = 20
- detection_alpha = 20
- epochs = 30000
- learning rate = 0.002
- batch_size = None
- train_size = 1.0
- posterior samples = 1000

COMPUTE:
- CPU threads = 8
- GPU = single GPU
- concurrency = 1

OUTPUT:
- /home/data/t070721/codex_workspace/Bladder Metabolism/cell2location_full_run/v2_30type
