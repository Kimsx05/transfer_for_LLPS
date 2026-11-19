from pathlib import Path
import hashlib,pandas as pd,json
O=Path('results/C5_Endothelial_NicheNet_spatial_v1');paths=[Path('results/llps_C5_distribution_audit_20261110_v1/v2/patient_sample_mapping_final.csv'),Path('results/llps_C5_distribution_audit_20261110_v1/v2/sample_celltype_counts.csv'),Path('results/fig2_like_C0/phase2_external/reference_validation/patient_celltype_pseudobulk.rds'),Path('results/fig2_like_C0/phase2_external/scripts/02_aggregate_reference.R'),Path('results/llps_C5_cellchat_v2/00_audit/sample_input_validation.csv'),Path('results/C5_spatial_exploratory_all63_v1/00_inputs_and_QC/spatial_tissue_patient_audit.tsv'),Path('results/C5_spatial_exploratory_all63_v1/00_inputs_and_QC/frozen_signatures.json')]
paths+=list(Path('results/C5_spatial_exploratory_all63_v1/01_spatial_scores').glob('*.parquet'))
rows=[]
for p in paths:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(2**20),b''):h.update(x)
 rows.append(dict(path=str(p.resolve()),bytes=p.stat().st_size,sha256=h.hexdigest()))
pd.DataFrame(rows).to_csv(O/'00_input_audit/source_checksums.tsv',sep='\t',index=False)
print('PROVENANCE_COMPLETE',len(rows))
