from pathlib import Path
import pandas as pd
O=Path('results/C5_Endothelial_GSEA_nominalNicheNet_v1');P=Path('results/C5_Endothelial_NicheNet_spatial_v1')
a=pd.read_csv(O/'03_nominal_NicheNet/checkpoints/Set_B.tsv',sep='\t');e=pd.read_csv(P/'06_CellChat_integration/complete_expression_activity_LR_target_evidence.tsv',sep='\t');old=pd.read_csv(P/'07_spatial_expression_validation/frozen_axes_subunits.tsv',sep='\t');e=e[(e.positive_patients>0)&~e.test_ligand.isin(old.ligand)].drop_duplicates(['test_ligand','receptor']);d=a[a['rank']<=20].merge(e[['test_ligand','receptor','required_subunits','Endo_min_subunit_pct','Endo_min_subunit_support_patients','patient_condition_median','positive_patients','eligible_patients']],on='test_ligand');d=d.sort_values(['rank','patient_condition_median'],ascending=[True,False]).drop_duplicates('test_ligand').head(3)
p=O/'06_CellChat_spatial_integration/frozen_new_spatial_axes.tsv'
if p.exists():
 prior=pd.read_csv(p,sep='\t');assert list(zip(d.test_ligand,d.receptor))==list(zip(prior.test_ligand,prior.receptor))
else:d.to_csv(p,sep='\t',index=False)
