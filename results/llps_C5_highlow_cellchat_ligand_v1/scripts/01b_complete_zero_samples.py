import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import pandas as pd,numpy as np
O=Path('/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1');S=O/'01_LLPS_stratification'
m=pd.read_csv(O/'00_input_audit/manifest.csv');s=pd.read_csv(S/'specimen_group_counts.csv');methods=s.method.unique().tolist();base=m[['specimen_id','patient_id','dataset','clinical_group']].merge(pd.DataFrame({'method':methods}),how='cross');s=base.merge(s,on=['specimen_id','patient_id','dataset','clinical_group','method'],how='left',validate='one_to_one')
for c in ['C5_cells','high','low','excluded','median_ties']:s[c]=s[c].fillna(0).astype(int)
s['group_status']=np.where(s.C5_cells==0,'NO_C5_CELLS',np.where((s.high==0)|(s.low==0),'EMPTY_HIGH_OR_LOW','GROUPED'));s.to_csv(S/'specimen_group_counts.csv',index=False,na_rep='NA')
e=s.merge(pd.DataFrame({'min_cells':[5,10,20,30]}),how='cross');e['eligible']=(e.high>e.min_cells)&(e.low>e.min_cells);e.assign(eligible=e.eligible.map({True:'TRUE',False:'FALSE'})).to_csv(S/'specimen_eligibility.csv',index=False,na_rep='NA')
ct=pd.read_csv(O.parent/'llps_C5_distribution_audit_20261110_v1/v2/sample_celltype_counts.csv');stromal=pd.read_csv(O/'00_input_audit/stromal_specimen_counts.csv');ct=ct.drop(columns=['Other_fibroblasts','Pericyte'],errors='ignore').merge(stromal,on='specimen_id',how='left',validate='one_to_one');ct[['Other_fibroblasts','Pericyte']]=ct[['Other_fibroblasts','Pericyte']].fillna(0).astype(int);partners=['Epithelial','Endothelial','Myeloid','T_NK','B','Plasma','Mast','Other_fibroblasts','Pericyte'];ct=ct.melt(id_vars='specimen_id',value_vars=partners,var_name='partner',value_name='partner_cells');ee=e.merge(ct,on='specimen_id');ee['edge_eligible']=ee.eligible&(ee.partner_cells>ee.min_cells);ee.assign(eligible=ee.eligible.map({True:'TRUE',False:'FALSE'}),edge_eligible=ee.edge_eligible.map({True:'TRUE',False:'FALSE'})).to_csv(S/'partner_eligibility.csv.gz',index=False,na_rep='NA')
rows=[]
for key,z in e.groupby(['method','min_cells','clinical_group']):rows.append(dict(zip(['method','min_cells','clinical_group'],key),specimens=z.loc[z.eligible,'specimen_id'].nunique(),patients=z.loc[z.eligible,'patient_id'].nunique()))
pd.DataFrame(rows).to_csv(S/'feasibility_summary.csv',index=False)
j=s.merge(m[['specimen_id','input_path','run_id']],on='specimen_id');j['job_id']=j.method+'__'+j.run_id;j['runnable']=(j.high>5)&(j.low>5);j.assign(runnable=j.runnable.map({True:'TRUE',False:'FALSE'})).to_csv(O/'02_CellChat_highlow/job_manifest.csv',index=False,na_rep='NA');assert len(j)==236 and j.runnable.sum()==55
print('All 59 specimens retained per method; no-C5 samples explicitly zero and ineligible')

ee.groupby(["method","min_cells","partner","clinical_group","dataset"]).apply(lambda z: pd.Series({"eligible_specimens":z.loc[z.edge_eligible,"specimen_id"].nunique(),"eligible_patients":z.loc[z.edge_eligible,"patient_id"].nunique()}),include_groups=False).reset_index().to_csv(S/"partner_coverage.csv",index=False)
