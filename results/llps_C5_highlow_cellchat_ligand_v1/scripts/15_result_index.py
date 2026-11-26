import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']: os.environ[k]='1'
from pathlib import Path
import pandas as pd,numpy as np
O=Path('/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1');T=O/'tables';T.mkdir(exist_ok=True)
network=[];focuslr=[];pathways=[];statuses=[];categories=[]
focus=pd.read_csv(O/'00_input_audit/candidate_overlap.csv').gene.tolist()
for p in sorted((O/'04_pathway_LR_analysis').glob('*_statistics.csv.gz')):
 r=pd.read_csv(p)
 network.append(r[r.kind.isin(['network','interaction_count'])])
 focuslr.append(r[(r.kind=='LR') & r.ligand.isin(focus) & r.partner.isin(['Endothelial','Epithelial','Myeloid','T_NK'])])
 pathways.append(r[r.kind=='pathway']);categories.append(r[r.kind=='category'])
for p in sorted((O/'03_patient_level_comparison').glob('*_DONE.csv')):statuses.append(pd.read_csv(p))
for name,rows in [('all_grouping_threshold_network_statistics',network),('focused_LR_all_grouping_threshold_statistics',focuslr),('all_pathway_grouping_threshold_statistics',pathways),('signal_category_grouping_threshold_statistics',categories)]:
 pd.concat(rows,ignore_index=True).to_csv(T/(name+'.csv.gz'),index=False)
pd.concat(statuses,ignore_index=True).to_csv(T/'analysis_completion_summary.csv',index=False)
main=pd.concat(network,ignore_index=True).query("method=='within_median' and min_cells==10")
main.to_csv(T/'primary_network_dataset_LODO.csv',index=False)
# Report each focus gene and every endothelial receptor, preserving null effects and denominators.
lr=pd.concat(focuslr,ignore_index=True).query("method=='within_median' and min_cells==10 and partner=='Endothelial'")
lr.to_csv(T/'primary_focused_Endothelial_LR_all_strata.csv.gz',index=False)
# Show whether primary effects change sign across grouping/threshold sensitivities, not as independent replication.
n=pd.concat(network,ignore_index=True);n=n[n.stratum=='Tumor_ALL'];n['effect_sign']=np.sign(n.difference_median)
n.to_csv(T/'network_sensitivity_overview.csv',index=False)
print('RESULT INDEX COMPLETE')
