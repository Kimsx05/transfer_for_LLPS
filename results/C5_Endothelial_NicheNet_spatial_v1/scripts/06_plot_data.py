from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
O=Path('results/C5_Endothelial_NicheNet_spatial_v1');S=Path('results/C5_spatial_exploratory_all63_v1');pilot=pd.read_csv(S/'00_inputs_and_QC/pilot_prespecified.tsv',sep='\t').section_id.tolist();axes=pd.read_csv(O/'07_spatial_expression_validation/frozen_axes_subunits.tsv',sep='\t');ranking=pd.read_csv(O/'04_NicheNet_ligand_activity/aggregate_priority.tsv',sep='\t');wanted=[g for g in ranking.test_ligand if g in set(axes.loc[axes.role=='candidate','ligand'])][:2];choice=[]
for sid in pilot:
 d=pd.read_parquet(O/'07_spatial_expression_validation'/sid/'expression.parquet')
 for l in wanted:
  ax=axes[axes.ligand==l].iloc[0];rec=d[ax.required_subunits.split(';')].min(1,skipna=False);sets={'Five-gene UCell':d.five_UCell,'Endothelial abundance':d.major_mean__Endothelial,l+' log1p(CP10k)':d[l],ax.receptor+' minimum subunit':rec};out=[]
  for label,v in sets.items():
   z=d[['spot_id','x','y']].copy();z['variable']=label;z['raw_value']=v;z['percentile']=rankdata(v)/len(v) if np.isfinite(v).all() else np.nan;out.append(z)
  pd.concat(out).to_csv(O/'figures'/f'D_{sid}_{l}_data.tsv.gz',sep='\t',index=False);choice.append(dict(section_id=sid,ligand=l,selection='Original seven prespecified pilot sections; first two ranked selected ligands, independent of spatial results'))
pd.DataFrame(choice).to_csv(O/'figures/D_selection.tsv',sep='\t',index=False)
