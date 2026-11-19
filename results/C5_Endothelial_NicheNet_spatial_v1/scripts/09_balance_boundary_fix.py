from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
O=Path('results/C5_Endothelial_NicheNet_spatial_v1');S=Path('results/C5_spatial_exploratory_all63_v1');changes=[]
for d in sorted((O/'07_spatial_expression_validation').glob('*/DONE')):
 folder=d.parent;sid=folder.name;b=pd.read_csv(folder/'balance.tsv',sep='\t');pairs=pd.read_csv(folder/'matched_spots.tsv',sep='\t');s=pd.read_parquet(S/'01_spatial_scores'/f'{sid}.parquet').set_index('spot_id');rc=s[['array_row','array_col']].to_numpy().astype(int);loc=set(map(tuple,rc));deg=np.array([sum((r+dr,c+dc) in loc for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)]) for r,c in rc])/6;s['degree']=deg
 for i,row in b[b.SMD_after.isna() & (b.n_pairs>=2)].iterrows():
  p=pairs[pairs.top_fraction==row.top_fraction];v=pd.Series(deg,index=s.index) if row.covariate=='degree' else pd.Series(rankdata(s[row.covariate])/len(s),index=s.index);a=v.loc[p.high_spot].to_numpy();c=v.loc[p.control_spot].to_numpy()
  if len(a)>=2 and np.var(a,ddof=1)+np.var(c,ddof=1)<1e-24 and abs(a.mean()-c.mean())<1e-12:b.loc[i,'SMD_after']=0;changes.append(dict(section_id=sid,top_fraction=row.top_fraction,covariate=row.covariate,reason='Both matched groups have identical constant covariate: SMD=0, not missing'))
 e=pd.read_csv(folder/'effects.tsv',sep='\t')
 for t,z in b.groupby('top_fraction'):
  balanced=z.n_pairs.iloc[0]>=10 and z.SMD_after.notna().all() and z.SMD_after.abs().max()<.2;e.loc[e.top_fraction==t,'balanced']=balanced
 b.to_csv(folder/'balance.tsv',sep='\t',index=False);e.to_csv(folder/'effects.tsv',sep='\t',index=False)
pd.DataFrame(changes).to_csv(O/'00_input_audit/constant_covariate_SMD_correction.tsv',sep='\t',index=False);print('CORRECTED',len(changes))
