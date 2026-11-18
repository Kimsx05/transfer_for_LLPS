from pathlib import Path
import pandas as pd,numpy as np
from scipy.stats import rankdata
O=Path('results/C5_fivegene_LLPS77_spatial_v1');P=Path('results/C5_spatial_exploratory_all63_v1');five=['MMP11','COL11A1','INHBA','COL5A1','CTHRC1']
def residual(y,cs):
 X=np.column_stack([np.ones(len(y))]+[rankdata(c)/len(y) for c in cs]);v=rankdata(y)/len(y);return v-X@np.linalg.lstsq(X,v,rcond=None)[0]
for sid in pd.read_csv(O/'00_input_audit/section_mapping.tsv',sep='\t').section_id:
 dest=O/'figures'/f'{sid}_plot_data.tsv.gz'
 if dest.exists():continue
 s=pd.read_parquet(P/'01_spatial_scores'/f'{sid}.parquet');d=pd.read_parquet(O/'01_LLPS77_spatial_scores'/f'{sid}.parquet');s=s.merge(d[['spot_id']+[c for c in d if c not in s]],on='spot_id',validate='one_to_one');cs=[s.major_mean__Fibroblast,s.nCount,s.nFeature,s.major_mean__Epithelial];s['five_adjusted_residual']=residual(s.five_UCell,cs);s['LLPS77_adjusted_residual']=residual(s.LLPS77_UCell,cs);prior=pd.read_csv(P/'02_spatial_maps'/f'{sid}_plot_data.tsv.gz',sep='\t').set_index('spot_id');s['epithelial_graph_distance']=prior.loc[s.spot_id,'epithelial_graph_distance'].to_numpy();s['HH_top10']=(s.five_UCell>s.five_UCell.quantile(.9))&(s.LLPS77_UCell>s.LLPS77_UCell.quantile(.9));cols=['spot_id','x','y','five_UCell','five_AUCell','A26_UCell','core15_UCell','C2C4_program_UCell','LLPS77_UCell','LLPS77_AUCell','LLPS77_common_UCell','major_mean__Fibroblast','major_mean__Epithelial','nCount','nFeature','five_adjusted_residual','LLPS77_adjusted_residual','epithelial_graph_distance','HH_top10']+[g+'_log1p_CP10k' for g in five];s[cols].to_csv(dest,sep='\t',index=False)
print('All plot inputs ready')
