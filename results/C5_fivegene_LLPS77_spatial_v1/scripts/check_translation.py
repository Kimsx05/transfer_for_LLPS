from pathlib import Path
import numpy as np,pandas as pd
from scipy.signal import fftconvolve
O=Path('results/C5_fivegene_LLPS77_spatial_v1');P=Path('results/C5_spatial_exploratory_all63_v1');rows=[]
for sid in pd.read_csv(O/'00_input_audit/section_mapping.tsv',sep='\t').section_id:
 s=pd.read_parquet(P/'01_spatial_scores'/f'{sid}.parquet',columns=['array_row','array_col']);coords=s.to_numpy(int);coords-=coords.min(0);grid=np.zeros(tuple(coords.max(0)+1));grid[coords[:,0],coords[:,1]]=1;corr=np.rint(fftconvolve(grid,grid[::-1,::-1],mode='full')).astype(int);lookup=set(map(tuple,coords))
 for dr,dc in [(0,2),(1,1),(5,11),(-9,15),(20,0)]:
  if abs(dr)>=grid.shape[0] or abs(dc)>=grid.shape[1]:continue
  actual=sum((r+dr,c+dc) in lookup for r,c in coords);fast=corr[dr+grid.shape[0]-1,dc+grid.shape[1]-1];assert actual==fast;rows.append(dict(section_id=sid,dr=dr,dc=dc,direct_count=actual,fft_count=fast))
pd.DataFrame(rows).to_csv(O/'00_input_audit/translation_overlap_validation.tsv',sep='\t',index=False);print('Exact FFT/direct overlap agreement',len(rows))
