from pathlib import Path
import json,sys,time
import numpy as np,pandas as pd,h5py
from scipy.stats import rankdata,spearmanr
O=Path('results/C5_spatial_exploratory_all63_v1');R=Path('/home/jinshengxi/LLPS/spatial_data/LLPS_spatial_transfer_manual');sys.path.insert(0,str(R/'scripts'));from data_io import index,blocks
qc=[]
sig=json.loads((O/'00_inputs_and_QC/frozen_signatures.json').read_text())['A26'];rows=[]
for r in pd.read_csv(O/'00_inputs_and_QC/section_mapping.tsv',sep='\t').itertuples():
 while not (O/'01_spatial_scores'/f'{r.section_id}.parquet').exists():time.sleep(5)
 s=pd.read_parquet(O/'01_spatial_scores'/f'{r.section_id}.parquet');score=s.A26_UCell.to_numpy();high=score>np.quantile(score,.9)
 with h5py.File(R/r.expression_relative_path,'r') as f:
  genes=index(f['var']);lookup={int(np.flatnonzero(genes==g)[0]):j for j,g in enumerate(sig) if g in genes};terms=np.zeros((len(s),26))
  for start,x in blocks(f):
   for i in range(x.shape[0]):
    a,b=x.indptr[i:i+2];rank=rankdata(-x.data[a:b],method='average');ids=x.indices[a:b]
    sel=np.isin(ids,list(lookup))
    for k,rr in zip(ids[sel],rank[sel]):terms[start+i,lookup[k]]=max(1500-rr,0)
 total=terms.sum(1);den=total[high].sum()
 n_present=len(lookup);reconstructed=total/(n_present*1500-n_present*(n_present+1)/2);error=float(np.max(abs(reconstructed-score)));qc.append(dict(section_id=r.section_id,max_UCell_reconstruction_error=error));assert error<1e-7,(r.section_id,error)
 for j,g in enumerate(sig):rows.append(dict(section_id=r.section_id,dataset=r.dataset,gene=g,high_rank_weight_share=terms[high,j].sum()/den if den else np.nan,spearman_with_A26=spearmanr(terms[:,j],score).statistic if terms[:,j].std()>0 else np.nan,mean_rank_weight=terms[:,j].mean()))
 print(r.section_id,flush=True)
pd.DataFrame(rows).to_csv(O/'04_fibroblast_and_depth_controls/A26_gene_rank_contributions.tsv',sep='\t',index=False)

pd.DataFrame(qc).to_csv(O/'00_inputs_and_QC/UCell_rank_reconstruction.tsv',sep='\t',index=False)
