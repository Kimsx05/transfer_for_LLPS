from pathlib import Path
import sys,os,json,tempfile,subprocess,shutil
import warnings
import pandas as pd,numpy as np,h5py,scipy.sparse as sp,pyarrow as pa
from scipy.stats import rankdata
from scipy.io import mmwrite,_fast_matrix_market
warnings.filterwarnings('ignore',category=pd.errors.PerformanceWarning)
_fast_matrix_market.PARALLELISM=1;pa.set_cpu_count(1);pa.set_io_thread_count(1)
O=Path('results/C5_fivegene_LLPS77_spatial_v1');P=Path('results/C5_spatial_exploratory_all63_v1');R=Path('/home/jinshengxi/LLPS/spatial_data/LLPS_spatial_transfer_manual');sys.path.insert(0,str(R/'scripts'));from data_io import index,blocks
sets=json.loads((O/'00_input_audit/scoring_sets.json').read_text());five=['MMP11','COL11A1','INHBA','COL5A1','CTHRC1'];mp=pd.read_csv(O/'00_input_audit/section_mapping.tsv',sep='\t').set_index('section_id')
for sid in sys.argv[1].split(','):
 dest=O/'01_LLPS77_spatial_scores'/f'{sid}.parquet'
 if dest.exists():print('REUSE',sid,flush=True);continue
 print('START',sid,flush=True);old=pd.read_parquet(P/'01_spatial_scores'/f'{sid}.parquet')
 with h5py.File(R/mp.loc[sid,'expression_relative_path']) as f:gn=index(f['var']);obs=index(f['obs']);x=sp.vstack([v for _,v in blocks(f)],format='csr')
 assert list(obs)==list(old.barcode);assert np.array_equal(np.asarray(x.sum(1)).ravel(),old.nCount)
 s=old[['spot_id','barcode']].copy();lookup={g:i for i,g in enumerate(gn)};allg=list(dict.fromkeys(sets['LLPS77']+five));ix=[lookup[g] for g in allg if g in lookup];present=[g for g in allg if g in lookup];terms=np.zeros((len(s),len(ix)))
 for j in range(len(s)):
  a,b=x.indptr[j:j+2];inds=x.indices[a:b];vals=x.data[a:b];rr=rankdata(-vals,method='average');dic=dict(zip(inds,np.maximum(1500-rr,0)));terms[j]=[dic.get(i,0) for i in ix]
 contributions=[]
 for k,g in enumerate(present):
  values=x[:,lookup[g]].toarray().ravel();s[g+'_log1p_CP10k']=np.log1p(values/old.nCount.to_numpy()*1e4);s[g+'_rank_contribution']=terms[:,k];contributions.append(dict(section_id=sid,gene=g,detection=(values>0).mean(),mean_rank_contribution=terms[:,k].mean()))
 base=sum(s[g+'_rank_contribution'].to_numpy() for g in five)/(5*1500-15);err=np.max(abs(base-old.five_UCell));assert err<1e-7,(sid,err)
 for g in five:s['five_without_'+g+'_UCell']=(base*(5*1500-15)-s[g+'_rank_contribution'])/(4*1500-10)
 for name,gs in sets.items():
  ids=[lookup[g] for g in gs if g in lookup];s[name+'_covered_genes']=len(ids);s[name+'_detected_genes']=np.asarray((x[:,ids]>0).sum(1)).ravel()
 with tempfile.TemporaryDirectory(prefix='five_llps_') as td:
  td=Path(td);mmwrite(td/'counts.mtx',x.T);(td/'genes.txt').write_text('\n'.join(gn)+'\n');(td/'barcodes.txt').write_text('\n'.join(obs)+'\n')
  with open(O/'00_input_audit'/f'{sid}_score.log','w') as log:subprocess.run(['Rscript',str(O/'scripts/score.R'),str(td),str(td/'scores.tsv'),str(O/'00_input_audit/scoring_sets.json')],stdout=log,stderr=log,check=True)
  sc=pd.read_csv(td/'scores.tsv',sep='\t');s=s.merge(sc,on='barcode',validate='one_to_one');shutil.copy2(td/'R_session.txt',O/'00_input_audit/R_session.txt')
 for name in sets:
  for method in ['UCell','AUCell']:
   col=name+'_'+method;assert s[col].between(0,1).all();s[col+'_percentile']=s[col].rank(pct=True)
 s.to_parquet(str(dest)+'.tmp',index=False);os.replace(str(dest)+'.tmp',dest)
 pd.DataFrame(contributions).to_csv(O/'00_input_audit'/f'{sid}_gene_detection.tsv',sep='\t',index=False)
 (O/'00_input_audit'/f'{sid}_technical.json').write_text(json.dumps(dict(n_spots=len(s),five_reconstruction_error=float(err),universe=len(gn))))
 print('DONE',sid,len(s),flush=True)
