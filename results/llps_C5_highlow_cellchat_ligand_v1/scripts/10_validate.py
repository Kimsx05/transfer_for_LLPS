import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import pandas as pd,numpy as np,hashlib,subprocess,struct
O=Path('/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1');checks=[]
def check(name,ok):
 checks.append(dict(check=name,pass_=bool(ok)))
 if not ok:raise AssertionError(name)
def bh(p):
 p=np.array(p,float);out=np.full(len(p),np.nan);ix=np.where(np.isfinite(p))[0];order=ix[np.argsort(p[ix])];out[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(ix)/np.arange(1,len(ix)+1))[::-1])[::-1]);return out
m=pd.read_csv(O/'00_input_audit/manifest.csv');c=pd.read_csv(O/'00_input_audit/C5_metadata.csv');check('59 specimens 42 patients',len(m)==59 and m.patient_id.nunique()==42);check('716 C5 37 specimens 31 patients',len(c)==716 and c.specimen_id.nunique()==37 and c.patient_id.nunique()==31)
for _,r in pd.read_csv(O/'00_input_audit/source_stat.csv').iterrows():
 st=Path(r.path).stat();check('source unchanged '+Path(r.path).name,st.st_size==r['size'] and abs(st.st_mtime-r.mtime)<1e-4)
for _,r in pd.read_csv(O/'00_input_audit/history_sources.csv').iterrows():check('history unchanged '+r.path,hashlib.sha256(Path(r.path).read_bytes()).hexdigest()==r.sha256)
s=pd.read_csv(O/'01_LLPS_stratification/cell_assignments.csv.gz');z=s[s.method=='within_median'];check('exact median rules',((z.state=='high')==(z.score>z['median'])).all())
jobs=pd.read_csv(O/'02_CellChat_highlow/job_manifest.csv')
for _,j in jobs[jobs.runnable].iterrows():
 d=O/'02_CellChat_highlow'/j.job_id;check(j.job_id+' complete',(d/'COMPLETE.csv').exists())
 for th in [5,10,20,30]:
  r=pd.read_csv(d/f'threshold{th}.csv.gz');good=(r.high_cells>th)&(r.low_cells>th)&(r.partner_cells>th);check(j.job_id+f' t{th} strict eligibility',(good==r.eligible).all());check(j.job_id+f' t{th} NA distinction',r.loc[~good,'strength'].isna().all() and np.isfinite(r.loc[good,'strength']).all());keys=['state','direction','partner'];lr=r[(r.kind=='LR')&good].groupby(keys).strength.sum();net=r[(r.kind=='network')&good].set_index(keys).strength;check(j.job_id+f' t{th} LR sum reconciles network',len(lr)==len(net) and np.allclose(lr.sort_index(),net.sort_index(),atol=1e-10))
for p in (O/'04_pathway_LR_analysis').glob('*_statistics.csv.gz'):
 r=pd.read_csv(p);check(p.name+' small n withheld',r.loc[r.n_patients<5,['p','FDR','CI_low','CI_high']].isna().all().all());check(p.name+' patient positivity bounded',(r.positive_high_patients<=r.n_patients).all() and (r.positive_low_patients<=r.n_patients).all())
 for _,z in r.groupby(['method','min_cells','kind','direction','stratum']):check(p.name+' BH '+str(_),np.allclose(bh(z.p),z.FDR,equal_nan=True))
for p in sorted((O/'figures').glob('*.pdf')):check('vector PDF '+p.name,len(subprocess.check_output(['pdfimages','-list',str(p)],text=True).splitlines())==2)
for p in sorted((O/'figures').glob('*.png')):
 b=p.read_bytes();i=b.index(b'pHYs');x,y,u=struct.unpack('>IIB',b[i+4:i+13]);check('PNG300 '+p.name,u==1 and abs(x*.0254-300)<1);check('plot data '+p.name,(p.parent/(p.stem+'_data.csv')).exists())
pd.DataFrame(checks).to_csv(O/'00_input_audit/final_QC.csv',index=False);print('FINAL QC PASS',len(checks))
