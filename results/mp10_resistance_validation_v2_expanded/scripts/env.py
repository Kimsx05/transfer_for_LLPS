from pathlib import Path
import os,sys,json,re,gzip,hashlib,csv,datetime
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');B=P/'results/mp10_resistance_validation_v2_expanded';V=P/'results/mp10_resistance_validation_v1';S=B/'00_plan_audit/sources'
sys.path=[p for p in sys.path if 'dist-packages' not in p];sys.path[:0]=[str(B/'runtime'),str(P/'spatial_ecology_top20_v2/runtime'),str(P/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
import numpy as np,pandas as pd
from scipy import stats
SEED=20261111
def save(x,p):
 p=B/p;p.parent.mkdir(parents=True,exist_ok=True);x.to_csv(p,sep='\t',index=False,na_rep='NA',compression={'method':'gzip','compresslevel':1,'mtime':1} if str(p).endswith('.gz') else None)
def read(p):return pd.read_csv(B/p,sep='\t',low_memory=False)
def soft(g):
 p=S/(g+'_metadata.soft')
 if not p.exists():p=S/(g+'_family.soft.gz')
 if not p.exists():return [],{}
 f=gzip.open(p,'rt') if p.suffix=='.gz' else p.open();rows=[];r=None;series={}
 for line in f:
  if line.startswith('^SAMPLE'):
   if r:rows.append(r)
   r={'GSM':line.strip().split(' = ')[1],'dataset':g}
  elif line.startswith('^PLATFORM'):
   if r:rows.append(r);r=None
  elif line.startswith('!Series_') and ' = ' in line:
   k,v=line.rstrip('\r\n').split(' = ',1);series.setdefault(k[8:],[]).append(v)
  elif r is not None and line.startswith('!Sample_') and ' = ' in line:
   k,v=line.rstrip('\r\n').split(' = ',1);k=k[8:]
   if k=='characteristics_ch1' and ': ' in v:k,v=v.split(': ',1)
   r[k]=r.get(k,'')+(' | ' if k in r else '')+v
 if r:rows.append(r)
 f.close();return rows,series
