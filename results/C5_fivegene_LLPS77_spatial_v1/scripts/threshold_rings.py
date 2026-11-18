from pathlib import Path
from collections import deque
import pandas as pd,numpy as np,scipy.sparse as sp,zlib
from scipy.stats import rankdata
from scipy.spatial import cKDTree
O=Path('results/C5_fivegene_LLPS77_spatial_v1');P=Path('results/C5_spatial_exploratory_all63_v1');labels=pd.read_csv(O/'00_input_audit/spatial_tissue_patient_audit.tsv',sep='\t');out=[];balances=[];pair_rows=[]
def distance(A,seed):
 d=np.full(A.shape[0],np.nan);q=deque(seed.tolist());d[seed]=0
 while q:
  i=q.popleft()
  for j in A.indices[A.indptr[i]:A.indptr[i+1]]:
   if np.isnan(d[j]):d[j]=d[i]+1;q.append(j)
 return d
def match(h,c,z,rng):
 ids=np.flatnonzero(c);tree=cKDTree(z[ids]);used=set();pairs=[]
 for i in rng.permutation(np.flatnonzero(h)):
  jj=np.asarray(tree.query_ball_point(z[i],np.sqrt(3)*.1),int)
  if not len(jj):continue
  jj=jj[np.argsort(((z[ids[jj]]-z[i])**2).sum(axis=1),kind='stable')]
  for j in jj:
   k=ids[j]
   if k not in used and max(abs(z[k]-z[i]))<=.1:pairs.append((i,k));used.add(k);break
 return tuple(np.array(pairs).T) if pairs else (np.array([],int),np.array([],int))
for row in labels.itertuples():
 sid=row.section_id;s=pd.read_parquet(P/'01_spatial_scores'/f'{sid}.parquet');new=pd.read_parquet(O/'01_LLPS77_spatial_scores'/f'{sid}.parquet');assert list(s.spot_id)==list(new.spot_id);info=dict(section_id=sid,dataset=row.dataset,patient_uid=row.patient_uid,tissue=row.tissue);n=len(s);v=[s.major_mean__Fibroblast.to_numpy(),s.nCount.to_numpy(),s.major_mean__Epithelial.to_numpy()];z=np.column_stack([rankdata(x)/n for x in v]);coords=s[['array_row','array_col']].to_numpy(int);lookup={tuple(x):i for i,x in enumerate(coords)};edges=np.array([(i,lookup[(r+dr,c+dc)]) for i,(r,c) in enumerate(coords) for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)] if (r+dr,c+dc) in lookup]);A=sp.csr_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(n,n));y=s.five_UCell.to_numpy();ll=new.LLPS77_UCell.to_numpy()
 for frac in [.05,.2]:
  for center in ['five','HH']:
   rng=np.random.default_rng(np.random.SeedSequence([123,zlib.crc32(sid.encode()),int(frac*100),int(center=='HH')]));h=y>np.quantile(y,1-frac);l=ll>np.quantile(ll,1-frac);c=~h if center=='five' else h&~l;h=h if center=='five' else h&l;a,b=match(h,c,z,rng);sm=[]
   for name,x in zip(['Fibroblast','UMI','Epithelial'],v):
    x=np.log1p(x);den=np.sqrt((np.var(x[a],ddof=1)+np.var(x[b],ddof=1))/2) if len(a)>1 else np.nan;d=(x[a].mean()-x[b].mean())/den if len(a)>1 and den>0 else np.nan;sm.append(d);balances.append(dict(**info,center=center,top_fraction=frac,covariate=name,n_pairs=len(a),n_high=h.sum(),match_rate=len(a)/max(h.sum(),1),SMD_after=d))
   for i,j in zip(a,b):pair_rows.append(dict(section_id=sid,center=center,top_fraction=frac,high_spot_id=s.spot_id.iloc[i],control_spot_id=s.spot_id.iloc[j]))
   balanced=len(a)>=10 and len(a)/max(h.sum(),1)>=.5 and np.all(np.abs(sm)<=.1)
   if not len(a):continue
   dh=distance(A,a);dc=distance(A,b)
   for ring in range(4):
    for cell in ['Endothelial','Myeloid']:
     x=np.log1p(s['major_mean__'+cell].to_numpy());ah=dh==ring;bc=dc==ring;out.append(dict(**info,center=center,matching='Fib_UMI_Epi',top_fraction=frac,variable=cell,ring=ring,effect=x[ah].mean()-x[bc].mean() if ah.any() and bc.any() else np.nan,balanced=balanced,n_pairs=len(a),n_high=ah.sum(),n_control=bc.sum(),overlap=(ah&bc).sum()))
base=pd.read_csv(O/'tables/rings.tsv',sep='\t')
if 'top_fraction' in base:base=base[base.top_fraction==.1]
base['top_fraction']=.1;pd.concat([base,pd.DataFrame(out)],ignore_index=True).to_csv(O/'tables/rings.tsv',sep='\t',index=False);pd.DataFrame(balances).to_csv(O/'tables/threshold_matching_balance.tsv',sep='\t',index=False)
print('Ring threshold sensitivity complete')

pd.DataFrame(pair_rows).to_parquet(O/'03_Fivegene_controls/threshold_matched_pairs.parquet',index=False)
