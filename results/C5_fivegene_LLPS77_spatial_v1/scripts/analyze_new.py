from pathlib import Path
import json,sys,gzip,time
import numpy as np,pandas as pd,scipy.sparse as sp
from scipy.stats import rankdata,spearmanr
from scipy.spatial import cKDTree
from scipy.signal import fftconvolve
from scipy.sparse.csgraph import connected_components
from collections import deque
O=Path('results/C5_fivegene_LLPS77_spatial_v1');P=Path('results/C5_spatial_exploratory_all63_v1');rng=np.random.default_rng(123)
labels=pd.read_csv(O/'00_input_audit/spatial_tissue_patient_audit.tsv',sep='\t');sets=json.loads((O/'00_input_audit/scoring_sets.json').read_text());five=['MMP11','COL11A1','INHBA','COL5A1','CTHRC1']
T=O/'tables'
if len(sys.argv)>1 and sys.argv[1]=='pilot':
 labels=labels[labels.section_id.isin(pd.read_csv(O/'00_input_audit/pilot_prespecified.tsv',sep='\t').section_id)];T=O/'00_input_audit/pilot_statistics';T.mkdir(exist_ok=True)
def rho(a,b):
 a=np.asarray(a);b=np.asarray(b);m=np.isfinite(a)&np.isfinite(b)
 return float(spearmanr(a[m],b[m]).statistic) if m.sum()>3 and np.std(a[m])>0 and np.std(b[m])>0 else np.nan
def resid(y,cs):
 z=np.column_stack([np.ones(len(y))]+[rankdata(c)/len(y) for c in cs]);v=rankdata(y)/len(y);return v-z@np.linalg.lstsq(z,v,rcond=None)[0]
def mi(A,v):
 v=np.asarray(v)-np.mean(v);return float(len(v)/A.sum()*(v@(A@v))/(v@v)) if v@v>0 and A.sum()>0 else np.nan
def dist(A,mask):
 d=np.full(A.shape[0],np.nan);q=deque(np.flatnonzero(mask).tolist());d[list(q)]=0
 while q:
  i=q.popleft()
  for j in A.indices[A.indptr[i]:A.indptr[i+1]]:
   if not np.isfinite(d[j]):d[j]=d[i]+1;q.append(j)
 return d
def match(h,c,cs):
 z=np.column_stack([rankdata(v)/len(v) for v in cs]);ids=np.flatnonzero(c);pairs=[];used=set()
 if not len(ids):return np.array([],int),np.array([],int)
 tree=cKDTree(z[ids])
 for i in rng.permutation(np.flatnonzero(h)):
  js=np.asarray(tree.query_ball_point(z[i],np.sqrt(len(cs))*.1),int);js=js[np.argsort(np.sum((z[ids[js]]-z[i])**2,axis=1),kind='stable')]
  for j in js:
   k=ids[j]
   if k not in used and max(abs(z[k]-z[i]))<=.1:pairs.append((i,k));used.add(k);break
 return tuple(np.array(pairs).T) if pairs else (np.array([],int),np.array([],int))
def smd(v,a,b):
 if len(a)<2 or len(b)<2:return np.nan
 den=np.sqrt((np.var(v[a],ddof=1)+np.var(v[b],ddof=1))/2)
 return (np.mean(v[a])-np.mean(v[b]))/den if den>0 else np.nan
def bh(x):
 x=np.asarray(x);r=np.full(len(x),np.nan);ii=np.flatnonzero(np.isfinite(x));ii=ii[np.argsort(x[ii])];r[ii]=np.minimum(1,np.minimum.accumulate((x[ii]*len(ii)/np.arange(1,len(ii)+1))[::-1])[::-1]);return r
ass=[];cl=[];bal=[];rr=[];qc=[];gene=[];co=[];nullrows=[];conf=[];ep=[];pairsout=[]
stores={'association':ass,'clustering':cl,'balance':bal,'rings':rr,'qc':qc,'gene_contribution':gene,'colocalization':co,'spatial_null':nullrows,'confounding':conf,'epithelial_distance':ep,'pairs':pairsout}
cache=O/'00_input_audit/statistic_checkpoints';cache.mkdir(exist_ok=True)
for row in labels.itertuples():
 sid=row.section_id
 import zlib
 rng=np.random.default_rng(np.random.SeedSequence([123,zlib.crc32(sid.encode())]))
 checkpoint=cache/(sid+'.json.gz')
 if checkpoint.exists():
  with gzip.open(checkpoint,'rt') as f:payload=json.load(f)
  for name,vals in payload.items():stores[name].extend(vals)
  print('REUSE statistics',sid,flush=True);continue
 starts={k:len(v) for k,v in stores.items()}
 dest=O/'01_LLPS77_spatial_scores'/f'{sid}.parquet'
 while not dest.exists():time.sleep(5)
 print(sid,flush=True);old=pd.read_parquet(P/'01_spatial_scores'/f'{sid}.parquet');new=pd.read_parquet(dest);s=old.merge(new[['spot_id']+[c for c in new if c not in old]],on='spot_id',validate='one_to_one');n=len(s);info=dict(section_id=sid,dataset=row.dataset,patient_uid=row.patient_uid,tissue=row.tissue)
 coords=s[['array_row','array_col']].to_numpy(int);lookup={tuple(v):i for i,v in enumerate(coords)};ee=np.array([(i,lookup[(r+dr,c+dc)]) for i,(r,c) in enumerate(coords) for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)] if (r+dr,c+dc) in lookup]);A=sp.csr_matrix((np.ones(len(ee)),(ee[:,0],ee[:,1])),shape=(n,n))
 fib=s.major_mean__Fibroblast.to_numpy();umi=s.nCount.to_numpy();nf=s.nFeature.to_numpy();epi=s.major_mean__Epithelial.to_numpy();cs=[fib,umi,nf,epi];y=s.five_UCell.to_numpy();ll=s.LLPS77_UCell.to_numpy();ry=resid(y,cs);rl=resid(ll,cs);s['five_adjusted_residual']=ry;s['LLPS77_adjusted_residual']=rl
 for source in ['five','A26','core15']:
  for method in ['UCell','AUCell']:
   v=s[source+'_'+method].to_numpy()
   for lname in sets:
    w=s[lname+'_'+method].to_numpy();cv=cs
    for adjustment,mask in [('raw',np.ones(n,bool)),('adjusted',np.ones(n,bool)),('fibroblast_high',fib>=np.quantile(fib,.75)),('positive_five',y>0)]:
     a=resid(v[mask],[c[mask] for c in cv]) if adjustment=='adjusted' else v[mask];b=resid(w[mask],[c[mask] for c in cv]) if adjustment=='adjusted' else w[mask]
     ass.append(dict(**info,signature=source,method=method,llps_set=lname,adjustment=adjustment,effect=rho(a,b),n_spots=mask.sum(),coverage=s.LLPS77_covered_genes.iloc[0]/77))
 for source in ['five','A26','core15','LLPS77','LLPS77_common']:
  for method in ['UCell','AUCell']:
   name=source+'_'+method;v=s[name].to_numpy();I=mi(A,v);p=(1+sum(mi(A,rng.permutation(v))>=I for _ in range(199)))/200
   qc.append(dict(**info,score=name,zero_fraction=(v==0).mean(),median=float(np.median(v)),moran_I=I,p=p,method_rho=rho(s[source+'_UCell'],s[source+'_AUCell']),detected_median=s[source+'_detected_genes'].median(),covered=s[source+'_covered_genes'].iloc[0],n_spots=n))
   for cov in ['major_mean__Fibroblast','major_mean__Mural','major_mean__Epithelial','nCount','nFeature','C2C4_program_UCell']:conf.append(dict(**info,score=name,covariate=cov,effect=rho(v,s[cov])))
   for frac in [.05,.1,.2]:
    high=v>np.quantile(v,1-frac);nc,lab=connected_components(A[high][:,high]) if high.any() else (0,np.array([],int));sizes=np.bincount(lab);cl.append(dict(**info,score=name,top_fraction=frac,n_high=int(high.sum()),actual_fraction=high.mean(),n_components=nc,largest_fraction=max(sizes,default=0)/max(high.sum(),1),largest_size=max(sizes,default=0)))
 for source in ['five','A26','core15']:
  v=s[source+'_UCell'].to_numpy();pos=v>0;qc.append(dict(**info,score=source+'_positive_only',moran_I=mi(A[pos][:,pos],v[pos]) if pos.sum()>3 else np.nan,n_spots=int(pos.sum())))
  qc.append(dict(**info,score=source+'_detection_adjusted',moran_I=mi(A,resid(v,[s[source+'_detected_genes'].to_numpy()]+cs)),n_spots=n))
 for name,v in [('five_adjusted',ry),('LLPS77_adjusted',rl)]:qc.append(dict(**info,score=name,moran_I=mi(A,v),n_spots=n))
 llgenes=[g for g in sets['LLPS77'] if g+'_rank_contribution' in s];lltotal=sum(s[g+'_rank_contribution'].to_numpy() for g in llgenes);lh=ll>np.quantile(ll,.9)
 for g in llgenes:
  gene.append(dict(**info,gene=g,program='LLPS77',detection=(s[g+'_log1p_CP10k']>0).mean(),rank_share_high=s.loc[lh,g+'_rank_contribution'].sum()/max(lltotal[lh].sum(),1)))
 for g in five:
  loo=s['five_without_'+g+'_UCell'].to_numpy();high=y>np.quantile(y,.9);lh=loo>np.quantile(loo,.9);con=s[g+'_rank_contribution'].to_numpy();total=sum(s[k+'_rank_contribution'].to_numpy() for k in five)
  gene.append(dict(**info,gene=g,program='five',detection=(s[g+'_log1p_CP10k']>0).mean(),rank_share_high=con[high].sum()/max(total[high].sum(),1),LOO_rho=rho(y,loo),LOO_high_jaccard=(high&lh).sum()/max((high|lh).sum(),1),LOO_LLPS_rho=rho(loo,ll),LOO_LLPS_adjusted=rho(resid(loo,cs),rl)))
 for frac in [.05,.1,.2]:
  h=y>np.quantile(y,1-frac)
  for target in ['LLPS77','A26','core15','C2C4_program']:
   l=s[target+'_UCell'].to_numpy()>np.quantile(s[target+'_UCell'],1-frac);dd=dist(A,l);co.append(dict(**info,target=target,top_fraction=frac,n_high_five=h.sum(),n_high_target=l.sum(),n_HH=(h&l).sum(),jaccard=(h&l).sum()/max((h|l).sum(),1),fraction_five_overlapping=(h&l).sum()/max(h.sum(),1),mean_graph_distance=np.nanmean(dd[h]) if h.any() else np.nan))
 degree=np.asarray(A.sum(1)).ravel();neighbor=np.divide(A@ll,degree,out=np.full(n,np.nan),where=degree>0);co.append(dict(**info,target='continuous_neighbor_LLPS',effect=rho(y,neighbor),adjusted_effect=rho(ry,resid(np.nan_to_num(neighbor,nan=np.nanmedian(neighbor)),cs))))
 # Spatial translation null, preserving local spatial pattern on observed overlap; no toroidal boundary wrapping.
 observed=rho(ry,rl);offsets=[(dr,dc) for dr in range(-int(np.ptp(coords[:,0])),int(np.ptp(coords[:,0]))+1) for dc in range(-int(np.ptp(coords[:,1])),int(np.ptp(coords[:,1]))+1) if (dr or dc) and (dr-dc)%2==0];rng.shuffle(offsets);accepted=[]
 grid=np.zeros(tuple(np.ptp(coords,axis=0)+1));origin=coords.min(axis=0);local=coords-origin;grid[local[:,0],local[:,1]]=1;overlaps=np.rint(fftconvolve(grid,grid[::-1,::-1],mode='full')).astype(int)
 for dr,dc in offsets:
  if overlaps[dr+grid.shape[0]-1,dc+grid.shape[1]-1]<.5*n:continue
  ij=[(i,lookup[(r+dr,c+dc)]) for i,(r,c) in enumerate(coords) if (r+dr,c+dc) in lookup]
  if len(ij)<.5*n:continue
  a,b=np.array(ij).T;nv=rho(ry[a],rl[b]);ov=rho(ry[a],rl[a]);accepted.append((nv,ov));nullrows.append(dict(**info,dr=dr,dc=dc,overlap_fraction=len(a)/n,null_rho=nv,observed_same_overlap_rho=ov))
  if len(accepted)==199:break
 vals=np.array(accepted);p=(1+np.sum(np.abs(vals[:,0])>=np.abs(vals[:,1])))/(len(vals)+1) if len(vals) else np.nan
 co.append(dict(**info,target='translation_null_adjusted',n_shifts=len(vals),effect=observed,p_spatial_sensitivity=p))
 for center in ['five','HH']:
  h=y>np.quantile(y,.9);l=ll>np.quantile(ll,.9);c=~h if center=='five' else h&~l;h=h if center=='five' else h&l
  for matchtype,mc in [('Fib_UMI',[fib,umi]),('Fib_UMI_Epi',[fib,umi,epi])]:
   a,b=match(h,c,mc);sm=[]
   for cov,v in [('Fibroblast',fib),('UMI',umi),('Epithelial',epi)]:
    v=np.log1p(v);before=smd(v,np.flatnonzero(h),np.flatnonzero(c));after=smd(v,a,b);sm.append(after);bal.append(dict(**info,center=center,matching=matchtype,covariate=cov,n_high=h.sum(),n_pairs=len(a),match_rate=len(a)/max(h.sum(),1),SMD_before=before,SMD_after=after))
   balanced=len(a)>=10 and len(a)/max(h.sum(),1)>=.5 and np.all(np.abs(sm)<=.1)
   if not len(a):continue
   for i,j in zip(a,b):pairsout.append(dict(section_id=sid,center=center,matching=matchtype,high_spot_id=s.spot_id.iloc[i],control_spot_id=s.spot_id.iloc[j],balanced_section=balanced))
   rr.append(dict(**info,center=center,matching=matchtype,variable='LLPS77_UCell',ring=0,effect=np.mean(ll[a]-ll[b]),balanced=balanced,n_pairs=len(a)))
   mh=np.zeros(n,bool);mc_=mh.copy();mh[a]=True;mc_[b]=True;dh=dist(A,mh);dc=dist(A,mc_)
   for ring in range(4):
    for cell in ['Endothelial','Myeloid']:
     v=np.log1p(s['major_mean__'+cell].to_numpy());ah=dh==ring;bc=dc==ring;rr.append(dict(**info,center=center,matching=matchtype,variable=cell,ring=ring,effect=v[ah].mean()-v[bc].mean() if ah.any() and bc.any() else np.nan,balanced=balanced,n_pairs=len(a),n_high=ah.sum(),n_control=bc.sum(),overlap=(ah&bc).sum()))
 ef=epi/s[[c for c in s if c.startswith('major_mean__')]].sum(axis=1).to_numpy()
 for threshold in [.2,.3,.4]:
  d=dist(A,ef>=threshold);mask=np.isfinite(d)
  if threshold==.3:s['epithelial_graph_distance']=d
  for adj,cc in [('raw',[]),('Fib_UMI_nFeature',[fib,umi,nf]),('Fib_UMI_nFeature_Epi',cs)]:
   v=resid(y[mask],[v[mask] for v in cc]) if cc else y[mask];w=resid(d[mask],[v[mask] for v in cc]) if cc else d[mask];ep.append(dict(**info,threshold=threshold,adjustment=adj,effect=rho(v,w),n_spots=mask.sum(),seed_fraction=(ef>=threshold).mean()))
 s['HH_top10']=(y>np.quantile(y,.9))&(ll>np.quantile(ll,.9))
 cols=['spot_id','x','y','five_UCell','five_AUCell','A26_UCell','core15_UCell','C2C4_program_UCell','LLPS77_UCell','LLPS77_AUCell','LLPS77_common_UCell','major_mean__Fibroblast','major_mean__Epithelial','nCount','nFeature','five_adjusted_residual','LLPS77_adjusted_residual','epithelial_graph_distance','HH_top10']+[g+'_log1p_CP10k' for g in five]
 s[cols].to_csv(O/'figures'/f'{sid}_plot_data.tsv.gz',sep='\t',index=False)
 with gzip.open(str(checkpoint)+'.tmp','wt') as f:json.dump({k:v[starts[k]:] for k,v in stores.items()},f,default=lambda x:x.item() if hasattr(x,'item') else str(x))
 Path(str(checkpoint)+'.tmp').replace(checkpoint)
for name,rows in [('association',ass),('clustering',cl),('balance',bal),('rings',rr),('qc',qc),('gene_contribution',gene),('colocalization',co),('spatial_null',nullrows),('confounding',conf),('epithelial_distance',ep)]:
 df=pd.DataFrame(rows)
 if name=='qc':
  df.loc[~np.isfinite(df.moran_I),'p']=np.nan;df['FDR_within_score']=df.groupby('score').p.transform(bh)
 if name=='colocalization':
  df.loc[(df.target=='translation_null_adjusted')&(df.n_shifts<19),'p_spatial_sensitivity']=np.nan;df['FDR_spatial_sensitivity']=bh(df.p_spatial_sensitivity)
 df.to_csv(T/f'{name}.tsv',sep='\t',index=False)
pd.DataFrame(pairsout).to_parquet(T/'matched_spot_pairs.parquet',index=False)
