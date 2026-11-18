import pyarrow as pa
pa.set_cpu_count(1);pa.set_io_thread_count(1)
from pathlib import Path
import json,time
import numpy as np,pandas as pd,scipy.sparse as sp
from scipy.stats import rankdata,spearmanr
from scipy.sparse.csgraph import connected_components,breadth_first_order
from scipy.spatial import cKDTree
O=Path('results/C5_spatial_exploratory_all63_v1');rng=np.random.default_rng(123)
labels=pd.read_csv(O/'00_inputs_and_QC/spatial_tissue_patient_audit.tsv',sep='\t');sigs=json.loads((O/'00_inputs_and_QC/frozen_signatures.json').read_text());stats=[];assoc=[];comp=[];rings=[];relations=[];balance=[];sections=[];limits={};scores=[k+'_'+m for k in sigs for m in ['UCell','AUCell']]
def rho(a,b):
 mask=np.isfinite(a)&np.isfinite(b)
 return spearmanr(np.asarray(a)[mask],np.asarray(b)[mask]).statistic if mask.sum()>3 and np.std(np.asarray(a)[mask])>0 and np.std(np.asarray(b)[mask])>0 else np.nan
def residual(y,c):
 X=np.column_stack([np.ones(len(y))]+[rankdata(v)/len(y) for v in c]);v=rankdata(y)/len(y);return v-X@np.linalg.lstsq(X,v,rcond=None)[0]
def dist(A,seeds):
 # BFS on observed graph: no pixel-unit assumption.
 from collections import deque
 d=np.full(A.shape[0],np.nan);q=deque(np.flatnonzero(seeds).tolist());d[list(q)]=0
 while q:
  i=q.popleft()
  for j in A.indices[A.indptr[i]:A.indptr[i+1]]:
   if not np.isfinite(d[j]):d[j]=d[i]+1;q.append(j)
 return d
def mi(A,y):
 v=y-np.mean(y);return float(len(v)/A.sum()*(v@(A@v))/(v@v)) if v@v>0 and A.sum()>0 else np.nan
def moran(A,y):
 val=mi(A,y);null=[mi(A,rng.permutation(y)) for _ in range(199)];return val,(1+sum(t>=val for t in null))/200 if np.isfinite(val) else np.nan
for row in labels.itertuples():
 sid=row.section_id
 while not (O/'01_spatial_scores'/f'{sid}.parquet').exists():time.sleep(5)
 time.sleep(1)
 s=pd.read_parquet(O/'01_spatial_scores'/f'{sid}.parquet');n=len(s);print(sid,flush=True)
 info=dict(section_id=sid,dataset=row.dataset,patient_uid=row.patient_uid,tissue=row.tissue)
 coords=s[['array_row','array_col']].to_numpy();lookup={tuple(v):i for i,v in enumerate(coords)};edges=[(i,lookup[(r+dr,c+dc)]) for i,(r,c) in enumerate(coords) for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)] if (r+dr,c+dc) in lookup]
 assert len(lookup)==n
 e=np.array(edges);A=sp.csr_matrix((np.ones(len(e)),(e[:,0],e[:,1])),shape=(n,n));fib=s.major_mean__Fibroblast.to_numpy();epi=s.major_mean__Epithelial.to_numpy();dep=s.nCount.to_numpy();nf=s.nFeature.to_numpy();y=s.A26_UCell.to_numpy();res=residual(y,[fib,dep,nf,epi]);s['A26_adjusted_rank_residual']=res
 for col in scores+['A26_adjusted_rank_residual']:
  val,p=moran(A,s[col].to_numpy());stats.append(dict(**info,score=col,moran_I=val,p_permutation=p,n_spots=n,graph_edges=len(e)//2,isolated_spots=int((np.diff(A.indptr)==0).sum())))
 for frac in [.05,.1,.2]:
  mask=y>np.quantile(y,1-frac);B=A[mask][:,mask];nc,lab=connected_components(B) if mask.any() else (0,np.array([],int));sizes=np.bincount(lab);comp.append(dict(**info,top_fraction=frac,actual_fraction=mask.mean(),n_high=mask.sum(),components=nc,largest_component=int(max(sizes,default=0)),largest_fraction_high=max(sizes,default=0)/max(mask.sum(),1)))
 for col in scores:
  v=s[col].to_numpy()
  for cov in ['major_mean__Fibroblast','major_mean__Mural','major_mean__Epithelial','nCount','nFeature','C2C4_program_UCell','core15_UCell','five_UCell','A26_AUCell']:
   assoc.append(dict(**info,score=col,variable=cov,spearman=rho(v,s[cov]),partial_spearman=rho(residual(v,[fib,dep,nf,epi]),residual(s[cov].to_numpy(),[fib,dep,nf,epi])) if cov not in ['major_mean__Fibroblast','major_mean__Epithelial','nCount','nFeature'] else np.nan))
 for cov,controls in [('Fibroblast',[dep,nf,epi]),('nCount',[fib,epi]),('nFeature',[fib,dep,epi])]:
  target={'Fibroblast':fib,'nCount':dep,'nFeature':nf}[cov];assoc.append(dict(**info,score='A26_conditional',variable=cov,spearman=rho(y,target),partial_spearman=rho(residual(y,controls),residual(target,controls))))
 high=y>np.quantile(y,.9);low=~high;ranks=np.column_stack([rankdata(fib)/n,rankdata(dep)/n]);loidx=np.flatnonzero(low);tree=cKDTree(ranks[low]);pairs=[];used=set()
 for i in rng.permutation(np.flatnonzero(high)):
  js=tree.query_ball_point(ranks[i],r=np.sqrt(2)*.1);js=sorted(js,key=lambda j:np.sum((ranks[loidx[j]]-ranks[i])**2))
  for j in js:
   k=loidx[j]
   if k not in used and np.max(abs(ranks[k]-ranks[i]))<=.1:pairs.append((i,k));used.add(k);break
 control=np.zeros(n,bool);matched=np.zeros(n,bool)
 if pairs:
  a,b=np.array(pairs).T;control[b]=True;matched[a]=True
  for name,v in [('Fibroblast',np.log1p(fib)),('nCount',np.log1p(dep)),('Epithelial',np.log1p(epi))]:balance.append(dict(**info,variable=name,n_matched=len(a),n_high=high.sum(),match_rate=len(a)/max(high.sum(),1),standardized_mean_difference=(v[a].mean()-v[b].mean())/max(v.std(),1e-12)))
  for name in ['A26_UCell','C2C4_program_UCell','major_mean__Endothelial','major_mean__Myeloid','major_mean__Epithelial']:
   v=s[name].to_numpy();v=np.log1p(v) if name.startswith('major') else v;relations.append(dict(**info,analysis='matched_same_spot',variable=name,effect=v[a].mean()-v[b].mean(),n_matched=len(a)))
  dh=dist(A,matched);dc=dist(A,control)
  for ring in range(4):
   for name in ['Endothelial','Myeloid']:
    v=np.log1p(s['major_mean__'+name].to_numpy());mh=dh==ring;mc=dc==ring;effect=v[mh].mean()-v[mc].mean() if mh.any() and mc.any() else np.nan
    rings.append(dict(**info,analysis='matched_seed_neighborhood',threshold=np.nan,ring=ring,variable=name,high_mean=v[mh].mean() if mh.any() else np.nan,control_mean=v[mc].mean() if mc.any() else np.nan,effect=effect,n_high=mh.sum(),n_control=mc.sum(),overlap=int((mh&mc).sum())))
 majorcols=[c for c in s if c.startswith('major_mean__')];ef=epi/s[majorcols].sum(axis=1).to_numpy()
 for threshold in [.2,.3,.4]:
  seeds=ef>=threshold;d=dist(A,seeds)
  if threshold==.3:s['epithelial_graph_distance']=d
  relations.append(dict(**info,analysis='epithelial_distance',variable=str(threshold),effect=rho(y,d),n_matched=np.isfinite(d).sum(),seed_fraction=seeds.mean(),q05_positive_fraction=(s.major_q05__Epithelial.to_numpy()[seeds]>0).mean() if seeds.any() else np.nan))
  for ring in range(5):
   mask=d==ring if ring<4 else d>=4
   rings.append(dict(**info,analysis='epithelial_distance',threshold=threshold,ring=ring,variable='A26_UCell',high_mean=y[mask].mean() if mask.any() else np.nan,control_mean=np.nan,effect=np.nan,n_high=mask.sum(),n_control=0,overlap=0))
 fh=fib>=np.quantile(fib,.75);fl=fib<=np.quantile(fib,.25);ch=s.C2C4_program_UCell.to_numpy()>s.C2C4_program_UCell.quantile(.9)
 sections.append(dict(**info,n_spots=n,A26_mean=y.mean(),A26_median=np.median(y),A26_fib_high_mean=y[fh].mean(),A26_fib_low_mean=y[fl].mean(),A26_high_in_fib_top_quartile=(high&fh).sum()/max(high.sum(),1),A26_high_in_fib_bottom_quartile=(high&fl).sum()/max(high.sum(),1),C2C4_high_overlap=(high&ch).sum()/max(high.sum(),1),C2C4_jaccard=(high&ch).sum()/max((high|ch).sum(),1),median_UMI=np.median(dep),median_features=np.median(nf),n_matched=len(pairs)))
 for col in ['major_mean__Fibroblast','major_mean__Epithelial','MMP11_log1p_CP10k','nCount']:limits.setdefault(col,[]).append(s[col].to_numpy())
 # Compact plotting table; no expression matrix copied.
 cols=['spot_id','x','y','A26_UCell','A26_AUCell','core15_UCell','five_UCell','MMP11_log1p_CP10k','major_mean__Fibroblast','major_mean__Epithelial','nCount','A26_UCell_percentile','C2C4_program_UCell','A26_adjusted_rank_residual','epithelial_graph_distance']
 s[cols].to_csv(O/'02_spatial_maps'/f'{sid}_plot_data.tsv.gz',sep='\t',index=False)
def bh(p):
 p=np.array(p,float);out=np.full(len(p),np.nan);ix=np.flatnonzero(np.isfinite(p));order=ix[np.argsort(p[ix])];out[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1]);return out
st=pd.DataFrame(stats);st['FDR_within_score_63']=st.groupby('score').p_permutation.transform(bh);st.to_csv(O/'03_spatial_clustering/moran.tsv',sep='\t',index=False)
pd.DataFrame(comp).to_csv(O/'03_spatial_clustering/components.tsv',sep='\t',index=False)
pd.DataFrame(assoc).to_csv(O/'04_fibroblast_and_depth_controls/correlations.tsv',sep='\t',index=False)
pd.DataFrame(balance).to_csv(O/'04_fibroblast_and_depth_controls/match_balance.tsv',sep='\t',index=False)
pd.DataFrame(rings).to_csv(O/'06_spatial_relationships/rings.tsv',sep='\t',index=False)
rr=pd.DataFrame(relations);rr.to_csv(O/'06_spatial_relationships/section_effects.tsv',sep='\t',index=False)
ss=pd.DataFrame(sections);ss=ss.merge(st[st.score=='A26_UCell'][['section_id','moran_I','FDR_within_score_63']],on='section_id');ss=ss.merge(pd.DataFrame(comp).query('top_fraction == 0.1')[['section_id','n_high']],on='section_id');ss['exploratory_clustering_flag']=(ss.moran_I>0)&(ss.FDR_within_score_63<.05)&(ss.n_high>=10);ss.to_csv(O/'tables/section_summary.tsv',sep='\t',index=False)
ss[['section_id','dataset','patient_uid','C2C4_high_overlap','C2C4_jaccard']].to_csv(O/'05_C4_C2_competition/overlap.tsv',sep='\t',index=False)
known=ss.dropna(subset=['patient_uid']).groupby(['dataset','patient_uid','tissue'],dropna=False).mean(numeric_only=True).reset_index();known.to_csv(O/'07_dataset_heterogeneity/patient_summary.tsv',sep='\t',index=False)
ss.groupby(['dataset','tissue']).agg(n_sections=('section_id','size'),n_spots=('n_spots','sum'),median_A26=('A26_median','median'),median_moran=('moran_I','median'),clustered=('exploratory_clustering_flag','sum')).reset_index().to_csv(O/'07_dataset_heterogeneity/dataset_tissue_summary.tsv',sep='\t',index=False)
summary=[]
for (analysis,var),group in rr.dropna(subset=['patient_uid','effect']).groupby(['analysis','variable']):
 pat=group.groupby(['dataset','patient_uid']).effect.mean().reset_index()
 for dataset in ['ALL_KNOWN_PATIENTS']+sorted(pat.dataset.unique().tolist())+['LEAVE_OUT:'+d for d in sorted(pat.dataset.unique())]:
  p=pat if dataset=='ALL_KNOWN_PATIENTS' else pat[pat.dataset!=dataset.split(':',1)[1]] if dataset.startswith('LEAVE_OUT:') else pat[pat.dataset==dataset]
  v=p.effect.to_numpy();n=len(v);boot=np.median(rng.choice(v,(1000,n),replace=True),axis=1) if n>=3 else np.array([np.nan]);summary.append(dict(analysis=analysis,variable=var,dataset=dataset,n_patients=n,median_effect=np.median(v) if n else np.nan,ci_low=np.nanquantile(boot,.025) if n>=3 else np.nan,ci_high=np.nanquantile(boot,.975) if n>=3 else np.nan,status='EXPLORATORY_PATIENT_BOOTSTRAP' if n>=3 else 'INSUFFICIENT_PATIENTS'))
pd.DataFrame(summary).to_csv(O/'07_dataset_heterogeneity/patient_effect_intervals.tsv',sep='\t',index=False)
(O/'00_inputs_and_QC/shared_color_limits.json').write_text(json.dumps({k:float(np.quantile(np.concatenate(v),.99)) for k,v in limits.items()},indent=2))
pd.concat([pd.read_csv(p,sep='\t') for p in sorted((O/'00_inputs_and_QC').glob('*_coverage.tsv'))]).to_csv(O/'00_inputs_and_QC/gene_coverage_detection.tsv',sep='\t',index=False)
