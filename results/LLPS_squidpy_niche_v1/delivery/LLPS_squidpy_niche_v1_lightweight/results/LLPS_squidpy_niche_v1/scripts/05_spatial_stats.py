from pathlib import Path
import itertools,json
import numpy as np,pandas as pd,scipy.sparse as sp
from scipy import stats
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
O=Path(__file__).resolve().parents[1];d=pd.read_parquet(O/'04_niche_characterization/spot_metadata.parquet');niches=sorted(d.niche.unique(),key=int);K=len(niches);ni={v:i for i,v in enumerate(niches)};rng=np.random.default_rng(123)
frag=[];pairs=[];corr=[];bias=[];components=[];models=[c for c in d if c.endswith('_niche')];native=pd.read_parquet(O/'03_clustering/native_labels.parquet');assert native.spot_id.tolist()==d.spot_id.tolist()
for c in [c for c in native if c.startswith('native_R')]:d[c]=native[c];models.append(c)
for model in models:
 for niche,g in d.groupby(model):
  counts=g.groupby('dataset').size();scount=g.groupby('section_id').size();bias.append(dict(model=model,niche=niche,n_spots=len(g),fraction=len(g)/len(d),n_sections=g.section_id.nunique(),n_datasets=g.dataset.nunique(),max_dataset_share=counts.max()/len(g),dominant_dataset=counts.idxmax(),max_section_share=scount.max()/len(g)))
for sid,g in d.groupby('section_id',sort=False):
 A=sp.load_npz(O/'01_spatial_graph'/f'{sid}_R1.npz');u,v=A.nonzero();info=dict(section_id=sid,dataset=g.dataset.iloc[0],patient_uid=g.patient_uid.iloc[0]);y=np.array([ni[c] for c in g.niche]);n=len(y)
 for model in models:
  labels=g[model].to_numpy();agree=np.mean(labels[u]==labels[v]) if len(u) else np.nan
  for niche in np.unique(labels):
   mask=labels==niche;nc,cl=connected_components(A[mask][:,mask]);sizes=np.bincount(cl);frag.append(dict(**info,model=model,niche=niche,n_spots=mask.sum(),n_components=nc,singleton_fraction=np.sum(sizes==1)/mask.sum(),median_component_size=np.median(sizes),largest_fraction=sizes.max()/mask.sum(),same_label_edge_fraction=agree))
   for size in sizes:components.append(dict(section_id=sid,model=model,niche=niche,component_size=size))
 obs=np.bincount(y[u]*K+y[v],minlength=K*K).reshape(K,K);sums=np.zeros((K,K));sq=np.zeros((K,K));rng=np.random.default_rng(123)
 for _ in range(1000):
  yy=rng.permutation(y);z=np.bincount(yy[u]*K+yy[v],minlength=K*K).reshape(K,K);sums+=z;sq+=z*z
 exp=sums/1000;sd=np.sqrt(np.maximum(sq/1000-exp*exp,0));present=np.bincount(y,minlength=K)>0
 for a,b in itertools.product(range(K),repeat=2):
  ok=present[a] and present[b];ratio=obs[a,b]/exp[a,b] if ok and exp[a,b]>0 else np.nan;zz=(obs[a,b]-exp[a,b])/sd[a,b] if ok and sd[a,b]>0 else np.nan
  pairs.append(dict(**info,niche_a=niches[a],niche_b=niches[b],observed=obs[a,b] if ok else np.nan,expected=exp[a,b] if ok else np.nan,sd=sd[a,b] if ok else np.nan,Z=zz,OE=ratio,direction='enrichment' if ratio>1 else 'depletion' if ratio<1 else 'NA',n_perms=1000))
 # Correlations: center vs neighborhood averages, within each section and niche.
 rec=pd.read_csv(O/'03_clustering/recommended_models.tsv',sep='\t');R=int(rec.loc[rec.role=='recommended','radius'].iloc[0]);C=sp.eye(n,format='csr')
 for _ in range(R):C=(C+C@A).astype(bool).astype('float32')
 den=np.asarray(C.sum(1)).ravel();features=['five_UCell','LLPS77_UCell','Endothelial_abundance'];vals=g[features].to_numpy()
 for scale,W in [('center',vals),('neighborhood',(C@vals)/den[:,None])]:
  fib=g.Fibroblast_abundance.to_numpy() if scale=='center' else (C@g.Fibroblast_abundance.to_numpy())/den
  for niche in ['ALL']+niches:
   mask=np.ones(n,bool) if niche=='ALL' else g.niche.eq(niche).to_numpy()
   for i,j in [(0,2),(1,2),(0,1)]:
    for method in ['pearson','spearman']:
     for adjustment in ['raw','fibroblast_adjusted']:
      a=W[mask,i];b=W[mask,j];f=fib[mask];N=len(a)
      if method=='spearman':a=stats.rankdata(a);b=stats.rankdata(b);f=stats.rankdata(f)
      if adjustment!='raw' and N:
       z=np.column_stack([np.ones(N),f]);a=a-z@np.linalg.lstsq(z,a,rcond=None)[0];b=b-z@np.linalg.lstsq(z,b,rcond=None)[0]
      r=float(stats.pearsonr(a,b).statistic) if N>=20 and np.std(a)>1e-12 and np.std(b)>1e-12 else np.nan
      corr.append(dict(**info,niche=niche,scale=scale,feature_a=features[i],feature_b=features[j],method=method,adjustment=adjustment,n_spots=N,r=r))
 print('spatial stats',sid,flush=True)
for file,rows in [('03_clustering/spatial_fragmentation_metrics.tsv',frag),('03_clustering/connected_component_sizes.tsv',components),('03_clustering/dataset_bias_audit.tsv',bias),('07_neighborhood_enrichment/section_pair_results.tsv',pairs),('06_C5_endothelial/section_correlations.tsv',corr)]:pd.DataFrame(rows).to_csv(O/file,sep='\t',index=False)
p=pd.DataFrame(pairs)
for sid,g in p.groupby('section_id'):
 for val in ['observed','expected','Z','OE']:g.pivot(index='niche_a',columns='niche_b',values=val).to_csv(O/'07_neighborhood_enrichment'/f'{sid}_{val}.tsv',sep='\t')
for balance in ['section','patient','dataset']:
 u=p
 if balance=='patient':u=p[p.patient_uid.notna()].groupby(['patient_uid','niche_a','niche_b'])[['OE','Z']].mean().reset_index()
 group=['niche_a','niche_b'] if balance!='dataset' else ['dataset','niche_a','niche_b']
 u.groupby(group).agg(mean_OE=('OE','mean'),median_OE=('OE','median'),mean_Z=('Z','mean'),n_units=('OE','count'),positive_fraction=('OE',lambda x:np.mean(x.dropna()>1) if x.notna().any() else np.nan)).reset_index().to_csv(O/'07_neighborhood_enrichment'/f'{balance}_balanced_summary.tsv',sep='\t',index=False)
c=pd.DataFrame(corr);out=[];keys=['niche','scale','feature_a','feature_b','method','adjustment'];rng=np.random.default_rng(123)
for key,g in c.groupby(keys):
 for unit in ['section','patient']:
  h=g if unit=='section' else g[g.patient_uid.notna()].groupby('patient_uid').r.mean().reset_index();x=h.r.dropna().to_numpy();n=len(x);z=np.arctanh(np.clip(x,-.999999,.999999));boot=np.tanh(rng.choice(z,(1000,n),replace=True).mean(1)) if n else np.array([np.nan]);lo,hi=np.quantile(boot,[.025,.975]);out.append(dict(zip(keys,key))|dict(unit=unit,n_units=n,r=np.tanh(z.mean()) if n else np.nan,ci_low=lo,ci_high=hi,positive_fraction=np.mean(x>0) if n else np.nan))
pd.DataFrame(out).to_csv(O/'06_C5_endothelial/correlation_summary.tsv',sep='\t',index=False)
c.groupby(['dataset']+keys).agg(mean_r=('r','mean'),n_sections=('r','count'),positive_fraction=('r',lambda x:np.mean(x.dropna()>0) if x.notna().any() else np.nan)).reset_index().to_csv(O/'06_C5_endothelial/dataset_correlations.tsv',sep='\t',index=False)
ag=[]
for a,b in itertools.combinations(models,2):ag.append(dict(model_a=a,model_b=b,ARI=adjusted_rand_score(d[a],d[b]),NMI=normalized_mutual_info_score(d[a],d[b])))
pd.DataFrame(ag).to_csv(O/'03_clustering/all_models_ARI_NMI.tsv',sep='\t',index=False)
(O/'07_neighborhood_enrichment/METHODS.md').write_text('Equivalent label-permutation adjacency count on symmetric directed entries of exact R1 graph (no diagonal). Within-section permutations1000 seed123, preserving niche counts. Expected and population SD from permutation moments. Absent pairs NA; zero-variance Z NA. O/E remains computable for present pairs with positive expectation even if Z undefined. Section/patient/dataset O/E averages and direction fractions are primary; average Z is descriptive, not pooled significance. Self-niche counts count each undirected edge twice, consistently in observed/expected. Niche definitions already incorporate spatial smoothing; adjacency is descriptive organization, not independent communication evidence.\n')
