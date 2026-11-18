"""Post-diagnostic sensitivity prompted by epithelial imbalance, not signature optimization."""
from pathlib import Path
import numpy as np,pandas as pd,pyarrow as pa
from scipy.stats import rankdata,spearmanr
from scipy.spatial import cKDTree
pa.set_cpu_count(1);pa.set_io_thread_count(1)
O=Path('results/C5_spatial_exploratory_all63_v1');labels=pd.read_csv(O/'00_inputs_and_QC/spatial_tissue_patient_audit.tsv',sep='\t');rows=[];distrows=[];rng=np.random.default_rng(123)
def residual(y,c):
 X=np.column_stack([np.ones(len(y))]+[rankdata(x)/len(y) for x in c]);v=rankdata(y)/len(y);return v-X@np.linalg.lstsq(X,v,rcond=None)[0]
for r in labels.itertuples():
 s=pd.read_parquet(O/'01_spatial_scores'/f'{r.section_id}.parquet');n=len(s);y=s.A26_UCell.to_numpy();high=y>np.quantile(y,.9);fib=s.major_mean__Fibroblast.to_numpy();dep=s.nCount.to_numpy();epi=s.major_mean__Epithelial.to_numpy();nf=s.nFeature.to_numpy();rank=np.column_stack([rankdata(x)/n for x in [fib,dep,epi]]);idx=np.flatnonzero(~high);tree=cKDTree(rank[idx]);used=set();pairs=[]
 for i in rng.permutation(np.flatnonzero(high)):
  js=tree.query_ball_point(rank[i],np.sqrt(3)*.1);js=sorted(js,key=lambda j:np.sum((rank[idx[j]]-rank[i])**2))
  for j in js:
   k=idx[j]
   if k not in used and max(abs(rank[k]-rank[i]))<=.1:pairs.append((i,k));used.add(k);break
 info=dict(section_id=r.section_id,dataset=r.dataset,patient_uid=r.patient_uid,n_matched=len(pairs),match_rate=len(pairs)/max(high.sum(),1))
 if pairs:
  a,b=np.array(pairs).T
  for name,v in [('Fibroblast',fib),('nCount',dep),('Epithelial',epi)]:
   v=np.log1p(v);info['SMD_'+name]=(v[a].mean()-v[b].mean())/max(v.std(),1e-12)
  for name in ['Endothelial','Myeloid']:
   v=np.log1p(s['major_mean__'+name].to_numpy());rows.append(dict(**info,variable=name,effect=v[a].mean()-v[b].mean(),status='EXPLORATORY' if len(a)>=10 else 'INSUFFICIENT_MATCHES'))
 d=pd.read_csv(O/'02_spatial_maps'/f'{r.section_id}_plot_data.tsv.gz',sep='\t').epithelial_graph_distance.to_numpy();ok=np.isfinite(d)
 for covname,covs in [('Fib_UMI_features',[fib,dep,nf]),('Fib_UMI_features_Epi',[fib,dep,nf,epi])]:
  controls=[x[ok] for x in covs]
  effect=spearmanr(residual(y[ok],controls),residual(d[ok],controls)).statistic if ok.sum()>10 and np.std(d[ok])>0 else np.nan
  distrows.append(dict(section_id=r.section_id,dataset=r.dataset,patient_uid=r.patient_uid,controls=covname,n_finite=ok.sum(),partial_spearman=effect))
q=pd.DataFrame(rows);q.to_csv(O/'06_spatial_relationships/epithelial_matched_sensitivity.tsv',sep='\t',index=False);pd.DataFrame(distrows).to_csv(O/'06_spatial_relationships/epithelial_distance_partial.tsv',sep='\t',index=False)
summary=[]
for name,z in q[q.status=='EXPLORATORY'].groupby('variable'):
 for subset,t in [('all_matched',z),('balanced_50pct',z[(z.match_rate>=.5)&(z[['SMD_Fibroblast','SMD_nCount','SMD_Epithelial']].abs().max(axis=1)<=.1)])]:
  for ds in ['ALL']+sorted(z.dataset.unique().tolist()):
   w=t if ds=='ALL' else t[t.dataset==ds];p=w.dropna(subset=['patient_uid']).groupby(['dataset','patient_uid']).effect.mean().to_numpy();n=len(p);boot=np.median(rng.choice(p,(1000,n),replace=True),axis=1) if n>=3 else np.array([np.nan]);summary.append(dict(variable=name,subset=subset,dataset=ds,n_sections=len(w),median_section_effect=w.effect.median(),n_patients=n,median_patient_effect=np.median(p) if n else np.nan,ci_low=np.quantile(boot,.025),ci_high=np.quantile(boot,.975)))
pd.DataFrame(summary).to_csv(O/'07_dataset_heterogeneity/epithelial_matching_sensitivity_summary.tsv',sep='\t',index=False)
