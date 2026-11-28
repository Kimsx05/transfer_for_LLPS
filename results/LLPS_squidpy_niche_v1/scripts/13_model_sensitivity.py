from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
O=Path(__file__).resolve().parents[1];d=pd.read_parquet(O/'04_niche_characterization/spot_metadata.parquet');models=['recommended_niche','subtype_best_niche','R2_representative_niche','R3_representative_niche','kmeans_niche','native_R1','native_R2','native_R3'];cols=[c for c in d if c.startswith('major_mean__')];X=d[cols].to_numpy();X=X/X.sum(1)[:,None];comp=[];contr=[]
for model in models:
 for kind in ['major','subtype']:
  cs=[c for c in d if c.startswith('major_mean__' if kind=='major' else 'sub_mean__')];values=d[cs].to_numpy();values=values/values.sum(1)[:,None]
  t=pd.DataFrame(values,columns=[c.split('__')[1] for c in cs]);t[['section_id','niche']]=d[['section_id',model]].rename(columns={model:'niche'});z=t.groupby(['section_id','niche']).mean().groupby('niche').mean()
  for niche,row in z.iterrows():
   for cell,value in row.items():comp.append(dict(model=model,kind=kind,niche=niche,celltype=cell,fraction=value))
 for sid,g in d.groupby('section_id'):
  info=dict(model=model,section_id=sid,dataset=g.dataset.iloc[0],patient_uid=g.patient_uid.iloc[0]);N=len(g)
  for niche,h in g.groupby(model):
   n=len(h)
   if n<20 or N-n<20:continue
   for f in ['five_UCell__z','LLPS77_UCell__z','Endothelial_abundance__z']:
    effect=h[f].mean()-(g[f].sum()-h[f].sum())/(N-n);contr.append(dict(**info,niche=niche,feature=f,effect=effect,n_niche=n,n_rest=N-n))
pd.DataFrame(comp).to_csv(O/'04_niche_characterization/model_cell_composition.tsv',sep='\t',index=False)
c=pd.DataFrame(contr);c.to_csv(O/'05_C5_LLPS_enrichment/model_section_contrasts.tsv',sep='\t',index=False);out=[];rng=np.random.default_rng(123)
for key,g in c.groupby(['model','niche','feature']):
 for unit in ['section','patient']:
  h=g if unit=='section' else g[g.patient_uid.notna()].groupby('patient_uid').effect.mean().reset_index();v=h.effect.dropna().to_numpy();n=len(v)
  lo,hi=np.quantile(rng.choice(v,(1000,n),replace=True).mean(1),[.025,.975]) if n>=3 else (np.nan,np.nan)
  p=stats.wilcoxon(v).pvalue if n>=3 and np.any(v!=0) else 1 if n>=3 else np.nan
  out.append(dict(zip(['model','niche','feature'],key))|dict(unit=unit,n_units=n,effect=v.mean() if n else np.nan,ci_low=lo,ci_high=hi,p=p,positive_fraction=np.mean(v>0) if n else np.nan))
z=pd.DataFrame(out);z['FDR']=np.nan
for _,g in z.groupby(['model','feature','unit']):
 idx=g.index[g.p.notna()];z.loc[idx,'FDR']=multipletests(z.loc[idx,'p'],method='fdr_bh')[1]
z.to_csv(O/'05_C5_LLPS_enrichment/model_enrichment_sensitivity.tsv',sep='\t',index=False)
inv=[]
for model in models:inv.append(dict(label_column=model,n_niches=d[model].nunique(),scope='all 265002 spots',composition='continuous' if not model.startswith('native') else 'dominant_cell_type',role='fixed primary' if model=='recommended_niche' else 'sensitivity/control'))
pd.DataFrame(inv).to_csv(O/'03_clustering/retained_model_inventory.tsv',sep='\t',index=False)
(O/'05_C5_LLPS_enrichment/MODEL_SENSITIVITY.md').write_text('Added to characterize model dependence after composition-only diagnostics showed partition uncertainty. Fixed labels are unchanged; no selection using score enrichment. Section and verified-patient unadjusted z-score contrasts, minimum20 per group, same bootstrap/signed-rank/BH families as main inference. Numeric niche IDs across models are not equivalent. This sensitivity does not replace the main fibroblast-adjusted analysis. All methods retained, including negative and low-support results.\n')
quant=[]
for radius in [1,2,3]:
 X=np.load(O/'02_neighborhood_composition'/f'dominant_R{radius}.npy');_,counts=np.unique(X,axis=0,return_counts=True);labels=d[f'native_R{radius}'];sizes=labels.value_counts();quant.append(dict(radius=radius,n_unique_profiles=len(counts),fraction_spots_in_profiles_with_at_least_30_copies=float(counts[counts>=30].sum()/len(X)),n_clusters=len(sizes),median_cluster_size=float(sizes.median()),n_clusters_below_100_spots=int((sizes<100).sum()),fraction_spots_in_clusters_below_100=float(sizes[sizes<100].sum()/len(X))))
pd.DataFrame(quant).to_csv(O/'03_clustering/native_quantization_audit.tsv',sep='\t',index=False)
