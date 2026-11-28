from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
O=Path(__file__).resolve().parents[1];d=pd.read_parquet(O/'04_niche_characterization/spot_metadata.parquet');out=[];models=['recommended_niche','subtype_best_niche','R2_representative_niche','R3_representative_niche','kmeans_niche']
for model in models:
 for sid,g in d.groupby('section_id'):
  info=dict(model=model,section_id=sid,dataset=g.dataset.iloc[0],patient_uid=g.patient_uid.iloc[0]);fib=g.Fibroblast_abundance.to_numpy();frac=g.Fibroblast_fraction.to_numpy()
  for niche in g[model].unique():
   a=g[model].eq(niche).to_numpy()
   if a.sum()<20 or (~a).sum()<20:continue
   lo=max(np.quantile(fib[a],.05),np.quantile(fib[~a],.05));hi=min(np.quantile(fib[a],.95),np.quantile(fib[~a],.95));keep=(fib>=lo)&(fib<=hi);support=min(np.mean(keep[a]),np.mean(keep[~a]));v=np.column_stack([np.log1p(fib[keep]),frac[keep]]);v=(v-v.mean(0))/np.maximum(v.std(0),1e-8);X=np.column_stack([np.ones(keep.sum()),a[keep].astype(float),v,v*v]);cond=np.linalg.cond(X) if len(X) else np.inf;valid=support>=.2 and (a&keep).sum()>=20 and (~a&keep).sum()>=20 and cond<1e6
   for f in ['five_UCell','LLPS77_UCell']:
    beta=np.linalg.lstsq(X,g.loc[keep,f+'__z'],rcond=None)[0][1] if valid else np.nan;out.append(dict(**info,niche=str(niche),feature=f,effect=beta,support_fraction=support,condition_number=cond,status='estimable' if valid else 'insufficient_overlap_or_collinearity'))
c=pd.DataFrame(out);c.to_csv(O/'05_C5_LLPS_enrichment/adjusted_model_section_effects.tsv',sep='\t',index=False)
old=pd.read_csv(O/'05_C5_LLPS_enrichment/niche_fibroblast_adjusted_results.tsv',sep='\t',dtype={'niche':str});old=old[old.feature.isin(['five_UCell','LLPS77_UCell'])];check=c[c.model=='recommended_niche'].merge(old,on=['section_id','niche','feature'],suffixes=('_new','_old'),validate='one_to_one');assert len(check)==len(old);assert np.allclose(check.effect_new,check.effect_old,equal_nan=True,atol=1e-10)
rows=[];rng=np.random.default_rng(123)
for key,g in c.groupby(['model','niche','feature']):
 for unit in ['section','patient']:
  h=g if unit=='section' else g[g.patient_uid.notna()].groupby('patient_uid').effect.mean().reset_index();x=h.effect.dropna().to_numpy();n=len(x);lo,hi=np.quantile(rng.choice(x,(1000,n),replace=True).mean(1),[.025,.975]) if n>=3 else (np.nan,np.nan);p=stats.wilcoxon(x).pvalue if n>=3 and np.any(x!=0) else 1 if n>=3 else np.nan
  rows.append(dict(zip(['model','niche','feature'],key))|dict(unit=unit,n_units=n,effect=x.mean() if n else np.nan,ci_low=lo,ci_high=hi,p=p,positive_fraction=np.mean(x>0) if n else np.nan))
s=pd.DataFrame(rows);s['FDR']=np.nan
for _,g in s.groupby(['model','feature','unit']):
 idx=g.index[g.p.notna()];s.loc[idx,'FDR']=multipletests(s.loc[idx,'p'],method='fdr_bh')[1]
s.to_csv(O/'05_C5_LLPS_enrichment/adjusted_model_sensitivity.tsv',sep='\t',index=False)
(O/'05_C5_LLPS_enrichment/adjusted_model_validation.txt').write_text(f'PASS: {len(check)} frozen-main section coefficients independently reproduced within1e-10, including NA patterns. Same support and regression design applied to all four continuous control models; no labels changed.\n')
