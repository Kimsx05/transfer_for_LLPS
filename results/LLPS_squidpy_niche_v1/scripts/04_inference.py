from pathlib import Path
import json,itertools
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
O=Path(__file__).resolve().parents[1];D=O/'05_C5_LLPS_enrichment';rng=np.random.default_rng(123)
d=pd.read_parquet(O/'00_input_audit/audited_spots.parquet');lab=pd.read_parquet(O/'03_clustering/spot_labels.parquet');d=d.merge(lab.drop(columns=['section_id','dataset']),on='spot_id',validate='one_to_one');cov=pd.read_parquet(O/'02_neighborhood_composition/coverage.parquet');d=d.merge(cov,on='spot_id',validate='one_to_one');d['niche']=d.recommended_niche.astype(str)
maj=[c for c in d if c.startswith('major_mean__')];sub=[c for c in d if c.startswith('sub_mean__')];tot=d[maj].sum(axis=1);d['Fibroblast_abundance']=d.major_mean__Fibroblast;d['Endothelial_abundance']=d.major_mean__Endothelial;d['Fibroblast_fraction']=d.major_mean__Fibroblast/tot;d['Endothelial_fraction']=d.major_mean__Endothelial/tot
features=['five_UCell','LLPS77_UCell','A26_UCell','Fibroblast_abundance','Endothelial_abundance'];overlap=[c for c in d if c.startswith('LLPS77_without_') and c.endswith('_UCell')];features+=overlap
for f in features:
 d[f+'__z']=d.groupby('section_id')[f].transform(lambda x:(x-x.mean())/x.std(ddof=0) if x.std(ddof=0)>0 else x*0)
 d[f+'__percentile']=d.groupby('section_id')[f].rank(pct=True,method='average')
for f in ['five_UCell','LLPS77_UCell']:
 d[f+'__adjusted']=np.nan
 for sid,g in d.groupby('section_id'):
  v=np.column_stack([np.log1p(g.Fibroblast_abundance),g.Fibroblast_fraction]);v=(v-v.mean(0))/np.maximum(v.std(0),1e-8);X=np.column_stack([np.ones(len(g)),v,v*v]);y=g[f+'__z'].to_numpy();d.loc[g.index,f+'__adjusted']=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
scorecols=[c for f in features for c in [f,f+'__z',f+'__percentile']]+['five_UCell__adjusted','LLPS77_UCell__adjusted']
summary=d.groupby(['section_id','dataset','patient_uid','niche'],dropna=False)[scorecols].agg(['mean','median','count']);summary.columns=['__'.join(c) for c in summary.columns];summary=summary.reset_index();summary.to_csv(D/'section_level_niche_scores.tsv',sep='\t',index=False)
summary[summary.patient_uid.notna()].groupby(['dataset','patient_uid','niche']).mean(numeric_only=True).reset_index().to_csv(D/'patient_level_niche_scores.tsv',sep='\t',index=False)
contr=[];adj=[]
for sid,g in d.groupby('section_id'):
 info=dict(section_id=sid,dataset=g.dataset.iloc[0],patient_uid=g.patient_uid.iloc[0]);fib=g.Fibroblast_abundance.to_numpy();frac=g.Fibroblast_fraction.to_numpy()
 for niche in sorted(d.niche.unique(),key=int):
  a=g.niche.eq(niche).to_numpy()
  for threshold in [None,.5,.75,.9]:
   mask=np.ones(len(g),bool) if threshold is None else fib>=np.quantile(fib,threshold)
   for minimum in [10,20,30]:
    aa=a&mask;bb=~a&mask
    for f in scorecols:
     ok=aa.sum()>=minimum and bb.sum()>=minimum
     contr.append(dict(**info,niche=niche,feature=f,model='unadjusted' if threshold is None else 'fibroblast_rich',threshold=threshold,min_spots=minimum,n_niche=int(aa.sum()),n_rest=int(bb.sum()),effect=float(g.loc[aa,f].mean()-g.loc[bb,f].mean()) if ok else np.nan))
  # Shared-support regression for each niche-versus-rest, both absolute and relative fibroblast.
  if a.sum()<20 or (~a).sum()<20:continue
  lo=max(np.quantile(fib[a],.05),np.quantile(fib[~a],.05));hi=min(np.quantile(fib[a],.95),np.quantile(fib[~a],.95));keep=(fib>=lo)&(fib<=hi);support=min(np.mean(keep[a]),np.mean(keep[~a]));v=np.column_stack([np.log1p(fib[keep]),frac[keep]]);v=(v-v.mean(0))/np.maximum(v.std(0),1e-8);X=np.column_stack([np.ones(keep.sum()),a[keep].astype(float),v,v*v]);cond=np.linalg.cond(X) if len(X) else np.inf
  valid=support>=.2 and (a&keep).sum()>=20 and (~a&keep).sum()>=20 and cond<1e6
  for f in ['five_UCell','LLPS77_UCell','A26_UCell']:
   beta=np.linalg.lstsq(X,g.loc[keep,f+'__z'],rcond=None)[0][1] if valid else np.nan
   adj.append(dict(**info,niche=niche,feature=f,effect=beta,support_fraction=support,common_low=lo,common_high=hi,design_condition_number=cond,status='estimable' if valid else 'insufficient_overlap_or_collinearity',n_niche=int((a&keep).sum()),n_rest=int((~a&keep).sum())))
 print('inference',sid,flush=True)
c=pd.DataFrame(contr);c.to_csv(D/'section_level_contrasts.tsv',sep='\t',index=False);ar=pd.DataFrame(adj);ar.to_csv(D/'niche_fibroblast_adjusted_results.tsv',sep='\t',index=False)
def inference(frame,keys):
 out=[]
 for key,g in frame.groupby(keys,dropna=False):
  info=dict(zip(keys,key if isinstance(key,tuple) else (key,)))
  for unit in ['section','patient']:
   h=g if unit=='section' else g[g.patient_uid.notna()].groupby(['dataset','patient_uid'],as_index=False).effect.mean()
   x=h.effect.dropna().to_numpy();n=len(x)
   if n<3:lo=hi=p=np.nan
   else:
    boot=np.mean(rng.choice(x,(1000,n),replace=True),axis=1);lo,hi=np.quantile(boot,[.025,.975]);p=stats.wilcoxon(x).pvalue if np.any(x!=0) else 1.
   out.append(info|dict(unit=unit,n_units=n,effect=np.mean(x) if n else np.nan,ci_low=lo,ci_high=hi,p=p,positive_fraction=np.mean(x>0) if n else np.nan))
 out=pd.DataFrame(out);out['FDR']=np.nan
 family=[k for k in keys if k!='niche']+['unit']
 for _,g in out.groupby(family,dropna=False):
  idx=g.index[g.p.notna()];out.loc[idx,'FDR']=multipletests(out.loc[idx,'p'],method='fdr_bh')[1]
 return out
s=inference(c,['niche','feature','model','threshold','min_spots']);s.to_csv(D/'niche_score_statistics.tsv',sep='\t',index=False)
for f,name in [('five_UCell','C5_fivegene'),('LLPS77_UCell','LLPS77'),('A26_UCell','C5_A26')]:s[s.feature.str.startswith(f)].to_csv(D/f'niche_{name}_enrichment.tsv',sep='\t',index=False)
inference(ar,['niche','feature']).to_csv(D/'fibroblast_adjusted_summary.tsv',sep='\t',index=False)
# Gene overlap is checked using frozen source sets, not labels.
sets=json.loads((O.parent/'C5_spatial_exploratory_all63_v1/00_inputs_and_QC/frozen_signatures.json').read_text());ll=json.loads((O.parent/'C5_fivegene_LLPS77_spatial_v1/00_input_audit/scoring_sets.json').read_text())['LLPS77']
pd.DataFrame([dict(signature=k,n_overlap=len(set(sets[k])&set(ll)),overlap_genes=';'.join(sorted(set(sets[k])&set(ll))),sensitivity='reused existing without-overlap score' if set(sets[k])&set(ll) else 'not applicable') for k in ['five','A26']]).to_csv(D/'niche_signature_overlap_audit.tsv',sep='\t',index=False)
d.to_parquet(O/'04_niche_characterization/spot_metadata.parquet',index=False,compression='zstd')
# Central and neighborhood composition, pooled and section/patient balanced.
rec=pd.read_csv(O/'03_clustering/recommended_models.tsv',sep='\t');radius=int(rec.loc[rec.role=='recommended','radius'].iloc[0]);comp=[]
for kind,cols in [('major',maj),('subtype',sub)]:
 for scale in ['center','neighborhood']:
  v=d[cols].to_numpy();v=v/v.sum(1)[:,None] if scale=='center' else np.load(O/'02_neighborhood_composition'/f'{kind}_R{radius}.npy')
  t=pd.DataFrame(v,columns=cols);t[['section_id','patient_uid','niche']]=d[['section_id','patient_uid','niche']]
  for balance in ['pooled','section','patient']:
   if balance=='pooled':z=t.groupby('niche')[cols].mean()
   else:
    u=t.groupby(['section_id','patient_uid','niche'],dropna=False)[cols].mean().reset_index()
    if balance=='patient':u=u[u.patient_uid.notna()].groupby(['patient_uid','niche'])[cols].mean().reset_index()
    z=u.groupby('niche')[cols].mean()
   for ni,row in z.iterrows():
    for cell,value in row.items():comp.append(dict(niche=ni,kind=kind,scale=scale,balance=balance,celltype=cell.split('__',1)[1],fraction=value))
pd.DataFrame(comp).to_csv(O/'04_niche_characterization/niche_composition.tsv',sep='\t',index=False)
prop=d.groupby(['section_id','dataset','patient_uid','tissue','source_sample_type','niche'],dropna=False).size().rename('n_spots').reset_index();prop['proportion']=prop.n_spots/prop.groupby('section_id').n_spots.transform('sum');prop.to_csv(O/'08_clinical_association/section_niche_proportions.tsv',sep='\t',index=False)
# Full zero-filled proportions needed for clinical absence, absent scores remain NA.
index=pd.MultiIndex.from_product([d.section_id.unique(),d.niche.unique()],names=['section_id','niche']);pr=prop.set_index(['section_id','niche'])[['n_spots','proportion']].reindex(index,fill_value=0).reset_index();meta=d[['section_id','dataset','patient_uid','tissue','source_sample_type']].drop_duplicates();pr=pr.merge(meta,on='section_id',validate='many_to_one');pr.to_csv(O/'08_clinical_association/section_niche_proportions_complete.tsv',sep='\t',index=False)
(O/'05_C5_LLPS_enrichment/STATISTICAL_METHODS.md').write_text('Section paired contrasts compare a niche to all other niches within the same section; minimum10/20/30 spots in both groups. Raw, within-section z and percentile reported. Fibroblast-rich thresholds are section P50/P75/P90 of posterior mean. Adjusted effects: section-specific niche indicator regression on within-section score z, log1p Fibroblast abundance and relative fraction (linear+quadratic), restricted to overlapping 5th–95th percentile abundance support; at least20 spots per group, at least20% support and condition number<1e6. Failed support returns NA. Section effects averaged within verified patients before patient inference. 1000 bootstrap unit-resampled mean95%CI; signed-rank P; BH across niches within each feature/model/threshold/minimum/unit family. Section inference may contain correlated repeated patients; patient results prioritized where known. CIs are nominal and not multiplicity adjusted. Section-specific regressions absorb section/dataset intercepts but do not remove all dataset-specific confounding.\n')
