from pathlib import Path
import pandas as pd,numpy as np,json,zlib
from scipy.stats import wilcoxon
O=Path('results/C5_fivegene_LLPS77_spatial_v1');rng=np.random.default_rng(123);out=[]
for table,keys in [('association',['signature','method','llps_set','adjustment']),('rings',['center','matching','variable','ring','top_fraction']),('rings_AUCell',['center','matching','variable','ring','top_fraction']),('epithelial_distance',['threshold','adjustment'])]:
 d=pd.read_csv(O/'tables'/f'{table}.tsv',sep='\t')
 if table.startswith('rings'):d=pd.concat([d.assign(subset='all_matched'),d[d.balanced==True].assign(subset='balanced')]);keys=keys+['subset']
 if table=='association':d=pd.concat([d.assign(subset='all_sections'),d[d.coverage>=.9].assign(subset='coverage_ge90')]);keys=keys+['subset']
 ds=d.groupby(keys+['dataset'],dropna=False).effect.agg(['median','count']).reset_index();ds.to_csv(O/'07_dataset_heterogeneity'/f'{table}_dataset.tsv',sep='\t',index=False)
 pats=d.dropna(subset=['patient_uid']).groupby(keys+['dataset','patient_uid'],dropna=False).effect.mean().reset_index();pats.to_csv(O/'07_dataset_heterogeneity'/f'{table}_patients.tsv',sep='\t',index=False)
 for group,g in pats.groupby(keys,dropna=False):
  info=dict(zip(keys,group if isinstance(group,tuple) else (group,)))
  for scope in ['ALL_KNOWN']+['DATASET:'+v for v in sorted(d.dataset.unique())]+['LEAVE_OUT:'+v for v in sorted(d.dataset.unique())]:
   pp=g if scope=='ALL_KNOWN' else g[g.dataset==scope.split(':',1)[1]] if scope.startswith('DATASET:') else g[g.dataset!=scope.split(':',1)[1]];pp=pp[pp.effect.notna()];v=pp.effect.to_numpy();n=len(v);brng=np.random.default_rng(np.random.SeedSequence([123,zlib.crc32('|'.join(pp.dataset.astype(str)+'::'+pp.patient_uid.astype(str)).encode())]));boot=np.median(brng.choice(v,(1000,n)),axis=1) if n>=3 else np.array([np.nan]);tautology=table=='rings' and info.get('center')=='HH' and info.get('variable')=='LLPS77_UCell';pval=float(wilcoxon(v).pvalue) if n>=5 and np.any(v!=0) and not tautology else np.nan;out.append(dict(table=table,**info,p_patient_wilcoxon=pval,scope=scope,n_patients=n,median=np.median(v) if n else np.nan,ci_low=np.quantile(boot,.025),ci_high=np.quantile(boot,.975),status='EXPLORATORY' if n>=3 else 'INSUFFICIENT_PATIENTS',omitted_dataset_has_known_patients=scope.split(':',1)[1] in set(g.dataset) if ':' in scope else True))
 # Unknown-patient sources remain section-level descriptions, not independent-patient CIs.
 sec=[]
 for group,g in d.groupby(keys,dropna=False):
  info=dict(zip(keys,group if isinstance(group,tuple) else (group,)))
  for excluded in ['NONE']+sorted(d.dataset.unique()):
   z=g if excluded=='NONE' else g[g.dataset!=excluded];sec.append(dict(**info,excluded_dataset=excluded,n_sections=z.section_id.nunique(),median_section_effect=z.effect.median(),median_of_dataset_medians=z.groupby('dataset').effect.median().median()))
 pd.DataFrame(sec).to_csv(O/'07_dataset_heterogeneity'/f'{table}_section_LODO_descriptive.tsv',sep='\t',index=False)
summary=pd.DataFrame(out)
summary['FDR_primary_scope_within_table']=np.nan
for tab,group in summary[summary.scope=='ALL_KNOWN'].groupby('table'):
 idx=group.index[group.p_patient_wilcoxon.notna()];idx=summary.loc[idx].sort_values('p_patient_wilcoxon').index;pv=summary.loc[idx,'p_patient_wilcoxon'].to_numpy();summary.loc[idx,'FDR_primary_scope_within_table']=np.minimum(1,np.minimum.accumulate((pv*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1])
summary.to_csv(O/'07_dataset_heterogeneity/patient_intervals.tsv',sep='\t',index=False)
q=pd.read_csv(O/'tables/qc.tsv',sep='\t');q.groupby(['dataset','tissue','score']).agg(n_sections=('section_id','size'),median_score=('median','median'),median_moran=('moran_I','median'),median_zero_fraction=('zero_fraction','median')).reset_index().to_csv(O/'07_dataset_heterogeneity/tissue_description.tsv',sep='\t',index=False)
# Broad comparison with prior A26 three-variable same-spot sensitivity, no causal interpretation.
p=Path('results/C5_spatial_exploratory_all63_v1/06_spatial_relationships/epithelial_matched_sensitivity.tsv');old=pd.read_csv(p,sep='\t');old.to_csv(O/'06_epithelial_and_TME_neighborhood/prior_A26_sensitivity_reference.tsv',sep='\t',index=False)
# organize tables without redundant copies
for folder,names in {'01_LLPS77_spatial_scores':['qc'],'02_Fivegene_spatial_distribution':['clustering','gene_contribution'],'03_Fivegene_controls':['balance','confounding'],'04_Fivegene_LLPS77_association':['association'],'05_spatial_colocalization':['colocalization','spatial_null'],'06_epithelial_and_TME_neighborhood':['rings','epithelial_distance']}.items():
 for name in names:
  link=O/folder/(name+'.tsv')
  if not link.exists():link.symlink_to('../tables/'+name+'.tsv')
print('Summaries done')

balance=pd.read_csv(O/'tables/balance.tsv',sep='\t')
wide=balance.pivot(index=['section_id','dataset','center','matching','n_high','n_pairs','match_rate'],columns='covariate',values='SMD_after').reset_index()
wide['balanced']=(wide.n_pairs>=10)&(wide.match_rate>=.5)&(wide[['Fibroblast','UMI','Epithelial']].abs()<=.1).all(axis=1)
wide.to_csv(O/'03_Fivegene_controls/section_balance_status.tsv',sep='\t',index=False)
wide[~wide.balanced].to_csv(O/'03_Fivegene_controls/insufficient_balance_sections.tsv',sep='\t',index=False)

summary[summary.scope=='ALL_KNOWN'].groupby('table').agg(n_contrasts=('scope','size'),n_tested=('p_patient_wilcoxon','count')).reset_index().to_csv(O/'00_input_audit/patient_BH_families.tsv',sep='\t',index=False)
