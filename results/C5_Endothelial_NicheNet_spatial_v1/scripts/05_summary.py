from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import wilcoxon
O=Path('results/C5_Endothelial_NicheNet_spatial_v1');D=O/'07_spatial_expression_validation'
def save(x,p):pd.DataFrame(x).to_csv(O/p,sep='\t',index=False)
def bh(p):
 p=np.array(p,float);out=np.full(len(p),np.nan);ok=np.isfinite(p);v=p[ok];ix=np.argsort(v);q=np.minimum.accumulate((v[ix]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1];a=np.empty(len(v));a[ix]=np.minimum(q,1);out[ok]=a;return out
fx=pd.concat([pd.read_csv(p,sep='\t') for p in sorted(D.glob('*/effects.tsv'))],ignore_index=True);fx['translation_FDR']=bh(fx.translation_p)
qkeys=['section_id','ligand','top_fraction','ring'];qc=fx[fx.role!='curated_program'][qkeys+['ligand_detection_sender','n_senders','receptor_joint_detection_receiver','n_receiver_spots']].copy();qc['technical_evaluable']=(qc.ligand_detection_sender*qc.n_senders>=9.999)&(qc.receptor_joint_detection_receiver*qc.n_receiver_spots>=9.999);fx=fx.merge(qc[qkeys+['technical_evaluable']],on=qkeys,how='left',validate='many_to_one');fx['technical_evaluable']=fx.technical_evaluable.fillna(False);save(fx,'tables/spatial_effects.tsv');balance=pd.concat([pd.read_csv(p,sep='\t') for p in D.glob('*/balance.tsv')],ignore_index=True);save(balance,'tables/matching_balance.tsv');coverage=pd.concat([pd.read_csv(p,sep='\t') for p in D.glob('*/coverage.tsv')],ignore_index=True);save(coverage,'tables/spatial_gene_coverage.tsv')
keys=['ligand','receptor','role','top_fraction','ring'];long=fx.melt(id_vars=keys+['section_id','dataset','patient_uid','balanced','technical_evaluable','n_pairs'],value_vars=['raw_rho','adjusted_rho','matched_product_effect'],var_name='metric',value_name='effect').dropna(subset=['effect']);save(long.groupby(keys+['dataset','metric'],dropna=False).agg(n_sections=('effect','size'),median_effect=('effect','median'),positive_fraction=('effect',lambda x:(x>0).mean())).reset_index(),'08_patient_dataset_sensitivity/dataset_summary.tsv')
save(long,'tables/all_section_effects_long.tsv');save(long[long.patient_uid.notna()].groupby(keys+['dataset','patient_uid','metric'],dropna=False).effect.mean().reset_index(),'08_patient_dataset_sensitivity/all_known_patient_descriptive_effects.tsv');long=long[long.technical_evaluable & ((long.metric!='matched_product_effect')|(long.n_pairs>=5))]
known=long[long.patient_uid.notna() & ~long.patient_uid.astype(str).isin(['UNRESOLVED','UNKNOWN'])];pp=known.groupby(keys+['dataset','patient_uid','metric'],dropna=False).effect.mean().reset_index();save(pp,'08_patient_dataset_sensitivity/patient_effects.tsv');out=[];rng=np.random.default_rng(123)
for key,g in pp.groupby(keys+['metric'],dropna=False):
 info=dict(zip(keys+['metric'],key));datasets=sorted(long.dataset.unique())
 for omit in ['NONE']+datasets:
  t=g if omit=='NONE' else g[g.dataset!=omit];v=t.effect.to_numpy();n=len(v)
  if not n:continue
  boots=[]
  if n>=3:
   boots=np.median(np.concatenate([rng.choice(z.effect.to_numpy(),(2000,len(z)),True) for _,z in t.groupby('dataset')],axis=1),axis=1).tolist()
  ci=np.quantile(boots,[.025,.975]) if boots else [np.nan,np.nan];pv=wilcoxon(v).pvalue if n>=5 and np.any(v!=0) else np.nan
  out.append(dict(**info,omitted_dataset=omit,n_patients=n,n_datasets=t.dataset.nunique(),median_effect=np.median(v),ci_low=ci[0],ci_high=ci[1],p=pv))
su=pd.DataFrame(out);su['BH_FDR']=np.nan
# Primary patient test family: all selected candidate/control axes Ring1 top10 adjusted and matched outcomes; sensitivities separate.
for (metric,omit),ix in su.groupby(['metric','omitted_dataset']).groups.items():su.loc[ix,'BH_FDR']=bh(su.loc[ix,'p'])
save(su,'08_patient_dataset_sensitivity/patient_bootstrap_LODO.tsv')
# Balanced subset is complete and separately reported, not selected by significance.
bb=known[(known.balanced==True)&(known.metric=='matched_product_effect')];save(bb.groupby(keys+['dataset','patient_uid','metric'],dropna=False).effect.mean().reset_index(),'08_patient_dataset_sensitivity/balanced_patient_effects.tsv')
print('SUMMARY_COMPLETE',len(fx),fx.section_id.nunique(),known.patient_uid.nunique())
