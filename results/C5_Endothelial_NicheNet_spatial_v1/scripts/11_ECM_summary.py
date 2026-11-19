from pathlib import Path
import numpy as np,pandas as pd
O=Path('results/C5_Endothelial_NicheNet_spatial_v1');e=pd.read_csv(O/'08_patient_dataset_sensitivity/ECM_depleted_program_effects.tsv',sep='\t');s=pd.read_csv(O/'tables/spatial_effects.tsv',sep='\t');q=s[(s.top_fraction==.1)&(s.ring==1)&(s.role!='curated_program')][['section_id','ligand','technical_evaluable']];e=e.merge(q,on=['section_id','ligand'],validate='many_to_one');e.to_csv(O/'08_patient_dataset_sensitivity/ECM_depleted_program_effects.tsv',sep='\t',index=False);p=e[e.technical_evaluable & e.patient_uid.notna()].groupby(['ligand','program','dataset','patient_uid'],as_index=False).adjusted_rho.mean().dropna();out=[]
for (l,g),t in p.groupby(['ligand','program']):
 for omit in ['NONE']+sorted(p.dataset.unique()):
  z=t if omit=='NONE' else t[t.dataset!=omit];v=z.adjusted_rho.to_numpy();n=len(v)
  if not n:continue
  rng=np.random.default_rng(123);bt=np.median(np.concatenate([rng.choice(k.adjusted_rho,(2000,len(k)),True) for _,k in z.groupby('dataset')],axis=1),axis=1) if n>=3 else np.array([np.nan]);out.append(dict(ligand=l,program=g,omitted_dataset=omit,n_patients=n,median_effect=np.median(v),ci_low=np.quantile(bt,.025),ci_high=np.quantile(bt,.975)))
r=pd.DataFrame(out);r.to_csv(O/'08_patient_dataset_sensitivity/ECM_depleted_patient_CI_LODO.tsv',sep='\t',index=False);print(r[(r.omitted_dataset=='NONE')&(r.ligand.isin(['TGFB1','FN1','HSPG2']))].to_string(index=False))
