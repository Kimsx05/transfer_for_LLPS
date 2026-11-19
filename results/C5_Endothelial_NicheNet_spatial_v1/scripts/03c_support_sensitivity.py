from pathlib import Path
import numpy as np,pandas as pd
O=Path('results/C5_Endothelial_NicheNet_spatial_v1');r=pd.read_csv(O/'03_candidate_ligands/expression_patient.tsv.gz',sep='\t');lr=pd.read_csv(O/'03_candidate_ligands/broad_LR_all.tsv',sep='\t')[['ligand','receptor']].drop_duplicates();a=pd.read_csv(O/'04_NicheNet_ligand_activity/all_program_activities.tsv',sep='\t');out=[]
for sender in ['Fibro_C5','Fibro_C0','Fibro_C4']:
 for threshold in [.05,.1,.2]:
  for omit in ['NONE']+sorted(r.dataset.unique()):
   t=r if omit=='NONE' else r[r.dataset!=omit];s=t[t.group==sender];e=t[t.group=='Endothelial']
   def summ(g,n):
    return g.groupby('gene').agg(total=('n_cells','sum'),detected=('n_detected','sum'),support=('pct',lambda x:np.nan)).drop(columns='support').assign(support=g[(g.pct>=threshold)&(g.n_cells>=n)].groupby('gene').size()).fillna({'support':0})
   ss=summ(s,5);ee=summ(e,20);lg=ss[(ss.detected/ss.total>=threshold)&(ss.support>=2)].index;rg=ee[(ee.detected/ee.total>=threshold)&(ee.support>=2)].index;valid=set(lr[lr.ligand.isin(lg)&lr.receptor.isin(rg)].ligand)
   for row in a.itertuples():out.append(dict(sender=sender,threshold=threshold,omitted_dataset=omit,program=row.program,ligand=row.test_ligand,fixed_prior_rank=row.rank,expression_gate_retained=row.test_ligand in valid))
pd.DataFrame(out).to_csv(O/'08_patient_dataset_sensitivity/ligand_expression_gate_sensitivity.tsv',sep='\t',index=False)
random=pd.read_csv(O/'04_NicheNet_ligand_activity/matched_random_programs.tsv.gz',sep='\t');summary=random.groupby(['program','test_ligand']).aupr_corrected.agg(['median','min','max']).reset_index();summary=summary.merge(a,on=['program','test_ligand']);summary.to_csv(O/'04_NicheNet_ligand_activity/random_program_comparison.tsv',sep='\t',index=False)
# Audit sample-stage confounding without forcing non-identifiable clinical terms.
s=pd.read_csv(O/'00_input_audit/specimen_sampling.tsv',sep='\t');pd.crosstab(s.dataset,s.clinical_group).to_csv(O/'00_input_audit/dataset_stage_crosstab.tsv',sep='\t');p=pd.read_csv(O/'01_receiver_pseudobulk/patient_sampling.tsv',sep='\t');p[p.eligible].groupby('dataset').fraction.agg(['count','min','max','std']).to_csv(O/'00_input_audit/within_dataset_fraction_variation.tsv',sep='\t')
print('SUPPORT_SENSITIVITY_COMPLETE')
