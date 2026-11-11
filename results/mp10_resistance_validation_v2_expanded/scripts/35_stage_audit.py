from env import *
f=pd.read_csv(B/'02_epi/cell_average_tie_scores.tsv.gz',sep='\t',usecols=['sample','dataset','T_stage','N_stage','M_stage'],low_memory=False);rows=[]
for sample,d in f.groupby('sample'):
 row=dict(sample=sample,dataset=d.dataset.iloc[0],n_cells=len(d),source='protected original cell_level_fidelity.parquet metadata',stage_type='recorded T/N/M; clinical versus pathological origin not independently resolved',treatment_timing='not present in available cell metadata')
 for c in ['T_stage','N_stage','M_stage']:
  v=d[c].dropna().unique();assert len(v)<=1,(sample,c,v);row[c]=v[0] if len(v)==1 else np.nan
 rows.append(row)
a=pd.DataFrame(rows);save(a,'00_plan_audit/epi_stage_timing_audit.tsv');sm=read('02_epi/pseudobulk_samples.tsv').drop(columns=['T_stage','N_stage','M_stage'],errors='ignore');sm=sm.merge(a[['sample','T_stage','N_stage','M_stage']],on='sample',how='left',validate='one_to_one');save(sm,'02_epi/pseudobulk_samples.tsv')
# No quantitative pathological tumor content column was found in public GEO case metadata for main cohorts.
save(pd.DataFrame([dict(dataset=g,analysis='pathological_tumor_content_adjustment',status='BLOCKED_DATA',needed='patient-linked quantitative pathology tumor-cell percentage; cannot substitute expression purity or assume microdissection100%') for g in ['GSE169455','GSE212810','GSE52219','GSE244266','GSE87304','GSE69795','GSE247185','GSE5287','Taber2020']]),'00_plan_audit/purity_adjustment_audit.tsv')
print('STAGE_AUDIT_COMPLETE',len(a))
