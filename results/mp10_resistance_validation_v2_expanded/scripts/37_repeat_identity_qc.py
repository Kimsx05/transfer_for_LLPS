from env import *
# Audit title-based repeat grouping without claiming a verified patient crosswalk.
s=read('00_plan_audit/prepared/GSE244266_samples.tsv')
s['base_title']=s.title.str.replace(r'_F2$','',regex=True)
x=pd.read_csv(B/'00_plan_audit/prepared/GSE244266_expression.tsv.gz',sep='\t',index_col=0)
rows=[]
for base,g in s.groupby('base_title',sort=True):
 if len(g)<2:continue
 ids=g.matrix_id.tolist();a=x[ids[0]].to_numpy();b=x[ids[1]].to_numpy();ok=np.isfinite(a)&np.isfinite(b)
 rows.append(dict(base_title=base,n_assays=len(g),GSMs=';'.join(g.GSM),titles=';'.join(g.title),full_background_genes=int(ok.sum()),expression_Pearson=float(np.corrcoef(a[ok],b[ok])[0,1]),regimen_consistent=g['treatment arm'].nunique(dropna=False)==1,cT_consistent=g['clinical_stage_strat_factor'].nunique(dropna=False)==1,identity_status='provisional title-based repeat group; patient/physical RNA identity not independently verified',needed='author patient-GSM crosswalk and F2 repeat definition'))
assert len(s)==189 and s.base_title.nunique()==182 and len(rows)==7
save(pd.DataFrame(rows),'00_plan_audit/GSE244266_repeat_QC.tsv')
c=read('05_clinical/patient_endpoint_audit.tsv')
if 'patient_identity_status' not in c:c['patient_identity_status']='see source-specific identity audit; no cross-study uniqueness assertion'
c.loc[c.dataset.eq('GSE244266'),'patient_identity_status']='provisional base-title group, not independently verified patient; no response inference'
c.loc[c.dataset.eq('GSE235067'),'patient_identity_status']='panel assay identifier; clinical patient crosswalk unresolved'
save(c,'05_clinical/patient_endpoint_audit.tsv')
print(pd.DataFrame(rows).to_string(index=False));print('REPEAT_IDENTITY_QC_COMPLETE')
