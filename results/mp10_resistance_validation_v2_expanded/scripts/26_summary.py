from env import *
from scipy.stats import false_discovery_control
import shutil
mf=read('01_signatures/signature_manifest.tsv');cov=read('00_plan_audit/coverage.tsv');parts=[]
def add(path,kind,estimate,low,high):
 p=B/path
 if not p.exists():return
 d=pd.read_csv(p,sep='\t');d['effect_type']=kind;d['effect']=d.get(estimate,np.nan);d['CI_low']=d.get(low,np.nan);d['CI_high']=d.get(high,np.nan)
 if 'module' not in d:d['module']='cell_models'
 if 'status' not in d:d['status']=np.where(d.P.notna(),'COMPLETE','BLOCKED_DESIGN')
 if 'n' not in d:d['n']=d.get('n_RNA',d.get('n_backgrounds',np.nan))
 if 'endpoint' not in d:d['endpoint']=d.get('state','association')
 parts.append(d)
add('04_measured_drug_response/paired_IC50_associations.tsv','correlation','r','CI_low','CI_high');add('03_cell_models/contrast_results.tsv','point_biserial_r','r','r_low','r_high');add('03_cell_models/GSE171537_background_summary.tsv','mean_delta_5background','mean_delta','CI_low','CI_high');add('04_measured_drug_response/GDSC_associations.tsv','correlation','r','CI_low','CI_high');add('05_clinical/model_results.tsv','OR_per1SD','OR','CI_low','CI_high');add('06_paired_states/paired_results.tsv','paired_mean_delta','mean_delta','CI_low','CI_high');add('02_epi/patient_summary.tsv','patient_Fisher_rho','rho','CI_low','CI_high');add('02_epi/pseudobulk_associations.tsv','pseudobulk_rho','rho','CI_low','CI_high')
e=pd.concat(parts,ignore_index=True);e['FDR']=np.nan
for fam,idx in e.groupby('family').groups.items():
 good=e.loc[idx,'P'].notna();ii=np.array(idx)[good]
 if len(ii):e.loc[ii,'FDR']=false_discovery_control(e.loc[ii,'P'].values,method='bh')
e['global_BH_all_estimable']=np.nan;ii=e.index[e.P.notna()];e.loc[ii,'global_BH_all_estimable']=false_discovery_control(e.loc[ii,'P'],method='bh');primary=e.module.isin(['cell_models'])|((e.module=='clinical')&e.main.fillna(False)&e.adjustment.eq('unadjusted'))|((e.module=='measured_drug')&e.main.fillna(False))|((e.module=='epi')&e.subset.eq('Tumor')&e.method.isin(['average_tie','pseudobulk_native'])&e.adjustment.eq('raw'))
e['is_primary_transparency_family']=primary;e['global_BH_primary']=np.nan;ii=e.index[primary&e.P.notna()];e.loc[ii,'global_BH_primary']=false_discovery_control(e.loc[ii,'P'],method='bh')
e=e.merge(mf[['signature','expected_resistance_direction','property','training','MP10_intersection']],on='signature',how='left',validate='many_to_one');signed=np.where(e.effect_type.eq('OR_per1SD'),np.log(e.effect),e.effect)*e.expected_resistance_direction;e['direction']=np.select([e.expected_resistance_direction.eq(0),signed>0,signed<0],['state_control_no_resistance_expectation','expected','opposite'],default='not_estimable');e.loc[e.module=='paired','direction']='post_state_no_causal_direction_claim';ok=e.status.eq('COMPLETE')&e.P.notna();e.loc[ok&(e.FDR>=.05),'status']='COMPLETE_NULL';e.loc[ok&(e.FDR<.05)&e.direction.eq('opposite'),'status']='COMPLETE_OPPOSITE'
cv=cov[['dataset','signature','coverage','n_present','n_original','version']].copy()
epi_cov=read('00_plan_audit/epi_coverage.tsv');epi_cov['version']=np.where(epi_cov.n_present==epi_cov.n_original,'full','partial_adapted')
pb_cov=read('00_plan_audit/epi_pseudobulk_coverage.tsv');pb_cov['version']=np.where(pb_cov.n_present==pb_cov.n_original,'full','partial_adapted')
cv=pd.concat([cv,epi_cov[cv.columns],pb_cov[cv.columns]],ignore_index=True).rename(columns={'dataset':'scoring_dataset'})
e['scoring_dataset']=e.dataset;e.loc[e.dataset.isin(['GDSC1','GDSC2']),'scoring_dataset']='GDSC_BLCA';e.loc[e.module.eq('epi'),'scoring_dataset']='Epi';pb=e.module.eq('epi')&e.method.eq('pseudobulk_native');e.loc[pb,'scoring_dataset']='EpiPB_'+e.loc[pb,'dataset']
e=e.merge(cv,on=['scoring_dataset','signature'],how='left',validate='many_to_one')
p=B/'01_signatures/implementation_by_dataset.tsv'
if p.exists():
 imp=pd.read_csv(p,sep='\t').rename(columns={'dataset':'scoring_dataset'});e=e.merge(imp[['scoring_dataset','signature','score_name','definition_status','native_exact']],on=['scoring_dataset','signature'],how='left',validate='many_to_one')
for _,missing_source in mf[mf.status!='READY'].iterrows():
 mask=e.signature.eq(missing_source.signature)&e.P.isna();e.loc[mask,'status']=missing_source.status;e.loc[mask,'reason']=missing_source.reason
e['independence_note']=e.training.fillna('');save(e,'07_summary/evidence_matrix.tsv.gz');save(e[e.signature=='MP10'],'07_summary/MP10_evidence.tsv');save(e[e.P.notna()&e.FDR.lt(.05)],'07_summary/all_FDR_lt005.tsv')
# Public scoring states and every unavailable branch are explicit, not missing from the deliverable.
qc=[]
for p in (B/'00_plan_audit/scores').glob('*_score_qc.tsv'):qc.append(pd.read_csv(p,sep='\t'))
scores=pd.concat(qc,ignore_index=True);scores['module']='scoring_QC';scores['endpoint']='score';scores['P']=np.nan
blocked=[]
for g in ['GSE169455','GSE212810','GSE52219','GSE244266','GSE87304','GSE69795','GSE247185','GSE5287','GSE124305','GSE235067','Taber2020_NAC','Taber2020_FL']:
 for _,m in mf.iterrows():
  if g=='GSE235067':blocked.append(dict(module='scoring_QC',dataset=g,signature=m.signature,endpoint='score',status='BLOCKED_COVERAGE',reason='41 endogenous genes; complete transcriptomic ranking background unavailable; no eligible complete native set',P=np.nan))
  if m.status!='READY':blocked.append(dict(module='signature_source',dataset=g,signature=m.signature,endpoint='all',status=m.status,reason=m.reason,P=np.nan))
# Score completion and unavailable same-direction sets remain explicit for every method.
qc_path=B/'02_epi/all_methods_score_QC.tsv'
full_qc=(B/'00_plan_audit/official_full_QC.json').exists() and qc_path.exists()
if full_qc:
 scores=pd.concat([scores,pd.read_csv(qc_path,sep='\t')],ignore_index=True)
ready=set(read('00_plan_audit/epi_coverage.tsv').query('status=="READY"').signature)
execution=json.loads((B/'logs/night_status.json').read_text()) if (B/'logs/night_status.json').exists() else {}
for ds in read('00_plan_audit/epi_patient_identity.tsv').dataset.unique():
 for _,m in mf.iterrows():
  for method in ['average_tie','official_AUCell','UCell']:
   if m.signature in ready:
    if full_qc:continue
    state='COMPLETE' if method=='average_tie' else ('RUNNING' if execution.get('state')=='running' else 'DEFERRED_RESOURCE')
    reason='full average-tie scoring complete' if method=='average_tie' else 'official full scoring execution state: '+str(execution.get('state','not_started'))
   else:
    state=m.status if m.status!='READY' else 'BLOCKED_DESIGN'
    reason=m.reason if m.status!='READY' else 'native signed/weighted/ratio definition is evaluated in pseudobulk; not a same-direction single-cell set'
   blocked.append(dict(module='scoring_QC',dataset=ds,signature=m.signature,endpoint='score',method=method,status=state,reason=reason,P=np.nan))
extra=[('GSE169455','strict_individual_cisplatin','patient-linked regimen/cycles absent; aggregate76GC/46ddMVAC/3other cannot identify patients'),('GSE247185','GC_specific','Gem annotated gentamicin in source; original regimen clarification needed'),('GSE87304','clinical_response','305 GSM/title to individual ypT/ypN/author_response/regimen table absent'),('GSE244266','clinical_response','patient/GSM/title arm to ypTN/response; slide_names.csv / slide_names_uncms.csv absent from public code'),('GSE5287','advanced_radiologic','30 patient-linked CR/PR/SD/PD or verified R/NR plus sample timing absent; survival not substituted'),('GSE124305','pre_post_pair','patient crosswalk to GSE87304, original116 paired identifiers absent'),('Taber2020_FL','strict_pre_FL_timing','individual expression sampling dates/tissue versus first-line treatment dates absent; author cohort association completed, contemporaneous pretreatment validation not asserted')]
for g,ep,reason in extra:
 for sn in mf.signature:blocked.append(dict(module='strict_branch',dataset=g,signature=sn,endpoint=ep,status='BLOCKED_LABELS',reason=reason,P=np.nan))
status=pd.concat([e,scores,pd.DataFrame(blocked)],ignore_index=True);save(status,'00_plan_audit/analysis_status.tsv.gz');save(status[status.status.str.startswith(('BLOCKED','DEFERRED'),na=False)],'00_plan_audit/blockers.tsv.gz')
plan=read('00_plan_audit/FDR_families.tsv');real=e.groupby('family').agg(n_result_rows=('P','size'),n_estimable=('P','count'),n_FDR_lt005=('FDR',lambda z:(z<.05).sum())).reset_index();save(plan.merge(real,on='family',how='outer'),'00_plan_audit/FDR_actual_vs_plan.tsv');save(status.groupby(['module','status'],dropna=False).size().reset_index(name='n_rows'),'07_summary/status_counts.tsv')
# Clinical assay / unique identity / qualifying population / evaluable endpoint counts directly from audit/results.
cl=read('05_clinical/patient_endpoint_audit.tsv');co=cl.groupby('dataset').agg(n_assays=('matrix_id','size'),n_patient_IDs=('patient_id','nunique')).reset_index();mm=e[(e.module=='clinical')&(e.signature=='MP10')&e.adjustment.eq('unadjusted')];save(co.merge(mm[['dataset','branch','endpoint','n_eligible','n','events','responders','n_missing','status']],on='dataset',how='left'),'05_clinical/cohort_counts.tsv')
# Code fixes versus protected v1; no source object changed.
bugs=[('HLA identity','v1 HLA-DRB3->HLA-DPB1 reverse alias confirmed','v2 HGNC exact identity; no HLA substitution; missing retained','epi_signature_feature_mapping.tsv'),('unmapped background collapse','v1 same unknown author name RF00019 merged distinct Ensembl','v2 preserves UNMAPPED|original_id; only same mapped gene merged','GSE171537_gene_mapping.tsv'),('Ensembl fallback','v1 configured but actual mapping omitted stable-ID fallback','v2 HGNC Ensembl/Entrez + version stripping + unique alias, ambiguity retained','all *_gene_mapping.tsv'),('adjusted_rho','v1 unconditional NA branch','v2 actual rank residual depth/cycle and nFeature sensitivity; pending is labelled separately','../02_epi/sample_effects.tsv'),('patient identity','v1 all patient placeholders unverified','same12 Inhouse specimens verified by user evidence; a/b combined; others remain descriptive','epi_patient_identity.tsv'),('report hardcodes','v1 hardcoded report n/status','v2 all tables/report generated from actual output files','../scripts/29_report.py'),('AUCell method name','v1 registered average-tie AUCell-like','preserved same definition with numerical error<1e-12; official AUCell/UCell separate versions','../02_epi/average_tie_execution.json')]
save(pd.DataFrame(bugs,columns=['issue','before','after','audit_file']),'00_plan_audit/bugfix_audit.tsv')
# Small reproducible source manifest (raw files not included in review archive).
manifest=[]
for p in sorted(S.iterdir()):
 if not p.is_file() or p.suffix=='.part' or 'INCOMPLETE' in p.name or 'range_part' in p.name:continue
 h=hashlib.sha256()
 with p.open('rb') as f:
  while z:=f.read(2**20):h.update(z)
 manifest.append(dict(file=str(p),size=p.stat().st_size,sha256=h.hexdigest(),source='see acquisition scripts + original GEO accession / article supplement filename',retrieved_local_time=datetime.datetime.fromtimestamp(p.stat().st_mtime).isoformat()))
save(pd.DataFrame(manifest),'00_plan_audit/source_file_manifest.tsv');print('SUMMARY',len(e),'effects',len(status),'states',len(manifest),'source files',flush=True)
