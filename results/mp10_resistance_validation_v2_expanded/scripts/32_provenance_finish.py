from env import *
import subprocess
# Numerical before/after effects of the corrected DN mapping, with original input preserved.
old=pd.read_csv(V/'01_epi/cell_scores.tsv.gz',sep='\t',usecols=['cell_id','MP10','DN']);new=pd.read_csv(B/'02_epi/cell_average_tie_scores.tsv.gz',sep='\t',usecols=['cell_id','sample','MP10','TSUNODA_DN']);a=old.merge(new,on='cell_id',suffixes=('_v1','_v2'),validate='one_to_one');a['DN_change']=a.TSUNODA_DN-a.DN;rows=[]
for sample,d in a.groupby('sample'):
 rows.append(dict(sample=sample,n_cells=len(d),MP10_max_abs_change=np.max(np.abs(d.MP10_v1-d.MP10_v2)),DN_mean_before=d.DN.mean(),DN_mean_after=d.TSUNODA_DN.mean(),DN_mean_abs_change=d.DN_change.abs().mean(),rho_before=stats.spearmanr(d.MP10_v1,d.DN).statistic if len(d)>=50 else np.nan,rho_after=stats.spearmanr(d.MP10_v2,d.TSUNODA_DN).statistic if len(d)>=50 else np.nan,reason='HLA-DRB3 true absence no HLA-DPB1 substitution; normalized gene coverage denominator changes'))
save(pd.DataFrame(rows),'00_plan_audit/v1_v2_epi_DN_before_after.tsv')
# Disclose implementation corrections and enumeration gaps without rewriting the original frozen plan.
amend=[dict(time=datetime.datetime.now().isoformat(),type='implementation_enumeration',change='Expand already specified Normal/CNV/recorded-Type and leave-one branches in result/status tables; original frozen hypothesis_table missed these sensitivity expansions. Original10985 rows retained; FDR_actual_vs_plan explicitly shows additional pre-specified sensitivity families.',reason='Task/config mandated sensitivity analyses, not selected by associations'),dict(time=datetime.datetime.now().isoformat(),type='implementation_alignment',change='GSE244266 expression-median merged182 IDs also deduplicated in clinical metadata before model join',reason='Frozen duplicate rule; no endpoint-guided choice'),dict(time=datetime.datetime.now().isoformat(),type='source_recovery',change='Public Taber96 expression/clinical recovered before freeze; actual IC50 table subsequently recovered and evaluated using frozen formula and design threshold',reason='Public source availability; no outcome-based model choice'),dict(time=datetime.datetime.now().isoformat(),type='random_control_precision',change='100 matched controls use pooled pseudobulk mean/detection quintiles, not single-cell detection. Labelled descriptive; no drug specificity claim.',reason='Matching basis clarified to match actual input aggregation')]
save(pd.DataFrame(amend),'00_plan_audit/amendments/implementation_history.tsv')
# Remaining gene/source blocks are different from a failure to execute available data.
manifest=read('00_plan_audit/dataset_manifest.tsv');clinical=read('05_clinical/patient_endpoint_audit.tsv');models=read('03_cell_models/model_identity.tsv');identity=[]
for _,r in manifest.iterrows():
 g=r.dataset;d=clinical[clinical.dataset==g];m=models[models.dataset==g];identity.append(dict(dataset=g,species='Homo sapiens',expression_scale=r.scale,n_assays=r.n_assays,n_analysis_samples=r.n_analysis_samples,independent_unit='patient within source cohort' if len(d) else 'cell background/study; RNA nested replicate' if len(m) else 'SIDM cell model',disease='urothelial source cohort; per-patient exclusions audited' if g=='GSE69795' else 'BLCA',timepoint=';'.join(d.pretreatment.astype(str).unique()) if len(d) else ';'.join(m.exposure.unique()) if len(m) else 'basal expression vs separate screening',treatment=';'.join(d.regimen.astype(str).unique()) if len(d) else ';'.join(m.drug.unique()) if len(m) else 'measured cisplatin',note=r.note))
save(pd.DataFrame(identity),'00_plan_audit/dataset_identity_manifest.tsv')
# Audit sample score finiteness with reasons; no NA is turned into zero.
rows=[]
for p in (B/'00_plan_audit/scores').glob('*_scores.tsv'):
 d=pd.read_csv(p,sep='\t')
 for (ds,sn),x in d.groupby(['dataset','signature']):rows.append(dict(dataset=ds,signature=sn,n=len(x),n_nonfinite=int((~np.isfinite(x.score)).sum()),n_nonfinite_z=int((~np.isfinite(x.score_z)).sum()),n_unique=x.score.nunique(),reason='constant score or denominator0: inspect score_qc' if (~np.isfinite(x.score)).any() or (~np.isfinite(x.score_z)).any() else 'OK'))
save(pd.DataFrame(rows),'00_plan_audit/nonfinite_score_audit.tsv')
print('PROVENANCE_COMPLETE',flush=True)
