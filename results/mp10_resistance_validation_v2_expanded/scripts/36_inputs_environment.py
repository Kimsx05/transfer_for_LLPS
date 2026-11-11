from env import *
import platform,importlib.metadata,shutil
C=P/'neighborhood_mp10_calibration_v5/01_score_fidelity/counts_cache';paths=[C/x for x in ['data.bin','indices.bin','indptr.bin','genes.txt','cells.txt']]+[P/'neighborhood_mp10_calibration_v5/01_score_fidelity/cell_level_fidelity.parquet',P/'mp_ucell_state_v1/00_input_audit/sample_patient_mapping.tsv'];rows=[]
for p in paths:
 h=hashlib.sha256()
 with p.open('rb') as f:
  while chunk:=f.read(2**20):h.update(chunk)
 rows.append(dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest(),role='protected existing RNA/Epi input',provenance='registered v1 frozen input paths and source workflow; final content hash recorded for replay; no write to input'))
save(pd.DataFrame(rows),'00_plan_audit/epi_input_manifest.tsv');vers={}
for package in ['numpy','pandas','scipy','pyarrow','openpyxl','pypdf']:
 try:vers[package]=importlib.metadata.version(package)
 except importlib.metadata.PackageNotFoundError:vers[package]='not_available'
(B/'00_plan_audit/python_environment.json').write_text(json.dumps(dict(executable=sys.executable,python=sys.version,platform=platform.platform(),packages=vers,BLAS_threads=1,only_new_runtime_packages='pypdf7.0.1 openpyxl3.1.5 et_xmlfile2.0.0 installed into v2/runtime; no global environment upgrade'),indent=2));shutil.copy2(B/'02_epi/pilot/sessionInfo.txt',B/'00_plan_audit/R_epi_sessionInfo.txt')
models=read('03_cell_models/model_identity.tsv');context=models[['dataset','study','background','state','drug','exposure','replicate_type']].drop_duplicates();context['Base47_subtype']='unknown: no authenticated per-line original classifier result recovered';context['driver_genotype']='not inferred; source blanks are not WT';context['culture_source']='original GEO growth/treatment protocol in prepared sample metadata';save(context,'03_cell_models/background_context.tsv')
scores=read('03_cell_models/sample_scores.tsv.gz');fields=['exposure','replicate_type','assay_treatment'];scores=scores.drop(columns=fields,errors='ignore').merge(models[['dataset','GSM']+fields],on=['dataset','GSM'],how='left',validate='many_to_one');save(scores,'03_cell_models/sample_scores.tsv.gz')
cl=read('05_clinical/patient_endpoint_audit.tsv');cl['sample_timing_evidence']=np.where(cl.dataset.str.startswith('Taber2020'),'author response-cohort tissue; individual sampling dates not independently recovered',np.where(cl.dataset=='GSE5287','individual pretreatment timing not established; response model blocked','see original GEO/source inclusion audit'));save(cl,'05_clinical/patient_endpoint_audit.tsv');print('INPUT_ENVIRONMENT_COMPLETE',sum(x['bytes'] for x in rows))
