from env import *
import shutil
mf=read('01_signatures/signature_manifest.tsv');members=read('01_signatures/members.tsv');cov=read('00_plan_audit/coverage.tsv');dm=read('00_plan_audit/dataset_manifest.tsv');qcs=pd.concat([pd.read_csv(p,sep='\t') for p in (B/'00_plan_audit/scores').glob('*_score_qc.tsv')],ignore_index=True);imp=[]
for _,r in cov.iterrows():
 m=mf.set_index('signature').loc[r.signature];q=qcs[(qcs.dataset==r.dataset)&(qcs.signature==r.signature)];used=q.n_members_used.iloc[0] if len(q) and 'n_members_used' in q and pd.notna(q.n_members_used.iloc[0]) else r.n_present;full=used==r.n_original;scale=dm.set_index('dataset').loc[r.dataset,'scale'];label='full' if full else f'partial_adapted_{int(used)}of{int(r.n_original)}';native=m.property=='native' and full
 if m.method=='linear_ratio':label='native_linear' if scale.startswith(('counts','linear')) else 'adapted_inverse_log2_intensity'
 status=q.status.iloc[0] if len(q) else r.status
 if status!='COMPLETE':label='not_scored_'+str(status);native=False;used=np.nan
 name=r.signature+'__'+label;imp.append(dict(dataset=r.dataset,signature=r.signature,score_name=name,definition_status=label,n_members_actually_used=used,native_exact=(status=='COMPLETE') and (native or (m.method=='linear_ratio' and label=='native_linear')),source_property=m.property,coverage=r.coverage,status=q.status.iloc[0] if len(q) else r.status))
implementation=pd.DataFrame(imp);save(implementation,'01_signatures/implementation_by_dataset.tsv')
for p in (B/'00_plan_audit/scores').glob('*_scores.tsv'):
 s=pd.read_csv(p,sep='\t');extra=implementation[['dataset','signature','score_name','definition_status']];s=s.drop(columns=['score_name','definition_status'],errors='ignore').merge(extra,on=['dataset','signature'],how='left',validate='many_to_one');save(s,p.relative_to(B))
# Definitions with identical members remain labelled, not counted as separate biological evidence.
redundancy=[]
for _,m in mf.iterrows():
 r=members[members.signature==m.signature];key=';'.join(sorted(r.gene.fillna('?')+':'+r.weight.astype(str)))+'|'+str(m.method);redundancy.append(dict(signature=m.signature,definition_sha256=hashlib.sha256(key.encode()).hexdigest(),status=m.status))
save(pd.DataFrame(redundancy),'01_signatures/definition_redundancy.tsv')
# Duplicate parent expressions in the two5637 derivation process series are explicitly measured.
a=read('00_plan_audit/prepared/GSE210921_expression.tsv.gz').set_index('gene');b=read('00_plan_audit/prepared/GSE210953_expression.tsv.gz').set_index('gene');mod=read('03_cell_models/model_identity.tsv');ai=mod[(mod.dataset=='GSE210921')&(mod.state=='parental')].matrix_id;bi=mod[(mod.dataset=='GSE210953')&(mod.state=='parental')].matrix_id;common=a.index.intersection(b.index);dups=[]
for x in ai:
 for y in bi:
  delta=np.abs(a.loc[common,x].values-b.loc[common,y].values);dups.append(dict(series1='GSE210921',matrix1=x,series2='GSE210953',matrix2=y,n_common_genes=len(common),fraction_equal=float(np.mean(delta==0)),max_absolute_difference=float(np.max(delta)),interpretation='identical source expression' if np.all(delta==0) else 'near-identical same-name P0; treat dependent; physical RNA identity not proven' if np.mean(delta==0)>.999 and x==y else 'different expression profiles'))
save(pd.DataFrame(dups),'00_plan_audit/5637_parent_expression_identity.tsv')
# Operational full20 sensitivity is a cohort restriction, not a newly optimized score.
e=read('07_summary/evidence_matrix.tsv.gz');full=set(cov[(cov.signature=='MP10')&(cov.n_present==20)].dataset);save(e[(e.signature=='MP10')&e.dataset.isin(full)],'07_summary/MP10_full20_coverage_sensitivity.tsv')
cl=read('05_clinical/patient_endpoint_audit.tsv');raw=read('00_plan_audit/GSE69795_samples.tsv');base=raw[raw['sample procured']=='TUR'].copy();base=base.merge(cl[cl.dataset=='GSE69795'][['GSM','strict_inclusion']],on='GSM');base['exclusion_reason']=np.select([base['source_name_ch1']=='UTUC',base['source_name_ch1']=='URETHRA',base.cstage.astype(str)=='1',base.cystectomy=='No'],['UTUC','urethral primary (not bladder)','bladder cT1','no cystectomy'],default='included');save(base[['GSM','title','source_name_ch1','cstage','cystectomy','pstage','strict_inclusion','exclusion_reason']],'00_plan_audit/GSE69795_inclusion_reconciliation.tsv')
# Protected version hashes are checked, not merely asserted.
snapshot=json.loads((B/'00_plan_audit/v1_snapshot.json').read_text());checks=[]
for r in snapshot['files']:
 p=Path(r['path']);actual=hashlib.sha256(p.read_bytes()).hexdigest();expected=r.get('sha256',r.get('hash'));checks.append(dict(path=str(p),expected=expected,actual=actual,unchanged=actual==expected))
save(pd.DataFrame(checks),'00_plan_audit/v1_protection_check.tsv');assert all(r['unchanged'] for r in checks),'protected v1 changed'
# Risk-focused integrity checks, using actual tables only.
flip=read('05_clinical/MP10_encoding_flip_check.tsv');assert flip.passed.all();assert np.allclose(flip.OR*flip.flipped_OR,1,atol=1e-4);assert len(read('02_epi/cell_cycle_scores.tsv.gz'))==123699
cy=read('02_epi/cell_cycle_scores.tsv.gz');assert not cy.cell_id.duplicated().any() and np.isfinite(cy[['S.Score','G2M.Score']]).all().all()
usage=read('05_clinical/model_patient_usage.tsv.gz');assert not usage.duplicated(['dataset','branch','signature','endpoint','adjustment','patient_id']).any();assert usage.nonresponse.isin([0,1]).all()
# Predeclared pathology risk cases: no simulated scientific result.
from importlib.util import spec_from_file_location,module_from_spec
import ast
src=(B/'scripts/12_design.py').read_text();tree=ast.parse(src);node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='endpoint');ns={'np':np,'re':re};exec(compile(ast.Module(body=[node],type_ignores=[]),'endpoint','exec'),ns);fn=ns['endpoint'];cases=[('0','0',True,(1,1)),('A','0',True,(1,0)),('1','X',True,(np.nan,0)),('2','X',True,(0,0)),('X','1',True,(0,0)),('0','X',True,(np.nan,np.nan)),('X','0',True,(np.nan,np.nan)),('0','0',False,(np.nan,np.nan))];test=[]
for t,n,s,expected in cases:
 got=fn(t,n,s);passed=all((pd.isna(a) and pd.isna(b)) or a==b for a,b in zip(got,expected));assert passed;test.append(dict(T=t,N=n,surgery=s,good=got[0],pCR=got[1],passed=passed))
save(pd.DataFrame(test),'00_plan_audit/endpoint_encoding_checks.tsv')
countids=read('03_cell_models/model_identity.tsv');deltas=read('03_cell_models/contrast_results.tsv').query('signature=="MP10"');states=deltas[['dataset','background','state']].copy();states['process']=states.state.str.replace(r'_P\d+$','',regex=True);resources=[]
for p in (B/'logs').glob('*.log'):
 text=p.read_text(errors='replace');mx=re.findall(r'Maximum resident set size \(kbytes\): (\d+)',text)
 if mx:resources.append(dict(log=p.name,peak_RSS_GB=max(map(int,mx))/1024**2))
save(pd.DataFrame(resources),'00_plan_audit/resource_measurements.tsv');notes=dict(timestamp=datetime.datetime.now().isoformat(),all_v1_hashes_unchanged=True,encoding_flip_checks=len(flip),encoding_flip_passed=int(flip.passed.sum()),cycle_cells=len(cy),max_single_process_RSS_GB=max([r['peak_RSS_GB'] for r in resources]),official_full_scoring_complete=(B/'02_epi/cell_official_scores.tsv.gz').exists() and (B/'00_plan_audit/official_full_QC.json').exists(),cell_model_assays=len(countids),baseline_model_assays=int(countids.baseline_include.sum()),paired_backgrounds=deltas.background.nunique(),studies=countids.study.nunique(),contrasts=len(deltas),study_background_drug_process_units=len(states[['dataset','background','process']].drop_duplicates()),signatures_registered=len(mf),signatures_members_ready=int(mf.status.eq('READY').sum()),clinical_unique_IDs_with_MP10_evaluable_primary=len(usage[(usage.signature=='MP10')&(usage.adjustment=='unadjusted')].merge(read('05_clinical/model_results.tsv').query('signature=="MP10" and main==True and adjustment=="unadjusted"')[['dataset','branch','endpoint']],on=['dataset','branch','endpoint'])[['dataset','patient_id']].drop_duplicates()),disk_bytes=sum(p.stat().st_size for p in B.rglob('*') if p.is_file()))
(B/'00_plan_audit/validation.json').write_text(json.dumps(notes,indent=2));shutil.copy2(B/'00_plan_audit/frozen_config.json',B/'00_plan_audit/analysis_plan.json');print(json.dumps(notes,indent=2),flush=True)
