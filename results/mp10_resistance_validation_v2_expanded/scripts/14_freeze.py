from env import *
assert not (B/'00_plan_audit/frozen_config.sha256').exists(), 'already frozen; use logged amendment'
mf=read('01_signatures/signature_manifest.tsv');members=read('01_signatures/members.tsv');cl=read('05_clinical/patient_endpoint_audit.tsv');model=read('03_cell_models/model_identity.tsv');prep=B/'00_plan_audit/prepared'
# Separate treatment populations before any native cross-sample z scoring.
t=cl[cl.dataset=='Taber2020'].copy();ex=read('00_plan_audit/prepared/Taber2020_expression.tsv.gz').set_index('gene')
for suffix,mask in [('NAC',t.therapy_intent.str.contains('NAC')),('FL',t.therapy_intent=='advanced_first_line')]:
 ids=t.loc[mask,'matrix_id'].tolist();name='Taber2020_'+suffix;save(ex[ids].reset_index(),f'00_plan_audit/prepared/{name}_expression.tsv.gz');save(pd.DataFrame({'matrix_id':ids}),f'00_plan_audit/prepared/{name}_samples.tsv');j=json.loads((prep/'Taber2020_input.json').read_text());j.update(dataset=name,n_assays=len(ids));(prep/f'{name}_input.json').write_text(json.dumps(j,indent=2));cl.loc[(cl.dataset=='Taber2020')&cl.matrix_id.isin(ids),'dataset']=name
save(cl,'05_clinical/patient_endpoint_audit.tsv')
# Clinical branches; all original labels remain separate from reconstructed pathology.
plans=[]
def cp(g,branch,mask,ends,main=True):
 d=cl[cl.dataset==g];ids=d.loc[mask(d),'matrix_id'].tolist()
 for ep in ends:
  plans.append(dict(dataset=g,branch=branch,endpoint=ep,matrix_ids=ids,main=main,family=('clinical_NAC_good' if ep=='good_response_ypT_le1_N0' else 'clinical_pCR' if ep=='pCR_ypT0N0' else 'clinical_advanced_response' if g=='Taber2020_FL' else 'clinical_author_response')+('' if main else '_sensitivity')))
ends=['good_response_ypT_le1_N0','pCR_ypT0N0','author_response']
for g in cl.dataset.unique():
 if g=='GSE169455':
  cp(g,'author_NAC_fallback',lambda d:d.fallback_inclusion==True,ends);cp(g,'induction_separate',lambda d:d.therapy_intent.str.lower().eq('induction'),ends,False)
 elif g=='GSE212810':
  cp(g,'strict_cisplatin',lambda d:d.strict_inclusion==True,ends);cp(g,'all_platinum',lambda d:d.fallback_inclusion==True,ends,False);cp(g,'cisplatin_cycles_ge3',lambda d:(d.strict_inclusion==True)&(pd.to_numeric(d.cycles,errors='coerce')>=3),ends,False)
 elif g=='GSE69795':
  cp(g,'BLCA_cT2_surgery_ddMVAC_bev',lambda d:d.strict_inclusion==True,ends);cp(g,'author_pre_cohort',lambda d:d.pretreatment==True,['author_response'],False)
 elif g=='GSE247185':
  cp(g,'author_cisplatin',lambda d:d.fallback_inclusion==True,ends);cp(g,'exclude_Patient12_pembro',lambda d:d.patient_id!='Patient12',['author_response'],False)
 elif g=='Taber2020_NAC':cp(g,'NAC_GC_author_sampling',lambda d:d.strict_inclusion==True,ends)
 elif g=='Taber2020_FL':
  cp(g,'FL_unambiguous_cisplatin',lambda d:d.regimen.isin(['GC','MVAC','Cisplatin + Etoposid']),['author_response']);cp(g,'FL_author_all_regimens',lambda d:d.fallback_inclusion==True,['author_response'],False)
 else:cp(g,'primary' if g!='GSE124305' else 'post_NAC_not_baseline',lambda d:d.pretreatment==True,ends)
(B/'00_plan_audit/clinical_plan.json').write_text(json.dumps(plans,indent=2))
# Reference sample population for gene z and score z, frozen independently of evaluability of endpoint.
references={}
for g in cl.dataset.unique():
 q=[p for p in plans if p['dataset']==g and p['main']];references[g]=q[0]['matrix_ids'] if q and q[0]['matrix_ids'] else cl.loc[cl.dataset==g,'matrix_id'].tolist()
# 244266 repeated titles: expression median merge before any scoring; metadata remains189 original assays.
for g in ['GSE244266']:
 d=cl[cl.dataset==g];x=read(f'00_plan_audit/prepared/{g}_expression.tsv.gz').set_index('gene');mp=dict(zip(d.matrix_id,d.patient_id));x=x.rename(columns=mp).T.groupby(level=0,sort=False).median().T;save(x.reset_index(),f'00_plan_audit/prepared/{g}_analysis_expression.tsv.gz');references[g]=x.columns.tolist();save(d[['GSM','matrix_id','patient_id']],f'00_plan_audit/{g}_repeat_merge.tsv')
# Coverage is assessed on actual complete expression rows; no missing entry is imputed as zero.
coverage=[];datasets=[];inputrows=[]
for p in sorted(prep.glob('*_input.json')):
 j=json.loads(p.read_text());g=j['dataset']
 if g=='Taber2020':continue
 ep=prep/f'{g}_analysis_expression.tsv.gz';ep=ep if ep.exists() else prep/f'{g}_expression.tsv.gz';x=pd.read_csv(ep,sep='\t').set_index('gene');good=np.isfinite(x.values).all(axis=1);available=set(x.index[good]);refs=[a for a in references.get(g,x.columns.tolist()) if a in x.columns];references[g]=refs
 j.update(expression_path=str(ep),n_complete_background=int(good.sum()),n_nonfinite_background_removed=int((~good).sum()),n_analysis_samples=len(x.columns),reference_n=len(refs));datasets.append(j)
 for _,m in mf.iterrows():
  r=members[members.signature==m.signature];present=r.gene.isin(available)&r.gene.notna();frac=present.sum()/len(r) if len(r) else 0;dir_ok=all(present[r.weight>0].mean()>=.8 if (r.weight>0).any() else True for _ in [0]) and (present[r.weight<0].mean()>=.8 if (r.weight<0).any() else True)
  status='READY';reason='';exact=bool(len(r) and present.all())
  if not len(r):status='BLOCKED_DATA';reason=str(m.reason)
  elif frac<.8 or not dir_ok or (m.method=='linear_ratio' and not exact):status='BLOCKED_COVERAGE';reason='requires >=80% overall and each signed part; ratio requires all6'
  elif j['scale']=='panel_counts':status='BLOCKED_COVERAGE';reason='41-gene panel lacks complete transcriptome ranking background; no eligible native set' if m.method=='ssGSEA' else reason
  elif m.method=='linear_ratio' and j['scale'] in ['SCAN','author_normalized_unknown_linear_basis','author_log2CPM']:status='BLOCKED_DESIGN';reason='linear scale/offset not sufficiently documented; no ratio of transformed values'
  if m.signature=='CisSig_gene_set_adapted':status='BLOCKED_DESIGN';reason='single-cell enrichment version only; use native CisSig19 for bulk'
  coverage.append(dict(dataset=g,signature=m.signature,n_original=len(r),n_mapped=int(r.gene.notna().sum()),n_present=int(present.sum()),coverage=frac,exact=exact,version='full' if exact else 'partial_adapted',missing=';'.join(r.loc[~present,'original_member'].astype(str)),status=status,reason=reason,n_background=len(available)))
 save(pd.DataFrame({'gene':x.index,'complete_across_population':good,'positive_rate':(x>0).mean(axis=1).values,'positive_definition':'count/linear >0 detection' if j['scale'].startswith(('counts','linear')) else 'normalized value>0 is not biological detection'}),f'00_plan_audit/{g}_background_qc.tsv.gz')
 for f in [ep,p]:inputrows.append(dict(dataset=g,path=str(f),size=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),species='human',scale=j['scale']))
save(pd.DataFrame(coverage),'00_plan_audit/coverage.tsv');save(pd.DataFrame(datasets),'00_plan_audit/dataset_manifest.tsv');save(pd.DataFrame(inputrows),'00_plan_audit/input_manifest.tsv');(B/'00_plan_audit/reference_samples.json').write_text(json.dumps(references,indent=2))
# Enumerate planned tests, including blocked tests (NA, excluded from realized BH; planned n disclosed).
h=[]
for p in plans:
 for sn in mf.signature:
  for adjustment in ['unadjusted']+(['cT','subtype','regimen','cycles'] if p['main'] else []):h.append(dict(module='clinical',dataset=p['dataset'],comparison=p['branch'],signature=sn,endpoint=p['endpoint'],method=adjustment,family=p['family']+('' if adjustment=='unadjusted' else '_adjusted')))
for (g,b),d in model[model.baseline_include].groupby(['dataset','background']):
 for state in d.state.unique():
  if state=='parental':continue
  for sn in mf.signature:h.append(dict(module='cell_models',dataset=g,comparison=b+'|'+state,signature=sn,endpoint='resistant_vs_parental',method='point_biserial',family='cell_'+str(d.loc[d.state==state,'drug'].iloc[0])))
for ds in ['GDSC2','GDSC1']:
 for sn in mf.signature:
  for ep in ['LN_IC50','AUC']:h.append(dict(module='measured_drug',dataset=ds,comparison='cisplatin',signature=sn,endpoint=ep,method='Spearman',family='measured_'+ds+'_'+ep))
ei=read('00_plan_audit/epi_patient_identity.tsv')
for ds in list(ei.dataset.unique())+['ALL_verified']:
 for sn in mf.loc[mf.signature!='MP10','signature']:
  for method in ['average_tie','official_AUCell','UCell','pseudobulk_native']:
   for adj in ['raw','depth_only','depth_cycle','nFeature_cycle'] if method!='pseudobulk_native' else ['raw']:
    h.append(dict(module='epi',dataset=ds,comparison='Tumor',signature=sn,endpoint='MP10_association',method=method+'|'+adj,family='epi_'+method+'_'+adj))
for g in ['GSE69795','GSE124305']:
 for sn in mf.signature:h.append(dict(module='paired',dataset=g,comparison='post_minus_pre',signature=sn,endpoint='paired_change',method='Wilcoxon',family='paired_state'))
save(pd.DataFrame(h),'00_plan_audit/hypothesis_table.tsv');save(pd.DataFrame(h).groupby('family',as_index=False).size().rename(columns={'size':'n_planned'}),'00_plan_audit/FDR_families.tsv')
config=dict(version='v2_expanded_post_v1_design_freeze',frozen_at=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),seed=SEED,resources=dict(total_threads=4,BLAS=1,preferred_RAM_GB=64,hard_RAM_GB=200,heavy_window='23:00-07:00 Asia/Shanghai'),genes=members[members.signature=='MP10'].gene.tolist(),coverage='>=80% total and each sign; MP10 >=16/20; CisSig/Xylinas full exact otherwise separately labelled partial; ratio all6',preprocessing='HGNC exact then unique stable IDs/unique alias; unmapped preserve original IDs; counts sum, normalized duplicate median; drop rows with nonfinite value across frozen population; retain zero-expressed genes',bulk=dict(counts='edgeR TMM logCPM prior.count1',linear='log2(x+1)',author_normalized='no re-log',ssGSEA='GSVA2.4.9 alpha0.25 normalizeFALSE then divide each signature by its score range across full scoring population; exactly matches v1 single-set normalizeTRUE; constant range NA',native_z='gene z ddof1 across frozen reference samples, zero SD excludes member then coverage rule; no outcome-specific scaling',ratio='CPM/linear normalized or inverse2^author_log2 normalized intensity adapted; zero denominator NA; no log ratio'),epi=dict(background=90011,aucMaxRank=4501,average_tie='preserve existing registered implementation',official_AUCell='normAUC TRUE; ties random seed20261111, block256, Serial',UCell='maxRank1500; ties average; package2.14.0',cycle='Seurat2019 S/G2M; AddModuleScore standard RNA normalized expression, ctrl=min(nS,nG2M), seed20261111',partial='rank x/y/covariates within sample; intercept OLS, Pearson residuals; n>=50, rank/full degrees checks',inference='patient equal Fisher-z; exact two-sided sign test H0 P(effect>0)=0.5, excludes zeros; bootstrap whole patients10000; unverified identities descriptive only',pseudobulk='sum raw counts sample then normalize/score within dataset; verified repeated samples average scores within patient, unverified descriptive rho only',random_sets='100 sets matched to MP10 expression/detection quintile bins, seed20261111, rank correlations; state control not drug-specific'),clinical=dict(primary='continuous score per1SD nonresponse; independent cohort branches',Firth='use when n<40 or min(outcome count)<10 or glm separation/nonconvergence/extremeSE>10; profile penalized LR CI/P; otherwise glm Wald CI/P',adjustment='each cT/subtype/regimen/cycles separately, baseline-only; n>=40 and >=10 minority outcomes per predictor df; missing complete cases; no post-outcome covariates',LOO='GSE212810 MP10 primary good response, same reference mean/SD',encoding_flip='MP10 reciprocal OR and unchanged P'),models=dict(delta='mean R-mean P; cohort-specific units not pooled',r='Pearson point-biserial from RNA replicate points; two-sided t test, Fisher-z CI; >=2 each group and positive variance',hedges='pooledSD Hedges small-sample J correction; noncentral t inversion CI; flag pooledSD<1e-6',biological_unit='background/study; RNA repeats never independent resistance induction',five_lines='GSE171537 5 line effects equal mean with t CI, exact two-sided sign test'),drug=dict(primary='GDSC2 LN_IC50 Spearman model-level;10000 whole-model bootstrap CI',secondary='GDSC1 and AUC separately; Pearson supplement; stableSIDM dedup median repeated assays',IC50fold='only actual paired linear units, no selection concentration substitute'),FDR='BH within frozen families across all estimable planned tests; blocked P stays NA; report planned and estimated n; additional BH across all primary unadjusted hypotheses. Sensitivity separate never substitutes primary',no_associations_inspected=True)
for f in ['clinical_plan.json','reference_samples.json','hypothesis_table.tsv','coverage.tsv','../01_signatures/members.tsv','../01_signatures/signature_manifest.tsv']:
 p=B/'00_plan_audit'/f;config.setdefault('locked_files',{})[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
p=B/'00_plan_audit/frozen_config.json';p.write_text(json.dumps(config,indent=2,ensure_ascii=False));(B/'00_plan_audit/frozen_config.sha256').write_text(hashlib.sha256(p.read_bytes()).hexdigest());print('FROZEN',config['frozen_at'],len(h),'planned tests',len(datasets),'expression populations')
