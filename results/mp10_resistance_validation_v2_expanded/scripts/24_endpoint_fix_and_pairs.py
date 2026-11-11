from env import *
import shutil
p=B/'05_clinical/patient_endpoint_audit.tsv';cl=pd.read_csv(p,sep='\t');original=hashlib.sha256(p.read_bytes()).hexdigest();changes=[]
rawmeta=read('00_plan_audit/GSE69795_samples.tsv').set_index('GSM')
for i,r in cl[cl.dataset=='GSE69795'].iterrows():
 raw=str(r.pathology_raw);s=re.sub(r'\s+','',raw).upper();mt=re.fullmatch(r'T(A|IS|[0-4][AB]?)N([0-3X])',s)
 if mt and r.pretreatment and (r.surgery or rawmeta.loc[r.GSM,'cystectomy']=='Upper tract'):
  t,n=mt.groups();good=1 if t in ['0','A','IS','1'] and n=='0' else 0 if n in ['1','2','3'] or re.match('[234]',t) else np.nan;pcr=1 if t=='0' and n=='0' else 0 if t in ['A','IS','1'] or n in ['1','2','3'] or re.match('[234]',t) else np.nan
  if pd.isna(r.author_response) or pd.isna(r.good_response_ypT_le1_N0):
   changes.append(dict(GSM=r.GSM,raw=raw,old_good=r.good_response_ypT_le1_N0,new_good=good,old_author=r.author_response,new_author=good,reason='literal TaN0 / Ta N0 should follow already frozen tri-state rules; numeric stage groups remain untouched'))
   cl.loc[i,'surgery']=True
   cl.loc[i,['ypT','ypN','good_response_ypT_le1_N0','pCR_ypT0N0','author_response']]=[t,n,good,pcr,good]
if changes:
 ad=B/'00_plan_audit/amendments';ad.mkdir(exist_ok=True)
 for src in ['05_clinical/patient_endpoint_audit.tsv','05_clinical/model_results.tsv']:
  dest=ad/('before_literal_pTN_fix_'+Path(src).name)
  if not dest.exists():shutil.copy2(B/src,dest)
 save(pd.DataFrame(changes),'00_plan_audit/amendments/literal_pTN_fix.tsv');save(cl,'05_clinical/patient_endpoint_audit.tsv');(ad/'literal_pTN_fix.json').write_text(json.dumps(dict(time=datetime.datetime.now().isoformat(),reason='implementation correction discovered by auditing raw pstage strings after initial run; no rule/threshold/sample selection change',before_sha256=original,after_sha256=hashlib.sha256(p.read_bytes()).hexdigest()),indent=2))
# Explicit GEO pre/post crosswalk; both samples scored on same cohort/reference scale.
m=read('00_plan_audit/GSE69795_samples.tsv');pre=m[m['sample procured']=='TUR'];post=m[m['sample procured']=='surgery'];pair=post[['GSM','title','pretreatment sample record-pair']].merge(pre[['GSM','title']],left_on='pretreatment sample record-pair',right_on='title',suffixes=('_post','_pre'),validate='one_to_one');assert len(pair)==23
pair=pair.rename(columns={'GSM_post':'post_id','GSM_pre':'pre_id','title_pre':'patient_id'});pair['source']='GEO pretreatment sample record-pair';save(pair,'06_paired_states/GSE69795_pair_identity.tsv');sc=read('00_plan_audit/scores/GSE69795_scores.tsv');rows=[];points=[];mf=read('01_signatures/signature_manifest.tsv')
for sn in mf.signature:
 ss=sc[sc.signature==sn].set_index('matrix_id').score;a=pair.copy();a['pre']=a.pre_id.map(ss);a['post']=a.post_id.map(ss);a['delta']=a.post-a.pre;a['signature']=sn;points.append(a);v=a.delta.dropna().values;n=len(v);rng=np.random.default_rng(SEED);lo,hi=np.quantile(v[rng.integers(n,size=(10000,n))].mean(axis=1),[.025,.975]) if n else [np.nan,np.nan];p=stats.wilcoxon(v,alternative='two-sided',method='auto').pvalue if n and np.any(v!=0) else np.nan
 rows.append(dict(module='paired',dataset='GSE69795',signature=sn,endpoint='post_minus_pre',n=n,mean_delta=v.mean() if n else np.nan,CI_low=lo,CI_high=hi,P=p,independent_unit='patient pair; all original source cohort, BLCA/UTUC mix retained in audit',family='paired_state',status='COMPLETE' if n else 'BLOCKED_COVERAGE',reason='residual tissue selection and composition confounding; not baseline prediction'))
for sn in mf.signature:rows.append(dict(module='paired',dataset='GSE124305',signature=sn,endpoint='post_minus_pre',n=0,P=np.nan,family='paired_state',status='BLOCKED_LABELS',reason='133 post expression scored; actual patient crosswalk to305 pre-GSE87304 absent; 116 reported pairs not reconstructable from numeric IDs'))
save(pd.concat(points),'06_paired_states/GSE69795_paired_scores.tsv');save(pd.DataFrame(rows),'06_paired_states/paired_results.tsv');print('endpoint changes',len(changes),'pairs',len(pair))
