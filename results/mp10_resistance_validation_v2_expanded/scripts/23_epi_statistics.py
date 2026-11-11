from env import *
from scipy.stats import rankdata
f=read('02_epi/cell_average_tie_scores.tsv.gz');mf=read('01_signatures/signature_manifest.tsv');cov=read('00_plan_audit/epi_coverage.tsv');signatures=[s for s in cov.loc[cov.status=='READY','signature'] if s in f.columns];cycles=B/'02_epi/cell_cycle_scores.tsv.gz'
if cycles.exists():f=f.merge(pd.read_csv(cycles,sep='\t'),on='cell_id',how='left',validate='one_to_one')
rows=[];method_tables=[('average_tie',f)]
p=B/'02_epi/cell_official_scores.tsv.gz'
if p.exists():
 qc=B/'00_plan_audit/official_full_QC.json'
 if not qc.exists() or qc.stat().st_mtime<p.stat().st_mtime:
  import subprocess
  subprocess.run([sys.executable,str(B/'scripts/38_official_method_qc.py')],check=True)
 off=pd.read_csv(p,sep='\t')
 for method,ss in off.groupby('method'):
  wide=ss.pivot(index='cell_id',columns='signature',values='score');meta=f.drop(columns=[s for s in signatures if s in f]);method_tables.append((method,meta.merge(wide,on='cell_id',how='left',validate='one_to_one')))
# Tumor follows original labels; CNV is an annotation sensitivity, not relabelling all Epi malignant.
cfg1=json.loads((V/'config/frozen_v1.json').read_text());tumor=f.Type.isin(cfg1['epi']['tumor'])
for method,data in method_tables:
 subsets={'Tumor':data.Type.isin(cfg1['epi']['tumor']),'Normal':data.Type.eq('Normal'),'CNV_supported_Tumor':data.Type.isin(cfg1['epi']['tumor'])&data.malignant_like_status.eq('malignant_like_main')}
 for tp in sorted(set(data.Type.dropna())-{'Tumor','Normal'}):subsets['recorded_Type_'+tp]=data.Type.eq(tp)
 if 'T_stage' in data:
  for stage in sorted(set(data.T_stage.dropna())-{'Normal'}):subsets['recorded_T_stage_'+str(stage)]=data.Type.isin(cfg1['epi']['tumor'])&data.T_stage.eq(stage)
 for subset,mask in subsets.items():
  for sample,d in data[mask].groupby('sample'):
   for adj,cs in [('raw',[]),('depth_only',['nCount_RNA']),('depth_cycle',['nCount_RNA','S.Score','G2M.Score']),('nFeature_cycle',['nFeature_RNA','S.Score','G2M.Score'])]:
    for sn in [s for s in signatures if s!='MP10' and s in d]:
     base=dict(dataset=d.dataset.iloc[0],sample=sample,patient=d.patient.iloc[0],verified=bool(d.verified.iloc[0]),subset=subset,signature=sn,method=method,adjustment=adj,n_cells=len(d),n_used=0,rho=np.nan,reason='',MP10_sd=d.MP10.std(),signature_sd=d[sn].std(),mean_detected_genes=d.get(sn+'_detected',pd.Series(dtype=float)).mean(),mean_detection_rate=d.get(sn+'_detected',pd.Series(dtype=float)).mean()/float(cov.loc[cov.signature==sn,'n_present'].iloc[0]))
     if any(c not in d for c in cs):base['reason']='DEFERRED_RESOURCE: cell-cycle computation pending';rows.append(base);continue
     a=d[['MP10',sn]+cs].replace([np.inf,-np.inf],np.nan).dropna();base['n_used']=len(a)
     if len(a)<50 or a.MP10.nunique()<2 or a[sn].nunique()<2:base['reason']='n<50 or constant score';rows.append(base);continue
     ranked=rankdata(a.values,axis=0);xy=ranked[:,:2]
     if cs:
      X=np.c_[np.ones(len(a)),ranked[:,2:]];rank=np.linalg.matrix_rank(X)
      if rank<len(cs)+1 or len(a)-rank<3:base['reason']='rank deficient covariates or insufficient residual degrees';rows.append(base);continue
      xy=xy-X@np.linalg.lstsq(X,xy,rcond=None)[0]
     if np.any(np.std(xy,axis=0)<1e-10):base['reason']='zero residual variance'
     else:base['rho']=np.corrcoef(xy.T)[0,1]
     rows.append(base)
eff=pd.DataFrame(rows);save(eff,'02_epi/sample_effects.tsv');pat=[];summ=[];loo=[];desc=[]
def summarize(d,ds,meta):
 vals=d.rho.dropna().values;n=len(vals);rng=np.random.default_rng(SEED)
 if n:
  z=np.arctanh(np.clip(vals,-.999999999,.999999999));point=np.tanh(z.mean());lo,hi=np.tanh(np.quantile(z[rng.integers(n,size=(10000,n))].mean(axis=1),[.025,.975])) if n>=2 else [np.nan,np.nan];nz=vals[vals!=0];p=stats.binomtest((nz>0).sum(),len(nz),.5).pvalue if len(nz)>=2 else np.nan
 else:point=lo=hi=p=np.nan
 return dict(module='epi',dataset=ds,**meta,n=n,rho=point,CI_low=lo,CI_high=hi,P=p,n_positive=int(sum(vals>0)),n_negative=int(sum(vals<0)),independent_unit='verified patient equal Fisher-z',family='epi_'+meta['method']+'_'+meta['adjustment']+('' if meta['subset']=='Tumor' else '_'+meta['subset']),status='COMPLETE' if n>=2 else 'BLOCKED_DESIGN',reason='' if n>=2 else 'fewer than2 verified eligible patients')
for key,d in eff.groupby(['dataset','subset','signature','method','adjustment']):
 ds,subset,sn,method,adj=key;meta=dict(subset=subset,signature=sn,method=method,adjustment=adj);desc.append(dict(dataset=ds,**meta,n_specimens=d.rho.notna().sum(),median_rho=d.rho.median(),min_rho=d.rho.min(),max_rho=d.rho.max(),reason='specimen descriptive; no patient inference for unverified identities'))
 dd=d[d.verified&d.rho.notna()];pr=[]
 for pid,h in dd.groupby('patient'):
  rho=np.tanh(np.mean(np.arctanh(np.clip(h.rho,-.999999999,.999999999))));r=dict(dataset=ds,patient=pid,**meta,n_samples=len(h),rho=rho);pr.append(r);pat.append(r)
 pt=pd.DataFrame(pr,columns=['patient','rho']);summ.append(summarize(pt,ds,meta))
 if len(pt)>=3 and subset=='Tumor':
  for pid in pt.patient:loo.append(dict(excluded_patient=pid,**summarize(pt[pt.patient!=pid],ds,meta)))
pp=pd.DataFrame(pat);save(pp,'02_epi/patient_effects.tsv');save(pd.DataFrame(desc),'02_epi/descriptive_dataset_effects.tsv')
if len(pp):
 for key,d in pp.groupby(['subset','signature','method','adjustment']):
  subset,sn,method,adj=key;meta=dict(subset=subset,signature=sn,method=method,adjustment=adj);summ.append(summarize(d,'ALL_verified',meta))
  for ds in sorted(f.dataset.unique()):loo.append(dict(excluded_dataset=ds,**summarize(d[d.dataset!=ds],'ALL_verified',meta)))
save(pd.DataFrame(summ),'02_epi/patient_summary.tsv');save(pd.DataFrame(loo),'02_epi/leave_one_diagnostics.tsv')
# Native sample scores: never duplicate one sample score across its cells.
sm=read('02_epi/pseudobulk_samples.tsv');scores=[]
for p in (B/'00_plan_audit/scores').glob('EpiPB_*_scores.tsv'):
 a=pd.read_csv(p,sep='\t');a['dataset']=a.dataset.str.replace('EpiPB_','',regex=False);scores.append(a)
if scores:
 sc=pd.concat(scores);sc=sc.rename(columns={'matrix_id':'sample'}).merge(sm,on=['dataset','sample'],validate='many_to_one');save(sc,'02_epi/pseudobulk_scores.tsv');pb=[];rand=[]
 random=read('02_epi/random_matched_scores.tsv') if (B/'02_epi/random_matched_scores.tsv').exists() else pd.DataFrame()
 groups=[]
 for ds,all_d in sc[sc.n_cells>=50].groupby('dataset'):
  groups.append((ds,'Tumor',all_d[all_d.Type.isin(cfg1['epi']['tumor'])]))
  for tp,dt in all_d.groupby('Type'):groups.append((ds,'Normal' if tp=='Normal' else 'recorded_Type_'+tp,dt))
  if 'T_stage' in all_d:
   for stage,dt in all_d[all_d.Type.isin(cfg1['epi']['tumor'])].groupby('T_stage'):groups.append((ds,'recorded_T_stage_'+str(stage),dt))
 for ds,subset,d in groups:
  if d.empty:continue
  wide=d.pivot(index='sample',columns='signature',values='score');mm=d[['sample','patient','verified']].drop_duplicates().set_index('sample');wd=wide.join(mm)
  for sn in [s for s in wide if s!='MP10']:
   a=wd[['MP10',sn,'patient','verified']].dropna(subset=['MP10',sn]);r,p=stats.spearmanr(a.MP10,a[sn]) if len(a)>=3 else (np.nan,np.nan);vp=a[a.verified].groupby('patient')[['MP10',sn]].mean();n=len(vp);rr=lo=hi=pv=np.nan
   if n>=4 and vp.MP10.nunique()>1 and vp[sn].nunique()>1:
    rr,pv=stats.spearmanr(vp.MP10,vp[sn]);rng=np.random.default_rng(SEED);ii=rng.integers(n,size=(10000,n));xx=rankdata(vp.MP10.values[ii],axis=1);yy=rankdata(vp[sn].values[ii],axis=1);xx-=xx.mean(axis=1)[:,None];yy-=yy.mean(axis=1)[:,None];den=np.sqrt((xx*xx).sum(1)*(yy*yy).sum(1));bs=np.divide((xx*yy).sum(1),den,out=np.full(10000,np.nan),where=den>0);lo,hi=np.nanquantile(bs,[.025,.975])
   pb.append(dict(module='epi',dataset=ds,subset=subset,signature=sn,method='pseudobulk_native',adjustment='raw',n=n,n_specimens=len(a),rho=rr,CI_low=lo,CI_high=hi,P=pv,descriptive_specimen_rho=r,independent_unit='verified patient; repeated samples equal mean score',family='epi_pseudobulk_native_raw'+('' if subset=='Tumor' else '_'+subset),status='COMPLETE' if np.isfinite(rr) else 'BLOCKED_DESIGN',reason='' if np.isfinite(rr) else 'fewer than4 verified patients or constant score'))
  if len(random):
   for ctrl,rc in random[(random.dataset==ds)&random['sample'].isin(wide.index)].groupby('control'):
    a=wide[['MP10']].join(rc.set_index('sample').score).dropna();rho=stats.spearmanr(a.MP10,a.score).statistic if len(a)>=3 else np.nan;rand.append(dict(dataset=ds,subset=subset,control=ctrl,n_specimens=len(a),descriptive_rho=rho,matching='pooled pseudobulk mean expression and pseudobulk detection quintiles; not per-cell detection'))
 save(pd.DataFrame(pb),'02_epi/pseudobulk_associations.tsv');save(pd.DataFrame(rand),'02_epi/random_control_associations.tsv')
print('COMPLETE',len(eff),'sample effects',len(pp),'verified patient effects',flush=True)
