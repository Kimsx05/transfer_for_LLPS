from env import *
from scipy.optimize import brentq
sc=pd.concat([pd.read_csv(p,sep='\t') for p in (B/'00_plan_audit/scores').glob('*_scores.tsv')],ignore_index=True);model=read('03_cell_models/model_identity.tsv');mf=read('01_signatures/signature_manifest.tsv');rows=[];heat=[];joined=model.merge(sc,on=['dataset','matrix_id'],how='inner',validate='one_to_many');save(joined,'03_cell_models/sample_scores.tsv.gz')
def corr_ci(x,y,method='spearman',nboot=10000):
 x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float);ok=np.isfinite(x)&np.isfinite(y);x=x[ok];y=y[ok];n=len(x)
 if n<4 or np.std(x)==0 or np.std(y)==0:return np.nan,np.nan,np.nan,np.nan,n
 r,p=(stats.spearmanr if method=='spearman' else stats.pearsonr)(x,y);rng=np.random.default_rng(SEED);ii=rng.integers(n,size=(nboot,n));xx=x[ii];yy=y[ii]
 if method=='spearman':xx=stats.rankdata(xx,axis=1);yy=stats.rankdata(yy,axis=1)
 xx-=xx.mean(axis=1)[:,None];yy-=yy.mean(axis=1)[:,None];den=np.sqrt((xx*xx).sum(axis=1)*(yy*yy).sum(axis=1));bs=np.divide((xx*yy).sum(axis=1),den,out=np.full(nboot,np.nan),where=den>0);lo,hi=np.nanquantile(bs,[.025,.975]);return r,lo,hi,p,n
for (g,b),md in model[model.baseline_include].groupby(['dataset','background']):
 for state in md.state.unique():
  if state=='parental':continue
  d=joined[(joined.dataset==g)&(joined.background==b)&joined.baseline_include];drug=md.loc[md.state==state,'drug'].iloc[0]
  for sn in mf.signature:
   h=d[d.signature==sn];a=h.loc[h.state=='parental','score'].dropna().values;c=h.loc[h.state==state,'score'].dropna().values;na,nc=len(a),len(c);row=dict(module='cell_models',dataset=g,study=md.study.iloc[0],background=b,state=state,signature=sn,drug=drug,endpoint='resistant_vs_parental',family='cell_'+drug,n_parental=na,n_resistant=nc,n_RNA=na+nc,n_background=1,n_derived_models=1,n_studies=1,independent_unit='RNA/passage replicates nested within one derived model',delta=np.mean(c)-np.mean(a) if na and nc else np.nan,r=np.nan,r_low=np.nan,r_high=np.nan,g=np.nan,g_low=np.nan,g_high=np.nan,P=np.nan,status='BLOCKED_COVERAGE',reason='score unavailable or parental expression absent')
   if na and nc:
    row.update(parental_mean=float(np.mean(a)),resistant_mean=float(np.mean(c)),parental_sd=float(np.std(a,ddof=1)) if na>1 else np.nan,resistant_sd=float(np.std(c,ddof=1)) if nc>1 else np.nan,status='COMPLETE',reason='')
    if min(na,nc)>=2 and np.std(np.r_[a,c])>0:
     y=np.r_[np.zeros(na),np.ones(nc)];r,p=stats.pearsonr(y,np.r_[a,c]);row.update(r=r,P=p);z=np.arctanh(np.clip(r,-.999999999,.999999999));row.update(r_low=np.tanh(z-1.96/np.sqrt(na+nc-3)),r_high=np.tanh(z+1.96/np.sqrt(na+nc-3)))
     df=na+nc-2;sp=np.sqrt(((na-1)*np.var(a,ddof=1)+(nc-1)*np.var(c,ddof=1))/df);J=np.exp(__import__('scipy').special.gammaln(df/2)-.5*np.log(df/2)-__import__('scipy').special.gammaln((df-1)/2));row['tiny_within_SD']=sp<1e-6
     if sp>0:
      dlt=row['delta']/sp;t=dlt/np.sqrt(1/na+1/nc);row['g']=J*dlt
      try:row['g_low']=J*brentq(lambda nc_:stats.nct.cdf(t,df,nc_)-.975,-1e4,1e4)*np.sqrt(1/na+1/nc);row['g_high']=J*brentq(lambda nc_:stats.nct.cdf(t,df,nc_)-.025,-1e4,1e4)*np.sqrt(1/na+1/nc)
      except ValueError:row['g_CI_reason']='noncentral t inversion failed'
     se=np.sqrt(np.var(a,ddof=1)/na+np.var(c,ddof=1)/nc);dfw=(np.var(a,ddof=1)/na+np.var(c,ddof=1)/nc)**2/((np.var(a,ddof=1)/na)**2/(na-1)+(np.var(c,ddof=1)/nc)**2/(nc-1)) if se>0 else np.nan;row.update(delta_low=row['delta']-stats.t.ppf(.975,dfw)*se,delta_high=row['delta']+stats.t.ppf(.975,dfw)*se)
    else:row['reason']='1+1 or no variance: descriptive delta only; r/g/CI/P not estimable'
   rows.append(row)
  ep=B/f'00_plan_audit/scores/{g}_signature_expression.tsv.gz'
  if ep.exists():
   x=pd.read_csv(ep,sep='\t').set_index('gene');a=[i for i in md.loc[md.state=='parental','matrix_id'] if i in x];c=[i for i in md.loc[md.state==state,'matrix_id'] if i in x]
   if a and c:
    for gene in read('01_signatures/members.tsv').query('signature=="MP10"').gene:
     heat.append(dict(dataset=g,background=b,state=state,gene=gene,delta_log_expression=float(x.loc[gene,c].mean()-x.loc[gene,a].mean()) if gene in x.index else np.nan))
r=pd.DataFrame(rows);r['large_standardized_effect']=r.g.abs()>5;save(r,'03_cell_models/contrast_results.tsv');save(pd.DataFrame(heat),'03_cell_models/MP10_gene_changes.tsv');summ=[]
for (state,sn),d in r[(r.dataset=='GSE171537')].groupby(['state','signature']):
 x=d.delta.dropna().values;n=len(x);nz=x[x!=0];se=np.std(x,ddof=1)/np.sqrt(n) if n>1 else np.nan;summ.append(dict(dataset='GSE171537',state=state,signature=sn,n_backgrounds=n,n_positive=int(sum(x>0)),n_negative=int(sum(x<0)),mean_delta=np.mean(x) if n else np.nan,CI_low=np.mean(x)-stats.t.ppf(.975,n-1)*se if n>1 else np.nan,CI_high=np.mean(x)+stats.t.ppf(.975,n-1)*se if n>1 else np.nan,P=stats.binomtest(sum(nz>0),len(nz),.5).pvalue if len(nz) else np.nan,method='equal5 backgrounds; two-sided exact sign',family='five_background_'+state))
save(pd.DataFrame(summ),'03_cell_models/GSE171537_background_summary.tsv')
drug=read('04_measured_drug_response/GDSC_cisplatin_measured.tsv');dr=[];mat=[]
for ds,dd in drug.groupby('DATASET'):
 d=dd.groupby('SANGER_MODEL_ID',as_index=False).agg(LN_IC50=('LN_IC50','median'),AUC=('AUC','median'),n_assay=('LN_IC50','size'),cell_line=('CELL_LINE_NAME','first'),dose_min=('MIN_CONC','min'),dose_max=('MAX_CONC','max'));d=d.merge(sc[sc.dataset=='GDSC_BLCA'],left_on='SANGER_MODEL_ID',right_on='matrix_id',validate='one_to_many');d['screen']=ds;mat.append(d)
 for sn in mf.signature:
  a=d[d.signature==sn]
  for endpoint in ['LN_IC50','AUC']:
   for method in ['spearman','pearson']:
    rho,lo,hi,p,n=corr_ci(a.score,a[endpoint],method);dr.append(dict(module='measured_drug',dataset=ds,signature=sn,endpoint=endpoint,method=method,r=rho,CI_low=lo,CI_high=hi,P=p,n=n,independent_unit='stable SIDM model',family='measured_'+ds+'_'+endpoint+('' if method=='spearman' else '_Pearson'),status='COMPLETE' if np.isfinite(rho) else 'BLOCKED_COVERAGE',reason='CisSig training-source overlap' if sn=='CisSig19' else '',main=ds=='GDSC2' and endpoint=='LN_IC50' and method=='spearman'))
save(pd.concat(mat),'04_measured_drug_response/GDSC_matched_scores.tsv.gz');save(pd.DataFrame(dr),'04_measured_drug_response/GDSC_associations.tsv');print('COMPLETE',len(r),'model contrasts',len(dr),'drug correlations')
