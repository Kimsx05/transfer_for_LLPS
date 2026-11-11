.libPaths(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/refdir/Rlib',.libPaths()))
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(GSVA);library(edgeR)})
setDTthreads(1);set.seed(20261111)
B='/home/data/t070721/codex_workspace/Bladder Metabolism/results/mp10_resistance_validation_v2_expanded'
rd=function(p)fread(file.path(B,p));wt=function(x,p)fwrite(x,file.path(B,p),sep='\t',na='NA')
cfg=fromJSON(file.path(B,'00_plan_audit/frozen_config.json'));stopifnot(digest::digest(file=file.path(B,'00_plan_audit/frozen_config.json'),algo='sha256')==readLines(file.path(B,'00_plan_audit/frozen_config.sha256')))
members=rd('01_signatures/members.tsv');mf=rd('01_signatures/signature_manifest.tsv');dm=rd('00_plan_audit/dataset_manifest.tsv');cov=rd('00_plan_audit/coverage.tsv');refs=fromJSON(file.path(B,'00_plan_audit/reference_samples.json'));dir.create(file.path(B,'00_plan_audit/scores'),showWarnings=FALSE)
score_one=function(g,path,scale,reference=NULL){
 xx=fread(path);genes=xx[[1]];x=as.matrix(xx[,-1]);rownames(x)=genes;storage.mode(x)='double';x=x[apply(x,1,function(z)all(is.finite(z))),,drop=FALSE];stopifnot(!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)))
 if(scale=='panel_counts')return(invisible(NULL))
 ids=intersect(if(is.null(reference))colnames(x) else reference,colnames(x));stopifnot(length(ids)>1)
 linear=NULL
 if(scale=='counts'){
  stopifnot(all(x>=0));y=calcNormFactors(DGEList(x),method='TMM');linear=cpm(y,log=FALSE);expr=cpm(y,log=TRUE,prior.count=1);wt(data.table(matrix_id=colnames(x),library_size=y$samples$lib.size,TMM=y$samples$norm.factors),paste0('00_plan_audit/scores/',g,'_normalization.tsv'))
 }else if(startsWith(scale,'linear')){stopifnot(all(x>=0));linear=x;expr=log2(x+1)}else{expr=x;if(scale=='author_log2')linear=2^x}
 means=rowMeans(expr[,ids,drop=FALSE]);sds=apply(expr[,ids,drop=FALSE],1,sd);z=sweep(sweep(expr,1,means,'-'),1,sds,'/');z[!is.finite(z)]=NA_real_
 cv=cov[dataset==g];sets=list();sc=list();qcs=list()
 for(k in seq_len(nrow(mf))){
  sn=mf$signature[k];m=mf[k];r=members[signature==sn];rr=r[!is.na(gene)&gene%in%rownames(expr)];method=m$method;status=cv[signature==sn,status];reason=cv[signature==sn,reason]
  if(!length(status))status='READY'
  if(status!='READY'){qcs[[sn]]=data.table(dataset=g,signature=sn,status=status,reason=reason);next}
  if(method=='ssGSEA'){sets[[sn]]=unique(rr$gene);next}
  if(method=='linear_ratio'){
   if(is.null(linear)){qcs[[sn]]=data.table(dataset=g,signature=sn,status='BLOCKED_DESIGN',reason='meaningful nonnegative linear scale unavailable');next}
   den=colMeans(linear[c('LRRC8A','LRRC8B','LRRC8C','LRRC8D','LRRC8E'),,drop=FALSE]);s=linear['NPEPPS',]/den;s[den<=0]=NA_real_;sc[[sn]]=s
  }else{
   rr=rr[is.finite(sds[gene])&sds[gene]>0];valid=nrow(rr)/nrow(r)>=.8 && all(vapply(c(-1,1),function(sign){a=r[sign(weight)==sign];!nrow(a)||sum(sign(rr$weight)==sign)/nrow(a)>=.8},logical(1)))
   if(!valid){qcs[[sn]]=data.table(dataset=g,signature=sn,status='BLOCKED_COVERAGE',reason='variable genes in fixed reference population below signed80%');next}
   zz=z[rr$gene,,drop=FALSE]
   if(method=='median_z')s=apply(zz,2,median)
   else if(method=='weighted_z')s=colMeans(zz*rr$weight)
   else {up=colMeans(zz[rr$weight>0,,drop=FALSE]);dn=colMeans(zz[rr$weight<0,,drop=FALSE]);s=up-dn;wt(data.table(matrix_id=colnames(x),UP_mean_z=up,DN_mean_z=dn),paste0('00_plan_audit/scores/',g,'_',sn,'_components.tsv'))}
   sc[[sn]]=s
  }
  qcs[[sn]]=data.table(dataset=g,signature=sn,status='COMPLETE',reason='',n_members_used=nrow(rr),partial_after_geneSD=nrow(rr)<nrow(r))
 }
 if(length(sets)){
  par=ssgseaParam(expr,sets,minSize=1,maxSize=100000,alpha=.25,normalize=FALSE,checkNA='yes');raw=gsva(par,verbose=FALSE,BPPARAM=BiocParallel::SerialParam())
  for(sn in rownames(raw)){range=diff(range(raw[sn,]));sc[[sn]]=if(range>0)raw[sn,]/range else rep(NA_real_,ncol(expr));qcs[[sn]]=data.table(dataset=g,signature=sn,status=if(range>0)'COMPLETE' else 'BLOCKED_DESIGN',reason=if(range>0)'' else 'constant ssGSEA range',ssGSEA_range=range,n_members_used=length(sets[[sn]]))}
  # Exact equality to old single-gene-set GSVA normalization is verified on each cohort MP10.
  if('MP10'%in%names(sets)){old=gsva(ssgseaParam(expr,list(MP10=sets$MP10),minSize=16,alpha=.25,normalize=TRUE),verbose=FALSE,BPPARAM=BiocParallel::SerialParam());stopifnot(max(abs(old[1,]-sc$MP10),na.rm=TRUE)<1e-8)}
 }
 out=rbindlist(lapply(names(sc),function(sn){v=sc[[sn]];mu=mean(v[match(ids,colnames(x))],na.rm=TRUE);ss=sd(v[match(ids,colnames(x))],na.rm=TRUE);data.table(dataset=g,matrix_id=colnames(x),signature=sn,score=as.numeric(v),score_z=if(is.finite(ss)&&ss>0)(as.numeric(v)-mu)/ss else NA_real_,reference_mean=mu,reference_sd=ss,reference_n=length(ids))}));wt(out,paste0('00_plan_audit/scores/',g,'_scores.tsv'));wt(rbindlist(qcs,fill=TRUE),paste0('00_plan_audit/scores/',g,'_score_qc.tsv'))
 selected=unique(members$gene);gg=intersect(selected,rownames(expr));wt(data.table(gene=gg,reference_mean=means[gg],reference_sd=sds[gg],expr[gg,,drop=FALSE],check.names=FALSE),paste0('00_plan_audit/scores/',g,'_signature_expression.tsv.gz'))
 cat('COMPLETE',g,ncol(expr),nrow(expr),'scores',length(sc),'\n');invisible(out)
}
amend=file.path(B,'00_plan_audit/amendments/GSE285648_DMSO_reference.json');if(file.exists(amend)){aa=fromJSON(amend);refs[['GSE285648']]=aa$reference_samples}
args=commandArgs(TRUE)
for(g in dm$dataset){if(length(args)&&!g%in%args)next;if(file.exists(file.path(B,paste0('00_plan_audit/scores/',g,'_scores.tsv')))&&Sys.getenv('RECOMPUTE')!='1')next;d=dm[dataset==g];score_one(g,d$expression_path,d$scale,refs[[g]]);gc()}
capture.output(sessionInfo(),file=file.path(B,'00_plan_audit/R_sessionInfo.txt'))
