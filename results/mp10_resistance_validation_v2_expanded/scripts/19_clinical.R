.libPaths(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/refdir/Rlib',.libPaths()));Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
suppressPackageStartupMessages({library(data.table);library(jsonlite)});setDTthreads(1);B='/home/data/t070721/codex_workspace/Bladder Metabolism/results/mp10_resistance_validation_v2_expanded';source(file.path(B,'scripts/firth.R'));rd=function(p)fread(file.path(B,p));wt=function(x,p)fwrite(x,file.path(B,p),sep='\t',na='NA');plans=fromJSON(file.path(B,'00_plan_audit/clinical_plan.json'),simplifyVector=FALSE);cl=rd('05_clinical/patient_endpoint_audit.tsv');mf=rd('01_signatures/signature_manifest.tsv');sc=rbindlist(lapply(list.files(file.path(B,'00_plan_audit/scores'),pattern='_scores.tsv$',full.names=TRUE),fread));res=list();used=list();dist=list();flips=list();loo=list()
fit1=function(d,adj='unadjusted'){
 x=d$score_z;y=1-d$response;n=nrow(d);if(length(unique(y))<2||n<5||!is.finite(sd(x))||sd(x)==0)stop('insufficient variable/evaluable independent patients')
 if(adj=='unadjusted')X=cbind('(Intercept)'=1,score=x) else X=model.matrix(~score_z+covariate,d)
 if(qr(X)$rank<ncol(X))stop('rank deficient design')
 mod=glm.fit(X,y,family=binomial());se=sqrt(diag(solve(crossprod(X*sqrt(mod$weights)))))[2];fir=n<40||min(table(y))<10||!mod$converged||any(mod$fitted.values<1e-8|mod$fitted.values>1-1e-8)||!is.finite(se)||se>10
 if(fir){a=firth_profile(X,y,2);method='Firth profile penalized LR'}else{b=mod$coefficients[2];a=list(beta=b,low=b-qnorm(.975)*se,high=b+qnorm(.975)*se,P=2*pnorm(-abs(b/se)));method='logistic MLE Wald'}
 list(beta=a$beta,OR=exp(a$beta),CI_low=exp(a$low),CI_high=exp(a$high),P=a$P,estimator=method,converged=TRUE)
}
for(p in plans){
 g=p$dataset;ids=unlist(p$matrix_ids);meta=cl[dataset==g&matrix_id%in%ids];if(g=='GSE244266'){meta[,matrix_id:=patient_id];meta=unique(meta,by='patient_id')};stopifnot(!anyDuplicated(meta$patient_id))
 for(sn in mf$signature){
  scores=sc[dataset==g&signature==sn];d=merge(meta,scores,by=c('dataset','matrix_id'),all.x=TRUE);d[,response:=as.numeric(get(p$endpoint))];eligible=nrow(meta)
  for(adj in c('unadjusted',if(isTRUE(p$main))c('cT','subtype','regimen','cycles') else NULL)){
   a=copy(d);base=list(module='clinical',dataset=g,branch=p$branch,signature=sn,endpoint=p$endpoint,adjustment=adj,main=p$main,family=paste0(p$family,if(adj!='unadjusted')'_adjusted' else ''),n_eligible=eligible,n=0L,events=NA_integer_,responders=NA_integer_,n_missing=eligible,OR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,status='BLOCKED_LABELS',reason='',independent_unit='patient')
   if(adj!='unadjusted'){
    if(adj=='cT'){a[,covariate:=as.numeric(sub('.*?([1-4]).*','\\1',as.character(cT)))];a[!covariate%in%1:4,covariate:=NA_real_]}
    else if(adj=='cycles')a[,covariate:=as.numeric(cycles)]
    else {a[,covariate:=as.character(get(adj))];a[is.na(covariate)|covariate%in%c('unknown','NA',''),covariate:=NA_character_];a[,covariate:=factor(covariate)]}
    a=a[!is.na(covariate)]
   }
   a=a[is.finite(score_z)&!is.na(response)];base$n=nrow(a);base$events=sum(a$response==0);base$responders=sum(a$response==1);base$n_missing=eligible-nrow(a)
   err=NULL
   if(!nrow(scores))err='score unavailable: inspect coverage and signature/linear-scale audit'
   else if(!nrow(a))err='no evaluable endpoint/complete baseline covariates'
   else if(adj!='unadjusted'&&(nrow(a)<40||length(unique(a$covariate))<2||min(table(a$response))<10*(if(is.factor(a$covariate))nlevels(droplevels(a$covariate)) else 2)))err='prespecified adjustment support threshold not met (n>=40, minority>=10 per predictor df)'
   if(is.null(err))ans=tryCatch(fit1(a,adj),error=function(e){err<<-conditionMessage(e);NULL}) else ans=NULL
   if(!is.null(ans)){base=modifyList(base,ans);base$status='COMPLETE';base$reason='';used[[length(used)+1]]=data.table(dataset=g,branch=p$branch,signature=sn,endpoint=p$endpoint,adjustment=adj,patient_id=a$patient_id,matrix_id=a$matrix_id,nonresponse=1-a$response,score=a$score,score_z=a$score_z)}else{base$reason=err;base$status=if(!nrow(scores))'BLOCKED_COVERAGE' else if(adj!='unadjusted')'BLOCKED_DESIGN' else 'BLOCKED_LABELS'}
   res[[length(res)+1]]=as.data.table(base)
   if(adj=='unadjusted'&&nrow(a)){
    z=a[,.(n=.N,mean=mean(score),median=median(score),sd=sd(score)),by=response];z[,`:=`(dataset=g,branch=p$branch,signature=sn,endpoint=p$endpoint,rank_P=if(length(unique(a$response))==2)wilcox.test(score~response,a,exact=FALSE)$p.value else NA_real_)];dist[[length(dist)+1]]=z
    if(sn=='MP10'&&!is.null(ans)){
     aa=copy(a);aa[,response:=1-response];flip=fit1(aa);flips[[length(flips)+1]]=data.table(dataset=g,branch=p$branch,endpoint=p$endpoint,OR=ans$OR,flipped_OR=flip$OR,product=ans$OR*flip$OR,P=ans$P,flipped_P=flip$P,passed=abs(log(ans$OR*flip$OR))<1e-4&&abs(ans$P-flip$P)<1e-5)
     if(g=='GSE212810'&&p$branch=='strict_cisplatin'&&p$endpoint=='good_response_ypT_le1_N0')for(pid in a$patient_id){ll=tryCatch(fit1(a[patient_id!=pid]),error=function(e)list(OR=NA_real_,P=NA_real_,reason=conditionMessage(e)));loo[[length(loo)+1]]=as.data.table(c(list(excluded_patient=pid,n=nrow(a)-1),ll))}
    }
   }
  }
 }
 cat('COMPLETE',g,p$branch,p$endpoint,'\n')
}
wt(rbindlist(res,fill=TRUE),'05_clinical/model_results.tsv');wt(rbindlist(used,fill=TRUE),'05_clinical/model_patient_usage.tsv.gz');wt(rbindlist(dist,fill=TRUE),'05_clinical/response_distributions.tsv');wt(rbindlist(flips,fill=TRUE),'05_clinical/MP10_encoding_flip_check.tsv');wt(rbindlist(loo,fill=TRUE),'05_clinical/GSE212810_MP10_leave_one_patient.tsv')
