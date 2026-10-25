source('results/fig2_like_C0/phase1B_review/scripts/common.R')
P2<-file.path(OUT,'phase2_external');RV<-file.path(P2,'reference_validation');RES<-file.path(P2,'resources')
.libPaths(c(file.path(P2,'Rlib'),.libPaths()))
suppressPackageStartupMessages({library(nnls);library(SingleCellExperiment)})
source(file.path(RES,'MuSiC_analysis.R'));source(file.path(RES,'MuSiC_construct.R'));source(file.path(RES,'MuSiC_utils.R'))
WT<-function(x,name)fwrite(x,file.path(RV,name),sep='\t',na='NA')
basis_aggregate<-function(counts,metadata,types,mincells=20){
 key<-paste(metadata$patient_key,metadata$cell_type,sep='::');lev<-unique(key)
 A<-sparseMatrix(i=seq_along(key),j=match(key,lev),x=1,dims=c(length(key),length(lev)))
 x<-as.matrix(counts%*%A)
 m<-metadata[,.(n_cells=sum(n_cells),patient_key=patient_key[1],cell_type=cell_type[1]),by=.(key=paste(patient_key,cell_type,sep='::'))]
 m<-m[match(lev,key)];ok<-m$n_cells>=mincells&colSums(x)>0;m<-m[ok];x<-x[,ok,drop=FALSE]
 coverage<-m[,.(n_patients=.N),by=cell_type];coverage<-merge(data.table(cell_type=types),coverage,all.x=TRUE);coverage[is.na(n_patients),n_patients:=0L]
 th<-sweep(x,2,colSums(x),'/');M<-V<-matrix(NA_real_,nrow(x),length(types),dimnames=list(rownames(x),types));S<-setNames(rep(NA_real_,length(types)),types)
 for(ct in types){j<-which(m$cell_type==ct);if(!length(j))next
  M[,ct]<-rowMeans(th[,j,drop=FALSE]);S[ct]<-mean(colSums(x[,j,drop=FALSE])/m$n_cells[j])
  if(length(j)>1)V[,ct]<-apply(th[,j,drop=FALSE],1,var)
 }
 list(D=sweep(M,2,S,'*'),M=M,Sigma=V,S=S,coverage=coverage)
}
# Exact sufficient-statistic equivalence check against official music_basis.
set.seed(123);fixture<-matrix(rpois(100*120,lambda=5),100,120,dimnames=list(paste0('g',1:100),paste0('c',1:120)))
fd<-data.table(patient_key=rep(paste0('P',1:4),each=30),cell_type=rep(rep(c('A','B','C'),each=10),4),n_cells=1)
sce<-SingleCellExperiment(list(counts=fixture),colData=as.data.frame(fd))
ref<-music_basis(sce,clusters='cell_type',samples='patient_key',verbose=FALSE)
agg<-basis_aggregate(fixture,fd,c('A','B','C'),1)
checks<-data.table(quantity=c('design','mean_relative_expression','variance','mean_library_size'),max_abs_difference=c(max(abs(ref$Disgn.mtx-agg$D)),max(abs(ref$M.theta-agg$M)),max(abs(ref$Sigma-agg$Sigma)),max(abs(ref$M.S-agg$S))))
stopifnot(all(is.finite(checks$max_abs_difference)),max(checks$max_abs_difference)<1e-10)
WT(checks,'MuSiC_aggregate_equivalence.tsv')
b<-readRDS(file.path(RV,'patient_celltype_pseudobulk.rds'));x<-b$counts;m<-as.data.table(b$metadata);types<-b$types
# One natural mixture per patient is selected before any deconvolution outcome is inspected.
sm<-m[,.(total_cells=sum(n_cells),n_stromal=sum(n_cells[cell_type%in%CLUSTERS]),n_C0=sum(n_cells[cell_type=='C0']),n_C5=sum(n_cells[cell_type=='C5'])),by=.(sample,patient_key,dataset,tissue_origin)]
setorder(sm,patient_key,-total_cells,sample);sm[,primary_patient_specimen:=seq_len(.N)==1,by=patient_key]
WT(sm,'natural_mixture_design.tsv')
platforms<-c('Full_RNA','GPL570','GPL6947','GPL6102','GPL6947_gene_efficiency_SD0.5')
coverage_sets<-lapply(c('GPL570','GPL6947','GPL6102'),function(p)readLines(file.path(P2,'tables',paste0(p,'_covered_genes.txt'))));names(coverage_sets)<-c('GPL570','GPL6947','GPL6102')
set.seed(123);efficiency<-setNames(2^rnorm(nrow(x),0,.5),rownames(x))
predictions<-controls<-audits<-training<-list();counter<-0;cc<-0
for(ds in sort(unique(m$dataset))){
 cat('Validation fold:',ds,'\n');flush.console()
 tr<-which(m$dataset!=ds);te<-which(m$dataset==ds)
 stopifnot(!length(intersect(m$patient_key[tr],m$patient_key[te])))
 bas<-basis_aggregate(x[,tr,drop=FALSE],m[tr],types)
 cv<-copy(bas$coverage);cv[,held_out_dataset:=ds];training[[ds]]<-cv
 if(any(cv$n_patients<3)){
  audits[[ds]]<-data.table(held_out_dataset=ds,scenario='ALL',status='FAIL_TRAINING_COVERAGE',reason=paste(cv[n_patients<3,paste0(cell_type,':',n_patients,' patients')],collapse='; '))
  cat('Fold not estimable under prespecified patient coverage:',audits[[ds]]$reason,'\n');next
 }
 genes<-rownames(x)[apply(bas$M,1,max)>=1e-6 & rowSums(is.finite(bas$Sigma))==length(types)]
 for(platform in platforms){
  g<-genes;if(platform!='Full_RNA')g<-intersect(g,coverage_sets[[sub('_gene_efficiency_SD0.5','',platform)]])
  D<-bas$D[g,,drop=FALSE];V<-bas$Sigma[g,,drop=FALSE]
  audits[[paste(ds,platform)]]<-data.table(held_out_dataset=ds,scenario=platform,status='ESTIMABLE',n_genes=length(g),rank=qr(D)$rank,condition_number=kappa(D),n_training_patients=uniqueN(m$patient_key[tr]),n_test_patients=uniqueN(m$patient_key[te]))
  stopifnot(length(g)>500,qr(D)$rank==length(types))
  fit_one<-function(y){
   y<-y[g];if(grepl('efficiency',platform))y<-y*efficiency[g]
   # Same bulk relative-abundance and nonzero-gene handling as music_prop.
   y<-y/sum(y);nz<-is.finite(y)&y>0
   fit<-music.iter(y[nz],D[nz,,drop=FALSE],bas$S,V[nz,,drop=FALSE],iter.max=1000,nu=.0001,eps=.01,centered=FALSE,normalize=FALSE)
   p<-setNames(fit$p.weight,types);den<-sum(p[CLUSTERS]);stopifnot(all(is.finite(p)),abs(sum(p)-1)<1e-8)
   list(p=p,den=den,C0=if(den>1e-8)p['C0']/den else NA_real_,C5=if(den>1e-8)p['C5']/den else NA_real_,convergence=fit$converge,R2=fit$R.squared)
  }
  for(s in sm[dataset==ds,sample]){
   ind<-which(m$sample==s);y<-Matrix::rowSums(x[,ind,drop=FALSE]);f<-fit_one(y);d<-sm[sample==s]
   counter<-counter+1;predictions[[counter]]<-cbind(d,data.table(scenario=platform,true_C0=if(d$n_stromal>0)d$n_C0/d$n_stromal else NA_real_,true_C5=if(d$n_stromal>0)d$n_C5/d$n_stromal else NA_real_,pred_C0=f$C0,pred_C5=f$C5,pred_stromal=f$den,convergence=f$convergence,R2=f$R2),as.data.table(as.list(setNames(f$p,paste0('estimate_',names(f$p))))))
  }
  # Held-out patient/type diagnostic pure profiles; no training cells reused.
  if(platform=='Full_RNA'){
   for(j in te[m$n_cells[te]>=20]){
    f<-fit_one(as.numeric(x[,j]) |> setNames(rownames(x)))
    cc<-cc+1;controls[[cc]]<-cbind(m[j],data.table(pred_C0=f$C0,pred_C5=f$C5,pred_stromal=f$den),as.data.table(as.list(setNames(f$p,paste0('estimate_',names(f$p))))))
   }
  }
 }
 WT(rbindlist(predictions,fill=TRUE),'natural_mixture_predictions.tsv')
 WT(rbindlist(controls,fill=TRUE),'held_out_pure_type_predictions.tsv')
}
WT(rbindlist(training,fill=TRUE),'training_patient_coverage.tsv');WT(rbindlist(audits,fill=TRUE),'validation_fold_audit.tsv')
p<-rbindlist(predictions,fill=TRUE);stopifnot(nrow(p)>0)
metric<-function(d){d<-d[is.finite(true_C0)&is.finite(pred_C0)];e<-d$pred_C0-d$true_C0;data.table(n_patients=nrow(d),pearson=if(nrow(d)>2)cor(d$true_C0,d$pred_C0)else NA_real_,spearman=if(nrow(d)>2)cor(d$true_C0,d$pred_C0,method='spearman')else NA_real_,MAE=mean(abs(e)),RMSE=sqrt(mean(e^2)),bias=mean(e))}
primary<-p[primary_patient_specimen==TRUE&n_stromal>=50]
metrics<-rbindlist(list(primary[,cbind(data.table(dataset='OVERALL'),metric(.SD)),by=scenario],primary[,metric(.SD),by=.(scenario,dataset)]),fill=TRUE)
WT(metrics,'C0_validation_metrics.tsv')
allmetrics<-p[primary_patient_specimen==TRUE&n_stromal>0,metric(.SD),by=scenario];WT(allmetrics,'C0_validation_metrics_any_stromal.tsv')
ctrl<-rbindlist(controls,fill=TRUE);cs<-ctrl[cell_type%in%c('C0','C5'),.(n_profiles=.N,median_pred_C0=median(pred_C0,na.rm=TRUE),median_pred_C5=median(pred_C5,na.rm=TRUE),median_pred_stromal=median(pred_stromal)),by=cell_type];WT(cs,'C0_C5_pure_profile_confusion.tsv')
# Gate is fixed in validation_plan.md; no downstream clinical fits if any critical criterion fails.
full<-metrics[dataset=='OVERALL'&scenario=='Full_RNA'];per<-metrics[dataset!='OVERALL'&scenario=='Full_RNA'&n_patients>=5]
gate<-data.table(criterion=c('complete_dataset_validation','overall_spearman_ge_0.70','overall_MAE_le_0.15','overall_RMSE_le_0.20','overall_abs_bias_le_0.10','dataset_MAE_le_0.20','C0_C5_identity_discrimination'),pass=c(all(rbindlist(audits,fill=TRUE)$status=='ESTIMABLE'),full$spearman>=.7,full$MAE<=.15,full$RMSE<=.20,abs(full$bias)<=.1,all(per$MAE<=.2),all(cs[cell_type=='C0',median_pred_C0]>=.7)&all(cs[cell_type=='C5',median_pred_C5]>=.7)&all(cs[cell_type=='C0',median_pred_C5]<=.2)&all(cs[cell_type=='C5',median_pred_C0]<=.2)))
WT(gate,'validation_gate.tsv');writeLines(if(all(gate$pass))'SCREEN_PASS_REQUIRES_PLATFORM_REVIEW'else'STOP_ABUNDANCE_CLINICAL_ANALYSIS',file.path(RV,'VALIDATION_STATUS.txt'))
# Small final reference statistics are retained for review, not used clinically unless validated.
final<-basis_aggregate(x,m,types);keep<-apply(final$M,1,max)>=1e-6
saveRDS(list(D=final$D[keep,],Sigma=final$Sigma[keep,],S=final$S,types=types,coverage=final$coverage),file.path(RV,'MuSiC_full_reference_statistics.rds'),compress='gzip')
capture.output(sessionInfo(),file=file.path(RV,'sessionInfo.txt'));print(metrics);print(gate)
