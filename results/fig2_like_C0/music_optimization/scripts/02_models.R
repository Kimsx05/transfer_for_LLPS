source('results/fig2_like_C0/music_optimization/scripts/common.R');CORE<-file.path(OPT,'02_model_comparison')
hg<-fread(file.path(FIBROOT,'beautified_SCP_protein_coding_top3/hgnc_complete_set.txt'))
uni<-Reduce(intersect,lapply(c('GPL570','GPL6947','GPL6102'),function(g)readLines(file.path(P2,'tables',paste0(g,'_covered_genes.txt')))))
uni<-intersect(intersect(uni,hg[locus_group=='protein-coding gene',symbol]),rownames(X));uni<-uni[!grepl('^(MT-|RPL|RPS)',uni)]
writeLines(sort(uni),file.path(CORE,'fixed_shared_gene_universe.txt'))
select_genes<-function(ba){
 g<-uni[apply(ba$M[uni,,drop=FALSE],1,max)>=5e-6]
 L<-log2(ba$M[g,,drop=FALSE]*1e6+1)
 # Specificity discounted by training-patient variability; no held-out profiles are consulted.
 scores<-sapply(colnames(L),function(ct){j<-which(ba$metadata$cell_type==ct);v<-apply(log2(ba$theta[g,j,drop=FALSE]*1e6+1),1,var);(L[,ct]-rowMeans(L[,setdiff(colnames(L),ct),drop=FALSE]))/(1+sqrt(pmax(v,0)))})
 unique(unlist(lapply(seq_len(ncol(scores)),function(j)g[head(order(scores[,j],decreasing=TRUE),100)])))
}
fit_hierarchy<-function(y,fine,coarse,genes){
 f1<-fit_music(y,coarse,genes);g<-f1$genes;nc<-setdiff(colnames(coarse$D),'Stroma')
 # q.weight multiplies D to fit 100*relative bulk expression in music.iter.
 nuisance<-as.vector(coarse$D[g,nc,drop=FALSE]%*%f1$q.weight[match(nc,colnames(coarse$D))])
 target<-100*f1$Y-nuisance;neg<-mean(target<0);target<-pmax(target,0);names(target)<-g
 keep<-target>0;sf<-list(D=fine$D[,CLUSTERS],S=fine$S[CLUSTERS],Sigma=fine$Sigma[,CLUSTERS])
 f2<-music.basic(target[keep],sf$D[g[keep],,drop=FALSE],sf$S,sf$Sigma[g[keep],,drop=FALSE],iter.max=1000,nu=.0001,eps=.01)
 p<-f1$p[nc];p<-c(setNames(f2$p.weight*unname(f1$p['Stroma']),CLUSTERS),p)
 list(p=p,negative_residual_fraction=neg,convergence=paste(f1$converge,f2$converge,sep='; '))
}
preds<-audit<-coverage<-genesets<-list();cnt<-0
sm<-fread(file.path(RV,'natural_mixture_design.tsv'))
# Reserve original 20 globally; additional available patients are also evaluated by leave-dataset-out.
EVAL<-sm[primary_patient_specimen==TRUE&n_stromal>=50];EVAL[,evaluation_set:=ifelse(patient_key%in%TESTIDS,'original20_development_evaluation','additional_cross_dataset_evaluation')]
WT(data.table(patient_key=unique(M$patient_key),reserved_original20=unique(M$patient_key)%in%TESTIDS),'patient_split_manifest.tsv')
for(ds in sort(unique(M$dataset))){
 cat('Dataset:',ds,'\n');flush.console()
 tr<-which(!M$patient_key%in%TESTIDS&M$dataset!=ds);te<-EVAL[dataset==ds]
 stopifnot(!length(intersect(M$patient_key[tr],TESTIDS)),!length(intersect(M$patient_key[tr],te$patient_key)))
 for(mode in c('Baseline_refit','Optimized_reference','Hierarchical_custom')){
  cat(mode,'\n');flush.console()
  opt<-mode!='Baseline_refit';ba<-basis(X[,tr,drop=FALSE],M[tr],if(opt)5 else 20,opt)
  cv<-merge(data.table(cell_type=TYPES),ba$coverage,all.x=TRUE);cv[is.na(n_patients),n_patients:=0];cv[,`:=`(method=mode,held_out_dataset=ds)];coverage[[paste(ds,mode)]]<-cv
  reason<-if(any(cv$n_patients<3))paste(cv[n_patients<3,paste0(cell_type,':',n_patients)],collapse=';')else''
  if(nzchar(reason)){
   audit[[paste(ds,mode)]]<-data.table(method=mode,dataset=ds,status='UNSUPPORTED_COVERAGE',reason=reason)
   for(i in seq_len(nrow(te))){cnt<-cnt+1;preds[[cnt]]<-cbind(te[i],data.table(method=mode,status='UNSUPPORTED_COVERAGE',pred_C0=NA_real_))};next
  }
  genes<-if(opt)select_genes(ba)else rownames(X)[apply(ba$M,1,max)>=1e-6]
  genesets[[paste(ds,mode)]]<-data.table(dataset=ds,method=mode,gene=genes)
  audit[[paste(ds,mode)]]<-data.table(method=mode,dataset=ds,status='ESTIMABLE',n_genes=length(genes),condition_number=kappa(ba$D[genes,]),n_training_patients=uniqueN(M$patient_key[tr]))
  if(mode=='Hierarchical_custom'){
   mc<-copy(M[tr]);mc[cell_type%in%CLUSTERS,cell_type:='Stroma'];coarse<-basis(X[,tr,drop=FALSE],mc,5,TRUE)
  }
  for(i in seq_len(nrow(te))){s<-te$sample[i];y<-Matrix::rowSums(X[,M$sample==s,drop=FALSE]);j<-which(M$sample==s)
   fit<-tryCatch(if(mode=='Hierarchical_custom')fit_hierarchy(y,ba,coarse,genes)else fit_music(y,ba,genes),error=function(e)e)
   cnt<-cnt+1
   if(inherits(fit,'error')){preds[[cnt]]<-cbind(te[i],data.table(method=mode,status=paste('FIT_ERROR',conditionMessage(fit)),pred_C0=NA_real_));next}
   pp<-fit$p[TYPES];den<-sum(pp[CLUSTERS]);pv<-if(den>1e-8)pp['C0']/den else NA_real_
   preds[[cnt]]<-cbind(te[i],data.table(method=mode,status=if(is.finite(pv))'OK'else'ZERO_ESTIMATED_STROMA',pred_C0=unname(pv),true_C0=te$n_C0[i]/te$n_stromal[i],pred_stromal=den,negative_residual_fraction=if(mode=='Hierarchical_custom')fit$negative_residual_fraction else NA_real_),as.data.table(as.list(setNames(pp,paste0('estimate_',TYPES)))))
  }
 }
 WT(rbindlist(preds,fill=TRUE),'optimized_model_predictions.tsv')
}
WT(rbindlist(audit,fill=TRUE),'model_fold_audit.tsv');WT(rbindlist(coverage,fill=TRUE),'model_training_coverage.tsv');WT(rbindlist(genesets,fill=TRUE),'training_selected_genes.tsv')
p<-rbindlist(preds,fill=TRUE);hist<-copy(TEST);hist[,`:=`(method='Historical_baseline',evaluation_set='original20_development_evaluation',status='OK')]
q<-rbindlist(list(p,hist),fill=TRUE);WT(q,'all_model_predictions.tsv')
summary<-q[,cbind(n_expected=.N,as.data.table(as.list(metric(true_C0,pred_C0)))),by=.(method,evaluation_set)]
summary[,coverage:=n/n_expected];WT(summary,'model_comparison.tsv')
WT(q[,as.data.table(as.list(metric(true_C0,pred_C0))),by=.(method,evaluation_set,dataset)],'model_comparison_by_dataset.tsv')
# Common subset is reported separately; never replace coverage with this selective denominator.
original<-q[evaluation_set=='original20_development_evaluation'];common<-original[status=='OK',.(nm=uniqueN(method)),by=sample][nm==4,sample]
WT(original[sample%in%common,as.data.table(as.list(metric(true_C0,pred_C0))),by=method],'common_evaluable_patient_comparison.tsv')
p1<-ggplot(original[is.finite(pred_C0)],aes(method,abs(pred_C0-true_C0),color=method))+geom_boxplot(outlier.shape=NA,width=.5)+geom_point(position=position_jitter(width=.1,seed=123),size=1.5)+theme_pub()+theme(legend.position='none',axis.text.x=element_text(angle=15,hjust=1))+labs(title='C0 estimation: bounded model comparison',x=NULL,y='Absolute error in C0 fraction',caption='Original 20 reserved patients; failures remain in the coverage table.\nHistorical baseline used different training pools; all new fits globally exclude the original 20.')
export_plot(p1,'FigM1_model_error_comparison',8.6,5.2)
print(summary)
