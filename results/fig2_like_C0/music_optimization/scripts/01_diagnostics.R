source('results/fig2_like_C0/music_optimization/scripts/common.R');CORE<-file.path(OPT,'01_diagnostics')
set.seed(123)
stat_ci<-function(d){a<-metric(d$true_C0,d$pred_C0);bs<-replicate(2000,{j<-sample(seq_len(nrow(d)),replace=TRUE);metric(d$true_C0[j],d$pred_C0[j])[-1]})
 z<-as.list(a);for(k in rownames(bs)){v<-bs[k,];v<-v[is.finite(v)];z[[paste0(k,'_low')]]<-if(length(v))quantile(v,.025)else NA_real_;z[[paste0(k,'_high')]]<-if(length(v))quantile(v,.975)else NA_real_};as.data.table(z)}
d<-copy(TEST);d[,error:=pred_C0-true_C0]
WT(rbindlist(list(cbind(dataset='OVERALL',stat_ci(d)),d[,stat_ci(.SD),by=dataset]),fill=TRUE),'dataset_metrics_bootstrap.tsv')
rt<-rank(d$true_C0);rp<-rank(d$pred_C0);adj<-cor(resid(lm(rt~factor(d$dataset))),resid(lm(rp~factor(d$dataset))))
perm<-replicate(5000,{q<-d[,.(tp=sample(pred_C0)),by=dataset];cor_safe(d[order(dataset),true_C0],q$tp)})
obs<-cor_safe(d$true_C0,d$pred_C0)
WT(data.table(overall_rho=obs,rank_partial_rho_dataset=adj,between_dataset_R2_true=summary(lm(rt~factor(d$dataset)))$r.squared,between_dataset_R2_pred=summary(lm(rp~factor(d$dataset)))$r.squared,within_dataset_permutation_null_median=median(perm),permutation_p=(1+sum(perm>=obs))/(1+length(perm))),'dataset_rank_decomposition.tsv')
WT(d[,.(sample,patient_key,dataset,true_C0,pred_C0,error,n_stromal)],'original_20_patient_predictions.tsv')
p<-ggplot(d,aes(true_C0,pred_C0,color=dataset))+geom_abline(slope=1,intercept=0,linetype='dashed',color='grey50')+geom_point(size=2.8)+scale_color_manual(values=DSCOL)+coord_equal(xlim=c(0,1),ylim=c(0,1))+theme_pub()+theme(legend.position='bottom')+labs(title='C0 recovery by dataset',subtitle=sprintf('n=20 | rho=%.3f | MAE=%.3f | RMSE=%.3f',obs,mean(abs(d$error)),sqrt(mean(d$error^2))),x='True C0 fraction',y='Predicted C0 fraction',color=NULL)
export_plot(p,'FigD1_C0_true_pred_by_dataset',7.3,6)
p<-ggplot(d,aes(dataset,error,fill=dataset))+geom_hline(yintercept=0,linetype='dashed',color='grey50')+geom_boxplot(outlier.shape=NA,width=.55)+geom_point(position=position_jitter(width=.08,seed=123),shape=21,size=2)+scale_fill_manual(values=DSCOL)+theme_pub()+theme(legend.position='none',axis.text.x=element_text(angle=20,hjust=1))+labs(title='C0 prediction error within each dataset',x=NULL,y='Predicted minus true fraction',caption='All 20 original evaluation patients retained.')
export_plot(p,'FigD2_C0_error_by_dataset',7.2,4.7)
truth<-M[cell_type%in%CLUSTERS,.(n_cells=sum(n_cells),umi=sum(umi)),by=.(sample,cell_type)]
truth<-merge(CJ(sample=d$sample,cell_type=CLUSTERS),truth,all.x=TRUE);truth[is.na(n_cells),`:=`(n_cells=0,umi=0)];truth[,`:=`(true_fraction=n_cells/sum(n_cells),true_RNA_fraction=umi/sum(umi)),by=sample]
pl<-melt(d,id.vars=c('sample','dataset','patient_key','pred_stromal'),measure.vars=paste0('estimate_',CLUSTERS),variable.name='cell_type',value.name='estimate');pl[,cell_type:=sub('estimate_','',cell_type)];pl[,pred_fraction:=estimate/pred_stromal]
z<-merge(truth,pl,by=c('sample','cell_type'));z[,error:=pred_fraction-true_fraction]
WT(z,'C0_C6_cell_and_RNA_truth_predictions.tsv');WT(z[,as.data.table(as.list(metric(true_fraction,pred_fraction))),by=cell_type],'C0_C6_error_metrics.tsv')
p<-ggplot(z,aes(factor(cell_type,levels=CLUSTERS),error,fill=cell_type))+geom_hline(yintercept=0,linetype='dashed',color='grey50')+geom_boxplot(outlier.shape=NA,width=.6)+geom_point(position=position_jitter(width=.12,seed=123),size=1.2,alpha=.7)+scale_fill_manual(values=COLORS)+theme_pub()+theme(legend.position='none')+labs(title='Prediction errors across stromal clusters',x=NULL,y='Predicted minus true within-stromal fraction')
export_plot(p,'FigD3_C0_C6_errors',7.5,4.7)
w<-dcast(z,sample~cell_type,value.var='error');ec<-cor(as.matrix(w[,-1]),method='spearman');WT(as.data.table(ec,keep.rownames='cluster'),'C0_C6_error_correlations.tsv')
ref<-readRDS(file.path(RV,'MuSiC_full_reference_statistics.rds'));v<-log2(sweep(ref$D[,CLUSTERS],2,ref$S[CLUSTERS],'/')*1e6+1);sim<-cor(v)
WT(as.data.table(sim,keep.rownames='cluster'),'reference_expression_similarity.tsv')
heat<-function(mat,title,stem){l<-as.data.table(as.table(mat));setnames(l,c('row','column','value'));p<-ggplot(l,aes(column,row,fill=value))+geom_tile(color='white')+geom_text(aes(label=sprintf('%.2f',value)),size=3)+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-1,1),name='Correlation')+theme_pub()+theme(axis.line=element_blank(),axis.ticks=element_blank())+labs(title=title,x=NULL,y=NULL);export_plot(p,stem,6,5.5)}
heat(sim,'Reference expression similarity (log CPM)','FigD4_reference_similarity');heat(ec,'Error correlation under a compositional constraint','FigD5_error_correlation')
pat<-M[,.(n_cells=sum(n_cells),umi=sum(umi)),by=.(dataset,patient_key,cell_type)]
WT(pat[cell_type%in%CLUSTERS,.(cells=sum(n_cells),patients_present=.N,patients_ge5=sum(n_cells>=5),patients_ge20=sum(n_cells>=20),median_UMI_per_cell=median(umi/n_cells)),by=.(dataset,cell_type)],'patient_coverage_and_RNA_content.tsv')
# Reconstruct baseline fold size factors, then compare cell vs RNA contribution explicitly.
rna<-rbindlist(lapply(unique(d$dataset),function(ds){bas<-basis(X[,M$dataset!=ds],M[dataset!=ds],20,FALSE);q<-copy(z[dataset==ds]);q[,estimated_RNA:=estimate*bas$S[cell_type]];q[,estimated_RNA_fraction:=estimated_RNA/sum(estimated_RNA),by=sample];q}))
WT(rna,'RNA_content_diagnostic.tsv')
WT(rna[cell_type=='C0',.(rho_cell=cor_safe(true_fraction,pred_fraction),rho_RNA=cor_safe(true_RNA_fraction,estimated_RNA_fraction),MAE_cell=mean(abs(pred_fraction-true_fraction)),MAE_RNA=mean(abs(estimated_RNA_fraction-true_RNA_fraction)),error_vs_RNA_cell_difference=cor_safe(error,true_RNA_fraction-true_fraction)),by=dataset],'RNA_vs_cell_error_summary.tsv')
# Dataset-specific reference centroid correlations are descriptive, not a new cell-level test.
cent<-list();for(ds in unique(M$dataset)){ba<-basis(X[,M$dataset==ds],M[dataset==ds],1,FALSE);cent[[ds]]<-ba$M}
cs<-rbindlist(lapply(CLUSTERS,function(ct){rbindlist(lapply(combn(names(cent),2,simplify=FALSE),function(pair){if(!ct%in%colnames(cent[[pair[1]]])||!ct%in%colnames(cent[[pair[2]]]))return(NULL);data.table(cluster=ct,dataset1=pair[1],dataset2=pair[2],cor_logCPM=cor_safe(log2(cent[[pair[1]]][,ct]*1e6+1),log2(cent[[pair[2]]][,ct]*1e6+1),'pearson'))}))}))
WT(cs,'cluster_dataset_expression_similarity.tsv');cat('Diagnostics complete\n')
