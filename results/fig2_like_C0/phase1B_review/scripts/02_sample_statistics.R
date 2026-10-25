source('results/fig2_like_C0/phase1B_review/scripts/common.R')
m<-fread(file.path(CORE,'sample_metadata_reviewed.tsv'));m[,tissue_origin:=factor(tissue_origin,levels=c('Normal','Tumor'))]
z<-m[analysis_included==TRUE]
WT(m,'C0_fraction_per_sample_reviewed.tsv')
coverage<-m[,.(n_source=.N,n_with_stroma=sum(n_all_Fibro_Pericyte>0),n_included=sum(analysis_included)),by=.(dataset,tissue_origin)]
WT(coverage,'dataset_tissue_sample_counts.tsv')

dataset_summary<-function(x,value,prefix) {
 out<-rbindlist(lapply(sort(unique(x$dataset)),function(ds){
  d<-x[dataset==ds];p<-paired(d,value);nN<-sum(d$tissue_origin=='Normal');nT<-sum(d$tissue_origin=='Tumor')
  method<-'DESCRIPTIVE_ONLY_INSUFFICIENT_INDEPENDENT_REPLICATES';pv<-NA_real_;lo<-hi<-effect<-NA_real_
  if(nN&&nT)effect<-mean(d[tissue_origin=='Tumor',get(value)])-mean(d[tissue_origin=='Normal',get(value)])
  if(nrow(p)>=3){
   delta<-p$Tumor-p$Normal;effect<-mean(delta);boot<-replicate(5000,mean(sample(delta,replace=TRUE)))
   ci<-quantile(boot,c(.025,.975));lo<-ci[1];hi<-ci[2]
   pv<-if(all(delta==0))1 else wilcox.test(p$Tumor,p$Normal,paired=TRUE,exact=FALSE)$p.value
   method<-'paired Wilcoxon; effect=mean paired difference; patient bootstrap CI'
  }else if(nrow(p)==0&&min(nN,nT)>=3){
   nn<-d[tissue_origin=='Normal',get(value)];tt<-d[tissue_origin=='Tumor',get(value)]
   boot<-replicate(5000,mean(sample(tt,replace=TRUE))-mean(sample(nn,replace=TRUE)));ci<-quantile(boot,c(.025,.975));lo<-ci[1];hi<-ci[2]
   pv<-wilcox.test(tt,nn,exact=FALSE)$p.value;method<-'unpaired Wilcoxon; mean difference; sample bootstrap CI'
  }
  data.table(dataset=ds,n_Normal=nN,n_Tumor=nT,n_verified_pairs=nrow(p),mean_TminusN=effect,
   lower=lo,upper=hi,p_value=pv,test=method,
   median_Normal=if(nN)median(d[tissue_origin=='Normal',get(value)])else NA_real_,
   median_Tumor=if(nT)median(d[tissue_origin=='Tumor',get(value)])else NA_real_)
 }));out[,FDR:=p.adjust(p_value,'BH')];WT(out,paste0(prefix,'_dataset_stats.tsv'));out
}
ds<-dataset_summary(z,'C0_proportion','C0_fraction')
plots_dataset<-function(x,value,title,ylab,filename,percent=FALSE) {
 x<-copy(x);x[,dataset:=factor(dataset,levels=sort(unique(dataset)))];x[,tissue_origin:=factor(tissue_origin,levels=c('Normal','Tumor'))]
 cnt<-x[,.N,by=.(dataset,tissue_origin)]
 p<-ggplot(x,aes(dataset,get(value),fill=tissue_origin))+
  geom_boxplot(position=position_dodge(.7),width=.55,outlier.shape=NA,linewidth=.4)+
  geom_point(position=position_jitterdodge(jitter.width=.12,dodge.width=.7,seed=123),shape=21,size=1.8,stroke=.25)+
  scale_fill_manual(values=GROUP_COLORS,name=NULL)+theme_pub()+
  theme(axis.text.x=element_text(angle=20,hjust=1),legend.position='top')+
  labs(title=title,x=NULL,y=ylab,caption=paste(vapply(levels(x$dataset),function(k){c<-cnt[dataset==k];paste0(k,': N=',sum(c[tissue_origin=='Normal',N]),', T=',sum(c[tissue_origin=='Tumor',N]))},character(1)),collapse=' | '))
 if(percent)p<-p+scale_y_continuous(labels=scales::percent)
 export_plot(p,filename,9,5.4)
}
plots_dataset(z,'C0_proportion','C0 proportion within each dataset','C0 / all Fibro + Pericyte','FigS_C0_fraction_dataset',TRUE)
ff<-ds[is.finite(mean_TminusN)];ff[,dataset:=factor(dataset,levels=rev(dataset))]
p<-ggplot(ff,aes(mean_TminusN,dataset))+geom_vline(xintercept=0,linetype='dashed',color='grey65')+
 geom_segment(aes(x=lower,xend=upper,yend=dataset),linewidth=.7,na.rm=TRUE)+geom_point(size=3,color=COLORS['C0'])+
 scale_x_continuous(labels=scales::percent)+theme_pub()+labs(title='Within-dataset change in C0 proportion',
 x='Tumor minus Normal (percentage points)',y=NULL,caption='Paired mean difference and patient-bootstrap 95% CI when >=3 pairs.\nGSE135337: descriptive difference only (one Normal); no CI or test.')
export_plot(p,'FigS_C0_fraction_forest',7.6,4.5)
pa<-paired(z,'C0_proportion');delta<-pa$Tumor-pa$Normal
pst<-data.table(n_pairs=nrow(pa),n_nonzero=sum(delta!=0),n_decreased=sum(delta<0),n_increased=sum(delta>0),
 median_paired_difference=median(delta),p_value=wilcox.test(pa$Tumor,pa$Normal,paired=TRUE,exact=FALSE)$p.value)
WT(pa,'C0_fraction_verified_pairs_updated.tsv');WT(pst,'C0_fraction_paired_stats_updated.tsv')
pl<-melt(pa,id.vars=c('patient_uid','dataset'),measure.vars=c('Normal','Tumor'),variable.name='tissue_origin',value.name='value')
p<-ggplot(pl,aes(tissue_origin,value,group=patient_uid))+geom_line(color='grey65',linewidth=.4)+geom_point(aes(fill=tissue_origin),shape=21,size=2.4)+
 scale_fill_manual(values=GROUP_COLORS)+scale_y_continuous(labels=scales::percent)+theme_pub()+theme(legend.position='none')+
 labs(title='C0 proportion: verified matched specimens',subtitle=sprintf('%d pairs | Paired Wilcoxon P = %s',nrow(pa),ptext(pst$p_value)),
 x=NULL,y='C0 / all Fibro + Pericyte',caption='Includes five newly verified HRA003620 pairs with nonzero denominators.')
export_plot(p,'FigS_C0_fraction_paired_updated',6.5,5)
models<-list(cluster_lm(z,C0_proportion~factor(dataset)+tumor,label='dataset_adjusted'))
for(ds0 in unique(z$dataset))models[[ds0]]<-cluster_lm(z[dataset!=ds0],C0_proportion~factor(dataset)+tumor,label=paste0('leave_out_',ds0))
models[['min50']]<-cluster_lm(z[n_all_Fibro_Pericyte>=50],C0_proportion~factor(dataset)+tumor,label='dataset_adjusted_min50')
WT(rbindlist(models,fill=TRUE),'C0_fraction_adjusted_and_leave_dataset_out.tsv')

# C0 state and C0 proportion are different outcomes. Predeclared median-score
# primary threshold is >=10 C0 cells; >=1 and >=20 are capture sensitivities.
q<-fread(file.path(CORE,'C0_LLPS_per_specimen.tsv'));q[,tissue_origin:=factor(tissue_origin,levels=c('Normal','Tumor'))]
l<-q[LLPS_primary==TRUE]
mods<-list(cluster_lm(l,median_LLPS~factor(dataset)+tumor,label='primary_min10_dataset_adjusted'),
 cluster_lm(l,median_LLPS~factor(dataset)+tumor+log10_features,label='plus_median_log10_nFeature'),
 cluster_lm(l,median_LLPS~factor(dataset)+tumor+log10_counts,label='plus_median_log10_nCount'),
 cluster_lm(q[LLPS_anycapture==TRUE],median_LLPS~factor(dataset)+tumor,label='any_C0_capture_sensitivity'),
 cluster_lm(q[primary_replication_allowed&n_C0>=20],median_LLPS~factor(dataset)+tumor,label='min20_C0_sensitivity'))
for(ds0 in unique(l$dataset))mods[[paste0('omit_',ds0)]]<-cluster_lm(l[dataset!=ds0],median_LLPS~factor(dataset)+tumor,label=paste0('leave_out_',ds0))
WT(rbindlist(mods,fill=TRUE),'C0_LLPS_models.tsv')
tech<-rbindlist(lapply(c('median_nFeature','median_nCount'),function(v){t<-cor.test(l$median_LLPS,l[[v]],method='spearman',exact=FALSE);data.table(covariate=v,n=nrow(l),rho=unname(t$estimate),p_value=t$p.value)}))
WT(tech,'C0_LLPS_technical_correlations.tsv')
p<-ggplot(l,aes(tissue_origin,median_LLPS,fill=tissue_origin))+geom_boxplot(width=.45,outlier.shape=NA,linewidth=.4)+
 geom_point(position=position_jitter(width=.10,height=0,seed=123),shape=21,size=2.2,stroke=.3)+scale_fill_manual(values=GROUP_COLORS)+
 theme_pub()+theme(legend.position='none')+labs(title='LLPS state within C0',
 subtitle=sprintf('Normal n = %d; Tumor n = %d | Adjusted P = %s',sum(l$tissue_origin=='Normal'),sum(l$tissue_origin=='Tumor'),ptext(mods[[1]]$p_value)),
 x=NULL,y='Specimen median LLPS77 AUCell',caption='>=10 C0 cells per specimen. Dataset-adjusted model; patient-cluster SE.')
export_plot(p,'FigS1_C0_LLPS_TumorNormal',6.6,4.9)
lp<-paired(l,'median_LLPS');WT(lp,'C0_LLPS_verified_pairs.tsv')
if(nrow(lp)>=3){
 pw<-wilcox.test(lp$Tumor,lp$Normal,paired=TRUE,exact=FALSE)
 WT(data.table(n_pairs=nrow(lp),n_decreased=sum(lp$Tumor<lp$Normal),median_paired_difference=median(lp$Tumor-lp$Normal),p_value=pw$p.value),'C0_LLPS_paired_stats.tsv')
 lp<-melt(lp,id.vars=c('patient_uid','dataset'),measure.vars=c('Normal','Tumor'),variable.name='tissue_origin',value.name='value')
 p<-ggplot(lp,aes(tissue_origin,value,group=patient_uid))+geom_line(color='grey65',linewidth=.45)+geom_point(aes(fill=tissue_origin),shape=21,size=2.6)+
 scale_fill_manual(values=GROUP_COLORS)+theme_pub()+theme(legend.position='none')+labs(title='C0 LLPS state in paired specimens',
 subtitle=sprintf('%d eligible pairs | Paired Wilcoxon P = %s',uniqueN(lp$patient_uid),ptext(pw$p.value)),x=NULL,y='Specimen median LLPS77 AUCell',caption='>=10 C0 cells in each specimen; verified patient pairs only.')
 export_plot(p,'FigS2_C0_LLPS_paired',6,5)
}
plots_dataset(l,'median_LLPS','C0 LLPS state across datasets','Specimen median LLPS77 AUCell','FigS3_C0_LLPS_dataset')
dataset_summary(l,'median_LLPS','C0_LLPS')
WT(l[,.(n=.N,median=median(median_LLPS),Q1=quantile(median_LLPS,.25),Q3=quantile(median_LLPS,.75)),by=tissue_origin],'C0_LLPS_group_summary.tsv')
cat('Sample-level review complete\n')
