source('results/llps_C5_parallel_v1/scripts/common.R')
D<-file.path(OUT,'03_Bulk_Survival');WEB<-file.path(BASE,'fig2_C0_multimethod_deconvolution/03_methods/CIBERSORTx_web')
res<-list()
for(co in c('TCGA_BLCA','UROMOL')){
 job<-if(co=='TCGA_BLCA')14 else 13
 f<-fread(file.path(WEB,'output',paste0('CIBERSORTx_Job',job,'_Adjusted.txt')))
 full<-fread(file.path(D,'by_dataset',co,'clinical_merged.tsv'))
 log<-readLines(file.path(WEB,'output',paste0('CIBERSORTx_Job',job,'_runtime_log.txt')))
 stopifnot(any(grepl('Reference_C0focused',log)),!('Fibro_C5'%in%names(f)),!anyDuplicated(f$Mixture),setequal(f$Mixture,full$CIBERSORTx_sample_id))
 f<-f[match(full$CIBERSORTx_sample_id,Mixture)]
 z<-data.table(sample=full$CIBERSORTx_sample_id,Full14_C0=full$Fibro_C0,focused_C0=f$Fibro_C0)
 wt(z,file.path(D,paste0(co,'_Full14_focused_C0.tsv')))
 res[[co]]<-data.table(dataset=co,N=nrow(z),Spearman=cor(z$Full14_C0,z$focused_C0,method='spearman'),Pearson=cor(z$Full14_C0,z$focused_C0),mean_absolute_difference=mean(abs(z$Full14_C0-z$focused_C0)),C5_comparison='Not possible: focused model has no independent C5',note='New descriptive comparison of existing fractions; not deconvolution accuracy validation')
}
wt(rbindlist(res),file.path(D,'Full14_focused_consistency.tsv'))
