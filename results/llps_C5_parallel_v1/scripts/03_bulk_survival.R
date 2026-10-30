# Thin driver around original C0 fit_cox, choose_covars, KM and forest plotting blocks.
BASE<-'/home/jinshengxi/bladder_metabolism';OUT<-file.path(BASE,'results/llps_C5_parallel_v1/03_Bulk_Survival')
.libPaths(c(file.path(BASE,'cellchat_nichechat_fibro7/Rlib'),file.path(BASE,'fig2_C0_multimethod_deconvolution/Rlib'),file.path(BASE,'results/fig2_like_C0/phase2_external/Rlib'),'/home/jinshengxi/Rlibrary/packages',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(survival);library(patchwork)})
set.seed(123);setDTthreads(1)
inv<-fread(file.path(OUT,'input_inventory.tsv'),na.strings=c('NA',''))
src<-paste(readLines('results/cibersortx_phase2_clinical_v1/scripts/02_analyze.R'),collapse='\n')
a<-regexpr('truthy<-function',src,fixed=TRUE)[1];b<-regexpr('all_summary<-',src,fixed=TRUE)[1];funcs<-substr(src,a,b-1)
a<-regexpr(' uni<-multi<-phs',src,fixed=TRUE)[1];b<-regexpr(' all_cox[[co]]',src,fixed=TRUE)[1];body0<-substr(src,a,b-1)
body0<-sub("eps<-setdiff(eps,'NONE')","eps<-intersect(eps,'OS')",body0,fixed=TRUE)
allcx<-allkm<-list()
for(target in c('C0','C5')){
 eval(parse(text=gsub('C0',target,funcs,fixed=TRUE)))
 # Original savepdf plus companion transparent 300 dpi PNG.
 savepdf<-function(p,path,w=7,h=5){if(inherits(p,'patchwork'))p<-p&theme(plot.background=element_rect(fill='transparent',colour=NA));ggsave(path,p,width=w,height=h,device=cairo_pdf,bg='transparent',limitsize=FALSE);ggsave(sub('[.]pdf$','.png',path),p,width=w,height=h,dpi=300,bg='transparent',limitsize=FALSE)}
 for(co in inv[grepl('(^|;)OS(;|$)',endpoint_list),dataset]){
  root<-file.path(OUT,'by_dataset',co,target);fig<-file.path(root,'figures');dir.create(fig,recursive=TRUE,showWarnings=FALSE)
  d<-fread(file.path(OUT,'by_dataset',co,'clinical_merged.tsv'),na.strings=c('NA','','#N/A'))
  val<-paste0('Fibro_',target);zval<-paste0('z_',target)
  d[,(zval):=(get(val)-mean(get(val)))/sd(get(val))];cohort_sd<-sd(d[[val]]);cutoff<-median(d[[val]])
  stopifnot(!anyDuplicated(d$analysis_unit_id),all(is.finite(d[[val]])),cohort_sd>0)
  eval(parse(text=gsub('C0',target,body0,fixed=TRUE)))
  cx<-rbind(u,m,fill=TRUE);cx[,`:=`(dataset=co,cluster=target)]
  # Generic columns for direct same-cohort comparison.
  setnames(cx,paste0('PH_',target,'_p'),'PH_abundance_p');setnames(cx,paste0(target,'_SD'),'abundance_SD')
  allcx[[paste(co,target)]]<-cx
  gd[,`:=`(cluster=target,zero_fraction=mean(d[[val]]==0),median_zero=median(d[[val]])==0,split_interpretation=if(median(d[[val]])==0)'Fixed median is zero; High/Low exactly equals nonzero/zero; not balanced halves'else'Fixed cohort-wide median; ties low')];allkm[[paste(co,target)]]<-gd
  wt(d[,c('original_sample_id','analysis_unit_id',val,zval,'OS_time_months','OS_event','OS_eligible'),with=FALSE],file.path(root,'OS_input.tsv'))
  cat(co,target,'OS complete\n');flush.console()
 }
}
cx<-rbindlist(allcx,fill=TRUE);cx[,BH_across_cohorts:=p.adjust(p_value,'BH'),by=.(cluster,model)]
km<-rbindlist(allkm,fill=TRUE);km[,BH_across_cohorts:=p.adjust(logrank_p,'BH'),by=cluster]
wt(cx,file.path(OUT,'C0_C5_Cox_comparison.tsv'));wt(km,file.path(OUT,'C0_C5_KM_comparison.tsv'))
f<-cx[cluster=='C5'&status=='OK'];f[,label:=paste0(dataset,' | ',model,' | n=',N,', events=',Event_N)]
p<-ggplot(f,aes(HR,factor(label,levels=rev(label)),colour=model))+geom_vline(xintercept=1,linetype=2,colour='grey60')+geom_errorbar(aes(xmin=lower95,xmax=upper95),orientation='y',width=.15)+geom_point(size=2)+scale_x_log10()+scale_colour_manual(values=c(Univariable='#1F78B4',Multivariable='#D87869'))+theme_p()+theme(legend.position='bottom')+labs(title='C5 estimated abundance and overall survival',x='HR per cohort-wide SD (95% CI)',y=NULL,caption='Full14 CIBERSORTx estimated abundance; no pooled effect.\nIdentical OS definitions, standardization and cutoff rules for C0 and C5.')
savepdf(p,file.path(OUT,'C5_Cox_forest.pdf'),9,6)
