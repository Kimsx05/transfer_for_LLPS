source('results/llps_C5_parallel_v1/scripts/common.R')
CORE<-file.path(OUT,'02_Clinical_Stage')
m<-fread(file.path(OUT,'01_Normal_Tumor/C5_fraction_per_sample_reviewed.tsv'))
w<-fread(file.path(BASE,'cellchat_nichechat_fibro7/11_C0_dual_geneset_preparation/00_preflight/user_workbook_relevant_scRNA_rows.tsv'))
stopifnot(!anyDuplicated(w$xlsx_row),all(m$workbook_row%in%w$xlsx_row))
m[,`:=`(TNM_raw=w$F[match(workbook_row,w$xlsx_row)],NMIBC_MIBC=fifelse(source_sample_type%in%c('NMIBC','MIBC'),source_sample_type,NA_character_))]
# Only explicit T-prefixed categories, never overall pathologic stage.
m[,T_stage:=ifelse(grepl('^(p|c|yp)?T(a|is|[1-4])',TNM_raw),sub('^(p|c|yp)?(T(a|is|[1-4])).*$','\\2',TNM_raw),NA_character_)]
m[tissue_origin!='Tumor',`:=`(T_stage=NA_character_,NMIBC_MIBC=NA_character_)]
m[,T_definition:=ifelse(grepl('^pT',TNM_raw),'pathologic T',ifelse(grepl('^cT',TNM_raw),'clinical T','reported T; modality unspecified'))]
m[,stage_eligible:=analysis_included&tissue_origin=='Tumor'];wt(m,file.path(CORE,'sample_stage_mapping.tsv'))
stage_cols<-c('Ta'='#D9EAF2','Tis'='#C9DAD0','T1'='#A6CEE3','T2'='#71AFCE','T3'='#4F98C5','T4'='#1F78B4','NMIBC'='#A6CEE3','MIBC'='#1F78B4')
stats<-list()
for(k in c('NMIBC_MIBC','T_stage')){
 lev<-if(k=='T_stage')c('Ta','T1','T2','T3','T4')else c('NMIBC','MIBC')
 v<-m[stage_eligible&get(k)%in%lev];v[,stage:=factor(get(k),levels=lev)]
 gs<-v[,.(N=.N,n_patients=uniqueN(patient_key),n_datasets=uniqueN(dataset),n_zero=sum(n_C5==0),mean=mean(C5_proportion),median=median(C5_proportion),Q1=quantile(C5_proportion,.25),Q3=quantile(C5_proportion,.75)),by=stage]
 wt(gs,file.path(CORE,paste0(k,'_group_summary.tsv')));wt(v[,.(N=.N),by=.(dataset,stage)],file.path(CORE,paste0(k,'_dataset_counts.tsv')))
 pv<-rho<-tp<-NA_real_;reason<-'Descriptive only: fewer than 3 samples in at least one group or fewer than two groups'
 if(nrow(gs)>=2&&all(gs$N>=3)){
 tst<-if(k=='NMIBC_MIBC')wilcox.test(C5_proportion~stage,v,exact=FALSE)else kruskal.test(C5_proportion~stage,v)
 pv<-tst$p.value;reason<-'Exploratory sample-level comparison; dataset composition may confound'
 if(k=='T_stage'){ct<-cor.test(as.numeric(v$stage),v$C5_proportion,method='spearman',exact=FALSE);rho<-unname(ct$estimate);tp<-ct$p.value}
 }
 stats[[k]]<-data.table(comparison=k,N=nrow(v),n_patients=uniqueN(v$patient_key),test=if(k=='NMIBC_MIBC')'Wilcoxon rank-sum'else'Kruskal-Wallis',P=pv,trend_rho=rho,trend_P=tp,note=reason)
 labsx<-setNames(paste0(gs$stage,'\n(n=',gs$N,')'),gs$stage)
 p<-ggplot(v,aes(stage,C5_proportion,fill=stage))+geom_boxplot(width=.55,outlier.shape=NA,linewidth=.4)+geom_point(position=position_jitter(width=.13,height=0,seed=123),size=1.1,alpha=.6)+scale_fill_manual(values=stage_cols)+scale_x_discrete(labels=labsx)+scale_y_continuous(labels=scales::percent,limits=c(0,1))+theme_pub()+theme(legend.position='none')+labs(title=paste('C5 abundance |',k),subtitle=paste0('N = ',nrow(v),' | P = ',ptxt(pv)),x=NULL,y='C5 / all Fibro + Pericyte',caption=reason)
 export(p,file.path(CORE,paste0(k,'_boxplot')),6.6,5.3)
 p<-ggplot(gs,aes(stage,mean,fill=stage))+geom_col(width=.55)+geom_point(data=v,aes(stage,C5_proportion),position=position_jitter(width=.13,height=0,seed=123),inherit.aes=FALSE,size=1.1,alpha=.6)+scale_fill_manual(values=stage_cols)+scale_x_discrete(labels=labsx)+scale_y_continuous(labels=scales::percent,limits=c(0,1))+theme_pub()+theme(legend.position='none')+labs(title=paste('C5 abundance |',k),subtitle='Bars: sample mean; points: individual specimens',x=NULL,y='C5 / all Fibro + Pericyte',caption=reason)
 export(p,file.path(CORE,paste0(k,'_barplot')),6.6,5.3)
}
wt(rbindlist(stats),file.path(CORE,'stage_statistics.tsv'))
