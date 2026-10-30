# C0 phase1B common statistics and style, with transparent 300 dpi export.
source('results/fig2_like_C0/phase1B_review/scripts/common.R')
OUT<-file.path(BASE,'results/llps_C5_parallel_v1');CORE<-file.path(OUT,'01_Normal_Tumor')
export_plot<-function(p,stem,w=7,h=5){p<-p+theme(plot.background=element_rect(fill='transparent',colour=NA),panel.background=element_rect(fill='transparent',colour=NA),legend.background=element_rect(fill='transparent',colour=NA),legend.key=element_rect(fill='transparent',colour=NA));ggsave(file.path(CORE,paste0(stem,'.pdf')),p,width=w,height=h,device=cairo_pdf,bg='transparent');ggsave(file.path(CORE,paste0(stem,'.png')),p,width=w,height=h,dpi=300,bg='transparent')}
m<-fread(file.path(BASE,'results/fig2_like_C0/phase1B_review/sample_metadata_reviewed.tsv'))
d<-fread(file.path(BASE,'results/fig2_like_C0/phase1_core/cell_metadata_umap_LLPSscore.tsv.gz'))
a<-fread(file.path(BASE,'cellchat_nichechat_fibro7/11_C0_dual_geneset_preparation/01_C0_abundance/A1_sample_abundance_manifest.tsv'))
ct<-d[,.(n_C5=sum(cluster=='C5'),denominator=.N),by=sample]
m<-merge(m,ct,by='sample',all.x=TRUE);m[is.na(n_C5),n_C5:=0L]
stopifnot(sum(m$n_C5)==716,all(m[!is.na(denominator),denominator==n_all_Fibro_Pericyte]),all(a$n_C5[match(m$sample,a$sample)]==m$n_C5))
m[,C5_proportion:=fifelse(n_all_Fibro_Pericyte>0,n_C5/n_all_Fibro_Pericyte,NA_real_)]
m<-merge(m,a[,.(sample,n_total_cells,source_sample_type,workbook_row)],by='sample')
m[,C5_fraction_all_cells:=n_C5/n_total_cells]
m[,tissue_origin:=factor(tissue_origin,levels=c('Normal','Tumor'))];z<-m[analysis_included==TRUE]
WT(m,'C5_fraction_per_sample_reviewed.tsv')
# Exact phase1B statistical and plotting block, terminating before LLPS analyses.
src<-paste(readLines('results/fig2_like_C0/phase1B_review/scripts/02_sample_statistics.R'),collapse='\n')
a0<-regexpr('coverage<-m',src,fixed=TRUE)[1];b0<-regexpr('# C0 state and C0 proportion',src,fixed=TRUE)[1]
code<-gsub('C0','C5',substr(src,a0,b0-1),fixed=TRUE);eval(parse(text=code))
# Main independent comparison excludes verified T/N pairs, using updated mapping.
paids<-paired(z,'C5_proportion')$patient_uid
zz<-copy(z[!patient_uid%in%paids]);tst<-wilcox.test(C5_proportion~tissue_origin,zz,exact=FALSE,correct=TRUE)
st<-zz[,.(N=.N,median=median(C5_proportion),n_zero=sum(n_C5==0)),by=tissue_origin];st[,p_value:=tst$p.value];WT(st,'C5_independent_TN_statistics.tsv')
src0<-paste(readLines('results/fig2_like_C0/scripts/01_core_plots_and_pseudobulk.R'),collapse='\n')
a0<-regexpr('box <- function',src0,fixed=TRUE)[1];b0<-regexpr('p <- box(z',src0,fixed=TRUE)[1]
eval(parse(text=gsub('C0','C5',substr(src0,a0,b0-1),fixed=TRUE)))
export_plot(box(zz,sprintf('Normal n = %d; Tumor n = %d | Wilcoxon P = %s',sum(zz$tissue_origin=='Normal'),sum(zz$tissue_origin=='Tumor'),ptext(tst$p.value)),'Independent specimens; verified T/N pairs analyzed separately.'),'C5_Normal_Tumor',5.6,4.8)
# All-atlas denominator sensitivity: existing exact sample/cluster count mapping.
allz<-m[primary_replication_allowed&n_total_cells>0]
WT(dataset_summary(allz,'C5_fraction_all_cells','C5_all_cells'),'C5_all_cells_dataset_statistics.tsv')
plots_dataset(allz,'C5_fraction_all_cells','C5 proportion among all atlas cells','C5 / all atlas cells','C5_all_cells_dataset',TRUE)
# Add explicit zero cluster entries before averaging specimen proportions.
grid<-CJ(sample=z$sample,cluster=CLUSTERS)
comp<-merge(grid,d[,.N,by=.(sample,cluster)],by=c('sample','cluster'),all.x=TRUE);comp[is.na(N),N:=0L]
comp<-merge(comp,z[,.(sample,dataset,tissue_origin,patient_key,n_all_Fibro_Pericyte)],by='sample');comp[,fraction:=N/n_all_Fibro_Pericyte]
WT(comp,'all_cluster_sample_composition.tsv')
for(kind in c('pooled_cells','sample_equal')){
 cs<-if(kind=='pooled_cells')comp[,.(value=sum(N)),by=.(tissue_origin,cluster)] else comp[,.(value=mean(fraction)),by=.(tissue_origin,cluster)]
 cs[,fraction:=value/sum(value),by=tissue_origin];cs[,cluster:=factor(cluster,levels=CLUSTERS)];WT(cs,paste0('stacked_',kind,'.tsv'))
 p<-ggplot(cs,aes(tissue_origin,fraction,fill=cluster))+geom_col(width=.65,position=position_stack(reverse=TRUE))+scale_fill_manual(values=COLORS,breaks=CLUSTERS,drop=FALSE)+scale_y_continuous(labels=scales::percent,limits=c(0,1),expand=c(0,0))+theme_pub()+labs(title='Fibroblast / pericyte composition',subtitle=if(kind=='pooled_cells')'Pooled cells; descriptive composition'else'Equal-weight mean of specimen proportions',x=NULL,y='Proportion',fill='Cluster')
 export_plot(p,paste0('Normal_Tumor_stacked_',kind),5.6,4.8)
}
WT(d[,.N,by=.(dataset,cluster)],'cluster_dataset_cell_counts.tsv')
