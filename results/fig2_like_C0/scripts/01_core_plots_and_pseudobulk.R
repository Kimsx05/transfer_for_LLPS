source('results/fig2_like_C0/scripts/common.R')
suppressPackageStartupMessages({library(Seurat);library(qs);library(ggrastr)})
future::plan(future::sequential)
status('A','Reading frozen object and validating labels, scores and palette')
obj <- qs::qread(OBJECT,nthreads=2)
old_colors <- fread(COLOR_SOURCE)
stopifnot(identical(old_colors$subcluster,CLUSTERS),identical(old_colors$hex,unname(COLORS)))
write_tsv(data.table(cluster=CLUSTERS,hex=unname(COLORS),source=COLOR_SOURCE),'cluster_color_mapping.tsv')
um <- Embeddings(obj,'umap')
d <- data.table(cell=colnames(obj),sample=as.character(obj$orig.ident),dataset=as.character(obj$dataset),
 cluster=as.character(obj$subcluster),LLPS77_AUC=as.numeric(obj$LLPS77_AUC),UMAP_1=um[,1],UMAP_2=um[,2])
frozen <- fread(file.path(FIBROOT,'global_res012_cell_labels_scores.csv'))
stopifnot(nrow(d)==20146L,setequal(d$cluster,CLUSTERS),!anyNA(d),
 identical(d$cluster,paste0('C',obj$clean_snn_res.0.12)))
# Independent barcode join to earlier frozen analysis metadata.
oldmd <- readRDS(file.path(BASE,'cellchat_nichechat_fibro7/01_metadata/frozen_analysis_inputs.rds'))$metadata
stopifnot(identical(d$cluster,as.character(oldmd$fibro_subtype[match(d$cell,oldmd$cell)])),
 identical(d$cluster,frozen$subcluster[match(d$cell,frozen$cell_id)]),
 isTRUE(all.equal(d$LLPS77_AUC,frozen$LLPS77_AUC[match(d$cell,frozen$cell_id)],tolerance=1e-14)))
d[,cluster:=factor(cluster,levels=CLUSTERS)]
write_tsv(d,'cell_metadata_umap_LLPSscore.tsv.gz')
centers <- d[,.(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),by=cluster]
ut <- theme_pub()+theme(axis.text=element_blank(),axis.ticks=element_blank())
p <- ggplot(d[sample(.N)],aes(UMAP_1,UMAP_2))+
 geom_point_rast(aes(color=cluster),size=.30,raster.dpi=450)+
 geom_label(data=centers,aes(label=cluster),size=3,label.padding=unit(.12,'lines'),fill='white',label.size=.15)+
 scale_color_manual(values=COLORS,breaks=CLUSTERS,name='Cluster',drop=FALSE)+
 guides(color=guide_legend(override.aes=list(size=3)))+coord_fixed()+ut+
 labs(title='Fibroblast / pericyte subclusters',x='UMAP 1',y='UMAP 2')
export_plot(p,'Fig2A_fibro_pericyte_umap',6.4,5.3)
status('A','COMPLETE')

settings <- readRDS(file.path(dirname(FIBROOT),'clean_v4_settings.rds'))
stopifnot(length(settings$global_limits)==2,all(d$LLPS77_AUC[d$LLPS77_AUC>0]>=settings$global_limits[1]),
 all(d$LLPS77_AUC<=settings$global_limits[2]))
p <- ggplot()+geom_point_rast(data=d[LLPS77_AUC>0][order(LLPS77_AUC)],aes(UMAP_1,UMAP_2,color=LLPS77_AUC),size=.30,raster.dpi=450)+
 geom_point_rast(data=d[LLPS77_AUC==0],aes(UMAP_1,UMAP_2),color='black',size=.30,raster.dpi=450)+
 scale_color_gradientn(colors=viridisLite::plasma(256),limits=settings$global_limits,name='LLPS score',
 guide=guide_colorbar(barheight=unit(32,'mm')))+coord_fixed()+ut+
 labs(title='LLPS score',x='UMAP 1',y='UMAP 2',caption='Inherited LLPS77 AUCell score; atlas-wide color scale. Zero = black.')
export_plot(p,'Fig2B1_LLPSscore_umap_continuous',6.4,5.3)
write_tsv(data.table(min=settings$global_limits[1],max=settings$global_limits[2],n_zero=sum(d$LLPS77_AUC==0),palette='plasma; exact zero black'),'LLPS_color_scale.tsv')
status('B1','COMPLETE')

s <- d[,.(n_cells=.N,n_samples=uniqueN(sample),mean=mean(LLPS77_AUC),median=median(LLPS77_AUC),
 Q1=unname(quantile(LLPS77_AUC,.25)),Q3=unname(quantile(LLPS77_AUC,.75)),SD=sd(LLPS77_AUC)),by=cluster]
ss <- d[,.(n_cells=.N,mean=mean(LLPS77_AUC),median=median(LLPS77_AUC)),by=.(sample,dataset,cluster)]
s <- merge(s,ss[,.(median_of_sample_means=median(mean),median_of_sample_medians=median(median)),by=cluster],by='cluster')
write_tsv(s,'cluster_LLPSscore_summary.tsv');write_tsv(ss,'sample_cluster_LLPSscore_summary.tsv')
p <- ggplot(d,aes(cluster,LLPS77_AUC,fill=cluster))+geom_violin(scale='width',trim=TRUE,linewidth=.35)+
 geom_boxplot(width=.13,fill='white',outlier.shape=NA,linewidth=.35)+
 scale_fill_manual(values=COLORS,drop=FALSE)+theme_pub()+theme(legend.position='none')+
 labs(title='LLPS score across stromal subclusters',x=NULL,y='LLPS77 AUCell score',caption='Cell distributions shown descriptively; box = median and IQR.')
export_plot(p,'Fig2B2_LLPSscore_by_cluster_violin',6.5,4.6)
status('B2','COMPLETE')

# All 59 source samples stay in the table, including undefined zero-denominator rows.
m <- fread(MANIFEST)
m <- m[,.(sample,dataset,tissue_origin,patient_uid,patient_mapping_status,paired_specimen_code,
 sample_accession,replicate_unit_status,primary_replication_allowed)]
ct <- d[,.(n_C0=sum(cluster=='C0'),n_all_Fibro_Pericyte=.N),by=sample]
m <- merge(m,ct,by='sample',all.x=TRUE,sort=FALSE)
m[is.na(n_C0),`:=`(n_C0=0L,n_all_Fibro_Pericyte=0L)]
m[,C0_proportion:=fifelse(n_all_Fibro_Pericyte>0,n_C0/n_all_Fibro_Pericyte,NA_real_)]
stopifnot(sum(m$n_C0)==10717,sum(m$n_all_Fibro_Pericyte)==20146,!anyDuplicated(m$sample),
 all(d$sample%in%m$sample),all(m$tissue_origin%in%c('Tumor','Normal')))
paired_ids <- m[!is.na(patient_uid)&n_all_Fibro_Pericyte>0&primary_replication_allowed==TRUE,
 .(nt=uniqueN(tissue_origin)),by=patient_uid][nt==2,patient_uid]
m[,known_TN_pair:=!is.na(patient_uid)&patient_uid%in%paired_ids]
m[,main_independent_included:=n_all_Fibro_Pericyte>0 & primary_replication_allowed==TRUE & !known_TN_pair]
m[,exclusion_reason:=fcase(n_all_Fibro_Pericyte==0,'ZERO_DENOMINATOR_UNDEFINED',
 primary_replication_allowed==FALSE,'KNOWN_SAME_PATIENT_ALTERNATE_SPECIMEN',
 known_TN_pair,'KNOWN_TN_PAIR_ANALYZED_SEPARATELY',default='INCLUDED')]
write_tsv(m,'C0_fraction_per_sample.tsv')
z <- m[main_independent_included==TRUE];z[,tissue_origin:=factor(tissue_origin,levels=c('Normal','Tumor'))]
wt <- wilcox.test(C0_proportion~tissue_origin,data=z,exact=FALSE,correct=TRUE)
stats <- data.table(analysis='independent_main',test='two-sided Wilcoxon rank-sum; asymptotic, continuity corrected',
 n_Normal=sum(z$tissue_origin=='Normal'),n_Tumor=sum(z$tissue_origin=='Tumor'),n_pairs=NA_integer_,
 median_Normal=median(z[tissue_origin=='Normal',C0_proportion]),median_Tumor=median(z[tissue_origin=='Tumor',C0_proportion]),p_value=wt$p.value)
box <- function(x,subtitle,caption) ggplot(x,aes(tissue_origin,C0_proportion,fill=tissue_origin))+
 geom_boxplot(width=.48,outlier.shape=NA,linewidth=.45)+
 geom_point(position=position_jitter(width=.10,height=0,seed=123),shape=21,size=2,stroke=.35)+
 scale_fill_manual(values=GROUP_COLORS)+scale_y_continuous(labels=scales::percent,limits=c(0,1),expand=expansion(mult=c(.02,.06)))+
 theme_pub()+theme(legend.position='none')+labs(title='C0 proportion',subtitle=subtitle,x=NULL,
 y='C0 / all Fibro + Pericyte',caption=caption)
p <- box(z,sprintf('Normal n = %d; Tumor n = %d | Wilcoxon P = %s',stats$n_Normal,stats$n_Tumor,ptext(wt$p.value)),
 'Independent-specimen comparison; known T/N pairs analyzed separately.')
export_plot(p,'Fig2C_C0_fraction_tumor_vs_normal',5.6,4.8)
# Threshold sensitivity changes sample eligibility only, never the denominator.
zs <- z[n_all_Fibro_Pericyte>=50]
ws <- wilcox.test(C0_proportion~tissue_origin,data=zs,exact=FALSE)
stats <- rbind(stats,data.table(analysis='independent_min50_sensitivity',test=stats$test,
 n_Normal=sum(zs$tissue_origin=='Normal'),n_Tumor=sum(zs$tissue_origin=='Tumor'),n_pairs=NA_integer_,
 median_Normal=median(zs[tissue_origin=='Normal',C0_proportion]),median_Tumor=median(zs[tissue_origin=='Tumor',C0_proportion]),p_value=ws$p.value))
if(length(paired_ids)) {
 pa <- m[known_TN_pair&primary_replication_allowed==TRUE];pa[,tissue_origin:=factor(tissue_origin,levels=c('Normal','Tumor'))]
 stopifnot(all(pa[,.N,by=.(patient_uid,tissue_origin)]$N==1))
 wide <- dcast(pa,patient_uid~tissue_origin,value.var='C0_proportion')
 wp <- wilcox.test(wide$Tumor,wide$Normal,paired=TRUE,exact=FALSE,correct=TRUE)
 p <- ggplot(pa,aes(tissue_origin,C0_proportion,group=patient_uid))+geom_line(color='grey65',linewidth=.4)+
 geom_point(aes(fill=tissue_origin),shape=21,size=2.6,stroke=.4)+scale_fill_manual(values=GROUP_COLORS)+
 scale_y_continuous(labels=scales::percent,limits=c(0,1))+theme_pub()+theme(legend.position='none')+
 labs(title='C0 proportion: matched specimens',subtitle=sprintf('%d verified pairs | Paired Wilcoxon P = %s',nrow(wide),ptext(wp$p.value)),
 x=NULL,y='C0 / all Fibro + Pericyte',caption='Pairs require verified study-patient IDs; no pairing inferred from names.')
 export_plot(p,'Fig2C_C0_fraction_tumor_vs_normal_paired',5.6,4.8)
 write_tsv(wide,'C0_fraction_verified_pairs.tsv')
 stats <- rbind(stats,data.table(analysis='verified_paired',test='two-sided paired Wilcoxon; asymptotic, continuity corrected',
 n_Normal=nrow(wide),n_Tumor=nrow(wide),n_pairs=nrow(wide),median_Normal=median(wide$Normal),median_Tumor=median(wide$Tumor),p_value=wp$p.value))
}
write_tsv(stats,'C0_fraction_statistics.tsv')
write_tsv(m[n_all_Fibro_Pericyte>0,.(n=.N,median=median(C0_proportion)),by=.(dataset,tissue_origin)],'C0_fraction_by_dataset.tsv')
status('C','COMPLETE')

# Counts only; no normalized or integrated expression for pseudobulk.
counts <- LayerData(obj,assay='RNA',layer='counts')
stopifnot(identical(colnames(counts),d$cell),all(counts@x>=0),all(counts@x==round(counts@x)),
 all(Matrix::colSums(counts)==obj$nCount_RNA))
groups <- paste(d$sample,ifelse(d$cluster=='C0','C0','Others'),sep='::')
lev <- sort(unique(groups)); agg <- sparseMatrix(i=seq_along(groups),j=match(groups,lev),x=1,
 dims=c(nrow(d),length(lev)))
pb <- counts%*%agg;colnames(pb)<-lev
pm <- d[,.(n_cells=.N),by=.(sample,group=ifelse(cluster=='C0','C0','Others'))]
pm[,library:=paste(sample,group,sep='::')]
pm <- merge(pm,m[,.(sample,dataset,tissue_origin,patient_uid,primary_replication_allowed)],by='sample')
pm <- pm[match(lev,library)]
stopifnot(identical(pm$library,colnames(pb)))
saveRDS(list(counts=pb,metadata=pm),file.path(CORE,'C0_Others_pseudobulk_counts.rds'),compress='gzip')
write_tsv(pm,'pseudobulk_library_metadata.tsv')
capture.output(sessionInfo(),file=file.path(OUT,'summary/sessionInfo_core.txt'))
status('PSEUDOBULK_INPUT','COMPLETE')
