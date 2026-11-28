Sys.setenv(OMP_NUM_THREADS=4,OPENBLAS_NUM_THREADS=4,MKL_NUM_THREADS=4)
library(data.table)
library(ggplot2)
setDTthreads(4)
s='/home/jinshengxi/.agents/skills/jins-routine-plotting'
for(f in c('common/theme.R','common/palettes.R','common/export.R','plot_types/Stacked_composition/plot.R','plot_types/Violin_distribution/plot.R','plot_types/UMAP_cluster/plot.R'))source(file.path(s,f))
jin_setup()
original_export=jin_export
jin_export=function(plot,path,...){
 args=list(...);required=paste0(path,c('.pdf','.png'))
 if(!is.null(args$legend_plot))required=c(required,paste0(path,c('_panel_only.pdf','_legend_only.pdf')))
 if(all(file.exists(required)) && all(file.info(required)$size>0))return(invisible(NULL))
 original_export(plot,path,...)
}
o='results/LLPS_squidpy_niche_v1';F=file.path(o,'09_figures');pal=fread(file.path(F,'niche_colors.tsv'),colClasses=c(niche='character'));colors=setNames(pal$color,pal$niche);ord=NULL
comp=fread(file.path(o,'04_niche_characterization/niche_composition.tsv'));comp[,niche:=as.character(niche)];ord=as.character(sort(as.integer(unique(comp$niche))))
old=fread('results/C5_spatial_colocalization_esda_v1/00_input_audit/frozen_cell_colors.tsv');ct=sort(unique(comp$celltype));cellcolors=SCP::palette_scp(factor(ct,levels=ct),type='discrete');names(cellcolors)=ct;cellcolors[old$celltype]=old$color;fwrite(data.table(celltype=ct,color=cellcolors[ct]),file.path(F,'celltype_colors.tsv'),sep='\t')
for(k in c('major','subtype'))for(sc in c('center','neighborhood'))for(b in c('pooled','section')){
 z=comp[kind==k & scale==sc & balance==b];z[,fraction:=fraction/sum(fraction),by=niche] # float32 closure correction only
 p=plot_stacked_composition(z,'niche','celltype','fraction',condition_order=ord,category_order=sort(unique(z$celltype)),colors=cellcolors,xlab='Niche')
 jin_export(p,file.path(F,'composition',paste(k,sc,b,'stacked',sep='_')),width=9,height=6,legend_plot=p,panel_width=6)
}
v=fread(file.path(F,'score_plot_data.tsv.gz'));v[,niche:=as.character(niche)]
for(f in c('five_UCell','LLPS77_UCell','A26_UCell','five_UCell__adjusted','LLPS77_UCell__adjusted')){
 label=switch(f,five_UCell='C5-associated score',LLPS77_UCell='LLPS score',A26_UCell='C5 A26 score',five_UCell__adjusted='Adjusted C5-associated score',LLPS77_UCell__adjusted='Adjusted LLPS score')
 p=plot_violin_distribution(v,'niche',f,method='',ylab=label,order=ord,colors=colors);jin_export(p,file.path(F,paste0(f,'_violin')),width=7,height=4.6)
}
ss=fread(file.path(o,'05_C5_LLPS_enrichment/section_level_niche_scores.tsv'));ss=ss[five_UCell__count>=20];ss[,niche:=factor(as.character(niche),levels=ord)]
for(f in c('five_UCell','LLPS77_UCell','A26_UCell')){
 col=paste0(f,'__mean');p=ggplot(ss,aes(x=niche,y=.data[[col]],fill=niche))+geom_boxplot(outlier.shape=NA,linewidth=.4)+geom_point(position=position_jitter(width=.15,seed=123),size=1)+scale_fill_manual(values=colors)+jin_theme()+theme(legend.position='none')+labs(x='Niche',y=switch(f,five_UCell='C5-associated score',LLPS77_UCell='LLPS score',A26_UCell='C5 A26 score'))
 jin_export(p,file.path(F,paste0(f,'_section_boxplot')),width=7,height=4.6)
}
for(role in c('major_best','subtype_best')){
 z=fread(file.path(F,paste0(role,'_umap_data.tsv.gz')));set.seed(123);if(nrow(z)>50000)z=z[sample(.N,50000)]
 for(group in c(paste0(role,'_niche'),'dataset','section_id','kmeans_niche')){
  z[,(group):=as.character(get(group))];lev=sort(unique(z[[group]]));cc=setNames(SCP::palette_scp(factor(lev,levels=lev),type='discrete'),lev)
  if(grepl('niche',group))cc=colors[lev]
  p=plot_umap_cluster(z,'x','y',group,order=lev,colors=cc,point_size=.2,label_size=if(group=='section_id')1 else 3,seed=123)
  jin_export(p,file.path(F,paste(role,group,'umap',sep='_')),width=if(group=='section_id')12 else 8,height=if(group=='section_id')12 else 5.4,legend_plot=p,panel_width=5.4)
 }
}
p=fread(file.path(o,'08_clinical_association/section_niche_proportions_complete.tsv'));p[,niche:=as.character(niche)]
q=plot_stacked_composition(p,'section_id','niche','proportion',condition_order=unique(p$section_id),category_order=ord,colors=colors)+theme(axis.text.x=element_text(angle=90,hjust=1,size=7))
jin_export(q,file.path(F,'clinical_section_composition'),width=18,height=6,legend_plot=q,panel_width=15)
writeLines(capture.output(sessionInfo()),file.path(o,'00_input_audit/R_session.txt'))

mc=fread(file.path(o,'04_niche_characterization/model_cell_composition.tsv'))
for(m in c('kmeans_niche','subtype_best_niche'))for(k in c('major','subtype')){
 z=mc[model==m & kind==k];z[,niche:=factor(niche,levels=sort(unique(niche)))];p=ggplot(z,aes(celltype,niche,fill=fraction))+geom_tile()+scale_fill_viridis_c(name='Cell fraction')+jin_theme()+theme(axis.text.x=element_text(angle=60,hjust=1))+labs(x=NULL,y='Niche')
 jin_export(p,file.path(F,paste(m,k,'composition_heatmap',sep='_')),width=if(k=='subtype')13 else 8,height=6,legend_plot=p,panel_width=if(k=='subtype')11 else 6)
}
