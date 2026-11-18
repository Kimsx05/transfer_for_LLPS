suppressPackageStartupMessages({library(ggplot2);library(data.table);library(jsonlite);library(patchwork);library(png);library(grid)})
setDTthreads(1);set.seed(123)
sk='/home/jinshengxi/.agents/skills/jins-routine-plotting';source(file.path(sk,'common/theme.R'));source(file.path(sk,'common/palettes.R'));source(file.path(sk,'common/export.R'));jin_setup();theme_set(theme_bw(base_family='Arial',base_size=12))
null_device_files=character();cowplot::set_null_device(function(width,height){f=tempfile(fileext='.pdf');null_device_files<<-c(null_device_files,f);cairo_pdf(f,width=width,height=height,family='Arial')})
force_overview=FALSE
o='results/C5_spatial_exploratory_all63_v1';root='/home/jinshengxi/LLPS/spatial_data/LLPS_spatial_transfer_manual'
m=fread(file.path(o,'00_inputs_and_QC/section_mapping.tsv'));labels=fread(file.path(o,'00_inputs_and_QC/spatial_tissue_patient_audit.tsv'));m=merge(m,labels[,.(section_id,tissue)],by='section_id');m=m[order(dataset,tissue,section_id)]
dataset_colors=setNames(JIN_CLUSTER_COLORS,sort(unique(m$dataset)))
fwrite(data.table(dataset=names(dataset_colors),color=unname(dataset_colors)),file.path(o,'tables/dataset_colors.tsv'),sep='\t')
lim=fromJSON(file.path(o,'00_inputs_and_QC/shared_color_limits.json'))
cols=c('A26_UCell','A26_AUCell','core15_UCell','five_UCell','MMP11_log1p_CP10k','major_mean__Fibroblast','major_mean__Epithelial','nCount','A26_UCell_percentile','C2C4_program_UCell')
titles=c('A26 UCell','A26 AUCell','Core15 UCell','Five-gene UCell','MMP11 log1p(CP10k)','Fibroblast mean','Epithelial mean','log10(1 + UMI)','A26 within-section percentile','C2/C4 program UCell')
export=function(p,path,w,h){ggsave(paste0(path,'.pdf'),p,device=cairo_pdf,width=w,height=h,bg='white',limitsize=FALSE);ggsave(paste0(path,'.png'),p,device=ragg::agg_png,width=w,height=h,dpi=300,bg='white',limitsize=FALSE)}
shard=as.integer(commandArgs(TRUE)[1]);m=m[(seq_len(nrow(m))-1) %% 3 == shard];force_overview=TRUE
thumbs=list();scale_rows=list()
for(i in seq_len(nrow(m))){
 sid=m$section_id[i];message(sid);while(!file.exists(file.path(o,'02_spatial_maps',paste0(sid,'_plot_data.tsv.gz'))))Sys.sleep(1);d=fread(file.path(o,'02_spatial_maps',paste0(sid,'_plot_data.tsv.gz')));plots=list();legs=list();iw=m$image_width[i];ih=m$image_height[i];dx=diff(range(d$x))*.03;dy=diff(range(d$y))*.03;xmin=min(d$x)-dx;xmax=max(d$x)+dx;ymin=min(d$y)-dy;ymax=max(d$y)+dy
 for(j in seq_along(cols)){
  col=cols[j];v=d[[col]];upper=if(col %in% names(lim))lim[[col]] else 1
  if(col=='nCount'){v=log10(1+v);upper=log10(1+upper)}
  dd=data.frame(x=d$x,y=d$y,value=v)
  p=ggplot(dd,aes(x,y,color=value))+geom_point(size=.65,shape=16)+scale_color_viridis_c(limits=c(0,upper),oob=scales::squish,name=NULL,guide=guide_colorbar(display='rectangles'))+scale_x_continuous(limits=c(xmin,xmax),expand=c(0,0))+scale_y_reverse(limits=c(ymax,ymin),expand=c(0,0))+coord_fixed()+labs(title=titles[j])+theme_void(base_family='Arial',base_size=12)+theme(plot.title=element_text(size=12),legend.position='bottom',plot.background=element_rect(fill='white',color=NA))
  plots[[j]]=p;legs[[j]]=cowplot::get_legend(p);scale_rows[[length(scale_rows)+1]]=data.frame(section_id=sid,panel=col,lower=0,upper=upper,n_clipped=sum(v>upper,na.rm=TRUE))
 }
 title=paste(sid,m$dataset[i],m$tissue[i],sep=' | ')
 p=wrap_plots(plots,ncol=5)+plot_annotation(title=title,subtitle='C5-associated signatures (exploratory; not C5 abundance). Shared raw scales; abundance/expression/UMI clipped at global 99th percentile.',caption='Original pixel coordinates, upper origin; window cropped to observed spots. No HE in primary maps. Percentile panel is within-section only.')
 path=file.path(o,'02_spatial_maps',paste0(sid,'_overview'))
 if(force_overview || !all(file.exists(paste0(path,c('.pdf','.png','_panel_only.pdf','_legend_only.pdf'))))){export(p,path,20,9);ggsave(paste0(path,'_panel_only.pdf'),p & theme(legend.position='none'),device=cairo_pdf,width=20,height=9,bg='white');lg=cowplot::plot_grid(plotlist=lapply(legs,cowplot::ggdraw),labels=titles,ncol=5,label_size=10,label_fontfamily='Arial');ggsave(paste0(path,'_legend_only.pdf'),lg,device=cairo_pdf,width=20,height=3,bg='white')}
 thumb=plots[[9]]+labs(title=sid,subtitle=paste(m$dataset[i],m$tissue[i]))+theme(legend.position='none',plot.title=element_text(size=9),plot.subtitle=element_text(size=8));thumbs[[sid]]=thumb
 # Provisional HE: use manifest geometric registration only, never independent validation.
 if(m$registration_status[i]=='PASS_ID_DIMENSIONS_BOUNDS'){
  im=readPNG(file.path(root,m$histology_relative_path[i]));dd=data.frame(x=d$x,y=d$y,value=d$A26_UCell_percentile)
  he=ggplot(dd,aes(x,y,color=value))+annotation_custom(rasterGrob(im,interpolate=FALSE),xmin=0,xmax=iw,ymin=-ih,ymax=0)+geom_point(size=.8,alpha=.7)+scale_color_viridis_c(limits=c(0,1),name='Within-section percentile',guide=guide_colorbar(display='rectangles'))+scale_x_continuous(limits=c(0,iw),expand=c(0,0))+scale_y_reverse(limits=c(ih,0),expand=c(0,0))+coord_fixed()+theme_void(base_family='Arial',base_size=12)+labs(title=sid,subtitle='Provisional HE overlay (geometry checked; histology unverified)',caption='A26-associated score; not validated C5 localization')+theme(plot.background=element_rect(fill='white',color=NA))
  hp=file.path(o,'02_spatial_maps',paste0(sid,'_HE_provisional'));if(!all(file.exists(paste0(hp,c('.pdf','.png','_panel_only.pdf','_legend_only.pdf')))))jin_export(he,hp,8,7)
 }
}
fwrite(rbindlist(scale_rows),file.path(o,paste0('tables/figure_color_scales_shard',shard,'.tsv')),sep='\t')

unlink(null_device_files)
