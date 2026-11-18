options(mc.cores=1L)
suppressPackageStartupMessages({library(ggplot2);library(data.table);library(patchwork)})
setDTthreads(1);theme_set(theme_bw(base_family='Arial',base_size=12));o='results/C5_fivegene_LLPS77_spatial_v1';d=fread(file.path(o,'figures/GSE171351_A1_plot_data.tsv.gz'));plots=list()
for(col in c('LLPS77_UCell','LLPS77_AUCell','five_UCell','five_AUCell')){dd=data.frame(x=d$x,y=d$y,score=d[[col]]);plots[[col]]=ggplot(dd,aes(x,y,color=score))+geom_point(size=1.2)+scale_color_viridis_c(limits=c(0,1),guide=guide_colorbar(display='rectangles'))+scale_y_reverse()+coord_fixed()+theme_void(base_family='Arial',base_size=12)+theme(legend.position='bottom')+labs(title=col,color=NULL)}
p=wrap_plots(plots,ncol=2)+plot_annotation(title='Technical pilot: GSE171351_A1');ggsave(file.path(o,'00_input_audit/pilot_preview.pdf'),p,device=cairo_pdf,width=10,height=9,bg='white');ggsave(file.path(o,'00_input_audit/pilot_preview.png'),p,device=ragg::agg_png,width=10,height=9,dpi=300,bg='white')
