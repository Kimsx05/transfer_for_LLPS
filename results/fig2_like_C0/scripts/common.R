# Frozen inputs; run scripts from /home/jinshengxi/bladder_metabolism.
BASE <- '/home/jinshengxi/bladder_metabolism'
OUT <- file.path(BASE, 'results/fig2_like_C0')
CORE <- file.path(OUT, 'phase1_core')
OLD <- file.path(BASE, 'cellchat_nichechat_fibro7/11_C0_dual_geneset_preparation')
ATLAS <- '/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes'
FIBROOT <- file.path(ATLAS, 'fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order')
OBJECT <- file.path(FIBROOT, 'global_res012_object.qs')
COLOR_SOURCE <- file.path(FIBROOT, 'beautified_SCP_protein_coding_top3/cluster_SCP_Paired_colors.csv')
MANIFEST <- file.path(OLD, '01_C0_abundance/A1_sample_abundance_manifest.tsv')
.libPaths(unique(c(file.path(BASE,'cellchat_nichechat_fibro7/Rlib'),
 '/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_atlas_stage1_v1/Rlib',
 '/home/jinshengxi/Rlibrary/LLPS_pipeline_v1_20260914_150459',
 '/home/jinshengxi/Rlibrary/packages', .libPaths())))
suppressPackageStartupMessages({library(data.table); library(Matrix); library(ggplot2)})
set.seed(20261025)
CLUSTERS <- paste0('C',0:6)
# Exact saved SCP Paired mapping; never regenerate palette from observed levels.
COLORS <- c(C0='#A6CEE3',C1='#1F78B4',C2='#B2DF8A',C3='#33A02C',
 C4='#FDBF6F',C5='#FF7F00',C6='#FB9A99')
GROUP_COLORS <- c(Normal='#8DB5C7',Tumor='#D87869')
write_tsv <- function(x,name) fwrite(x,file.path(CORE,name),sep='\t',na='NA')
theme_pub <- function() theme_classic(base_size=11,base_family='sans')+
 theme(plot.title=element_text(size=13,face='bold'),plot.subtitle=element_text(size=9),
 plot.caption=element_text(size=8,hjust=0),axis.text=element_text(color='black'),
 legend.title=element_text(size=10),plot.margin=margin(10,12,10,10))
export_plot <- function(p,stem,w=6,h=4.6) {
 ggsave(file.path(CORE,paste0(stem,'.pdf')),p,width=w,height=h,device=cairo_pdf,bg='white')
 ggsave(file.path(CORE,paste0(stem,'.png')),p,width=w,height=h,dpi=350,device=ragg::agg_png,bg='white')
}
ptext <- function(p) if(is.na(p)) 'NA' else formatC(p,format='g',digits=3)
status <- function(step,message) {
 cat(format(Sys.time()),step,message,'\n')
 cat(paste(step,message,sep='\t'),'\n',file=file.path(OUT,'summary/execution_status.tsv'),append=TRUE)
}
