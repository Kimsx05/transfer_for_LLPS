set.seed(123)
ROOT <- '/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7'
.libPaths(unique(c(file.path(ROOT,'Rlib'), '/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_atlas_stage1_v1/Rlib', '/home/jinshengxi/Rlibrary/LLPS_pipeline_v1_20260914_150459', '/home/jinshengxi/Rlibrary/packages', .libPaths())))
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table);library(qs);library(qs2)})
future::plan(future::sequential)
GROUPS <- c(paste0('Fibro_C',0:6),'Epithelial','Endothelial','Myeloid','T_NK','B','Plasma','Mast')
MAJOR_OBJECT <- '/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/results/objects/major_annotated_with_LLPS77_AUC.qs'
FIBRO_OBJECT <- '/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order/global_res012_object.qs'
FIBRO_ROOT <- dirname(FIBRO_OBJECT)
write_tsv <- function(x,path) fwrite(x,file.path(ROOT,path),sep='\t',na='NA',quote='auto')
sha256 <- function(path) strsplit(system2('sha256sum',shQuote(path),stdout=TRUE),' ')[[1]][1]
blank_matrix <- function(value=NA_real_) matrix(value,14,14,dimnames=list(GROUPS,GROUPS))
