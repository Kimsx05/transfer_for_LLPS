set.seed(123)
ROOT<-'/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7'
NN_LIB<-file.path(ROOT,'Rlib_nichenet')
.libPaths(unique(c(NN_LIB,file.path(ROOT,'Rlib'),'/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_atlas_stage1_v1/Rlib','/home/jinshengxi/Rlibrary/LLPS_pipeline_v1_20260914_150459','/home/jinshengxi/Rlibrary/packages',.libPaths())))
suppressPackageStartupMessages({library(Matrix);library(data.table);library(qs2)})
future::plan(future::sequential)
GROUPS<-readRDS(file.path(ROOT,'01_metadata/frozen_analysis_inputs.rds'))$groups
write_tsv<-function(x,path)fwrite(x,file.path(ROOT,path),sep='\t',na='NA',quote='auto')
sha256<-function(path)strsplit(system2('sha256sum',shQuote(path),stdout=TRUE),' ')[[1]][1]
NN_OUT<-file.path(ROOT,'05_nichenet_preflight')
NN_RESOURCES<-c(ligand_target='ligand_target_matrix_nsga2r_final.rds',ligand_receptor='lr_network_human_21122021.rds',weighted_networks='weighted_networks_nsga2r_final.rds')
