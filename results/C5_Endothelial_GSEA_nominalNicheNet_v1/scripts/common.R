options(mc.cores=1L)
source('cellchat_nichechat_fibro7/scripts/nichenet_common.R');setDTthreads(1)
O<-'results/C5_Endothelial_GSEA_nominalNicheNet_v1';OLD<-'results/C5_Endothelial_NicheNet_spatial_v1'
w<-function(x,p)fwrite(x,file.path(O,p),sep='\t',na='NA')
rd<-function(p)fread(file.path(OLD,p))
cfg<-jsonlite::fromJSON(file.path(O,'00_input_audit/prespecified_parameters.json'))
sets<-readRDS('results/llps_C0_functional_signature_v2/resources/gene_sets.rds')$H
bg<-rd('02_receiver_geneset/expressed_background.tsv')
ligs<-sort(unique(rd('04_NicheNet_ligand_activity/all_program_activities.tsv')$test_ligand))
stopifnot(nrow(bg)==7070,length(ligs)==246,length(sets)==50)
models<-c('primary','GSE222315','InhouseData','stage_adjusted','sampling_adjusted')
DE<-setNames(lapply(models,function(m)rd(paste0('01_receiver_pseudobulk/DE_',m,'.tsv'))),models)
selectgenes<-function(d,ih,gs,p=.05){Reduce(intersect,list(d[logFC>0 & P.Value<p,gene],ih[logFC>0,gene],gs[logFC>0,gene],bg$gene))}
