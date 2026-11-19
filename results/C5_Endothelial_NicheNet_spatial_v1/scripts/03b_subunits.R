options(mc.cores=1L);source('cellchat_nichechat_fibro7/scripts/nichenet_common.R');setDTthreads(1);library(CellChat)
O<-'results/C5_Endothelial_NicheNet_spatial_v1';a<-fread(file.path(O,'07_spatial_expression_validation/frozen_axes.tsv'));db<-CellChatDB.human
subunits<-function(g){if(g%in%rownames(db$complex)){v<-unlist(db$complex[g,,drop=FALSE],use.names=FALSE);paste(v[!is.na(v)&nzchar(v)],collapse=';')}else g}
a[,required_subunits:=vapply(receptor,subunits,character(1))];su<-fread(file.path(O,'03_candidate_ligands/expression_summary.tsv'))[group=='Endothelial'];lr<-as.data.table(readRDS(file.path(NN_OUT,'lr_network_human_21122021.rds')))
a[,all_subunits_expression_gate:=vapply(strsplit(required_subunits,';',fixed=TRUE),function(gs)all(gs%in%su$gene)&&all(su[match(gs,gene),pooled_pct]>=.1 & su[match(gs,gene),support_patients]>=2),logical(1))]
a[,NicheNet_direct_receptor_subunits:=mapply(function(l,ss)paste(intersect(strsplit(ss,';',fixed=TRUE)[[1]],lr[from==l,to]),collapse=';'),ligand,required_subunits)]
a[,receptor_definition:='CellChatDB.human 2.1.2 complex mandatory subunits; NicheNet pairs do not replace this definition'];fwrite(a,file.path(O,'07_spatial_expression_validation/frozen_axes_subunits.tsv'),sep='\t')
# All requested axes retained for expression QC even if not selected for spatial testing.
q<-as.data.table(db$interaction,keep.rownames='interaction');q<-q[ligand%in%c('ADM','MDK','SEMA3C','WNT2','FN1','COL1A1','COL1A2')];q[,required_subunits:=vapply(receptor,subunits,character(1))];fwrite(q,file.path(O,'06_CellChat_integration/requested_axes_receptor_definitions.tsv'),sep='\t')
lt<-readRDS(file.path(NN_OUT,'ligand_target_matrix_nsga2r_final.rds'));defs<-fread(file.path(O,'02_receiver_geneset/frozen_program_genes.tsv'));acts<-fread(file.path(O,'04_NicheNet_ligand_activity/all_program_activities.tsv'));out<-list()
for(nm in unique(acts$program))for(l in unique(a$ligand)){g<-defs[program==nm & in_background==TRUE,gene];q<-as.data.table(nichenetr::get_weighted_ligand_target_links(l,g,lt,n=250));q[,program:=nm];out[[length(out)+1]]<-q}
fwrite(rbindlist(out,fill=TRUE),file.path(O,'05_ligand_target_links/selected_axis_links.tsv'),sep='\t');cat('SUBUNIT_COMPLETE\n')
