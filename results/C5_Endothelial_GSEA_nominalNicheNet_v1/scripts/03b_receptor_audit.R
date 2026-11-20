source('results/C5_Endothelial_GSEA_nominalNicheNet_v1/scripts/common.R')
lr<-as.data.table(readRDS(file.path(NN_OUT,'lr_network_human_21122021.rds')));w(unique(lr[,.(ligand=from,receptor=to)]),'00_input_audit/official_human_v2_LR_pairs.tsv')
capture.output(sessionInfo(),file=file.path(O,'00_input_audit/sessionInfo_expression.txt'))
