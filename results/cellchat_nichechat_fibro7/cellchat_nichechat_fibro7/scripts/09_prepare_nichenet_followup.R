source('/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/scripts/nichenet_common.R')
stats<-fread(file.path(ROOT,'03_cellchat_consensus/consensus_edge_statistics.tsv'))
nn<-readRDS(file.path(NN_OUT,'frozen_nichenet_preflight_inputs.rds'))
pairs<-stats[source%in%GROUPS[1:7] & target%in%GROUPS[8:14]]
pairs[,primary_receiver:=target%in%c('Epithelial','Myeloid','Endothelial','T_NK')]
pairs[,followup_candidate:=retained & consensus_strength>0]
setorder(pairs,-followup_candidate,-primary_receiver,-prevalence,-consensus_strength,na.last=TRUE)
write_tsv(pairs,'07_nichenet_templates/CellChat_fibro_receiver_pair_priority.tsv')
design<-pairs[followup_candidate==TRUE,.(pair_id=paste(source,target,sep='__'),source,receiver=target,CellChat_n_eligible=n_eligible,CellChat_prevalence=prevalence,CellChat_consensus_strength=consensus_strength,receiver_program_label='TO_BE_DEFINED',geneset_oi_path=NA_character_,program_source_reference=NA_character_,selection_rule=NA_character_,sample_replication_design=NA_character_,program_frozen=FALSE,seed=123)]
write_tsv(design,'07_nichenet_templates/receiver_program_design_after_CellChat.tsv')
lr<-fread(file.path(ROOT,'08_cross_method/CellChat_fibro_LR_for_NicheNet_followup.tsv'))
lr<-lr[pairs[followup_candidate==TRUE,.(source,target)],on=.(source,target),nomatch=0L]
lr<-lr[CellChat_support==TRUE]
official<-nn$lr_network
lr[,NicheNet_ligand_in_model:=ligand%in%nn$modeled_ligands]
lr[,NicheNet_exact_LR_link:=paste(ligand,receptor,sep='::')%in%paste(official$ligand,official$receptor,sep='::')]
complexes<-fread(file.path(ROOT,'00_audit/CellChatDB_human_complex_definitions.tsv'))
subunits<-function(g){
 if(!g%in%complexes$complex_id)return(g)
 z<-unlist(complexes[complex_id==g,!c('complex_id'),with=FALSE],use.names=FALSE)
 paste(unique(z[!is.na(z)&nzchar(z)]),collapse=';')
}
lr[,ligand_subunits:=vapply(ligand,subunits,character(1))]
lr[,receptor_subunits:=vapply(receptor,subunits,character(1))]
lr[,complex_match_note:=ifelse(ligand%in%complexes$complex_id | receptor%in%complexes$complex_id,'COMPLEX: exact identifier mismatch is not negative evidence; requires explicit subunit-aware reconciliation','SINGLE_GENE: exact official symbol match auditable')]
send<-fread(file.path(ROOT,'06_nichenet_expression_audit/sender_ligand_expression_summary.tsv'))
receive<-fread(file.path(ROOT,'06_nichenet_expression_audit/receiver_receptor_expression_summary.tsv'))
lr<-merge(lr,send[,.(source=group,ligand=gene,sender_ligand_n_evaluable=n_evaluable,sender_ligand_median_fraction=median_pct_expressing,sender_ligand_ge10pct_prevalence=prevalence_ge10pct)],by=c('source','ligand'),all.x=TRUE,sort=FALSE)
lr<-merge(lr,receive[,.(target=group,receptor=gene,receiver_receptor_n_evaluable=n_evaluable,receiver_receptor_median_fraction=median_pct_expressing,receiver_receptor_ge10pct_prevalence=prevalence_ge10pct)],by=c('target','receptor'),all.x=TRUE,sort=FALSE)
lr[,NicheNet_ligand_activity:=NA_real_];lr[,NicheNet_target_genes:=as.character(NicheNet_target_genes)]
write_tsv(lr,'08_cross_method/CellChat_supported_LR_with_NicheNet_preflight_evidence.tsv')
writeLines(c('# NicheNet follow-up after CellChat','',paste0('- Stable fibro-to-major pairs with positive primary median: ',nrow(design),'. Full candidate receiver universe retains all seven major groups.'),'- Priority table first marks positive stable CellChat pairs, then the four requested primary receivers, then prevalence and strength. It is a follow-up queue, not ligand activity.', '- receiver_program_design_after_CellChat.tsv lists exact pairs and leaves receiver program, gene file, provenance and selection rule undefined. No program is silently selected.', '- LR evidence joins CellChat feature support and audited sender/receptor expression where exact single-gene identifiers match. Missing complex-identifier joins remain NA; CellChat complex subunits are retained for later explicit reconciliation.', '- Expression prevalence across all eligible sender/receiver groups is contextual. Mechanism validation must also examine joint eligible samples for that specific source-target pair; separate cohort marginals do not prove same-sample co-expression.', '- Ligand-target prior membership and exact LR membership are resource compatibility evidence, not inferred activity or mechanism confidence. No BOTH_SUPPORTED label is assigned at this stage.', '- Next required design: an interpretable receiver transcriptional program with independent provenance; expressed background and sample replication plan; explicit handling of receptor/ligand complexes and resource mismatches. Preserve activity metrics and sample-level coverage, without disease-stage comparisons.'),file.path(ROOT,'07_nichenet_templates/NicheNet_followup_design_notes.md'))
cat('NICHENET FOLLOWUP TEMPLATE COMPLETE; LIGAND ACTIVITY NOT RUN\n')
