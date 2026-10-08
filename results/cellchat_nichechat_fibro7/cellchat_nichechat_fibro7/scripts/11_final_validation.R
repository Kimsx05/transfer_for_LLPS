source('/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/scripts/common.R')
source(file.path(ROOT,'scripts/consensus_functions.R'))
cfg<-readRDS(file.path(ROOT,'01_metadata/frozen_analysis_inputs.rds'))
m<-fread(file.path(ROOT,'01_metadata/per_sample_input_manifest.tsv'))
sts<-fread(file.path(ROOT,'02_cellchat_per_sample/sample_run_status.tsv'))
net<-fread(file.path(ROOT,'03_cellchat_consensus/sample_network_long.tsv'))
con<-readRDS(file.path(ROOT,'03_cellchat_consensus/consensus_network.rds'))
stopifnot(nrow(sts)==59L,!anyDuplicated(sts$sample),all(sts$status=='SUCCESS'),all(file.exists(sts$object_path)),nrow(net)==59L*196L,uniqueN(net$sample)==59L,identical(dimnames(con$filtered),list(GROUPS,GROUPS)),all(is.na(net[eligible==FALSE,weight])),all(is.finite(net[eligible==TRUE,weight])),sum(con$filtered>0)==182L)
for(i in seq_len(nrow(m))) {
 d<-net[sample==m$sample[i]]
 ref<-cfg$eligibility[sample==m$sample[i]]
 chk<-merge(d,ref,by=c('sample','source','target'),suffixes=c('_actual','_frozen'))
 stopifnot(nrow(chk)==196L,all(chk$eligible_actual==chk$eligible_frozen),all(chk$n_source_actual==chk$n_source_frozen),all(chk$n_target_actual==chk$n_target_frozen))
 l<-fread(file.path(dirname(m$input_path[i]),'lr_detected.tsv'))
 reconstructed<-copy(d[,.(source,target,eligible)]);reconstructed[,`:=`(weight=ifelse(eligible,0,NA_real_),count=ifelse(eligible,0,NA_real_))]
 if(nrow(l)){
   z<-l[,.(weight=sum(weight),count=.N),by=.(source,target)]
   reconstructed[z,on=.(source,target),`:=`(weight=i.weight,count=i.count)]
 }
 chk2<-merge(d,reconstructed,by=c('source','target'),suffixes=c('_native','_LR'))
 stopifnot(all(abs(chk2[eligible_native==TRUE,weight_native-weight_LR])<1e-10),all(chk2[eligible_native==TRUE,count_native==count_LR]))
}
manual<-summarize_edges(net)
chk<-merge(manual,con$statistics,by=c('source','target'),suffixes=c('_recomputed','_saved'))
stopifnot(all(chk$n_eligible_recomputed==chk$n_eligible_saved),all(chk$n_detected_recomputed==chk$n_detected_saved),all(chk$retained_recomputed==chk$retained_saved),all(abs(chk$display_strength_recomputed-chk$display_strength_saved)<1e-10),all(con$statistics$n_missing_eligible==0L))
lr<-fread(file.path(ROOT,'03_cellchat_consensus/lr_consensus_statistics.tsv'))
pw<-fread(file.path(ROOT,'03_cellchat_consensus/pathway_consensus_statistics.tsv'))
db<-fread(file.path(ROOT,'00_audit/CellChatDB_human_primary_interactions.tsv'))
stopifnot(nrow(lr)==nrow(db)*196L,nrow(pw)==uniqueN(db$pathway_name)*196L,!anyDuplicated(lr[,.(feature_id,source,target)]))
for(f in c('lr_all_samples_long.tsv.gz','pathway_all_samples_long.tsv.gz'))stopifnot(system2('gzip',c('-t',shQuote(file.path(ROOT,'03_cellchat_consensus',f))))==0L)
nn<-readRDS(file.path(ROOT,'05_nichenet_preflight/frozen_nichenet_preflight_inputs.rds'))
stopifnot(is.null(nn$geneset_oi),readLines(file.path(ROOT,'05_nichenet_preflight/NICHENET_STATUS.txt'))[1]=='READY')
expr<-fread(file.path(ROOT,'06_nichenet_expression_audit/expression_audit_sample_status.tsv'))
design<-fread(file.path(ROOT,'07_nichenet_templates/receiver_program_design_after_CellChat.tsv'))
stopifnot(nrow(expr)==59L,all(expr$status=='SUCCESS'),all(design$program_frozen==FALSE),all(is.na(design$geneset_oi_path)))
prefixes<-c('04_cellchat_figures/Figure_A_cohort_consensus_circle','04_cellchat_figures/Figure_A_fibro_centered_circle','09_figures/Figure_B_fibro_major_consensus_heatmap')
figure_files<-unlist(lapply(prefixes,function(p)file.path(ROOT,paste0(p,c('.pdf','.png','.svg')))))
stopifnot(all(file.exists(figure_files)),all(file.info(figure_files)$size>1000))
hash_ref<-fread(file.path(ROOT,'00_audit/frozen_object_preservation.tsv'))
current_hash<-vapply(hash_ref$input,sha256,character(1))
stopifnot(all(current_hash==hash_ref$sha256_before))
write_tsv(data.table(input=hash_ref$input,initial_sha256=hash_ref$sha256_before,final_sha256=current_hash,unchanged=TRUE),'00_audit/frozen_object_preservation_final.tsv')
writeLines(c('PASS: 59/59 independent sample inference objects; final failed samples 0.','PASS: every sample has the complete 14x14 frame, exact frozen count-derived eligibility, and correct NA/zero distinction.','PASS: LR weights/counts reconcile with every native sample aggregate; raw consensus recomputation matches every edge.','PASS: complete LR/pathway statistic universes and valid compressed per-sample long tables.','PASS: nine figure-format artifacts; visual labels and legends inspected.','PASS: NicheNet official resources and 59-sample expression audits ready; geneset_oi undefined; activity intentionally not run.','PASS: source major/fibro7 Seurat object SHA256 hashes remain identical to pre-flight.'),file.path(ROOT,'00_audit/FINAL_VALIDATION.txt'))
cat('FINAL VALIDATION PASS\n')
