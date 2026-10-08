source('/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/scripts/common.R')
source(file.path(ROOT,'scripts/consensus_functions.R'))
# All sample objects are complete before the cohort runner reaches this stage.
# Repair known export-only compatibility failures without rerunning inference.
source(file.path(ROOT,'scripts/repair_cellchat_exports.R'))
cfg<-readRDS(file.path(ROOT,'01_metadata/frozen_analysis_inputs.rds'))
manifest<-fread(file.path(ROOT,'01_metadata/per_sample_input_manifest.tsv'))
db<-fread(file.path(ROOT,'00_audit/CellChatDB_human_primary_interactions.tsv'))
statuses<-rbindlist(lapply(seq_len(nrow(manifest)),function(i){
 p<-file.path(dirname(manifest$input_path[i]),'status.tsv')
 if(file.exists(p))return(fread(p))
 data.table(sample=manifest$sample[i],sample_id=manifest$sample_id[i],patient_uid=NA_character_,dataset=manifest$dataset[i],status='NOT_ATTEMPTED',error_message='No per-sample status file',runtime_seconds=NA_real_,cell_number=manifest$included_cells[i],eligible_groups=sum(cfg$counts[sample==manifest$sample[i],n_cells]>=30),object_path=NA_character_,warnings='')
}),fill=TRUE)
stopifnot(!anyDuplicated(statuses$sample))
write_tsv(statuses,'02_cellchat_per_sample/sample_run_status.tsv')
if(!any(statuses$status=='SUCCESS'))stop('No successful sample inference; consensus is not estimable')
network<-rbindlist(lapply(seq_len(nrow(manifest)),function(i){
 p<-file.path(dirname(manifest$input_path[i]),'network_long.tsv')
 if(file.exists(p))return(fread(p))
 d<-copy(cfg$eligibility[sample==manifest$sample[i]])
 d[,`:=`(patient_uid=NA_character_,dataset=manifest$dataset[i],inference_status='NOT_ATTEMPTED',weight=NA_real_,count=NA_real_)];d
}),fill=TRUE)
stopifnot(nrow(network)==nrow(manifest)*length(GROUPS)^2,!anyDuplicated(network[,.(sample,source,target)]),all(is.na(network[eligible==FALSE,weight])))
write_tsv(network,'03_cellchat_consensus/sample_network_long.tsv')
stat<-summarize_edges(network)
write_tsv(stat,'03_cellchat_consensus/consensus_edge_statistics.tsv')
raw<-matrix_from_statistics(stat,'consensus_strength');prev<-matrix_from_statistics(stat,'prevalence');display<-matrix_from_statistics(stat,'display_strength')
write_matrix(raw,'03_cellchat_consensus/consensus_strength_raw.tsv')
write_matrix(prev,'03_cellchat_consensus/consensus_prevalence.tsv')
write_matrix(display,'03_cellchat_consensus/consensus_strength_filtered.tsv')
saveRDS(list(raw=raw,prevalence=prev,filtered=display,statistics=stat),file.path(ROOT,'03_cellchat_consensus/consensus_network.rds'))
outgoing<-stat[grepl('^Fibro_C',source)];incoming<-stat[grepl('^Fibro_C',target)]
setorder(outgoing,-prevalence,-consensus_strength,na.last=TRUE);setorder(incoming,-prevalence,-consensus_strength,na.last=TRUE)
write_tsv(outgoing,'03_cellchat_consensus/fibro_outgoing_consensus.tsv');write_tsv(incoming,'03_cellchat_consensus/fibro_incoming_consensus.tsv')

aggregate_feature_level<-function(kind,universe,chunk_size=50L) {
  detected<-rbindlist(lapply(which(statuses$status=='SUCCESS'),function(i){
    p<-file.path(dirname(manifest$input_path[match(statuses$sample[i],manifest$sample)]),paste0(kind,'_detected.tsv'))
    d<-fread(p)
    if(kind=='lr')d<-d[,.(sample,source,target,feature_id=interaction_id,weight)] else d<-d[,.(sample,source,target,feature_id=pathway_name,weight)]
    d
  }),fill=TRUE)
  stopifnot(!anyDuplicated(detected[,.(sample,source,target,feature_id)]))
  base<-copy(network[,.(sample,patient_uid,dataset,source,target,eligible,inference_status)])
  output<-file.path(ROOT,'03_cellchat_consensus',paste0(kind,'_all_samples_long.tsv.gz'))
  if(file.exists(output))stop('Refusing feature-long output overwrite: ',output)
  result<-list();blocks<-split(universe,ceiling(seq_along(universe)/chunk_size))
  for(b in seq_along(blocks)) {
    features<-blocks[[b]]
    long<-base[rep(seq_len(.N),each=length(features))];long[,feature_id:=rep(features,nrow(base))]
    long[,weight:=ifelse(eligible & inference_status=='SUCCESS',0,NA_real_)]
    long[detected[feature_id%in%features],on=.(sample,source,target,feature_id),weight:=i.weight]
    long[,detected:=ifelse(eligible & inference_status=='SUCCESS',weight>0,NA)]
    if(kind=='lr')long<-merge(long,db[,.(feature_id=interaction_id,ligand,receptor,pathway_name,annotation)],by='feature_id',sort=FALSE) else long[,pathway_name:=feature_id]
    fwrite(long,output,sep='\t',na='NA',compress='gzip',append=b>1L,col.names=b==1L)
    ss<-summarize_edges(long,by=c('source','target','feature_id'))
    result[[b]]<-ss
    cat(kind,' consensus block ',b,'/',length(blocks),'\n',sep='');flush.console()
  }
  ans<-rbindlist(result)
  if(kind=='lr')ans<-merge(ans,db[,.(feature_id=interaction_id,ligand,receptor,pathway_name,annotation)],by='feature_id',sort=FALSE) else ans[,pathway_name:=feature_id]
  write_tsv(ans,paste0('03_cellchat_consensus/',kind,'_consensus_statistics.tsv'))
  fibro<-ans[grepl('^Fibro_C',source)|grepl('^Fibro_C',target)]
  setorder(fibro,-prevalence,-consensus_strength,na.last=TRUE)
  write_tsv(fibro,paste0('03_cellchat_consensus/fibro_',kind,'_consensus.tsv'))
  ans
}
pathway_stats<-aggregate_feature_level('pathway',sort(unique(db$pathway_name)))
lr_stats<-aggregate_feature_level('lr',db$interaction_id)

# Sensitivity: normalize within each successful biological sample.
totals<-network[,.(sum_eligible_weight=if(all(inference_status=='SUCCESS'))sum(weight[eligible],na.rm=TRUE) else NA_real_),by=sample]
normalized<-merge(copy(network),totals,by='sample',sort=FALSE)
normalized[,raw_weight:=weight]
normalized[,weight:=ifelse(eligible & inference_status=='SUCCESS' & sum_eligible_weight>0,raw_weight/sum_eligible_weight,NA_real_)]
write_tsv(totals,'03_cellchat_consensus/sensitivity_sample_total_weights.tsv')
write_tsv(normalized,'03_cellchat_consensus/sensitivity_normalized_sample_network_long.tsv')
ns<-summarize_edges(normalized)
write_tsv(ns,'03_cellchat_consensus/sensitivity_normalized_edge_statistics.tsv')
nm<-matrix_from_statistics(ns,'display_strength');write_matrix(nm,'03_cellchat_consensus/sensitivity_normalized_filtered.tsv')
joined<-merge(stat,ns,by=c('source','target'),suffixes=c('_raw','_normalized'))
write_tsv(joined,'03_cellchat_consensus/sensitivity_edge_comparison.tsv')
good<-joined[is.finite(consensus_strength_raw)&is.finite(consensus_strength_normalized)]
rho<-if(nrow(good)>2 && sd(good$consensus_strength_raw)>0 && sd(good$consensus_strength_normalized)>0)cor(good$consensus_strength_raw,good$consensus_strength_normalized,method='spearman') else NA_real_
top_targets<-function(d){
 d<-d[grepl('^Fibro_C',source)&!grepl('^Fibro_C',target)]
 setorder(d,source,-display_strength,-prevalence)
 d[,.(target=paste(head(target[display_strength>0],3),collapse=';')),by=source]
}
tt<-merge(top_targets(stat),top_targets(ns),by='source',all=TRUE,suffixes=c('_raw','_normalized'))
tt[,same_top3_targets:=target_raw==target_normalized];write_tsv(tt,'03_cellchat_consensus/sensitivity_fibro_top_targets.tsv')
patterns<-data.table(fibro_subtype=paste0('Fibro_C',0:6),outgoing_raw=rowSums(display)[1:7],incoming_raw=colSums(display)[1:7],outgoing_normalized=rowSums(nm)[1:7],incoming_normalized=colSums(nm)[1:7])
patterns[,`:=`(outgoing_rank_raw=rank(-outgoing_raw,ties.method='min'),outgoing_rank_normalized=rank(-outgoing_normalized,ties.method='min'),incoming_rank_raw=rank(-incoming_raw,ties.method='min'),incoming_rank_normalized=rank(-incoming_normalized,ties.method='min'))]
write_tsv(patterns,'03_cellchat_consensus/sensitivity_fibro_relative_patterns.tsv')
writeLines(c('# CellChat consensus sensitivity','', '- Primary remains unnormalized per-sample net$weight -> median over all eligible samples including zeros -> n_eligible>=3 and prevalence>=0.30.', '- Sensitivity divides each sample weight by the sum of all eligible weights in that same sample; no pooling.',paste0('- Samples with undefined/zero network total: ',sum(!is.finite(totals$sum_eligible_weight)|totals$sum_eligible_weight<=0),'. Such samples supply NA, not invented normalized zeros.'),paste0('- Spearman edge-strength ranking correlation: ',format(rho,digits=4),'.'),paste0('- Raw positive displayed edges: ',sum(display>0),'; normalized: ',sum(nm>0),'; topology mismatches: ',sum((display>0)!=(nm>0)),'.'),paste0('- Same fibro top-3 external targets: ',sum(tt$same_top3_targets,na.rm=TRUE),'/',nrow(tt),'.'),'- Retained-by-rule edges with median=0 remain in the audit but cannot have visible positive edge width. This follows the requested zero-inclusive median definition.', '', '## Fibro target comparison', '', paste(capture.output(print(tt)),collapse='\n'), '', '## Relative C0–C6 pattern', '', paste(capture.output(print(patterns)),collapse='\n'), '', '- Rank/top-target discrepancies are descriptive checks, not statistical significance tests. Incomplete edges caused by failed inference are explicitly unavailable; no failure is converted to a zero.'),file.path(ROOT,'03_cellchat_consensus/CellChat_consensus_sensitivity.md'))

# Descriptive LLPS integration happens only after inference and consensus.
llps<-fread(file.path(FIBRO_ROOT,'llps_global_res012_cluster_summary.csv'))
llps[,fibro_subtype:=paste0('Fibro_',subcluster)]
integration<-llps[,.(fibro_subtype,LLPS_summary=paste0('frozen LLPS77_AUC median=',format(median,digits=6),'; IQR=',format(IQR,digits=6)))]
integration[,CellChat_outgoing_strength:=rowSums(display)[match(fibro_subtype,GROUPS)]]
integration[,CellChat_incoming_strength:=colSums(display)[match(fibro_subtype,GROUPS)]]
integration[,top_CellChat_targets:=vapply(fibro_subtype,function(f){z<-outgoing[source==f & !grepl('^Fibro_C',target) & display_strength>0];paste(head(z$target,3),collapse=';')},character(1))]
integration[,top_CellChat_pathways:=vapply(fibro_subtype,function(f){z<-pathway_stats[source==f & display_strength>0];z<-z[,.(strength=sum(display_strength)),by=pathway_name];setorder(z,-strength);paste(head(z$pathway_name,5),collapse=';')},character(1))]
integration[,`:=`(top_NicheNet_ligands=NA_character_,top_NicheNet_targets=NA_character_,NicheNet_status=if(file.exists(file.path(ROOT,'05_nichenet_preflight/NICHENET_STATUS.txt')))trimws(readLines(file.path(ROOT,'05_nichenet_preflight/NICHENET_STATUS.txt'))[1]) else 'PREFLIGHT',NicheNet_activity_status='NOT_RUN_RECEIVER_PROGRAM_REQUIRED')]
write_tsv(integration,'08_cross_method/fibro_LLPS_descriptive_integration.tsv')

# CellChat establishes pair support only. NicheNet activity remains uncomputed.
cross<-lr_stats[grepl('^Fibro_C',source)&!grepl('^Fibro_C',target),.(source,target,ligand,receptor,pathway_if_available=pathway_name,CellChat_support=retained,CellChat_n_eligible=n_eligible,CellChat_n_detected=n_detected,CellChat_prevalence=prevalence,CellChat_consensus_strength=consensus_strength,NicheNet_ligand_activity=NA_real_,NicheNet_target_genes=NA_character_,comparison_status='DEFERRED_RECEIVER_GENESET_OI_REQUIRED')]
write_tsv(cross,'08_cross_method/CellChat_fibro_LR_for_NicheNet_followup.tsv')
cat('CONSENSUS COMPLETE\n')
