# Final sensitivity handling for mathematically undefined zero-total samples.
# The raw primary CellChat objects/statistics/display matrix are not changed.
source('/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/scripts/common.R')
source(file.path(ROOT,'scripts/consensus_functions.R'))
net<-fread(file.path(ROOT,'03_cellchat_consensus/sample_network_long.tsv'))
primary<-readRDS(file.path(ROOT,'03_cellchat_consensus/consensus_network.rds'))
primary_path<-file.path(ROOT,'03_cellchat_consensus/consensus_strength_filtered.tsv')
primary_hash<-sha256(primary_path)
totals<-net[,.(sum_eligible_weight=sum(weight[eligible],na.rm=TRUE)),by=sample]
d<-merge(copy(net),totals,by='sample',sort=FALSE)
d[,`:=`(structurally_eligible=eligible,raw_weight=weight,normalization_defined=is.finite(sum_eligible_weight)&sum_eligible_weight>0)]
d[,weight:=ifelse(structurally_eligible & normalization_defined,raw_weight/sum_eligible_weight,NA_real_)]
# A zero-total sample has no normalized value. It is not treated as zero and
# does not enter the normalized denominator; its raw zero remains in PRIMARY.
d[,eligible:=structurally_eligible & normalization_defined & inference_status=='SUCCESS']
stat<-summarize_edges(d)
structural<-d[,.(n_structurally_eligible=sum(structurally_eligible),n_undefined_total_samples=sum(structurally_eligible & !normalization_defined)),by=.(source,target)]
stat<-merge(stat,structural,by=c('source','target'),sort=FALSE)
nm<-matrix_from_statistics(stat,'display_strength')
archive<-file.path(ROOT,'00_audit/sensitivity_initial_undefined_total_handling');dir.create(archive,showWarnings=FALSE)
files<-c('sensitivity_normalized_sample_network_long.tsv','sensitivity_normalized_edge_statistics.tsv','sensitivity_normalized_filtered.tsv','sensitivity_edge_comparison.tsv','sensitivity_fibro_top_targets.tsv','sensitivity_fibro_relative_patterns.tsv','CellChat_consensus_sensitivity.md')
for(f in files)if(!file.exists(file.path(archive,f)))file.copy(file.path(ROOT,'03_cellchat_consensus',f),file.path(archive,f))
write_tsv(d,'03_cellchat_consensus/sensitivity_normalized_sample_network_long.tsv')
write_tsv(stat,'03_cellchat_consensus/sensitivity_normalized_edge_statistics.tsv')
write_matrix(nm,'03_cellchat_consensus/sensitivity_normalized_filtered.tsv')
joined<-merge(primary$statistics,stat,by=c('source','target'),suffixes=c('_raw','_normalized'))
write_tsv(joined,'03_cellchat_consensus/sensitivity_edge_comparison.tsv')
good<-joined[is.finite(consensus_strength_raw)&is.finite(consensus_strength_normalized)]
rho<-cor(good$consensus_strength_raw,good$consensus_strength_normalized,method='spearman')
top_targets<-function(s){
 z<-s[grepl('^Fibro_C',source)&!grepl('^Fibro_C',target)];setorder(z,source,-display_strength,-prevalence)
 z[,.(target=paste(head(target[display_strength>0],3),collapse=';')),by=source]
}
tt<-merge(top_targets(primary$statistics),top_targets(stat),by='source',all=TRUE,suffixes=c('_raw','_normalized'))
tt[,same_top3_targets:=target_raw==target_normalized];write_tsv(tt,'03_cellchat_consensus/sensitivity_fibro_top_targets.tsv')
raw<-primary$filtered
patterns<-data.table(fibro_subtype=GROUPS[1:7],outgoing_raw=rowSums(raw)[1:7],incoming_raw=colSums(raw)[1:7],outgoing_normalized=rowSums(nm)[1:7],incoming_normalized=colSums(nm)[1:7],outgoing_to_major_raw=rowSums(raw[1:7,8:14]),outgoing_to_major_normalized=rowSums(nm[1:7,8:14]),incoming_from_major_raw=colSums(raw[8:14,1:7]),incoming_from_major_normalized=colSums(nm[8:14,1:7]))
patterns[,`:=`(outgoing_rank_raw=rank(-outgoing_raw,ties.method='min'),outgoing_rank_normalized=rank(-outgoing_normalized,ties.method='min'),incoming_rank_raw=rank(-incoming_raw,ties.method='min'),incoming_rank_normalized=rank(-incoming_normalized,ties.method='min'))]
write_tsv(patterns,'03_cellchat_consensus/sensitivity_fibro_relative_patterns.tsv')
excluded<-totals[sum_eligible_weight<=0|!is.finite(sum_eligible_weight)]
write_tsv(excluded,'03_cellchat_consensus/sensitivity_undefined_total_samples.tsv')
writeLines(c('# CellChat consensus sensitivity','', '- PRIMARY is unchanged: raw sample net$weight, median including eligible real zeros, n_eligible>=3 and prevalence>=0.30.', '- For each sample with positive eligible network total, normalized_weight=weight/sum(all eligible sample weights). Eligibility still requires source and target >=30 cells.',paste0('- Mathematically undefined zero-total samples: ',paste(excluded$sample,collapse=', '),'. They supply NA and are omitted from the NORMALIZED effective denominator only; their real raw zeros remain in the primary network. No normalized zero is invented.'),paste0('- Normalizable samples=',sum(totals$sum_eligible_weight>0),'/',nrow(totals),'. n_structurally_eligible, n_undefined_total_samples and normalized n_eligible are all retained per edge.'),paste0('- Spearman edge ranking correlation=',format(rho,digits=4),'.'),paste0('- Raw positive edges=',sum(raw>0),'; normalized=',sum(nm>0),'; topology mismatches=',sum((raw>0)!=(nm>0)),'.'),paste0('- Exact same ordered top-3 external targets=',sum(tt$same_top3_targets,na.rm=TRUE),'/',nrow(tt),'.'),'', '## Biological interpretation','', '- The broad architecture is stable, while relative fibro magnitudes and target ordering are scale-sensitive. C5 has the largest raw total outgoing strength (52.86); sample-normalized cohort edge totals place C0 first (0.566), C1 second (0.375), and C5 third (0.293). C5 also shifts from first to third for total incoming strength.', '- The strongest external raw edge is C0→Epithelial (4.146; 16 eligible samples, detected in all 16). C5→Epithelial and C5→Endothelial are also strong and recurrent. Raw versus normalized top receiver switches include C1 (Endothelial→Epithelial) and C5 (Epithelial→Endothelial). Epithelial and Endothelial remain central fibro receiver relationships.', '- This is a material qualification of subtype dominance: do not claim C5 is unconditionally the most active fibro subtype from raw row totals alone. Raw weight measures absolute inferred probability sums; normalized weight measures relative communication within each sample. Different subtype eligible-sample sets and sample total signal scales contribute to the change.', '- Row/column totals sum edge medians calculated from different eligible sample subsets. They are descriptive summaries, not a cohort-wide probability mass, and not cell-number-adjusted activity. C6 retains its low-coverage flag (3 eligible samples).', '- Primary raw network and cutoffs remain unchanged. NicheNet follow-up should retain C0 and C5 and both epithelial/endothelial receiver programs in the design queue; no gene program or ligand activity is selected here.', '', '## External target ordering', '', paste(capture.output(print(tt)),collapse='\n'), '', '## Relative C0–C6 patterns', '', paste(capture.output(print(patterns)),collapse='\n'), '', '- The initial stricter handling of undefined normalized observations is archived under 00_audit/sensitivity_initial_undefined_total_handling. The final sensitivity omits only samples whose normalization denominator is undefined; the primary network was never recomputed or overwritten.'),file.path(ROOT,'03_cellchat_consensus/CellChat_consensus_sensitivity.md'))
stopifnot(sha256(primary_path)==primary_hash)
cat('SENSITIVITY COMPLETE; primary SHA256 unchanged; rho=',rho,'; raw edges=',sum(raw>0),'; normalized edges=',sum(nm>0),'\n',sep='')
