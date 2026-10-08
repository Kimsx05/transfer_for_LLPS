source('/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/scripts/nichenet_common.R')
cfg<-readRDS(file.path(ROOT,'01_metadata/frozen_analysis_inputs.rds'))
nn<-readRDS(file.path(NN_OUT,'frozen_nichenet_preflight_inputs.rds'))
manifest<-fread(file.path(ROOT,'01_metadata/per_sample_input_manifest.tsv'))
out<-file.path(ROOT,'06_nichenet_expression_audit')
dir.create(out,showWarnings=FALSE)
sender_file<-file.path(out,'sender_ligand_expression_per_sample.tsv.gz')
receiver_file<-file.path(out,'receiver_receptor_expression_per_sample.tsv.gz')
background_file<-file.path(out,'receiver_expressed_background_candidates.tsv.gz')
if(any(file.exists(c(sender_file,receiver_file,background_file))))stop('Expression audit output exists; refusing accidental overwrite')
all_genes<-sort(unique(c(nn$ligands,nn$receptors,nn$target_genes)))
long_parts<-list(sender=list(),receiver=list());sample_status<-list()
for(i in seq_len(nrow(manifest))) {
  start<-Sys.time();s<-manifest$sample[i]
  input<-qs2::qs_read(manifest$input_path[i],nthreads=1);x<-input$expression;sm<-input$metadata
  stopifnot(identical(colnames(x),sm$cell),all(sm$sample==s))
  present_genes<-intersect(all_genes,rownames(x))
  groups<-factor(sm$cellchat_group,levels=GROUPS);n<-as.numeric(table(groups))
  G<-sparseMatrix(i=seq_len(ncol(x)),j=as.integer(groups),x=1,dims=c(ncol(x),length(GROUPS)),dimnames=list(colnames(x),GROUPS))
  use<-x[present_genes,,drop=FALSE]
  means<-sweep(as.matrix(use%*%G),2,n,'/');fraction<-sweep(as.matrix((use>0)%*%G),2,n,'/')
  means[,n==0]<-NA_real_;fraction[,n==0]<-NA_real_
  make_long<-function(role,genes,group_names) {
    d<-CJ(gene=genes,group=group_names)
    gi<-match(d$gene,present_genes);gj<-match(d$group,GROUPS)
    d[,`:=`(sample=s,patient_uid=NA_character_,dataset=manifest$dataset[i],role=role,n_cells=n[gj],gene_in_RNA=!is.na(gi))]
    d[,eligible:=n_cells>=30L]
    d[,mean_normalized_expression:=means[cbind(gi,gj)]];d[,pct_expressing:=fraction[cbind(gi,gj)]]
    d[,expressed_ge10pct:=ifelse(eligible & gene_in_RNA,pct_expressing>=.1,NA)]
    d[,any_expression:=ifelse(eligible & gene_in_RNA,pct_expressing>0,NA)]
    d[eligible==FALSE,`:=`(mean_normalized_expression=NA_real_,pct_expressing=NA_real_)]
    d
  }
  sender<-make_long('sender_ligand',nn$ligands,GROUPS[1:7]);receiver<-make_long('receiver_receptor',nn$receptors,GROUPS[8:14])
  fwrite(sender,sender_file,sep='\t',na='NA',compress='gzip',append=i>1L,col.names=i==1L)
  fwrite(receiver,receiver_file,sep='\t',na='NA',compress='gzip',append=i>1L,col.names=i==1L)
  background<-make_long('receiver_background_candidate',nn$target_genes,GROUPS[8:14])
  background<-background[eligible==TRUE & gene_in_RNA==TRUE & pct_expressing>=.10]
  fwrite(background,background_file,sep='\t',na='NA',compress='gzip',append=i>1L,col.names=i==1L)
  long_parts$sender[[i]]<-sender;long_parts$receiver[[i]]<-receiver
  sample_status[[i]]<-data.table(sample=s,sample_id=manifest$sample_id[i],status='SUCCESS',included_cells=ncol(x),runtime_seconds=as.numeric(difftime(Sys.time(),start,units='secs')))
  cat('EXPRESSION AUDIT ',i,'/',nrow(manifest),' ',s,'\n',sep='');flush.console()
  rm(input,x,use,means,fraction,background);gc(FALSE)
}
summarize_expression<-function(d) {
  d[,{
    valid<-eligible & gene_in_RNA;v<-pct_expressing[valid];m<-mean_normalized_expression[valid]
    list(n_eligible=sum(eligible),n_evaluable=sum(valid),gene_in_RNA=all(gene_in_RNA),n_samples_any_expression=sum(v>0),n_samples_ge10pct=sum(v>=.1),prevalence_any_expression=if(length(v))mean(v>0) else NA_real_,prevalence_ge10pct=if(length(v))mean(v>=.1) else NA_real_,median_pct_expressing=if(length(v))median(v) else NA_real_,median_mean_normalized_expression=if(length(m))median(m) else NA_real_,max_pct_expressing=if(length(v))max(v) else NA_real_)
  },by=.(group,gene,role)]
}
sender_summary<-summarize_expression(rbindlist(long_parts$sender));receiver_summary<-summarize_expression(rbindlist(long_parts$receiver))
write_tsv(sender_summary,'06_nichenet_expression_audit/sender_ligand_expression_summary.tsv')
write_tsv(receiver_summary,'06_nichenet_expression_audit/receiver_receptor_expression_summary.tsv')
write_tsv(rbindlist(sample_status),'06_nichenet_expression_audit/expression_audit_sample_status.tsv')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo_expression_audit.txt'))
writeLines(c('# NicheNet expression audit','', '- All 59 biological samples are read independently from the frozen RNA/data inputs; no pooling or resampling.', '- Sender groups: Fibro_C0–C6; receiver groups: all seven major populations. CellChat will determine follow-up priority.', '- Ligand/receptor universe is the verified official human NicheNet v2 LR resource. Missing RNA feature -> NA, not assumed zero. Absent/<30-cell groups -> NA. Evaluable but unexpressed genes -> true zero fraction and mean.', '- Detection means RNA/data>0; mean is arithmetic mean of existing log-normalized RNA/data values. 10% expression is an explicit diagnostic flag and retains sample prevalence, not a hard biological significance cutoff.', '- Cohort expression medians include eligible observed zeros. Sender/receptor summary n_eligible and n_evaluable are distinct for missing genes.', '- receiver_expressed_background_candidates.tsv.gz contains per-sample expressed target-universe genes for a future background. It is NOT geneset_oi, and does not freeze the eventual program/background replication design.', '- No ligand activity, DE test, receiver program selection, score input, or annotation/reduction mutation was performed.'),file.path(out,'NicheNet_expression_audit.md'))
if(requireNamespace('nichenetr',quietly=TRUE) && file.exists(file.path(ROOT,'07_nichenet_templates/receiver_program_design_template.tsv')))writeLines('READY',file.path(NN_OUT,'NICHENET_STATUS.txt'))
cat('NICHENET EXPRESSION AUDIT COMPLETE\n')
