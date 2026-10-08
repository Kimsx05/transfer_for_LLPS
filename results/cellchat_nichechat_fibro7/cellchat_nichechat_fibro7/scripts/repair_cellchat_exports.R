# Recover export-only failures from complete native CellChat objects; no inference rerun.
source('/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/scripts/common.R')
suppressPackageStartupMessages(library(CellChat))
cfg<-readRDS(file.path(ROOT,'01_metadata/frozen_analysis_inputs.rds'))
manifest<-fread(file.path(ROOT,'01_metadata/per_sample_input_manifest.tsv'))
db<-fread(file.path(ROOT,'00_audit/CellChatDB_human_primary_interactions.tsv'))
code<-readLines(file.path(ROOT,'scripts/03_cellchat_per_sample.R'))
a<-which(grepl('    full_w<-blank_matrix',code,fixed=TRUE));b<-which(grepl('    TRUE',code,fixed=TRUE))
stopifnot(length(a)==1L,length(b)==1L,a<b)
export_code<-parse(text=code[a:(b-1L)])
for(i in seq_len(nrow(manifest))) {
  row<-manifest[i];s<-row$sample;dir<-dirname(row$input_path)
  status_path<-file.path(dir,'status.tsv');object_path<-file.path(dir,'CellChat_object.rds')
  if(!file.exists(status_path)||!file.exists(object_path))next
  old<-fread(status_path)
  if(old$status=='SUCCESS')next
  if(!grepl("eligible.*calling scope",old$error_message))next
  started<-Sys.time();cc<-readRDS(object_path)
  # Verify a complete native inference with the exact frozen analysis parameters.
  p<-cc@options$parameter
  stopifnot(p$type.mean=='truncatedMean',p$trim==.1,identical(p$population.size,FALSE),p$nboot==100,p$seed.use==123,
            identical(cc@options$cohort_analysis$sample,s),cc@options$cohort_analysis$min.cells==30,
            identical(colnames(cc@data),cc@meta$cell),all(cc@meta$sample==s),all(is.finite(cc@net$weight)),
            all(c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')%in%unique(cc@DB$interaction$annotation)))
  elig<-copy(cfg$eligibility[sample==s])
  archive<-file.path(dir,'export_failure_before_repair');dir.create(archive,showWarnings=FALSE)
  if(!file.exists(file.path(archive,'status.tsv')))file.copy(status_path,file.path(archive,'status.tsv'))
  eval(export_code)
  # LR retention must reproduce native aggregateNet weight/count exactly.
  ll<-fread(file.path(dir,'lr_detected.tsv'))
  totals<-ll[,.(lr_weight=sum(weight),lr_count=.N),by=.(source,target)]
  for(k in seq_len(nrow(totals)))stopifnot(abs(full_w[totals$source[k],totals$target[k]]-totals$lr_weight[k])<1e-10,full_c[totals$source[k],totals$target[k]]==totals$lr_count[k])
  recovered_object_path<-object_path
  old[,`:=`(status='SUCCESS',error_message='',object_path=recovered_object_path,export_repaired=TRUE,original_inference_runtime_seconds=runtime_seconds,export_repair_runtime_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))]
  fwrite(old,status_path,sep='\t',na='NA')
  writeLines(c('Complete native inference object validated. Export-only data.table compatibility failure repaired.','No createCellChat/computeCommunProb call, reseeding or recomputation was performed.','Exact CellChat parameters, biological sample identity, DB categories and LR-to-net weight/count reconciliation passed.'),file.path(dir,'export_repair_audit.txt'))
  cat('EXPORT REPAIRED ',s,'\n',sep='');flush.console();rm(cc);gc(FALSE)
}
statuses<-rbindlist(lapply(file.path(dirname(manifest$input_path),'status.tsv'),function(p)if(file.exists(p))fread(p) else NULL),fill=TRUE)
write_tsv(statuses,'02_cellchat_per_sample/sample_run_status.tsv')
