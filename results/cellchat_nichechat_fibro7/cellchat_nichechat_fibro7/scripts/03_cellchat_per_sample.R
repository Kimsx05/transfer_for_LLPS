source('/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/scripts/common.R')
suppressPackageStartupMessages(library(CellChat))
stopifnot(packageVersion('CellChat')=='2.1.2')
cfg<-readRDS(file.path(ROOT,'01_metadata/frozen_analysis_inputs.rds'))
manifest<-fread(file.path(ROOT,'01_metadata/per_sample_input_manifest.tsv'))
DB<-CellChatDB.human
categories<-c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact')
DB$interaction<-DB$interaction[DB$interaction$annotation%in%categories,,drop=FALSE]
stopifnot(setequal(unique(DB$interaction$annotation),categories))
db<-as.data.table(DB$interaction,keep.rownames='interaction_id')
if(!file.exists(file.path(ROOT,'00_audit/CellChatDB_human_primary.rds'))) {
write_tsv(db,'00_audit/CellChatDB_human_primary_interactions.tsv')
write_tsv(as.data.table(DB$complex,keep.rownames='complex_id'),'00_audit/CellChatDB_human_complex_definitions.tsv')
write_tsv(db[,.N,by=annotation],'00_audit/CellChatDB_primary_category_counts.tsv')
saveRDS(DB,file.path(ROOT,'00_audit/CellChatDB_human_primary.rds'))
}
parameters<-list(min.cells=30L,type='truncatedMean',trim=.1,population.size=FALSE,nboot=100L,seed.use=123L,raw.use=TRUE,significance_threshold=.05,categories=categories,identifyOverExpressedGenes_do.fast=FALSE,workers_within_sample=1L)
if(!file.exists(file.path(ROOT,'00_audit/CellChat_exact_parameters.rds'))) {
api<-c('createCellChat','subsetData','identifyOverExpressedGenes','identifyOverExpressedInteractions','computeCommunProb','filterCommunication','computeCommunProbPathway','aggregateNet')
capture.output(lapply(api,function(n)list(function_name=n,formals=formals(get(n,asNamespace('CellChat'))))),file=file.path(ROOT,'00_audit/CellChat_actual_API.txt'))
saveRDS(parameters,file.path(ROOT,'00_audit/CellChat_exact_parameters.rds'))
writeLines(c('# Exact CellChat method','',paste0('CellChat 2.1.2; release v2.1.2; commit a0d3b2d231d46c8787177fffeac908270c253747. R ',getRversion(),'.'),'- Human database restricted to Secreted Signaling, ECM-Receptor and Cell-Cell Contact; no category is omitted.','- One object per orig.ident biological sample. All included cells of present groups enter inference; group levels absent from that sample are dropped inside CellChat and restored to the 14x14 output afterward.','- createCellChat(expression, meta, group.by="cellchat_group", datatype="RNA"); subsetData; identifyOverExpressedGenes(do.fast=FALSE); identifyOverExpressedInteractions; computeCommunProb(type="truncatedMean", trim=0.1, raw.use=TRUE, population.size=FALSE, nboot=100, seed.use=123); filterCommunication(min.cells=30); computeCommunProbPathway(thresh=0.05); aggregateNet(thresh=0.05, remove.isolate=FALSE).','- do.fast=FALSE uses the package Wilcoxon implementation without an optional presto dependency. Other gene-selection defaults remain unchanged (thresh.pc=0, thresh.fc=0, thresh.p=0.05).','- aggregateNet and LR retained weights use bootstrap p<0.05, matching installed API. Native pathway output uses p<=0.05 in computeCommunProbPathway; this API boundary difference is retained and disclosed. Raw candidate LR probabilities/p-values are also saved.','- No projection, downsampling, pooling, stage comparison, score input, or upstream mutation. Global and per-sample seed reset to 123; sequential future plan within each sample.','- Full 14x14 outputs: source/target <30 -> NA; successful eligible undetected -> 0; failed sample -> NA with FAILED status. Failed samples never supply zeros.','- The only adaptive failure handling is catch-and-record per sample; no parameter changes or invented empty objects.', '', 'Executable source: scripts/03_cellchat_per_sample.R'),file.path(ROOT,'00_audit/CellChat_method_and_actual_code.md'))
audit<-fread(file.path(ROOT,'00_audit/preflight_audit.tsv'));audit[metric=='CellChat_version',value:=as.character(packageVersion('CellChat'))];write_tsv(audit,'00_audit/preflight_audit.tsv')
}

run_one<-function(i) {
  set.seed(123);future::plan(future::sequential)
  row<-manifest[i];s<-row$sample;dir<-dirname(row$input_path);status_path<-file.path(dir,'status.tsv')
  if(file.exists(status_path)) {
    previous<-fread(status_path)
    if(previous$status=='SUCCESS' && file.exists(previous$object_path))return(previous)
    # A rerun is explicit: archive any previous attempt within the new result tree.
    stop('Prior non-success attempt exists for ',s,'; use a separately named retry script/directory')
  }
  started<-Sys.time();error<-'';warnings<-character();object_path<-file.path(dir,'CellChat_object.rds')
  elig<-copy(cfg$eligibility[sample==s]);smcounts<-cfg$counts[sample==s]
  ng<-sum(smcounts$n_cells>=30)
  log_con<-file(file.path(dir,'CellChat.log'),open='wt');sink(log_con);sink(log_con,type='message')
  ok<-tryCatch(withCallingHandlers({
    input<-qs2::qs_read(row$input_path,nthreads=1)
    sm<-as.data.frame(input$metadata);rownames(sm)<-sm$cell
    sm$samples<-sm$sample
    sm$cellchat_group<-droplevels(factor(sm$cellchat_group,levels=GROUPS))
    stopifnot(ncol(input$expression)==nrow(sm),identical(colnames(input$expression),sm$cell),length(unique(sm$sample))==1L,unique(sm$sample)==s)
    if(ng==0L)stop('No group has >=30 cells; no eligible edges')
    if(nlevels(sm$cellchat_group)<2L)stop('CellChat overexpression inference needs at least two present groups')
    cc<-createCellChat(object=input$expression,meta=sm,group.by='cellchat_group',datatype='RNA')
    cc@DB<-DB
    cc<-subsetData(cc)
    set.seed(123);cc<-identifyOverExpressedGenes(cc,do.fast=FALSE)
    cc<-identifyOverExpressedInteractions(cc)
    if(nrow(cc@LR$LRsig)==0L)stop('No overexpressed LR candidates; native inference unavailable')
    cc<-computeCommunProb(cc,type='truncatedMean',trim=.1,population.size=FALSE,nboot=100,seed.use=123L,raw.use=TRUE)
    cc<-filterCommunication(cc,min.cells=30)
    cc<-computeCommunProbPathway(cc,thresh=.05)
    cc<-aggregateNet(cc,thresh=.05,remove.isolate=FALSE)
    stopifnot(!is.null(cc@net$weight),all(is.finite(cc@net$weight)),all(cc@net$weight>=0))
    cc@options$cohort_analysis<-list(sample=s,seed=123,groups=GROUPS,min.cells=30,parameters=parameters)
    saveRDS(cc,object_path,compress=FALSE)
    full_w<-blank_matrix();full_c<-blank_matrix()
    ij<-cbind(match(elig$source,GROUPS),match(elig$target,GROUPS))
    full_w[ij[elig$eligible,,drop=FALSE]]<-0;full_c[ij[elig$eligible,,drop=FALSE]]<-0
    present<-rownames(cc@net$weight)
    for(src in present)for(tgt in present)if(elig[source==src & target==tgt,eligible]) {
      full_w[src,tgt]<-cc@net$weight[src,tgt];full_c[src,tgt]<-cc@net$count[src,tgt]
    }
    edge<-merge(elig,row[,.(sample,patient_uid,dataset)],by='sample',sort=FALSE)
    edge[,inference_status:='SUCCESS'];edge[,weight:=full_w[cbind(match(source,GROUPS),match(target,GROUPS))]];edge[,count:=full_c[cbind(match(source,GROUPS),match(target,GROUPS))]]
    stopifnot(all(is.na(edge[eligible==FALSE,weight])),all(!is.na(edge[eligible==TRUE,weight])))
    fwrite(edge,file.path(dir,'network_long.tsv'),sep='\t',na='NA')
    saveRDS(list(weight=full_w,count=full_c),file.path(dir,'network_14x14.rds'))
    fwrite(data.table(source=GROUPS,as.data.table(full_w)),file.path(dir,'net_weight_14x14.tsv'),sep='\t',na='NA')
    fwrite(data.table(source=GROUPS,as.data.table(full_c)),file.path(dir,'net_count_14x14.tsv'),sep='\t',na='NA')
    # Candidate LR arrays retain non-significant entries and native bootstrap p-values.
    lr<-as.data.table(as.table(cc@net$prob));setnames(lr,c('source','target','interaction_id','raw_probability'))
    lr[,pval:=as.vector(cc@net$pval)];lr[,c('source','target','interaction_id'):=lapply(.SD,as.character),.SDcols=c('source','target','interaction_id')]
    lr<-merge(lr,elig[,.(source,target,eligible)],by=c('source','target'),sort=FALSE)
    lr<-merge(lr,db[,.(interaction_id,ligand,receptor,pathway_name,annotation)],by='interaction_id',sort=FALSE)
    lr[,weight:=ifelse(eligible,ifelse(pval<.05,raw_probability,0),NA_real_)]
    lr[,detected:=eligible & !is.na(weight) & weight>0];lr[,`:=`(sample=s,patient_uid=row$patient_uid,dataset=row$dataset)]
    fwrite(lr,file.path(dir,'lr_candidates.tsv.gz'),sep='\t',na='NA',compress='gzip')
    fwrite(lr[detected==TRUE],file.path(dir,'lr_detected.tsv'),sep='\t',na='NA')
    p<-cc@netP$prob
    if(length(p)) {
      pw<-as.data.table(as.table(p));setnames(pw,c('source','target','pathway_name','weight'));pw[,c('source','target','pathway_name'):=lapply(.SD,as.character),.SDcols=c('source','target','pathway_name')]
      pw<-merge(pw,elig[,.(source,target,eligible)],by=c('source','target'),sort=FALSE);pw[eligible==FALSE,weight:=NA_real_]
      pw<-pw[eligible & weight>0];pw[,`:=`(sample=s,patient_uid=row$patient_uid,dataset=row$dataset)]
    } else pw<-data.table(source=character(),target=character(),pathway_name=character(),weight=numeric(),eligible=logical(),sample=character(),patient_uid=character(),dataset=character())
    fwrite(pw,file.path(dir,'pathway_detected.tsv'),sep='\t',na='NA')
    TRUE
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){error<<-conditionMessage(e);FALSE})
  sink(type='message');sink();close(log_con)
  if(!ok) {
    edge<-merge(elig,row[,.(sample,patient_uid,dataset)],by='sample',sort=FALSE)
    edge[,`:=`(inference_status='FAILED',weight=NA_real_,count=NA_real_)]
    fwrite(edge,file.path(dir,'network_long.tsv'),sep='\t',na='NA')
    saveRDS(list(weight=blank_matrix(),count=blank_matrix()),file.path(dir,'network_14x14.rds'))
  }
  status<-data.table(sample=s,sample_id=row$sample_id,patient_uid=row$patient_uid,dataset=row$dataset,status=if(ok)'SUCCESS' else 'FAILED',error_message=error,runtime_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),cell_number=row$included_cells,eligible_groups=ng,object_path=if(ok)object_path else NA_character_,warnings=paste(unique(warnings),collapse=' | '))
  fwrite(status,status_path,sep='\t',na='NA');cat('FINISHED ',s,' ',status$status,' seconds=',round(status$runtime_seconds,1),'\n',sep='');flush.console()
  status
}
args<-commandArgs(trailingOnly=TRUE)
indices<-if(length(args))as.integer(args) else seq_len(nrow(manifest))
stopifnot(all(indices%in%seq_len(nrow(manifest))))
# Separate R processes may execute different manifest rows; each uses one sample.
for(i in indices)run_one(i)
status_files<-file.path(dirname(manifest$input_path),'status.tsv');status_files<-status_files[file.exists(status_files)]
if(length(indices)==nrow(manifest) && length(status_files))write_tsv(rbindlist(lapply(status_files,fread),fill=TRUE),'02_cellchat_per_sample/sample_run_status.tsv')
if(!file.exists(file.path(ROOT,'00_audit/sessionInfo_CellChat_inference.txt')))capture.output(sessionInfo(),file=file.path(ROOT,'00_audit/sessionInfo_CellChat_inference.txt'))
