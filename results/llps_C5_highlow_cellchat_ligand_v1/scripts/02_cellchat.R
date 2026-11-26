source('/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1/scripts/common.R')
args<-commandArgs(TRUE);worker<-as.integer(args[1]);st<-fread(file.path(OUT,'01_LLPS_stratification/cell_assignments.csv.gz'));ss<-fread(file.path(OUT,'01_LLPS_stratification/specimen_group_counts.csv'));man<-fread(file.path(OUT,'00_input_audit/manifest.csv'));jobs<-merge(ss,man[,.(specimen_id,input_path,run_id)],by='specimen_id');jobs[,job_id:=paste(method,run_id,sep='__')];jobs[,runnable:=high>5&low>5];writecsv(jobs,file.path(OUT,'02_CellChat_highlow/job_manifest.csv'));jobs<-jobs[runnable==TRUE];jobs<-jobs[(seq_len(.N)-1)%%3==worker];db<-as.data.table(CellChatDB.human$interaction,keep.rownames='feature');identity<-fread(file.path(OUT,'00_input_audit/stromal_identity_final.csv.gz'));partners<-c(GROUPS[8:14],'Other_fibroblasts','Pericyte');groups<-c('C5_high','C5_low',partners)
for(i in seq_len(nrow(jobs))){j<-jobs[i];dest<-file.path(OUT,'02_CellChat_highlow',j$job_id);dir.create(dest,showWarnings=FALSE);if(file.exists(file.path(dest,'COMPLETE.csv')))next
tryCatch({cat('START',j$job_id,'\n');flush.console();cache<-file.path(dest,'unfiltered.rds')
if(!file.exists(cache)){
 same<-fread(file.path(OUT,'00_input_audit/identical_assignment_audit.csv'));ids<-strsplit(same[job_id==j$job_id,identical_jobs],';',fixed=TRUE)[[1]];sources<-file.path(OUT,'02_CellChat_highlow',setdiff(ids,j$job_id),'unfiltered.rds');sources<-sources[file.exists(sources)];if(length(sources)){stopifnot(file.copy(sources[1],paste0(cache,'.tmp')),file.rename(paste0(cache,'.tmp'),cache));writeLines(paste('Exact barcode/state identity; source',sources[1]),file.path(dest,'identical_assignment_reuse.txt'))}
}
if(!file.exists(cache)){
 inp<-qs2::qs_read(j$input_path,nthreads=1);md<-as.data.table(inp$metadata);if(!'cell'%in%names(md))stop('Missing cell barcode');stopifnot(identical(md$cell,colnames(inp$expression)));z<-st[method==j$method&specimen_id==j$specimen_id];orig<-as.character(md$cellchat_group);stopifnot(setequal(md$cell[orig=='Fibro_C5'],z$cell));labels<-orig;fi<-which(orig%in%paste0('Fibro_C',c(0:4,6)));labels[fi]<-identity$analysis_group[match(md$cell[fi],identity$cell)];stopifnot(!anyNA(labels[fi]));ix<-which(orig=='Fibro_C5');state<-z$state[match(md$cell[ix],z$cell)];labels[ix]<-ifelse(is.na(state),NA,paste0('C5_',state));keep<-!is.na(labels);meta<-data.frame(group=factor(labels[keep],levels=groups),row.names=md$cell[keep]);meta$group<-droplevels(meta$group);counts<-table(factor(labels[keep],levels=groups));set.seed(123)
 cc<-createCellChat(inp$expression[,keep,drop=FALSE],meta=meta,group.by='group');cc@DB<-CellChatDB.human;cc<-subsetData(cc);cc<-identifyOverExpressedGenes(cc,do.fast=FALSE);cc<-identifyOverExpressedInteractions(cc);cc<-computeCommunProb(cc,type='truncatedMean',trim=.1,population.size=FALSE,nboot=100,seed.use=123L,raw.use=TRUE)
 u<-list(net=cc@net,LR=cc@LR,idents=cc@idents,meta=cc@meta,options=cc@options,counts=counts,job=j);atomicRDS(u,cache);rm(inp,cc);gc(FALSE)
}else u<-readRDS(cache)
for(th in c(5,10,20,30)){
 cc<-new('CellChat');cc@net<-u$net;cc@LR<-u$LR;cc@idents<-u$idents;cc@meta<-u$meta;cc@options<-u$options;cc@DB<-CellChatDB.human;cc<-filterCommunication(cc,min.cells=th);cc<-computeCommunProbPathway(cc,thresh=.05);cc<-aggregateNet(cc,thresh=.05,remove.isolate=FALSE)
 out<-list();k<-0L
 for(part in partners)for(direction in c('outgoing','incoming')){
 eligible<-u$counts['C5_high']>th&&u$counts['C5_low']>th&&u$counts[part]>th
 for(state in c('high','low')){focal<-paste0('C5_',state);se<-if(direction=='outgoing')focal else part;re<-if(direction=='outgoing')part else focal
 # Database-wide common LR universe: non-candidate entries are undetected zeros, not invented p-values.
 ids<-db$feature;v<-pv<-rep(NA_real_,length(ids));if(eligible){v[]<-0;hit<-match(dimnames(u$net$prob)[[3]],ids);v[hit]<-u$net$prob[se,re,];pv[hit]<-u$net$pval[se,re,];v[hit]<-ifelse(pv[hit]<.05,v[hit],0)}
 k<-k+1L;out[[k]]<-data.table(kind='LR',feature=ids,state=state,direction=direction,partner=part,eligible=eligible,strength=v,cellchat_p=pv)
 pnames<-unique(db$pathway_name);v<-rep(if(eligible)0 else NA_real_,length(pnames));if(eligible&&length(dim(cc@netP$prob))==3){v[match(dimnames(cc@netP$prob)[[3]],pnames)]<-cc@netP$prob[se,re,]};k<-k+1L;out[[k]]<-data.table(kind='pathway',feature=pnames,state=state,direction=direction,partner=part,eligible=eligible,strength=v,cellchat_p=NA_real_)
 k<-k+1L;out[[k]]<-data.table(kind=c('network','interaction_count'),feature=c('total','detected_LR'),state=state,direction=direction,partner=part,eligible=eligible,strength=if(eligible)c(cc@net$weight[se,re],cc@net$count[se,re])else c(NA_real_,NA_real_),cellchat_p=NA_real_)
 }}
 tab<-rbindlist(out);tab[,`:=`(method=j$method,specimen_id=j$specimen_id,patient_id=j$patient_id,dataset=j$dataset,clinical_group=j$clinical_group,min_cells=th,high_cells=u$counts['C5_high'],low_cells=u$counts['C5_low'])];tab[,partner_cells:=as.numeric(u$counts[partner])];fwrite(tab,file.path(dest,paste0('threshold',th,'.csv.gz')),na='NA',compress='gzip');atomicRDS(list(net=cc@net,netP=cc@netP,LR=cc@LR,counts=u$counts,job=j,min_cells=th),file.path(dest,paste0('threshold',th,'_net.rds')))
}
writecsv(data.table(job_id=j$job_id,status='SUCCESS'),file.path(dest,'COMPLETE.csv'));cat('DONE',j$job_id,'\n');flush.console();rm(u,cc,tab,out);gc(FALSE)
},error=function(e){writecsv(data.table(job_id=j$job_id,error=conditionMessage(e)),file.path(dest,'ERROR.csv'));cat('FAILED',j$job_id,conditionMessage(e),'\n');flush.console()})
}
