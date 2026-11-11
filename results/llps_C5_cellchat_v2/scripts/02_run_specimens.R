source('/home/jinshengxi/bladder_metabolism/results/llps_C5_cellchat_v2/scripts/common.R')
cfg<-readRDS(file.path(OUT,'00_audit/config.rds'));manifest<-fread(file.path(OUT,'00_audit/manifest.csv'));cnt<-fread(file.path(AUD,'sample_celltype_counts.csv'));db<-as.data.table(cfg$DB$interaction,keep.rownames='interaction_id')
argv<-commandArgs(TRUE);export_only<-'--export-only'%in%argv;args<-as.integer(argv[argv!='--export-only']);indices<-seq_len(nrow(manifest));if(length(args))indices<-indices[(indices-1L)%%3L==args[1]]
basecache<-file.path(OUT,'00_audit/unfiltered');dir.create(basecache,showWarnings=FALSE)
if(export_only)indices<-indices[file.exists(file.path(basecache,paste0(manifest$run_id[indices],'_status.csv'))) & file.exists(file.path(basecache,paste0(manifest$run_id[indices],'.rds')))];for(i in indices){row<-manifest[i];s<-row$specimen_id;cache<-file.path(basecache,paste0(row$run_id,'.rds'));statusfile<-file.path(basecache,paste0(row$run_id,'_status.csv'));if(file.exists(statusfile)&&all(file.exists(vapply(cfg$thresholds,function(t)file.path(TD(t),row$run_id,'SUCCESS.csv'),character(1))))&&fread(statusfile)$status=='SUCCESS')next;started<-Sys.time();if(file.exists(statusfile)){prev<-fread(statusfile);if(prev$status=='FAILED'){archive<-paste0(statusfile,'.attempt_',format(Sys.time(),'%Y%m%d%H%M%S'));file.copy(statusfile,archive,overwrite=FALSE)}};cat('START ',s,'\n');flush.console()
 tryCatch({
 if(!file.exists(cache)){
 input<-qs2::qs_read(row$input_path,nthreads=1);sm<-as.data.frame(input$metadata);rownames(sm)<-sm$cell;sm$samples<-factor(s);sm$specimen_id<-s;sm$patient_id<-row$patient_id;sm$clinical_group<-row$clinical_group;sm$dataset<-row$dataset;sm$cellchat_group<-droplevels(factor(sm$cellchat_group,levels=GROUPS));stopifnot(identical(colnames(input$expression),sm$cell))
 set.seed(123);cc<-createCellChat(input$expression,meta=sm,group.by='cellchat_group',datatype='RNA');cc@DB<-cfg$DB;cc<-subsetData(cc);cc<-identifyOverExpressedGenes(cc,do.fast=FALSE);cc<-identifyOverExpressedInteractions(cc)
 if(nrow(cc@LR$LRsig)==0)stop('No overexpressed LR candidates')
 cc<-computeCommunProb(cc,type='truncatedMean',trim=.1,population.size=FALSE,nboot=100,seed.use=123L,raw.use=TRUE)
 saved<-list(net=cc@net,LR=cc@LR,idents=cc@idents,meta=cc@meta,options=cc@options,specimen_id=s,parameters=cfg$parameters,DB_n=nrow(cfg$DB$interaction));atomicRDS(saved,cache);rm(cc,input);gc(FALSE)
 }else{saved<-readRDS(cache);stopifnot(saved$specimen_id==s,identical(saved$parameters,cfg$parameters),saved$DB_n==3233)}
 if(!'patient_id'%in%names(saved$meta)){saved$meta$patient_id<-row$patient_id;saved$meta$specimen_id<-s;saved$meta$clinical_group<-row$clinical_group;saved$meta$dataset<-row$dataset;atomicRDS(saved,cache)}
 counts<-setNames(as.numeric(cnt[specimen_id==s,c(paste0('C',0:6),GROUPS[8:14]),with=FALSE]),GROUPS)
 for(t in cfg$thresholds){dest<-file.path(TD(t),row$run_id);dir.create(dest,recursive=TRUE,showWarnings=FALSE);if(file.exists(file.path(dest,'SUCCESS.csv')))next
 cc<-new('CellChat');cc@net<-saved$net;cc@LR<-saved$LR;cc@idents<-saved$idents;cc@meta<-saved$meta;cc@options<-saved$options;cc@DB<-cfg$DB
 cc<-filterCommunication(cc,min.cells=t);cc<-computeCommunProbPathway(cc,thresh=.05);cc<-aggregateNet(cc,thresh=.05,remove.isolate=FALSE)
 edge<-CJ(sender=GROUPS,receiver=GROUPS);edge[,`:=`(sender_cells=counts[sender],receiver_cells=counts[receiver])];edge[,eligible:=sender_cells>t & receiver_cells>t];edge[,`:=`(weight=NA_real_,LR_count=NA_real_)];ok<-which(edge$eligible);ij<-cbind(match(edge$sender[ok],rownames(cc@net$weight)),match(edge$receiver[ok],colnames(cc@net$weight)));edge[ok,`:=`(weight=cc@net$weight[ij],LR_count=cc@net$count[ij])]
 edge[,`:=`(specimen_id=s,patient_id=row$patient_id,dataset=row$dataset,clinical_group=row$clinical_group,min_cells=t,C5_cells=counts['Fibro_C5'],status='SUCCESS')];writecsv(edge,file.path(dest,'network_all_14_groups.csv'))
 # Complete candidate arrays preserve native probability/p-values; eligibility separately masks excluded entries.
 atomicRDS(list(net=cc@net,netP=cc@netP,LR=cc@LR,eligible_groups=counts>t,counts=counts,metadata=row,threshold=t),file.path(dest,'compact_networks.rds'))
 lr<-as.data.table(as.table(saved$net$prob));setnames(lr,c('sender','receiver','interaction_id','raw_probability'));lr[,p_value:=as.vector(saved$net$pval)];lr[,c('sender','receiver','interaction_id'):=lapply(.SD,as.character),.SDcols=c('sender','receiver','interaction_id')]
 lr<-merge(lr,db[,.(interaction_id,ligand,receptor,pathway=pathway_name,annotation)],by='interaction_id');lr<-merge(lr,edge[,.(sender,receiver,sender_cells,receiver_cells,eligible)],by=c('sender','receiver'))
 lr[,`:=`(probability=ifelse(eligible,raw_probability,NA_real_),cellchat_p_value=ifelse(eligible,p_value,NA_real_),weight=ifelse(eligible,ifelse(p_value<.05,raw_probability,0),NA_real_))];lr[,detected:=!is.na(weight)&weight>0];lr[,p_BH_within_specimen_direction:=if(all(!eligible)) NA_real_ else p.adjust(p_value,'BH'),by=.(sender,receiver)]
 lr[,`:=`(specimen_id=s,patient_id=row$patient_id,dataset=row$dataset,clinical_group=row$clinical_group,min_cells=t,C5_cells=counts['Fibro_C5'])]
 writecsv(lr[detected==TRUE],file.path(dest,'LR_detected_all_groups.csv'));fwrite(lr[sender=='Fibro_C5'|receiver=='Fibro_C5'],file.path(dest,'C5_LR_candidates.csv.gz'),na='NA',compress='gzip')
 # LR native aggregate boundary p<.05; pathway API p<=.05, deliberately retain package semantics.
 pw<-as.data.table(as.table(cc@netP$prob));if(ncol(pw)==4){setnames(pw,c('sender','receiver','pathway','weight'));pw[,c('sender','receiver','pathway'):=lapply(.SD,as.character),.SDcols=c('sender','receiver','pathway')];pw<-merge(pw,edge[,.(sender,receiver,eligible)],by=c('sender','receiver'));pw[eligible==FALSE,weight:=NA_real_];pw<-pw[eligible&weight>0]}else pw<-data.table(sender=character(),receiver=character(),pathway=character(),weight=numeric(),eligible=logical())
 pw[,`:=`(specimen_id=rep(s,.N),patient_id=rep(row$patient_id,.N),dataset=rep(row$dataset,.N),clinical_group=rep(row$clinical_group,.N),min_cells=rep(t,.N))];writecsv(pw,file.path(dest,'pathway_detected_all_groups.csv'))
 recon<-lr[detected==TRUE,.(w=sum(weight),n=.N),by=.(sender,receiver)];check<-merge(recon,edge,by=c('sender','receiver'));stopifnot(all(abs(check$w-check$weight)<1e-9),all(check$n==check$LR_count),all(is.na(edge[eligible==FALSE,weight])))
 writecsv(data.table(specimen_id=s,min_cells=t,status='SUCCESS',eligible_C5=counts['Fibro_C5']>t,n_candidates=length(dimnames(saved$net$prob)[[3]]),runtime_seconds=as.numeric(difftime(Sys.time(),started,units='secs'))),file.path(dest,'SUCCESS.csv'));rm(cc,lr,pw);gc(FALSE)
 }
 writecsv(data.table(specimen_id=s,status='SUCCESS',error='',runtime_seconds=as.numeric(difftime(Sys.time(),started,units='secs'))),statusfile);cat('DONE ',s,'\n');flush.console()
 },error=function(e){writecsv(data.table(specimen_id=s,status='FAILED',error=conditionMessage(e),runtime_seconds=as.numeric(difftime(Sys.time(),started,units='secs'))),statusfile);cat('FAILED ',s,': ',conditionMessage(e),'\n');flush.console()})
 rm(list=intersect(c('saved','cc','input','lr','pw'),ls()));gc(FALSE)
}
