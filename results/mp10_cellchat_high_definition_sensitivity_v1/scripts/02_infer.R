source('mp10_cellchat_high_definition_sensitivity_v1/scripts/common.R')
suppressPackageStartupMessages(library(CellChat));future::plan('sequential');options(warn=1)
cfg<-fromJSON(file.path(root,'config/run_config.json'));stopifnot(as.character(packageVersion('CellChat'))=='1.6.1')
DB<-CellChatDB.human;stopifnot(nrow(DB$interaction)==1939)
writeLines(capture.output(sessionInfo()),file.path(root,'logs/inference_sessionInfo.txt'))
for(fn in c('computeCommunProb','filterCommunication','computeCommunProbPathway','aggregateNet','computeExpr_LR','computeExpr_complex','computeExpr_coreceptor','computeExpr_agonist','computeExpr_antagonist','extractGene'))writeLines(capture.output(get(fn,asNamespace('CellChat'))),file.path(root,'config',paste0('CellChat_1.6.1_',fn,'.txt')))
x<-qread(cfg$expression,nthreads=1);stopifnot(all(is.finite(x@x)),all(x@x>=0),!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)))
full<-fread(file.path(project,'subtype/v2/final_annotation_pre_reference/05_gene_symbol_audit_all_genes.csv'))$original_symbol
oldcov<-fread(file.path(project,'mp10_subtype_screen_round1/tables/CellChat_gene_coverage.tsv'))
required<-unique(c(extractGene(DB),unlist(DB$complex),unlist(DB$cofactor)));required<-required[!is.na(required)&required!='']; missing_cache<-intersect(setdiff(required,rownames(x)),full)
# Only if required genes were omitted from cache, read original same-normalization RNA:data once; never renormalize.
if(length(missing_cache)){
 suppressPackageStartupMessages(library(Seurat));obj<-qread(file.path(project,'subtype/v2/final_annotation_pre_reference/final_annotation_v2_pre_reference.qs'),nthreads=1)
 full<-rownames(obj[['RNA']]);add<-intersect(setdiff(required,rownames(x)),full)
 if(length(add)){xa<-LayerData(obj,assay='RNA',layer='data',features=add);x<-rbind(x,xa[,colnames(x),drop=FALSE])};rm(obj);gc()
}
gcov<-data.table(gene=sort(unique(c(full,required))),in_full_feature_universe=sort(unique(c(full,required)))%in%full)
gcov[,in_expression_cache:=gene%in%rownames(x)];gzero<-setNames(Matrix::rowSums(x!=0)==0,rownames(x));gcov[,all_cells_zero:=unname(gzero[gene])];gcov[,required_by_DB:=gene%in%required];put(gcov,'gene_universe_coverage.tsv.gz')
expand<-function(n,kind){if(is.na(n)||n=='')return(character());if(kind=='LR'){if(n%in%rownames(DB$complex)){v<-unlist(DB$complex[n,,drop=FALSE],use.names=FALSE);return(v[v!=''&!is.na(v)])};return(n)};stopifnot(n%in%rownames(DB$cofactor));v<-unlist(DB$cofactor[n,,drop=FALSE],use.names=FALSE);v[v!=''&!is.na(v)]}
records<-list();for(i in seq_len(nrow(DB$interaction))){z<-DB$interaction[i,];req<-unique(c(expand(z$ligand,'LR'),expand(z$receptor,'LR')));aux<-unique(unlist(lapply(as.character(z[1,c('agonist','antagonist','co_A_receptor','co_I_receptor')]),expand,kind='aux')));mis<-setdiff(req,rownames(x));ma<-setdiff(aux,rownames(x));records[[i]]<-data.table(interaction_id=rownames(DB$interaction)[i],required_genes=paste(req,collapse=';'),auxiliary_genes=paste(aux,collapse=';'),missing_required=paste(mis,collapse=';'),missing_auxiliary=paste(ma,collapse=';'),computable=length(mis)==0,reason=if(length(mis))'required_gene_absent_from_original_RNA'else if(length(ma))'computable_native_1.6.1_auxiliary_intersection_omits_absent_genes'else'complete_gene_coverage')}
audit<-rbindlist(records);put(audit,'LR_coverage_audit.tsv');LR<-DB$interaction[audit$computable,,drop=FALSE];old<-fread(file.path(baseline,'tables/cellchat_lr_domain.tsv'));oldids<-old[scope=='ALL'&network_level=='fine',interaction_id];stopifnot(length(oldids)==1875)
domain<-as.data.table(LR,keep.rownames='interaction_id');domain[,`:=`(LR_order=seq_len(.N),in_legacy_1875=interaction_id%in%oldids)];put(domain,'fixed_LR_domain.tsv');put(data.table(interaction_id=oldids,in_fixed=oldids%in%rownames(LR)),'legacy_1875_domain.tsv')
for(n in c('interaction','complex','cofactor','geneInfo'))put(as.data.table(DB[[n]],keep.rownames='definition_id'),paste0('DB_',n,'.tsv'))
qsave(DB,file.path(root,'objects/CellChatDB.human_1.6.1.qs'),nthreads=1)
provenance<-list(db_hash=digest(DB,algo='sha256'),LR_hash=digest(LR,algo='sha256'),gene_order_hash=seqhash(rownames(x)),actual_expression_content_sha256=digest(x,algo='sha256'),supplemented_expression_sha256=if(length(missing_cache))digest(x[intersect(missing_cache,rownames(x)),,drop=FALSE],algo='sha256')else NA_character_,original_RNA_file_sha256=digest(file=file.path(project,'subtype/v2/final_annotation_pre_reference/final_annotation_v2_pre_reference.qs'),algo='sha256'),source_expression_file_sha256=digest(file=cfg$expression,algo='sha256'),full_feature_hash=seqhash(full),n_features=nrow(x),n_DB=1939,n_computable=nrow(LR),n_legacy_intersection=sum(domain$in_legacy_1875),supplemented_genes=missing_cache,expression='existing nonnegative LogNormalize RNA:data; no normalization on signaling subset',auxiliary_handling='CellChat 1.6.1 native intersect with features; absent co-receptors default 1; missing required complex subunits excluded; zero expression retained')
if(length(missing_cache))qsave(x[intersect(missing_cache,rownames(x)),,drop=FALSE],file.path(root,'objects/supplemental_two_genes.qs'),nthreads=1)
write_json(provenance,file.path(root,'config/domain_provenance.json'),pretty=TRUE,auto_unbox=TRUE)
Rcpp::sourceCpp(file.path(root,'scripts/fast_trimean.cpp'),cacheDir=file.path(root,'objects/rcpp_cache'))
compute_fast<-CellChat::computeCommunProb;env<-new.env(parent=environment(compute_fast));env$aggregate<-function(x,by,FUN,...){if(!is.matrix(x)||length(by)!=1||!identical(FUN,CellChat::triMean))return(stats::aggregate(x,by,FUN,...));f<-droplevels(factor(by[[1]]));a<-grouped_trimean_cpp(x,as.integer(f),nlevels(f));colnames(a)<-colnames(x);cbind(data.frame(Group.1=factor(levels(f),levels=levels(f))),as.data.frame(a,check.names=FALSE))};environment(compute_fast)<-env
set.seed(123);test<-matrix(sample(c(rep(0,20),runif(20)),503*8,TRUE),503,8);grp<-rep(letters[1:9],length.out=503);err<-max(abs(as.matrix(env$aggregate(test,list(grp),CellChat::triMean)[,-1])-as.matrix(stats::aggregate(test,list(grp),CellChat::triMean)[,-1])));stopifnot(err<1e-12);put(data.table(check='accelerated_trimean_type7_zero_ties',max_abs_error=err,pass=TRUE),'trimean_validation.tsv')
p<-fread(tb('network_plan.tsv'));params<-cfg[c('CellChat_version','type','raw.use','population.size','nboot','seed.use','model_nodes')];codehash<-digest(file=file.path(root,'scripts/fast_trimean.cpp'),algo='sha256')
worker<-function(i){z<-p[i];lf<-file.path(root,'logs',paste0(z$network_id,'.log'));sink(lf);on.exit(sink());tryCatch({
 if(z$status=='NE')return(z)
 md<-fread(tb(z$membership_file));stopifnot(seqhash(md$barcode)==z$input_hash);md<-md[modeled==TRUE];stopifnot(seqhash(md$barcode)==z$modeled_hash,all(md$barcode%in%colnames(x)))
 key<-digest(list(barcodes=md$barcode,nodes=md$node,params=params,provenance=provenance,accelerator=codehash),algo='sha256');fp<-file.path(root,'objects/networks',paste0(z$network_id,'.qs'))
 if(file.exists(fp)){o<-qread(fp,nthreads=1);if(identical(o$cache_key,key)){z[,`:=`(status='completed_reused',cache_key=key,data_max=o$data_max,n_LR=nrow(LR))];return(z)}}
 cat('START',z$network_id,format(Sys.time()),nrow(md),'cells\n')
 cc<-createCellChat(x[,md$barcode,drop=FALSE],meta=data.frame(labels=factor(md$node),row.names=md$barcode),group.by='labels');cc@DB<-DB;cc<-subsetData(cc);cc@LR$LRsig<-LR
 stopifnot(all(unique(unlist(strsplit(audit[computable==TRUE,required_genes],';')))%in%rownames(cc@data.signaling)))
 dm<-max(cc@data.signaling);stopifnot(is.finite(dm),dm>0)
 cc<-compute_fast(cc,LR.use=LR,type='triMean',raw.use=TRUE,population.size=FALSE,nboot=100,seed.use=123)
 stopifnot(identical(dimnames(cc@net$prob)[[3]],rownames(LR)),identical(rownames(cc@LR$LRsig),rownames(LR)),all(is.finite(cc@net$prob)),all(cc@net$prob>=0))
 counts<-table(cc@idents);before<-cc@net$prob;checks<-list();for(th in c(10,30)){filtered<-filterCommunication(cc,min.cells=th);expected<-before;small<-which(counts<=th);if(length(small)){expected[small,,]<-0;expected[,small,]<-0};stopifnot(identical(filtered@net$prob,expected),identical(cc@net$prob,before));checks[[as.character(th)]]<-data.table(network_id=z$network_id,min_cells=th,filter_boundary_verified=TRUE,n_removed=sum(counts<=th))};put(rbindlist(checks),paste0('membership/',z$network_id,'__filter_QC.tsv'))
 o<-list(prob=cc@net$prob,pval=cc@net$pval,counts=counts,LR=LR,cache_key=key,input_hash=z$modeled_hash,provenance=provenance,params=params,data_max=dm,network_id=z$network_id,cohort=z$cohort,scope=z$scope,high=z$high,n_cells=nrow(md))
 qsave(o,paste0(fp,'.tmp'),preset='high',nthreads=1);stopifnot(file.rename(paste0(fp,'.tmp'),fp));cat('DONE',format(Sys.time()),'\n');z[,`:=`(status='completed',cache_key=key,data_max=dm,n_LR=nrow(LR))];put(z,paste0('membership/',z$network_id,'__status.tsv'));z
 },error=function(e){cat('ERROR',conditionMessage(e),'\n');z[,`:=`(status='failed',reason=conditionMessage(e))];put(z,paste0('membership/',z$network_id,'__status.tsv'));z})}
r<-parallel::mclapply(seq_len(nrow(p)),worker,mc.cores=cfg$workers,mc.preschedule=FALSE);put(rbindlist(r,fill=TRUE),'network_status.tsv');stopifnot(all(vapply(r,function(z)z$status%in%c('completed','completed_reused','NE'),logical(1))));cat('INFERENCE DONE\n')
