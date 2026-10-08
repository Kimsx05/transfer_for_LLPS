source('/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/scripts/common.R')
cfg<-readRDS(file.path(ROOT,'01_metadata/frozen_analysis_inputs.rds'))
new_object<-file.path(ROOT,'01_metadata/major_frozen_with_cellchat_group.qs2')
if(file.exists(new_object))stop('Analysis-copy object already exists; refusing overwrite')
hash_before<-c(major=sha256(MAJOR_OBJECT),fibro=sha256(FIBRO_OBJECT))
o<-qs2::qs_read(MAJOR_OBJECT,nthreads=2)
md<-o[[]];ids<-Idents(o);reductions<-o@reductions;graphs<-o@graphs;commands<-o@commands
stopifnot(identical(rownames(md),cfg$metadata$cell),!('cellchat_group'%in%names(md)))
o$cellchat_group<-factor(cfg$metadata$cellchat_group,levels=GROUPS)
stopifnot(identical(o[[]][,names(md),drop=FALSE],md),identical(Idents(o),ids),identical(o@reductions,reductions),identical(o@graphs,graphs),identical(o@commands,commands))
x<-LayerData(o,assay='RNA',layer='data')
stopifnot(inherits(x,'dgCMatrix'),identical(colnames(x),colnames(o)),!anyDuplicated(rownames(x)),all(is.finite(x@x)),min(x@x)>=0)
normalization<-o@commands[['NormalizeData.RNA']]@params
stopifnot(normalization$normalization.method=='LogNormalize',normalization$scale.factor==10000)
# Verify existing normalization on the first cell of every biological sample.
# This calculates a check only; it never writes expression back to the object.
check_cells<-vapply(sort(unique(md$orig.ident)),function(s)which(md$orig.ident==s)[1],integer(1))
c<-LayerData(o,assay='RNA',layer='counts')[,check_cells,drop=FALSE]
expected<-log1p(t(t(c)/Matrix::colSums(c))*10000)
normalization_error<-max(abs(expected-x[,check_cells,drop=FALSE]))
stopifnot(normalization_error<1e-10)
write_tsv(data.table(metric=c('check_rule','checked_cells','max_absolute_error','normalized_expression_overwritten'),value=c('first cell per biological sample; full gene universe',length(check_cells),format(normalization_error,digits=17),'FALSE')),'00_audit/normalized_layer_verification.tsv')
qs2::qs_save(o,new_object,nthreads=2)
manifest<-copy(cfg$samples);manifest[,sample_id:=sprintf('S%03d',seq_len(.N))]
manifest[,input_path:=file.path(ROOT,'02_cellchat_per_sample',sample_id,'RNA_data_input.qs2')]
for(i in seq_len(nrow(manifest))) {
  dir.create(dirname(manifest$input_path[i]),recursive=TRUE,showWarnings=FALSE)
  if(file.exists(manifest$input_path[i]))stop('Sample input already exists: ',manifest$sample[i])
  cells<-cfg$metadata[sample==manifest$sample[i] & !is.na(cellchat_group),cell]
  sm<-cfg$metadata[match(cells,cell),.(cell,sample,patient_uid,dataset,cellchat_group)]
  stopifnot(identical(sm$cell,cells),all(sm$sample==manifest$sample[i]))
  qs2::qs_save(list(expression=x[,cells,drop=FALSE],metadata=sm),manifest$input_path[i],nthreads=2)
  cat('PREPARED ',manifest$sample_id[i],' ',manifest$sample[i],' n=',length(cells),'\n',sep='');flush.console()
}
write_tsv(manifest,'01_metadata/per_sample_input_manifest.tsv')
hash_after<-c(major=sha256(MAJOR_OBJECT),fibro=sha256(FIBRO_OBJECT))
stopifnot(identical(hash_before,hash_after))
write_tsv(data.table(input=c(MAJOR_OBJECT,FIBRO_OBJECT),sha256_before=hash_before,sha256_after=hash_after,unchanged=hash_before==hash_after),'00_audit/frozen_object_preservation.tsv')
writeLines(c('PASS: Added only cellchat_group to a new full-atlas copy.','All original metadata, Idents, reductions, graphs and commands are identical.','RNA assay content was read without mutation. Input source SHA256 values remain identical.',paste0('New object: ',new_object)),file.path(ROOT,'00_audit/input_preservation_audit.txt'))
cat('INPUT PREPARATION PASS\n')
