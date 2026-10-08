set.seed(123)
.libPaths(unique(c('/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_atlas_stage1_v1/Rlib', '/home/jinshengxi/Rlibrary/LLPS_pipeline_v1_20260914_150459', '/home/jinshengxi/Rlibrary/packages', .libPaths())))
suppressPackageStartupMessages({library(Seurat); library(qs); library(data.table)})
out <- '/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7'
paths <- c(major='/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/results/objects/major_annotated_with_LLPS77_AUC.qs', fibro='/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order/global_res012_object.qs')
for (nm in names(paths)) {
  cat('\nREADING ', nm, ': ', paths[[nm]], '\n', sep=''); flush.console()
  o <- if(nm=='major') qs2::qs_read(paths[[nm]], nthreads=2) else qs::qread(paths[[nm]], nthreads=2)
  md <- o[[]]
  fwrite(data.table(cell=rownames(md), md), file.path(out,'00_audit',paste0(nm,'_actual_metadata.tsv')), sep='\t', na='NA')
  capture.output({
    cat('Object path:', paths[[nm]], '\n')
    print(dim(o)); print(o@version); print(colnames(md)); print(names(o@assays)); print(names(o@reductions)); print(o@misc)
    for(a in names(o@assays)) {cat('\nASSAY:',a,'\n'); print(class(o[[a]])); print(Layers(o[[a]])); for(l in Layers(o[[a]])) {x<-LayerData(o,assay=a,layer=l); print(list(layer=l,dim=dim(x),class=class(x),range=range(x),n_empty_features=sum(Matrix::rowSums(x != 0)==0))); rm(x)}}
    cat('\nCOMMANDS\n'); print(lapply(o@commands,function(x) list(call=x@call.string,params=x@params)))
    for (f in colnames(md)) if(!is.numeric(md[[f]]) || grepl('cluster',f,ignore.case=TRUE)) {cat('\nFIELD: ',f,'\n',sep=''); print(head(sort(table(md[[f]],useNA='ifany'),decreasing=TRUE),100))}
    cat('\nIDENTITY\n'); print(table(Idents(o)))
  },file=file.path(out,'00_audit',paste0(nm,'_object_inspection.txt')))
  saveRDS(list(path=paths[[nm]],metadata=md,identity=Idents(o),assays=names(o@assays),reductions=names(o@reductions),version=o@version,commands=lapply(o@commands,function(x) list(call=x@call.string,params=x@params))),file.path(out,'00_audit',paste0(nm,'_inspection.rds')))
  cat('DONE ',nm,' cells=',ncol(o),'\n'); flush.console(); rm(o,md); gc()
}
capture.output(sessionInfo(),file=file.path(out,'00_audit','sessionInfo_preflight.txt'))
ip<-installed.packages(); fwrite(as.data.table(ip,keep.rownames='installed_name'),file.path(out,'00_audit','installed_R_packages.tsv'),sep='\t')
