source("/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/code/standardMAD_20260914_191846/00_config.R")
start_time <- Sys.time(); stamp <- format(start_time,"%Y%m%d_%H%M%S")
input_path <- "/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/checkpoints/02_harmony_clustered_res0.8_standardMAD_20260915_092711.qs2"
out_dir <- file.path(paths$harmony,paste0("SCP_UMAP_standardMAD_",stamp))
if (file.exists(out_dir)||dir.exists(out_dir)) stop("Refusing existing output directory")
dir.create(out_dir,recursive=FALSE,showWarnings=FALSE)
obj <- qs2::qs_read(input_path,nthreads=config$max_threads)
if (!"umap" %in% names(obj@reductions)) stop("UMAP missing")
cluster_col <- grep("RNA_snn_res\\.0\\.8$",colnames(obj[[]]),value=TRUE)
if(length(cluster_col)!=1L) stop("Resolution 0.8 cluster column not unique")
obj$cluster_res0.8 <- obj[[cluster_col,drop=TRUE]]
plots <- list(
  UMAP_cluster=SCP::CellDimPlot(obj,group.by="cluster_res0.8",reduction="umap",label=TRUE,raster=TRUE,theme_use="theme_blank",title="Clusters",seed=config$seed),
  UMAP_dataset=SCP::CellDimPlot(obj,group.by="dataset",reduction="umap",raster=TRUE,theme_use="theme_blank",title="Dataset",seed=config$seed),
  UMAP_orig_ident=SCP::CellDimPlot(obj,group.by="orig.ident",reduction="umap",raster=TRUE,theme_use="theme_blank",title="orig.ident",seed=config$seed),
  UMAP_split_by_dataset=SCP::CellDimPlot(obj,group.by="cluster_res0.8",reduction="umap",split.by="dataset",raster=TRUE,theme_use="theme_blank",title="Clusters by dataset",seed=config$seed)
)
for(nm in names(plots)){
  w <- if(nm=="UMAP_split_by_dataset") 18 else 10; h <- if(nm=="UMAP_split_by_dataset") 10 else 8
  ggplot2::ggsave(file.path(out_dir,paste0(nm,".pdf")),plots[[nm]],width=w,height=h,bg="white",limitsize=FALSE)
  ggplot2::ggsave(file.path(out_dir,paste0(nm,".png")),plots[[nm]],width=w,height=h,dpi=200,bg="white",limitsize=FALSE)
}
write_lines_new(c(paste("input",input_path,sep="\t"),paste("SCP_version",as.character(utils::packageVersion("SCP")),sep="\t"),"function\tSCP::CellDimPlot","theme_use\ttheme_blank",paste("elapsed_seconds",round(as.numeric(difftime(Sys.time(),start_time,units="secs")),3),sep="\t")),file.path(out_dir,"SCP_UMAP_summary.txt"))
message(format(Sys.time())," END output=",out_dir)
