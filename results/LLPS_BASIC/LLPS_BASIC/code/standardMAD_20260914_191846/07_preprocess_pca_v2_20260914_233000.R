source("/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/code/standardMAD_20260914_191846/00_config.R")
start_time <- Sys.time(); stamp <- format(start_time, "%Y%m%d_%H%M%S")
input_path <- "/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/checkpoints/01d_final_QC_DF_singlets_20260914_232441.qs2"
out_dir <- file.path(paths$pca, paste0("standardMAD_", stamp))
if (file.exists(out_dir) || dir.exists(out_dir)) stop("Refusing existing output directory: ", out_dir)
dir.create(out_dir, recursive=FALSE, showWarnings=FALSE)
message(format(Sys.time()), " START cells checkpoint=", input_path)
obj <- qs2::qs_read(input_path, nthreads=config$max_threads)
if (ncol(obj) != 372323L) stop("Unexpected final cell count: ", ncol(obj))
Seurat::DefaultAssay(obj) <- "RNA"
count_layers <- grep("^counts", SeuratObject::Layers(obj[["RNA"]]), value=TRUE)
if (!length(count_layers)) stop("No RNA counts layer")
join_applied <- length(count_layers) > 1L || !identical(count_layers, "counts")
if (join_applied) obj <- SeuratObject::JoinLayers(obj, assay="RNA")
set.seed(config$seed)
obj <- Seurat::NormalizeData(obj, assay="RNA", normalization.method="LogNormalize", scale.factor=10000, margin=1, verbose=TRUE)
obj <- Seurat::FindVariableFeatures(obj, assay="RNA", selection.method="vst", nfeatures=3000, verbose=TRUE)
if (length(SeuratObject::VariableFeatures(obj)) != 3000L) stop("Did not obtain exactly 3000 HVGs")
obj <- Seurat::ScaleData(obj, assay="RNA", features=SeuratObject::VariableFeatures(obj), do.scale=TRUE, do.center=TRUE, verbose=TRUE)
obj <- Seurat::RunPCA(obj, assay="RNA", features=SeuratObject::VariableFeatures(obj), npcs=30, seed.use=config$seed, approx=TRUE, verbose=TRUE)
if (ncol(Seurat::Embeddings(obj, "pca")) < 30L) stop("Fewer than 30 PCs produced")
blank <- ggplot2::theme_void() + ggplot2::theme(plot.background=ggplot2::element_rect(fill="white", colour=NA))
p1 <- Seurat::ElbowPlot(obj, ndims=30) + ggplot2::theme_bw()
p2 <- Seurat::DimPlot(obj, reduction="pca", group.by="dataset", raster=TRUE) + blank
p3 <- Seurat::DimPlot(obj, reduction="pca", group.by="orig.ident", raster=TRUE) + blank
for (nm in c("PCA_elbow","PCA_by_dataset","PCA_by_orig_ident")) {
  p <- get(c(PCA_elbow="p1",PCA_by_dataset="p2",PCA_by_orig_ident="p3")[[nm]])
  ggplot2::ggsave(file.path(out_dir,paste0(nm,".pdf")),p,width=10,height=7,bg="white")
  ggplot2::ggsave(file.path(out_dir,paste0(nm,".png")),p,width=10,height=7,dpi=200,bg="white")
}
checkpoint <- timestamped_path(paths$checkpoints, "01e_preprocessed_PCA30_standardMAD", "qs2")
save_qs2_new(obj, checkpoint)
write_lines_new(c(paste("input",input_path,sep="\t"),paste("cells",ncol(obj),sep="\t"),paste("HVG",length(SeuratObject::VariableFeatures(obj)),sep="\t"),paste("PCs",30,sep="\t"),paste("JoinLayers_applied",join_applied,sep="\t"),paste("checkpoint",checkpoint,sep="\t")), file.path(out_dir,"preprocess_pca_summary.txt"))
message(format(Sys.time()), " END checkpoint=", checkpoint, " elapsed_seconds=", round(as.numeric(difftime(Sys.time(),start_time,units="secs")),3))
