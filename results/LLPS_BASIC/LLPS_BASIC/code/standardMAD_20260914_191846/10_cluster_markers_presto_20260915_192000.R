source("/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/code/standardMAD_20260914_191846/00_config.R")
start_time <- Sys.time(); stamp <- format(start_time, "%Y%m%d_%H%M%S")
if (!requireNamespace("presto", quietly=TRUE)) stop("presto is required for this accelerated marker run")
input_path <- "/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/checkpoints/02_harmony_clustered_res0.8_standardMAD_20260915_092711.qs2"
out_dir <- file.path(paths$markers, paste0("standardMAD_res0.8_", stamp))
if (file.exists(out_dir) || dir.exists(out_dir)) stop("Refusing existing marker output directory: ", out_dir)
dir.create(out_dir, recursive=FALSE, showWarnings=FALSE)
if (!dir.exists(out_dir)) stop("Could not create marker output directory")
message(format(Sys.time()), " START marker checkpoint=", input_path)
obj <- qs2::qs_read(input_path, nthreads=config$max_threads)
if (ncol(obj) != 372323L) stop("Unexpected cell count: ", ncol(obj))
Seurat::DefaultAssay(obj) <- "RNA"
cluster_col <- grep("RNA_snn_res\\.0\\.8$", colnames(obj[[]]), value=TRUE)
if (length(cluster_col) != 1L) stop("Resolution 0.8 cluster column is not unique")
Seurat::Idents(obj) <- factor(obj[[cluster_col, drop=TRUE]])
if (length(levels(Seurat::Idents(obj))) != 48L) stop("Expected 48 clusters, found ", length(levels(Seurat::Idents(obj))))
data_layers <- grep("^data", SeuratObject::Layers(obj[["RNA"]]), value=TRUE)
join_applied <- length(data_layers) > 1L || (length(data_layers)==1L && !identical(data_layers,"data"))
if (join_applied) obj <- SeuratObject::JoinLayers(obj, assay="RNA")
if (!"data" %in% SeuratObject::Layers(obj[["RNA"]])) stop("Joined RNA data layer unavailable")
future::plan(future::sequential)
set.seed(config$seed)
warnings_seen <- character()
markers <- withCallingHandlers(
  Seurat::FindAllMarkers(obj, assay="RNA", slot="data", only.pos=TRUE,
                         min.pct=0.1, logfc.threshold=0.25,
                         test.use="wilcox", random.seed=config$seed,
                         verbose=TRUE, return.thresh=0.01),
  warning=function(w) { warnings_seen <<- c(warnings_seen, conditionMessage(w)); invokeRestart("muffleWarning") }
)
if (!is.data.frame(markers) || !nrow(markers)) stop("FindAllMarkers returned no markers")
fc_candidates <- c("avg_log2FC", "avg_logFC", "avg_diff")
fc_col <- fc_candidates[fc_candidates %in% colnames(markers)][1]
if (is.na(fc_col)) stop("No recognized logFC column: ", paste(colnames(markers),collapse=", "))
if (!"cluster" %in% colnames(markers)) stop("Marker result lacks cluster column")
all_path <- file.path(out_dir, "marker_all.csv")
top_path <- file.path(out_dir, "marker_top50.csv")
assert_new_file(all_path); utils::write.csv(markers, all_path, row.names=FALSE, quote=TRUE)
ord <- order(as.character(markers$cluster), -markers[[fc_col]], na.last=TRUE)
ordered <- markers[ord,,drop=FALSE]
top50 <- do.call(rbind, lapply(split(ordered, as.character(ordered$cluster)), function(z) utils::head(z,50L)))
rownames(top50) <- NULL
assert_new_file(top_path); utils::write.csv(top50, top_path, row.names=FALSE, quote=TRUE)
write_lines_new(c("assay\tRNA","slot\tdata","identity\tRNA_snn_res.0.8","only.pos\tTRUE","min.pct\t0.1","logfc.threshold\t0.25","test.use\twilcox","return.thresh\t0.01",paste("presto_version",as.character(utils::packageVersion("presto")),sep="\t"),paste("logFC_column",fc_col,sep="\t"),paste("clusters",length(unique(markers$cluster)),sep="\t"),paste("all_marker_rows",nrow(markers),sep="\t"),paste("top50_rows",nrow(top50),sep="\t"),paste("JoinLayers_applied",join_applied,sep="\t"),paste("input",input_path,sep="\t")), file.path(out_dir,"marker_parameters.txt"))
write_lines_new(if(length(warnings_seen)) unique(warnings_seen) else "NONE", file.path(out_dir,"marker_warnings.txt"))
message(format(Sys.time()), " END all_rows=",nrow(markers)," top50_rows=",nrow(top50)," output=",out_dir," elapsed_seconds=",round(as.numeric(difftime(Sys.time(),start_time,units="secs")),3))
