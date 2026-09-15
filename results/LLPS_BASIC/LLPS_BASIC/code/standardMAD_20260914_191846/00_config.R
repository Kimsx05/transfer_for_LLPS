# Phase 4 configuration extension. The immutable Phase 1/2 configuration remains
# in the parent code directory for auditability.
source("/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/code/00_config.R")

config$code_version <- "20260914_152000"
config$qs_compress_level <- 4L
config$normalization_method <- "LogNormalize"
config$normalization_scale_factor <- 10000
config$hvg_method <- "vst"
config$scale_features <- "all_RNA_features"
config$pca_npcs <- 30L
config$neighbors_dims <- 1:30
config$umap_dims <- 1:30
config$harmony_group <- "orig.ident"
config$cluster_algorithm <- 1L
config$marker_test <- "wilcox"

timestamped_path <- function(directory, stem, extension) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  candidate <- file.path(directory, paste0(stem, "_", stamp, ".", extension))
  if (file.exists(candidate) || dir.exists(candidate)) {
    candidate <- file.path(directory, paste0(stem, "_", stamp, "_", Sys.getpid(), ".", extension))
  }
  assert_new_file(candidate)
  candidate
}

read_manifest_object <- function(file_path, object_name) {
  extension <- tolower(tools::file_ext(file_path))
  if (extension == "rdata") {
    holder <- new.env(parent = emptyenv())
    loaded <- load(file_path, envir = holder, verbose = FALSE)
    if (!object_name %in% loaded) {
      stop("Manifest object not found in RData: ", object_name, call. = FALSE)
    }
    return(holder[[object_name]])
  }
  if (extension == "qs") return(qs::qread(file_path, nthreads = config$max_threads))
  if (extension == "qs2") return(qs2::qs_read(file_path, nthreads = config$max_threads))
  stop("Unsupported input type: ", extension, call. = FALSE)
}

latest_checkpoint <- function(pattern) {
  candidates <- list.files(paths$checkpoints, pattern = pattern, full.names = TRUE)
  if (!length(candidates)) stop("No checkpoint matching: ", pattern, call. = FALSE)
  candidates[[which.max(file.info(candidates)$mtime)]]
}

save_qs2_new <- function(object, path) {
  assert_new_file(path)
  qs2::qs_save(
    object = object,
    file = path,
    compress_level = config$qs_compress_level,
    shuffle = TRUE,
    nthreads = config$max_threads
  )
  if (!file.exists(path) || file.info(path)$size <= 0) stop("Checkpoint was not written: ", path, call. = FALSE)
  invisible(path)
}
