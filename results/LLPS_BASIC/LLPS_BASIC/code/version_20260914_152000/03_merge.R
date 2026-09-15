source(file.path(dirname(normalizePath(sys.frame(1)$ofile)), "00_config.R"))
start_time <- Sys.time()
log_message("03_merge", "START")

manifest <- read.delim(file.path(paths$audit, "input_manifest.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
objects <- vector("list", nrow(manifest))
names(objects) <- manifest$dataset_name

for (i in seq_len(nrow(manifest))) {
  log_message("03_merge", "READ_START", manifest$file_path[[i]])
  object <- read_manifest_object(manifest$file_path[[i]], manifest$object_name[[i]])
  if (!inherits(object, "Seurat")) stop("Manifest object is no longer Seurat: ", manifest$dataset_name[[i]], call. = FALSE)
  if (!"RNA" %in% names(object@assays)) stop("RNA assay missing: ", manifest$dataset_name[[i]], call. = FALSE)
  object$dataset <- manifest$dataset_name[[i]]
  objects[[i]] <- object
  log_message("03_merge", "READ_END", manifest$dataset_name[[i]], "cells=", ncol(object))
}

merged <- Seurat::merge(
  x = objects[[1L]],
  y = objects[-1L],
  add.cell.ids = names(objects),
  merge.data = FALSE,
  merge.dr = FALSE,
  project = "LLPS_scRNA_standard_v1"
)

if (ncol(merged) != sum(manifest$n_cells)) stop("Merged cell count mismatch", call. = FALSE)
if (anyDuplicated(colnames(merged))) stop("Merged cell barcodes are not unique", call. = FALSE)
if (!all(c("dataset", "orig.ident") %in% colnames(merged[[]]))) stop("Required metadata missing after merge", call. = FALSE)

dataset_counts <- as.data.frame(table(merged$dataset), stringsAsFactors = FALSE)
colnames(dataset_counts) <- c("group", "n_cells")
dataset_counts$group_type <- "dataset"
orig_counts <- as.data.frame(table(merged$orig.ident), stringsAsFactors = FALSE)
colnames(orig_counts) <- c("group", "n_cells")
orig_counts$group_type <- "orig.ident"
summary_table <- rbind(
  data.frame(group = "ALL", n_cells = ncol(merged), group_type = "total"),
  dataset_counts[, c("group", "n_cells", "group_type")],
  orig_counts[, c("group", "n_cells", "group_type")]
)
summary_table$n_features_union <- nrow(merged[["RNA"]])
write_tsv_new(summary_table, file.path(paths$audit, "merge_summary.tsv"))

checkpoint <- timestamped_path(paths$checkpoints, "01_merged_full", "qs2")
save_qs2_new(merged, checkpoint)
write_lines_new(checkpoint, file.path(paths$logs, "01_merged_checkpoint_path.txt"))
log_message("03_merge", "END", "cells=", ncol(merged), "features=", nrow(merged[["RNA"]]),
            "checkpoint=", checkpoint,
            "elapsed_seconds=", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 3))
