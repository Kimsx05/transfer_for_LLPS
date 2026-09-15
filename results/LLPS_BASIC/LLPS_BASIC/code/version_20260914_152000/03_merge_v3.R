source(file.path(dirname(normalizePath(sys.frame(1)$ofile)), "00_config.R"))
start_time <- Sys.time()
log_message("03_merge_v3", "START", "reason=checkpoint_before_summary_and_unnamed_group_vectors")

manifest <- read.delim(file.path(paths$audit, "input_manifest.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
objects <- vector("list", nrow(manifest))
names(objects) <- manifest$dataset_name
for (i in seq_len(nrow(manifest))) {
  log_message("03_merge_v3", "READ_START", manifest$file_path[[i]])
  object <- read_manifest_object(manifest$file_path[[i]], manifest$object_name[[i]])
  if (!inherits(object, "Seurat") || !"RNA" %in% names(object@assays)) {
    stop("Invalid manifest object: ", manifest$dataset_name[[i]], call. = FALSE)
  }
  object$dataset <- manifest$dataset_name[[i]]
  objects[[i]] <- object
  log_message("03_merge_v3", "READ_END", manifest$dataset_name[[i]], "cells=", ncol(object))
}

merged <- merge(
  x = objects[[1L]], y = objects[-1L], add.cell.ids = names(objects),
  merge.data = FALSE, merge.dr = FALSE, project = "LLPS_scRNA_standard_v1"
)
log_message("03_merge_v3", "MERGE_RETURNED")
if (ncol(merged) != sum(manifest$n_cells)) stop("Merged cell count mismatch", call. = FALSE)
if (anyDuplicated(colnames(merged))) stop("Merged cell barcodes are not unique", call. = FALSE)
if (!all(c("dataset", "orig.ident") %in% colnames(merged[[]]))) stop("Required metadata missing after merge", call. = FALSE)

checkpoint <- timestamped_path(paths$checkpoints, "01_merged_full", "qs2")
log_message("03_merge_v3", "CHECKPOINT_WRITE_START", checkpoint)
save_qs2_new(merged, checkpoint)
log_message("03_merge_v3", "CHECKPOINT_WRITE_END", checkpoint)
write_lines_new(checkpoint, file.path(paths$logs, "01_merged_checkpoint_path_v3.txt"))

make_count_rows <- function(values, group_type) {
  counts <- table(unname(as.character(values)), useNA = "ifany")
  data.frame(
    group = unname(names(counts)), n_cells = unname(as.integer(counts)),
    group_type = rep(group_type, length(counts)), stringsAsFactors = FALSE
  )
}
summary_table <- rbind(
  data.frame(group = "ALL", n_cells = ncol(merged), group_type = "total", stringsAsFactors = FALSE),
  make_count_rows(merged$dataset, "dataset"),
  make_count_rows(merged$orig.ident, "orig.ident")
)
summary_table$n_features_union <- rep(nrow(merged[["RNA"]]), nrow(summary_table))
merge_summary_path <- timestamped_path(paths$audit, "merge_summary", "tsv")
write_tsv_new(summary_table, merge_summary_path)
log_message("03_merge_v3", "END", "cells=", ncol(merged), "features=", nrow(merged[["RNA"]]),
            "checkpoint=", checkpoint,
            "elapsed_seconds=", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 3))
