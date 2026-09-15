source(file.path(dirname(normalizePath(sys.frame(1)$ofile)), "00_config.R"))
start_time <- Sys.time()
log_message("03_merge_v5", "START", "strategy=repair_malformed_feature_metadata_then_multi_merge")

manifest <- read.delim(file.path(paths$audit, "input_manifest.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
objects <- vector("list", nrow(manifest))
repairs <- vector("list", nrow(manifest))
names(objects) <- manifest$dataset_name

for (i in seq_len(nrow(manifest))) {
  object <- read_manifest_object(manifest$file_path[[i]], manifest$object_name[[i]])
  if (!inherits(object, "Seurat") || !"RNA" %in% names(object@assays)) stop("Invalid input object", call. = FALSE)
  assay <- object[["RNA"]]
  feature_meta <- slot(assay, "meta.data")
  needs_repair <- nrow(feature_meta) != nrow(assay) || length(rownames(feature_meta)) != nrow(assay)
  repairs[[i]] <- data.frame(
    dataset = manifest$dataset_name[[i]], n_features = nrow(assay),
    old_feature_metadata_rows = nrow(feature_meta), old_feature_metadata_columns = ncol(feature_meta),
    repair_applied = needs_repair,
    repair_action = if (needs_repair) "Reset malformed, recomputable RNA feature metadata in memory; counts/features unchanged" else "None",
    stringsAsFactors = FALSE
  )
  if (needs_repair) {
    slot(assay, "meta.data") <- data.frame(row.names = rownames(assay))
    object[["RNA"]] <- assay
    validity <- validObject(object, test = TRUE)
    if (!isTRUE(validity)) stop("In-memory repair failed validity for ", manifest$dataset_name[[i]], ": ", validity, call. = FALSE)
    log_message("03_merge_v5", "FEATURE_METADATA_REPAIRED", manifest$dataset_name[[i]],
                "features=", nrow(assay), "old_meta_rows=", nrow(feature_meta))
  }
  object$dataset <- manifest$dataset_name[[i]]
  objects[[i]] <- object
  log_message("03_merge_v5", "READ_READY", manifest$dataset_name[[i]], "cells=", ncol(object))
}

repair_path <- timestamped_path(paths$audit, "feature_metadata_repairs", "tsv")
write_tsv_new(do.call(rbind, repairs), repair_path)

merged <- merge(
  x = objects[[1L]], y = objects[-1L], add.cell.ids = names(objects),
  merge.data = FALSE, merge.dr = FALSE, project = "LLPS_scRNA_standard_v1"
)
log_message("03_merge_v5", "MERGE_RETURNED", "cells=", ncol(merged), "features=", nrow(merged[["RNA"]]))
if (ncol(merged) != sum(manifest$n_cells)) stop("Merged cell count mismatch", call. = FALSE)
if (anyDuplicated(colnames(merged))) stop("Merged cell barcodes are not unique", call. = FALSE)
if (!all(c("dataset", "orig.ident") %in% colnames(merged[[]]))) stop("Required metadata missing after merge", call. = FALSE)

checkpoint <- timestamped_path(paths$checkpoints, "01_merged_full", "qs2")
save_qs2_new(merged, checkpoint)
write_lines_new(checkpoint, file.path(paths$logs, "01_merged_checkpoint_path_v5.txt"))

make_count_rows <- function(values, group_type) {
  counts <- table(unname(as.character(values)), useNA = "ifany")
  data.frame(group = unname(names(counts)), n_cells = unname(as.integer(counts)),
             group_type = rep(group_type, length(counts)), stringsAsFactors = FALSE)
}
summary_table <- rbind(
  data.frame(group = "ALL", n_cells = ncol(merged), group_type = "total", stringsAsFactors = FALSE),
  make_count_rows(merged$dataset, "dataset"), make_count_rows(merged$orig.ident, "orig.ident")
)
summary_table$n_features_union <- rep(nrow(merged[["RNA"]]), nrow(summary_table))
merge_summary_path <- timestamped_path(paths$audit, "merge_summary", "tsv")
write_tsv_new(summary_table, merge_summary_path)
log_message("03_merge_v5", "END", "checkpoint=", checkpoint,
            "elapsed_seconds=", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 3))
