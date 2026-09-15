source(file.path(dirname(normalizePath(sys.frame(1)$ofile)), "00_config.R"))
start_time <- Sys.time()
log_message("04_qc_v2", "START", "reason=Seurat_5.5.1_S3_subset")

merged_path <- latest_checkpoint("^01_merged_full_.*\\.qs2$")
merged <- qs2::qs_read(merged_path, nthreads = config$max_threads)
Seurat::DefaultAssay(merged) <- "RNA"

count_layers <- grep("^counts", SeuratObject::Layers(merged[["RNA"]]), value = TRUE)
if (!length(count_layers)) stop("No RNA counts layers found", call. = FALSE)

compute_feature_percent <- function(object, feature_names) {
  totals <- setNames(numeric(ncol(object)), colnames(object))
  selected <- character()
  for (layer in count_layers) {
    matrix <- SeuratObject::LayerData(object[["RNA"]], layer = layer)
    present <- intersect(feature_names, rownames(matrix))
    if (length(present)) {
      totals[colnames(matrix)] <- totals[colnames(matrix)] + Matrix::colSums(matrix[present, , drop = FALSE])
      selected <- union(selected, present)
    }
  }
  denominator <- as.numeric(object$nCount_RNA)
  percentage <- ifelse(denominator > 0, 100 * totals[colnames(object)] / denominator, NA_real_)
  list(percent = percentage, genes_present = selected)
}

all_features <- rownames(merged[["RNA"]])
mt_genes <- grep("^MT-", all_features, value = TRUE)
hb_requested <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
hb_genes <- intersect(hb_requested, all_features)
if (!length(mt_genes)) stop("No human mitochondrial genes matching ^MT- were found", call. = FALSE)
if (!length(hb_genes)) stop("None of the configured hemoglobin genes were found", call. = FALSE)

mt_result <- compute_feature_percent(merged, mt_genes)
hb_result <- compute_feature_percent(merged, hb_genes)
merged$percent.mt <- mt_result$percent
merged$percent.hb <- hb_result$percent

meta <- merged[[]]
meta$cell_id <- rownames(meta)
orig_levels <- sort(unique(as.character(meta$orig.ident)))

threshold_rows <- lapply(orig_levels, function(group) {
  idx <- which(as.character(meta$orig.ident) == group)
  log_count <- log1p(meta$nCount_RNA[idx])
  log_feature <- log1p(meta$nFeature_RNA[idx])
  count_median <- median(log_count, na.rm = TRUE)
  count_mad <- stats::mad(log_count, center = count_median, constant = 1, na.rm = TRUE)
  feature_median <- median(log_feature, na.rm = TRUE)
  feature_mad <- stats::mad(log_feature, center = feature_median, constant = 1, na.rm = TRUE)
  mt_median <- median(meta$percent.mt[idx], na.rm = TRUE)
  mt_mad <- stats::mad(meta$percent.mt[idx], center = mt_median, constant = 1, na.rm = TRUE)
  data.frame(
    orig.ident = group, initial_cells = length(idx),
    log1p_nCount_median = count_median, log1p_nCount_MAD = count_mad,
    log1p_nCount_lower = count_median - config$qc_mad_multiplier * count_mad,
    log1p_nCount_upper = count_median + config$qc_mad_multiplier * count_mad,
    nCount_lower = pmax(0, expm1(count_median - config$qc_mad_multiplier * count_mad)),
    nCount_upper = expm1(count_median + config$qc_mad_multiplier * count_mad),
    log1p_nFeature_median = feature_median, log1p_nFeature_MAD = feature_mad,
    log1p_nFeature_lower = feature_median - config$qc_mad_multiplier * feature_mad,
    log1p_nFeature_upper = feature_median + config$qc_mad_multiplier * feature_mad,
    nFeature_lower = pmax(0, expm1(feature_median - config$qc_mad_multiplier * feature_mad)),
    nFeature_upper = expm1(feature_median + config$qc_mad_multiplier * feature_mad),
    percent_mt_median = mt_median, percent_mt_MAD = mt_mad,
    percent_mt_MAD_upper = mt_median + config$qc_mad_multiplier * mt_mad,
    mt_cutoff = pmin(mt_median + config$qc_mad_multiplier * mt_mad, config$mt_hard_ceiling),
    hb_cutoff = config$hb_cutoff,
    stringsAsFactors = FALSE
  )
})
thresholds <- do.call(rbind, threshold_rows)
write_tsv_new(thresholds, timestamped_path(paths$qc, "qc_thresholds_by_orig_ident", "tsv"))

lookup <- match(as.character(meta$orig.ident), thresholds$orig.ident)
log_count_all <- log1p(meta$nCount_RNA)
log_feature_all <- log1p(meta$nFeature_RNA)
merged$fail_nCount <- log_count_all < thresholds$log1p_nCount_lower[lookup] |
  log_count_all > thresholds$log1p_nCount_upper[lookup] | is.na(log_count_all)
merged$fail_nFeature <- log_feature_all < thresholds$log1p_nFeature_lower[lookup] |
  log_feature_all > thresholds$log1p_nFeature_upper[lookup] | is.na(log_feature_all)
merged$fail_mt <- merged$percent.mt > thresholds$mt_cutoff[lookup] | is.na(merged$percent.mt)
merged$fail_hb <- merged$percent.hb > config$hb_cutoff | is.na(merged$percent.hb)
merged$QC_pass <- !merged$fail_nCount & !merged$fail_nFeature & !merged$fail_mt & !merged$fail_hb

qc_meta <- merged[[]]
counts_rows <- lapply(orig_levels, function(group) {
  x <- qc_meta[as.character(qc_meta$orig.ident) == group, , drop = FALSE]
  n_fail_total <- rowSums(x[, c("fail_nCount", "fail_nFeature", "fail_mt", "fail_hb")])
  data.frame(
    orig.ident = group, initial_cells = nrow(x),
    fail_nCount = sum(x$fail_nCount), fail_nFeature = sum(x$fail_nFeature),
    fail_mt = sum(x$fail_mt), fail_hb = sum(x$fail_hb),
    overlap = sum(n_fail_total > 1L), final_QC_pass = sum(x$QC_pass),
    retention_rate = sum(x$QC_pass) / nrow(x), stringsAsFactors = FALSE
  )
})
qc_counts <- do.call(rbind, counts_rows)
write_tsv_new(qc_counts, timestamped_path(paths$qc, "qc_cell_counts", "tsv"))
write_tsv_new(data.frame(
  marker_group = c(rep("mitochondrial", length(mt_genes)), rep("hemoglobin", length(hb_requested))),
  gene = c(mt_genes, hb_requested),
  present = c(rep(TRUE, length(mt_genes)), hb_requested %in% hb_genes)
), timestamped_path(paths$qc, "qc_gene_presence", "tsv"))

plot_data <- qc_meta[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.hb")]
metric_names <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.hb")
plot_long <- data.frame(
  orig.ident = rep(plot_data$orig.ident, times = length(metric_names)),
  metric = rep(metric_names, each = nrow(plot_data)),
  value = unlist(plot_data[metric_names], use.names = FALSE),
  stringsAsFactors = FALSE
)
plot_long$metric <- factor(plot_long$metric, levels = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.hb"))
qc_plot <- ggplot2::ggplot(plot_long, ggplot2::aes(x = orig.ident, y = value)) +
  ggplot2::geom_violin(scale = "width", fill = "#4C78A8", color = NA, trim = TRUE) +
  ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 1) +
  ggplot2::labs(x = "orig.ident", y = NULL, title = "QC distributions before filtering") +
  ggplot2::theme_bw(base_size = 8) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1), panel.grid = ggplot2::element_blank())
pdf_path <- timestamped_path(paths$qc, "qc_distributions_by_orig_ident", "pdf")
png_path <- timestamped_path(paths$qc, "qc_distributions_by_orig_ident", "png")
assert_new_file(pdf_path); assert_new_file(png_path)
ggplot2::ggsave(pdf_path, qc_plot, width = 22, height = 18, limitsize = FALSE)
ggplot2::ggsave(png_path, qc_plot, width = 22, height = 18, dpi = 180, limitsize = FALSE)

qc_filtered <- subset(merged, cells = colnames(merged)[merged$QC_pass])
if (ncol(qc_filtered) != sum(merged$QC_pass)) stop("QC subset cell count mismatch", call. = FALSE)
checkpoint <- timestamped_path(paths$checkpoints, "01b_qc_filtered", "qs2")
save_qs2_new(qc_filtered, checkpoint)
write_lines_new(checkpoint, file.path(paths$logs, "01b_qc_checkpoint_path_v2.txt"))
log_message("04_qc_v2", "END", "initial=", ncol(merged), "QC_pass=", ncol(qc_filtered),
            "checkpoint=", checkpoint,
            "elapsed_seconds=", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 3))
