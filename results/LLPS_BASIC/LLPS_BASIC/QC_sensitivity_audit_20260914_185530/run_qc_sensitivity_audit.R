# Read-only QC sensitivity audit for LLPS scRNA-seq pipeline v1.
# This script never subsets or saves a Seurat object.

audit_start <- Sys.time()
audit_id <- "20260914_185530"
checkpoint <- "/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/checkpoints/01_merged_full_20260914_153425.qs2"
audit_dir <- "/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/QC_sensitivity_audit_20260914_185530"
plot_dir <- file.path(audit_dir, "plots")
isolated_library <- "/home/jinshengxi/Rlibrary/LLPS_pipeline_v1_20260914_150459"
shared_library <- "/home/jinshengxi/Rlibrary/packages"
.libPaths(unique(c(isolated_library, shared_library, .libPaths())))

max_threads <- 8L
Sys.setenv(OMP_NUM_THREADS = max_threads, OPENBLAS_NUM_THREADS = max_threads,
           MKL_NUM_THREADS = max_threads, RCPP_PARALLEL_NUM_THREADS = max_threads)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(max_threads)
  RhpcBLASctl::omp_set_num_threads(max_threads)
}
if (requireNamespace("RcppParallel", quietly = TRUE)) RcppParallel::setThreadOptions(numThreads = max_threads)
set.seed(123L)

required <- c("qs2", "SeuratObject", "Matrix", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
if (!file.exists(checkpoint)) stop("Checkpoint not found: ", checkpoint, call. = FALSE)
if (file.exists(audit_dir) || dir.exists(audit_dir)) stop("Refusing existing audit directory: ", audit_dir, call. = FALSE)
dir.create(audit_dir, recursive = FALSE, showWarnings = FALSE)
dir.create(plot_dir, recursive = FALSE, showWarnings = FALSE)

new_path <- function(name) {
  path <- file.path(audit_dir, name)
  if (file.exists(path) || dir.exists(path)) stop("Refusing overwrite: ", path, call. = FALSE)
  path
}
write_tsv <- function(x, name) write.table(x, new_path(name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
quant <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 7))
safe_rate <- function(n, d) ifelse(d > 0, 100 * n / d, NA_real_)

message("Reading immutable merge checkpoint: ", checkpoint)
object <- qs2::qs_read(checkpoint, nthreads = max_threads)
if (!inherits(object, "Seurat")) stop("Checkpoint is not a Seurat object", call. = FALSE)
meta <- object[[]]
if (!all(c("orig.ident", "dataset", "nCount_RNA", "nFeature_RNA") %in% colnames(meta))) {
  stop("Required metadata missing from merge checkpoint", call. = FALSE)
}
meta$cell_id <- rownames(meta)
meta$orig.ident <- as.character(meta$orig.ident)
meta$dataset <- as.character(meta$dataset)

# Compute MT/HB percentages directly from immutable counts layers. No object metadata is changed.
hb_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
count_layers <- grep("^counts", SeuratObject::Layers(object[["RNA"]]), value = TRUE)
if (!length(count_layers)) stop("No counts layers found", call. = FALSE)
mt_counts <- setNames(numeric(nrow(meta)), meta$cell_id)
hb_counts <- setNames(numeric(nrow(meta)), meta$cell_id)
orig_mt_features <- setNames(vector("list", length(unique(meta$orig.ident))), sort(unique(meta$orig.ident)))
orig_hb_features <- setNames(vector("list", length(unique(meta$orig.ident))), sort(unique(meta$orig.ident)))
dataset_mt_features <- setNames(vector("list", length(unique(meta$dataset))), sort(unique(meta$dataset)))
dataset_hb_features <- setNames(vector("list", length(unique(meta$dataset))), sort(unique(meta$dataset)))
inhouse_features <- character()

for (layer in count_layers) {
  matrix <- SeuratObject::LayerData(object[["RNA"]], layer = layer, fast = FALSE)
  cells <- colnames(matrix)
  cell_index <- match(cells, meta$cell_id)
  if (anyNA(cell_index)) stop("Layer contains cells absent from metadata: ", layer, call. = FALSE)
  features <- rownames(matrix)
  mt <- grep("^MT-", features, value = TRUE)
  hb <- intersect(hb_genes, features)
  if (length(mt)) mt_counts[cells] <- mt_counts[cells] + Matrix::colSums(matrix[mt, , drop = FALSE])
  if (length(hb)) hb_counts[cells] <- hb_counts[cells] + Matrix::colSums(matrix[hb, , drop = FALSE])
  orig_here <- unique(meta$orig.ident[cell_index])
  dataset_here <- unique(meta$dataset[cell_index])
  for (group in orig_here) {
    orig_mt_features[[group]] <- union(orig_mt_features[[group]], mt)
    orig_hb_features[[group]] <- union(orig_hb_features[[group]], hb)
  }
  for (group in dataset_here) {
    dataset_mt_features[[group]] <- union(dataset_mt_features[[group]], mt)
    dataset_hb_features[[group]] <- union(dataset_hb_features[[group]], hb)
  }
  if ("InhouseData" %in% dataset_here) inhouse_features <- union(inhouse_features, features)
}

denominator <- as.numeric(meta$nCount_RNA)
meta$percent.mt.audit <- ifelse(denominator > 0, 100 * mt_counts[meta$cell_id] / denominator, NA_real_)
meta$percent.hb.audit <- ifelse(denominator > 0, 100 * hb_counts[meta$cell_id] / denominator, NA_real_)
meta$log1p_nCount <- log1p(meta$nCount_RNA)
meta$log1p_nFeature <- log1p(meta$nFeature_RNA)

# Reproduce formal QC exactly: R mad with constant=1 was explicitly used in the formal script.
orig_levels <- sort(unique(meta$orig.ident))
thresholds <- do.call(rbind, lapply(orig_levels, function(group) {
  x <- meta[meta$orig.ident == group, , drop = FALSE]
  count_med <- median(x$log1p_nCount, na.rm = TRUE)
  count_mad <- stats::mad(x$log1p_nCount, center = count_med, constant = 1, na.rm = TRUE)
  feature_med <- median(x$log1p_nFeature, na.rm = TRUE)
  feature_mad <- stats::mad(x$log1p_nFeature, center = feature_med, constant = 1, na.rm = TRUE)
  mt_med <- median(x$percent.mt.audit, na.rm = TRUE)
  mt_mad <- stats::mad(x$percent.mt.audit, center = mt_med, constant = 1, na.rm = TRUE)
  data.frame(
    orig.ident = group, dataset = unique(x$dataset), total_cells = nrow(x),
    count_median_log = count_med, count_MAD_log = count_mad,
    count_lower_log = count_med - 3 * count_mad, count_upper3_log = count_med + 3 * count_mad,
    count_upper5_log = count_med + 5 * count_mad,
    nCount_lower_cutoff = pmax(0, expm1(count_med - 3 * count_mad)),
    nCount_upper_cutoff = expm1(count_med + 3 * count_mad),
    feature_median_log = feature_med, feature_MAD_log = feature_mad,
    feature_lower_log = feature_med - 3 * feature_mad, feature_upper3_log = feature_med + 3 * feature_mad,
    feature_upper5_log = feature_med + 5 * feature_mad,
    nFeature_lower_cutoff = pmax(0, expm1(feature_med - 3 * feature_mad)),
    nFeature_upper_cutoff = expm1(feature_med + 3 * feature_mad),
    mt_median = mt_med, mt_MAD = mt_mad, mt_MAD_upper = mt_med + 3 * mt_mad,
    mt_final_cutoff = pmin(mt_med + 3 * mt_mad, 25), stringsAsFactors = FALSE
  )
}))

idx <- match(meta$orig.ident, thresholds$orig.ident)
meta$fail_nCount_low <- is.na(meta$log1p_nCount) | meta$log1p_nCount < thresholds$count_lower_log[idx]
meta$fail_nCount_high3 <- !is.na(meta$log1p_nCount) & meta$log1p_nCount > thresholds$count_upper3_log[idx]
meta$fail_nCount_high5 <- !is.na(meta$log1p_nCount) & meta$log1p_nCount > thresholds$count_upper5_log[idx]
meta$fail_nFeature_low <- is.na(meta$log1p_nFeature) | meta$log1p_nFeature < thresholds$feature_lower_log[idx]
meta$fail_nFeature_high3 <- !is.na(meta$log1p_nFeature) & meta$log1p_nFeature > thresholds$feature_upper3_log[idx]
meta$fail_nFeature_high5 <- !is.na(meta$log1p_nFeature) & meta$log1p_nFeature > thresholds$feature_upper5_log[idx]
meta$fail_mt <- is.na(meta$percent.mt.audit) | meta$percent.mt.audit > thresholds$mt_final_cutoff[idx]
meta$fail_hb <- is.na(meta$percent.hb.audit) | meta$percent.hb.audit > 5
meta$retain_A <- !(meta$fail_nCount_low | meta$fail_nCount_high3 | meta$fail_nFeature_low |
                     meta$fail_nFeature_high3 | meta$fail_mt | meta$fail_hb)
meta$retain_B <- !(meta$fail_nCount_low | meta$fail_nCount_high5 | meta$fail_nFeature_low |
                     meta$fail_nFeature_high5 | meta$fail_mt | meta$fail_hb)
meta$retain_C <- !(meta$fail_nCount_low | meta$fail_nFeature_low | meta$fail_mt | meta$fail_hb)
meta$rescue_B <- !meta$retain_A & meta$retain_B
meta$rescue_C <- !meta$retain_A & meta$retain_C

tail_by_orig <- do.call(rbind, lapply(orig_levels, function(group) {
  x <- meta[meta$orig.ident == group, , drop = FALSE]
  th <- thresholds[thresholds$orig.ident == group, , drop = FALSE]
  data.frame(
    orig.ident = group, dataset = unique(x$dataset), total_cells = nrow(x),
    nCount_lower_cutoff = th$nCount_lower_cutoff, nCount_upper_cutoff = th$nCount_upper_cutoff,
    fail_nCount_low_count = sum(x$fail_nCount_low), fail_nCount_high_count = sum(x$fail_nCount_high3),
    nFeature_lower_cutoff = th$nFeature_lower_cutoff, nFeature_upper_cutoff = th$nFeature_upper_cutoff,
    fail_nFeature_low_count = sum(x$fail_nFeature_low), fail_nFeature_high_count = sum(x$fail_nFeature_high3),
    stringsAsFactors = FALSE
  )
}))
write_tsv(tail_by_orig, "01_qc_tail_decomposition_by_orig_ident.tsv")
tail_global <- data.frame(
  total_cells = nrow(meta),
  fail_nCount_low_count = sum(meta$fail_nCount_low), fail_nCount_high_count = sum(meta$fail_nCount_high3),
  fail_nFeature_low_count = sum(meta$fail_nFeature_low), fail_nFeature_high_count = sum(meta$fail_nFeature_high3)
)
write_tsv(tail_global, "02_qc_tail_decomposition_global.tsv")

sensitivity_row <- function(x, group_name, group_value) {
  initial <- nrow(x); a <- sum(x$retain_A); b <- sum(x$retain_B); c <- sum(x$retain_C)
  data.frame(
    group_name = group_name, group = group_value, initial_cells = initial,
    Rule_A_retained = a, Rule_A_removed = initial - a, Rule_A_retention_percent = safe_rate(a, initial),
    Rule_B_retained = b, Rule_B_removed = initial - b, Rule_B_retention_percent = safe_rate(b, initial),
    Rule_C_retained = c, Rule_C_removed = initial - c, Rule_C_retention_percent = safe_rate(c, initial),
    B_minus_A_cells = b - a, C_minus_A_cells = c - a, stringsAsFactors = FALSE
  )
}
sens_orig <- do.call(rbind, lapply(orig_levels, function(group) sensitivity_row(meta[meta$orig.ident == group, ], "orig.ident", group)))
sens_orig$dataset <- thresholds$dataset[match(sens_orig$group, thresholds$orig.ident)]
sens_orig <- sens_orig[, c("dataset", setdiff(names(sens_orig), "dataset"))]
write_tsv(sens_orig, "03_qc_sensitivity_by_orig_ident.tsv")
dataset_levels <- sort(unique(meta$dataset))
sens_dataset <- do.call(rbind, lapply(dataset_levels, function(group) sensitivity_row(meta[meta$dataset == group, ], "dataset", group)))
write_tsv(sens_dataset, "04_qc_sensitivity_by_dataset.tsv")

rescue_reason <- function(x, rescue_col) {
  rescued <- x[[rescue_col]]
  count_high <- x$fail_nCount_high3
  feature_high <- x$fail_nFeature_high3
  c(
    nCount_high_only = sum(rescued & count_high & !feature_high),
    nFeature_high_only = sum(rescued & !count_high & feature_high),
    both_high = sum(rescued & count_high & feature_high),
    mixed_reasons = sum(rescued & !(xor(count_high, feature_high) | (count_high & feature_high)))
  )
}
global <- sensitivity_row(meta, "global", "ALL")
b_reasons <- rescue_reason(meta, "rescue_B"); c_reasons <- rescue_reason(meta, "rescue_C")
for (nm in names(b_reasons)) global[[paste0("A_removed_B_retained_", nm)]] <- b_reasons[[nm]]
for (nm in names(c_reasons)) global[[paste0("A_removed_C_retained_", nm)]] <- c_reasons[[nm]]
global$A_removed_B_retained_total <- sum(meta$rescue_B)
global$A_removed_C_retained_total <- sum(meta$rescue_C)
write_tsv(global, "05_qc_sensitivity_global.tsv")

# Inhouse mitochondrial naming audit.
inhouse <- meta[meta$dataset == "InhouseData", , drop = FALSE]
common_mt <- c("MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8", "MT-CYB")
ens_mask <- grepl("^ENS[A-Z]*G[0-9]+", inhouse_features)
candidate_mask <- grepl("(^MT-|^mt-|^Mt-|mitochond|(^|[-_.])(ND[1-6]|COX?[1-3]|ATP6|ATP8|CYTB)([-_.]|$))",
                        inhouse_features, ignore.case = TRUE, perl = TRUE)
candidates <- sort(unique(union(inhouse_features[candidate_mask], intersect(common_mt, inhouse_features))))
candidate_table <- data.frame(
  gene = head(candidates, 50),
  match_MT_upper = grepl("^MT-", head(candidates, 50)),
  match_mt_lower = grepl("^mt-", head(candidates, 50)),
  match_Mt_mixed = grepl("^Mt-", head(candidates, 50)),
  common_mt_gene = head(candidates, 50) %in% common_mt,
  ensembl_id = grepl("^ENS[A-Z]*G[0-9]+", head(candidates, 50)), stringsAsFactors = FALSE
)
write_tsv(candidate_table, "07_inhouse_mt_candidate_genes.tsv")
inhouse_stats <- c(
  min = min(inhouse$percent.mt.audit, na.rm = TRUE), median = median(inhouse$percent.mt.audit, na.rm = TRUE),
  mean = mean(inhouse$percent.mt.audit, na.rm = TRUE), p90 = quant(inhouse$percent.mt.audit, .90),
  p95 = quant(inhouse$percent.mt.audit, .95), p99 = quant(inhouse$percent.mt.audit, .99),
  max = max(inhouse$percent.mt.audit, na.rm = TRUE), unique_values = length(unique(inhouse$percent.mt.audit)),
  all_zero = all(inhouse$percent.mt.audit == 0, na.rm = TRUE)
)
inhouse_lines <- c(
  paste("dataset", "InhouseData", sep = "\t"),
  paste("n_features_in_inhouse_layers", length(inhouse_features), sep = "\t"),
  paste("genes_matching_^MT-", sum(grepl("^MT-", inhouse_features)), sep = "\t"),
  paste("genes_matching_^mt-", sum(grepl("^mt-", inhouse_features)), sep = "\t"),
  paste("genes_matching_^Mt-", sum(grepl("^Mt-", inhouse_features)), sep = "\t"),
  paste("ensembl_gene_id_count", sum(ens_mask), sep = "\t"),
  paste("ensembl_gene_id_percent", safe_rate(sum(ens_mask), length(inhouse_features)), sep = "\t"),
  paste("largely_ensembl_ids", mean(ens_mask) >= 0.5, sep = "\t"),
  paste("common_mt_gene", common_mt, common_mt %in% inhouse_features, sep = "\t"),
  paste("percent.mt", names(inhouse_stats), inhouse_stats, sep = "\t"),
  "first_50_suspected_mitochondrial_genes",
  if (nrow(candidate_table)) candidate_table$gene else "NONE"
)
writeLines(inhouse_lines, new_path("06_inhouse_mt_audit.txt"), useBytes = TRUE)

mt_orig <- do.call(rbind, lapply(orig_levels, function(group) {
  x <- meta[meta$orig.ident == group, ]; th <- thresholds[thresholds$orig.ident == group, ]
  data.frame(dataset = unique(x$dataset), orig.ident = group,
             n_MT_genes_matched = length(orig_mt_features[[group]]),
             median_percent_mt = median(x$percent.mt.audit), p95_percent_mt = quant(x$percent.mt.audit, .95),
             p99_percent_mt = quant(x$percent.mt.audit, .99), max_percent_mt = max(x$percent.mt.audit),
             MAD_upper_cutoff = th$mt_MAD_upper, final_cutoff = th$mt_final_cutoff,
             fail_mt_cells = sum(x$fail_mt), fail_mt_percent = safe_rate(sum(x$fail_mt), nrow(x)))
}))
write_tsv(mt_orig, "08_mt_qc_audit_by_orig_ident.tsv")
mt_dataset <- do.call(rbind, lapply(dataset_levels, function(group) {
  x <- meta[meta$dataset == group, ]; groups <- unique(x$orig.ident)
  data.frame(dataset = group, n_MT_genes_matched = length(dataset_mt_features[[group]]),
             median_percent_mt = median(x$percent.mt.audit), p95_percent_mt = quant(x$percent.mt.audit, .95),
             p99_percent_mt = quant(x$percent.mt.audit, .99), max_percent_mt = max(x$percent.mt.audit),
             MAD_upper_cutoff_min = min(thresholds$mt_MAD_upper[thresholds$orig.ident %in% groups]),
             MAD_upper_cutoff_max = max(thresholds$mt_MAD_upper[thresholds$orig.ident %in% groups]),
             final_cutoff_min = min(thresholds$mt_final_cutoff[thresholds$orig.ident %in% groups]),
             final_cutoff_max = max(thresholds$mt_final_cutoff[thresholds$orig.ident %in% groups]),
             fail_mt_cells = sum(x$fail_mt), fail_mt_percent = safe_rate(sum(x$fail_mt), nrow(x)))
}))
write_tsv(mt_dataset, "09_mt_qc_audit_by_dataset.tsv")

hb_rows <- c(
  lapply(dataset_levels, function(group) {
    x <- meta[meta$dataset == group, ]
    data.frame(level = "dataset", dataset = group, orig.ident = "ALL",
               n_HB_genes_present = length(dataset_hb_features[[group]]),
               median_percent_hb = median(x$percent.hb.audit), p95_percent_hb = quant(x$percent.hb.audit, .95),
               p99_percent_hb = quant(x$percent.hb.audit, .99), max_percent_hb = max(x$percent.hb.audit),
               above_5pct_cells = sum(x$percent.hb.audit > 5), above_5pct_percent = safe_rate(sum(x$percent.hb.audit > 5), nrow(x)))
  }),
  lapply(orig_levels, function(group) {
    x <- meta[meta$orig.ident == group, ]
    data.frame(level = "orig.ident", dataset = unique(x$dataset), orig.ident = group,
               n_HB_genes_present = length(orig_hb_features[[group]]),
               median_percent_hb = median(x$percent.hb.audit), p95_percent_hb = quant(x$percent.hb.audit, .95),
               p99_percent_hb = quant(x$percent.hb.audit, .99), max_percent_hb = max(x$percent.hb.audit),
               above_5pct_cells = sum(x$percent.hb.audit > 5), above_5pct_percent = safe_rate(sum(x$percent.hb.audit > 5), nrow(x)))
  })
)
hb_audit <- do.call(rbind, hb_rows)
write_tsv(hb_audit, "10_hb_qc_audit.tsv")

# Simple audit plots. Per-orig.ident cutoffs are overlaid within dataset facets.
theme_audit <- ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid = ggplot2::element_blank())
save_plot <- function(plot, stem, width = 11, height = 7) {
  pdf <- file.path(plot_dir, paste0(stem, ".pdf")); png <- file.path(plot_dir, paste0(stem, ".png"))
  if (file.exists(pdf) || file.exists(png)) stop("Refusing plot overwrite: ", stem, call. = FALSE)
  ggplot2::ggsave(pdf, plot, width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(png, plot, width = width, height = height, dpi = 180, limitsize = FALSE)
}
p_count <- ggplot2::ggplot(meta, ggplot2::aes(log1p_nCount, color = dataset)) +
  ggplot2::geom_density(linewidth = 0.5, show.legend = FALSE) +
  ggplot2::geom_vline(data = thresholds, ggplot2::aes(xintercept = count_lower_log), color = "#2166AC", alpha = .35) +
  ggplot2::geom_vline(data = thresholds, ggplot2::aes(xintercept = count_upper3_log), color = "#B2182B", alpha = .35) +
  ggplot2::geom_vline(data = thresholds, ggplot2::aes(xintercept = count_upper5_log), color = "#EF8A62", alpha = .35, linetype = 2) +
  ggplot2::facet_wrap(~dataset, scales = "free_y") + ggplot2::labs(title = "log1p(nCount_RNA): lower 3 MAD, upper 3/5 MAD", y = "Density") + theme_audit
save_plot(p_count, "01_log1p_nCount_density_thresholds")
p_feature <- ggplot2::ggplot(meta, ggplot2::aes(log1p_nFeature, color = dataset)) +
  ggplot2::geom_density(linewidth = 0.5, show.legend = FALSE) +
  ggplot2::geom_vline(data = thresholds, ggplot2::aes(xintercept = feature_lower_log), color = "#2166AC", alpha = .35) +
  ggplot2::geom_vline(data = thresholds, ggplot2::aes(xintercept = feature_upper3_log), color = "#B2182B", alpha = .35) +
  ggplot2::geom_vline(data = thresholds, ggplot2::aes(xintercept = feature_upper5_log), color = "#EF8A62", alpha = .35, linetype = 2) +
  ggplot2::facet_wrap(~dataset, scales = "free_y") + ggplot2::labs(title = "log1p(nFeature_RNA): lower 3 MAD, upper 3/5 MAD", y = "Density") + theme_audit
save_plot(p_feature, "02_log1p_nFeature_density_thresholds")
p_mt <- ggplot2::ggplot(meta, ggplot2::aes(dataset, percent.mt.audit, fill = dataset)) +
  ggplot2::geom_violin(scale = "width", trim = TRUE, show.legend = FALSE) +
  ggplot2::labs(title = "percent.mt by dataset", x = NULL, y = "percent.mt") + theme_audit
save_plot(p_mt, "03_percent_mt_by_dataset")
p_hb <- ggplot2::ggplot(meta, ggplot2::aes(dataset, percent.hb.audit, fill = dataset)) +
  ggplot2::geom_violin(scale = "width", trim = TRUE, show.legend = FALSE) +
  ggplot2::geom_hline(yintercept = 5, color = "#B2182B", linetype = 2) +
  ggplot2::labs(title = "percent.hb by dataset", x = NULL, y = "percent.hb") + theme_audit
save_plot(p_hb, "04_percent_hb_by_dataset")
retention_plot <- rbind(
  data.frame(dataset = sens_dataset$group, rule = "Rule A", retention = sens_dataset$Rule_A_retention_percent),
  data.frame(dataset = sens_dataset$group, rule = "Rule B", retention = sens_dataset$Rule_B_retention_percent),
  data.frame(dataset = sens_dataset$group, rule = "Rule C", retention = sens_dataset$Rule_C_retention_percent)
)
p_ret <- ggplot2::ggplot(retention_plot, ggplot2::aes(dataset, retention, fill = rule)) +
  ggplot2::geom_col(position = "dodge") + ggplot2::labs(title = "QC sensitivity", x = NULL, y = "Retention (%)") +
  ggplot2::coord_cartesian(ylim = c(0, 100)) + theme_audit
save_plot(p_ret, "05_rule_ABC_retention_comparison")

# Evidence-based recommendation without changing any rule.
a_rate <- global$Rule_A_retention_percent
b_gain <- global$B_minus_A_cells; c_gain <- global$C_minus_A_cells
b_gain_pct <- safe_rate(b_gain, nrow(meta)); c_gain_pct <- safe_rate(c_gain, nrow(meta))
upper_total <- tail_global$fail_nCount_high_count + tail_global$fail_nFeature_high_count
lower_total <- tail_global$fail_nCount_low_count + tail_global$fail_nFeature_low_count
most_sensitive <- head(sens_orig[order(-sens_orig$B_minus_A_cells), c("dataset", "group", "B_minus_A_cells", "C_minus_A_cells")], 10)
inhouse_n_mt <- length(dataset_mt_features[["InhouseData"]])
inhouse_all_zero <- isTRUE(as.logical(inhouse_stats[["all_zero"]]))
recommendation <- if (b_gain_pct >= 2) "Rule B" else "Rule A"
strict_answer <- if (b_gain_pct >= 5) "当前 Rule A 对 upper tail 明显偏严。" else if (b_gain_pct >= 2) "当前 Rule A 对 upper tail 有中度敏感性。" else "没有证据表明当前 Rule A 总体明显过严。"
mt_answer <- if (inhouse_n_mt == 0) {
  "是：InhouseData counts layers 中没有任何 ^MT- feature，MT=0 主要是 feature naming/mapping/feature omission 问题，不能解释为真实线粒体转录为零。"
} else if (inhouse_all_zero) {
  "很可能是 feature mapping 或 counts 内容问题：存在 ^MT- feature，但所有计算值为 0。"
} else {
  "否：存在 ^MT- features 且 percent.mt 并非全零；无 MT failures 可由分布与样本阈值解释。"
}
report <- c(
  "# QC sensitivity audit summary", "",
  paste0("- Audit ID: `", audit_id, "`"),
  paste0("- Immutable source checkpoint: `", checkpoint, "`"),
  paste0("- Cells audited: ", format(nrow(meta), big.mark = ",")),
  "- No subset, filtered object, checkpoint modification, DoubletFinder, PCA, Harmony, clustering, or UMAP was performed.", "",
  "## Direct answers", "",
  paste0("1. **Is current QC obviously too strict?** ", strict_answer, " Rule A retention = ", round(a_rate, 2), "%; Rule B adds ", format(b_gain, big.mark = ","), " cells (", round(b_gain_pct, 2), "% of input)."),
  paste0("2. **Lower versus upper tails:** summed flag counts (overlap possible) are lower = ", format(lower_total, big.mark = ","), " and upper = ", format(upper_total, big.mark = ","), ". See tail tables for non-overlapping interpretation."),
  paste0("3. **Rule B versus A:** +", format(b_gain, big.mark = ","), " retained cells."),
  paste0("4. **Rule C versus A:** +", format(c_gain, big.mark = ","), " retained cells (", round(c_gain_pct, 2), "% of input)."),
  "5. **Most upper-tail-sensitive groups:**", paste0("   - ", most_sensitive$dataset, " / ", most_sensitive$group, ": B+A gain ", most_sensitive$B_minus_A_cells),
  paste0("6. **Why InhouseData MT is zero/no-failure:** ", mt_answer),
  paste0("7. **Does HB 5% have little impact?** Yes. It flags ", sum(meta$fail_hb), " cells (", round(safe_rate(sum(meta$fail_hb), nrow(meta)), 3), "%)."),
  paste0("8. **Recommendation for future consideration:** ", recommendation, ". This audit does not apply the recommendation. Rule C is not preferred by default because it removes all upper-tail safeguards before downstream doublet assessment."), "",
  "## Important implementation note", "",
  "The formal v1 QC code explicitly used `stats::mad(..., constant = 1)`. This audit reproduces that exact implementation. It is stricter than R's default consistency-scaled MAD (`constant = 1.4826`) and should be considered when interpreting sensitivity.", "",
  paste0("Audit elapsed seconds: ", round(as.numeric(difftime(Sys.time(), audit_start, units = "secs")), 2))
)
writeLines(report, new_path("QC_sensitivity_audit_summary.md"), useBytes = TRUE)
writeLines(c(capture.output(sessionInfo()), "", paste("Completed", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))),
           new_path("sessionInfo_audit.txt"), useBytes = TRUE)
message("AUDIT COMPLETE: ", audit_dir)
