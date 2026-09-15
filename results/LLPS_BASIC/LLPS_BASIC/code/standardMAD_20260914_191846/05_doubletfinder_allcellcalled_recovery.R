source(file.path(dirname(normalizePath(sys.frame(1)$ofile)), "00_config.R"))
start_time <- Sys.time()
log_message("05_doubletfinder_allCellCalled_recovery", "START", "strategy=validate_54_then_run_missing")

qc_path <- "/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/checkpoints/01b_QCflags_only_noFiltering_20260914_195631.qs2"
qc_object <- qs2::qs_read(qc_path, nthreads = config$max_threads)
preqc <- read.delim(file.path(paths$audit, "orig_ident_cell_counts.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
if (anyDuplicated(preqc$orig_ident)) stop("orig.ident is not globally unique across datasets; cannot map pre-QC counts safely", call. = FALSE)

df_output_dir <- file.path(paths$doubletfinder, "allCellCalled_20260914_195247")
if (!dir.exists(df_output_dir)) stop("Expected new all-cell-called DoubletFinder directory is absent: ", df_output_dir, call. = FALSE)
cell_output_dir <- file.path(df_output_dir, "cell_classifications")
if (file.exists(cell_output_dir) && !dir.exists(cell_output_dir)) stop("Classification directory path is a file", call. = FALSE)
if (!dir.exists(cell_output_dir)) dir.create(cell_output_dir, recursive = FALSE, showWarnings = FALSE)

sanitize <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
orig_levels <- sort(unique(as.character(qc_object$orig.ident)))
parameter_rows <- vector("list", length(orig_levels))
summary_rows <- vector("list", length(orig_levels))
classification_paths <- character(length(orig_levels))

process_sample <- function(orig_ident, index) {
  sample_start <- Sys.time()
  log_message("05_doubletfinder", "SAMPLE_START", index, "of", length(orig_levels), orig_ident)
  warning_messages <- character()
  result <- tryCatch(
    withCallingHandlers({
      cells <- colnames(qc_object)[as.character(qc_object$orig.ident) == orig_ident]
      current_n <- length(cells)
      preqc_n <- preqc$n_cells[match(orig_ident, preqc$orig_ident)]
      if (length(preqc_n) != 1L || is.na(preqc_n)) stop("Missing pre-QC cell count")
      if (current_n < 500L) stop("Too few QC-pass cells for reliable DoubletFinder paramSweep: ", current_n)

      sample_object <- subset(qc_object, cells = cells)
      Seurat::DefaultAssay(sample_object) <- "RNA"
      count_layers <- grep("^counts", SeuratObject::Layers(sample_object[["RNA"]]), value = TRUE)
      if (!length(count_layers)) stop("No counts layer after sample subset")
      join_applied <- length(count_layers) > 1L || !identical(count_layers, "counts")
      if (join_applied) {
        sample_object <- SeuratObject::JoinLayers(
          object = sample_object, assay = "RNA", layers = count_layers, new = "counts"
        )
      }
      log_message("05_doubletfinder", "JOINLAYERS", orig_ident, "applied=", join_applied,
                  "step=DoubletFinder_sample_preparation")

      set.seed(config$seed + index)
      sample_object <- Seurat::NormalizeData(
        object = sample_object, assay = "RNA", normalization.method = config$normalization_method,
        scale.factor = config$normalization_scale_factor, margin = 1, verbose = FALSE
      )
      sample_object <- Seurat::FindVariableFeatures(
        object = sample_object, assay = "RNA", selection.method = config$hvg_method,
        nfeatures = min(config$hvg_n, nrow(sample_object[["RNA"]])), verbose = FALSE
      )
      hvg <- SeuratObject::VariableFeatures(sample_object)
      if (length(hvg) < 100L) stop("Fewer than 100 variable features found")
      sample_object <- Seurat::ScaleData(
        object = sample_object, assay = "RNA", features = hvg,
        do.scale = TRUE, do.center = TRUE, verbose = FALSE
      )
      max_pcs <- min(config$pca_npcs, length(hvg) - 1L, current_n - 1L)
      if (max_pcs < max(config$pcs)) stop("Insufficient PCs available: ", max_pcs)
      sample_object <- Seurat::RunPCA(
        object = sample_object, assay = "RNA", features = hvg,
        npcs = config$pca_npcs, seed.use = config$seed + index,
        approx = TRUE, verbose = FALSE
      )
      sample_object <- Seurat::FindNeighbors(
        object = sample_object, reduction = "pca", dims = config$pcs,
        k.param = 20, compute.SNN = TRUE, prune.SNN = 1 / 15,
        nn.method = "annoy", annoy.metric = "euclidean",
        n.trees = 50, verbose = FALSE
      )
      sample_object <- Seurat::FindClusters(
        object = sample_object, resolution = config$cluster_resolution,
        algorithm = config$cluster_algorithm, random.seed = config$seed + index,
        verbose = FALSE
      )

      sweep <- DoubletFinder::paramSweep(
        seu = sample_object, PCs = config$pcs, sct = FALSE, num.cores = 1
      )
      sweep_stats <- DoubletFinder::summarizeSweep(sweep.list = sweep, GT = FALSE, GT.calls = NULL)
      bcmvn <- DoubletFinder::find.pK(sweep_stats)
      if (!nrow(bcmvn) || all(is.na(bcmvn$BCmetric))) stop("paramSweep produced no valid BCmetric")
      best <- which.max(bcmvn$BCmetric)
      chosen_pK <- suppressWarnings(as.numeric(as.character(bcmvn$pK[[best]])))
      if (!is.finite(chosen_pK)) stop("Could not parse selected pK")

      expected_rate <- min(0.008 * preqc_n / 1000, 0.08)
      nExp_poi <- round(current_n * expected_rate)
      clusters <- as.character(Seurat::Idents(sample_object))
      homotypic_prop <- DoubletFinder::modelHomotypic(clusters)
      nExp_adj <- round(nExp_poi * (1 - homotypic_prop))
      if (nExp_adj < 1L) stop("Homotypic-adjusted expected doublets is below 1")

      before_cols <- colnames(sample_object[[]])
      sample_object <- DoubletFinder::doubletFinder(
        seu = sample_object, PCs = config$pcs, pN = config$doubletfinder_pN,
        pK = chosen_pK, nExp = nExp_adj, reuse.pANN = NULL,
        sct = FALSE, annotations = NULL
      )
      added_cols <- setdiff(colnames(sample_object[[]]), before_cols)
      class_col <- grep("^DF.classifications_", added_cols, value = TRUE)
      pann_col <- grep("^pANN_", added_cols, value = TRUE)
      if (length(class_col) != 1L) stop("Expected one new DF.classifications column; found: ", paste(class_col, collapse = ","))
      if (length(pann_col) != 1L) stop("Expected one new pANN column; found: ", paste(pann_col, collapse = ","))
      classification <- as.character(sample_object[[class_col, drop = TRUE]])
      pann <- as.numeric(sample_object[[pann_col, drop = TRUE]])
      if (!all(classification %in% c("Singlet", "Doublet"))) stop("Unexpected DoubletFinder classification values")

      cell_table <- data.frame(
        cell_id = colnames(sample_object), orig.ident = orig_ident,
        dataset = as.character(sample_object$dataset),
        df_classification_column = class_col,
        df_classification = classification,
        df_pANN_column = pann_col, df_pANN = pann,
        stringsAsFactors = FALSE
      )
      cell_path <- file.path(cell_output_dir, paste0(sanitize(orig_ident), "_", run_stamp, ".tsv"))
      write_tsv_new(cell_table, cell_path)
      if (length(warning_messages)) {
        warning_path <- file.path(cell_output_dir, paste0(sanitize(orig_ident), "_warnings_", run_stamp, ".txt"))
        write_lines_new(unique(warning_messages), warning_path)
      }
      list(
        parameter = data.frame(
          orig.ident = orig_ident, dataset = unique(as.character(sample_object$dataset)),
          n_cells_preQC = preqc_n, n_cells = current_n,
          expected_rate = expected_rate, pN = config$doubletfinder_pN,
          chosen_pK = chosen_pK, nExp_poi = nExp_poi,
          homotypic.prop = homotypic_prop, nExp_adj = nExp_adj,
          homotypic_source = "direct_modelHomotypic",
          preliminary_HVG = length(hvg), preliminary_PCs = paste(range(config$pcs), collapse = ":"),
          preliminary_resolution = config$cluster_resolution,
          JoinLayers_applied = join_applied, warning_count = length(unique(warning_messages)),
          elapsed_seconds = as.numeric(difftime(Sys.time(), sample_start, units = "secs")),
          stringsAsFactors = FALSE
        ),
        summary = data.frame(
          orig.ident = orig_ident, doublet_number = sum(classification == "Doublet"),
          singlet_number = sum(classification == "Singlet"), stringsAsFactors = FALSE
        ),
        cell_path = cell_path
      )
    }, warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      error_path <- file.path(cell_output_dir, paste0(sanitize(orig_ident), "_ERROR_", run_stamp, ".txt"))
      write_lines_new(c(
        paste("orig.ident", orig_ident, sep = "\t"),
        paste("time", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), sep = "\t"),
        paste("error", conditionMessage(e), sep = "\t"),
        if (length(warning_messages)) paste("warnings", paste(unique(warning_messages), collapse = " | "), sep = "\t")
      ), error_path)
      stop("DoubletFinder failed for orig.ident ", orig_ident, ": ", conditionMessage(e), call. = FALSE)
    }
  )
  log_message("05_doubletfinder", "SAMPLE_END", index, "of", length(orig_levels), orig_ident,
              "doublets=", result$summary$doublet_number, "singlets=", result$summary$singlet_number,
              "elapsed_seconds=", round(result$parameter$elapsed_seconds, 3))
  result
}

recover_existing <- function(orig_ident, index) {
  safe_name <- sanitize(orig_ident)
  existing <- file.path(cell_output_dir, paste0(safe_name, "_20260914_195744.tsv"))
  if (!file.exists(existing)) return(NULL)
  tab <- read.delim(existing, stringsAsFactors = FALSE, check.names = FALSE)
  expected_cells <- colnames(qc_object)[as.character(qc_object$orig.ident) == orig_ident]
  if (nrow(tab) != length(expected_cells) || !setequal(tab$cell_id, expected_cells)) {
    stop("Existing classification validation failed for ", orig_ident, call. = FALSE)
  }
  if (!all(tab$df_classification %in% c("Singlet", "Doublet"))) {
    stop("Invalid existing classifications for ", orig_ident, call. = FALSE)
  }
  class_names <- unique(tab$df_classification_column)
  if (length(class_names) != 1L) stop("Non-unique classification column for ", orig_ident, call. = FALSE)
  parts <- strsplit(class_names, "_", fixed = TRUE)[[1L]]
  if (length(parts) < 4L) stop("Cannot parse DoubletFinder column for ", orig_ident, call. = FALSE)
  chosen_pK <- as.numeric(parts[[length(parts) - 1L]])
  nExp_adj <- as.integer(parts[[length(parts)]])
  preqc_n <- preqc$n_cells[match(orig_ident, preqc$orig_ident)]
  current_n <- nrow(tab)
  expected_rate <- min(0.008 * preqc_n / 1000, 0.08)
  nExp_poi <- round(current_n * expected_rate)
  homotypic_reconstructed <- 1 - nExp_adj / nExp_poi
  log_message("05_doubletfinder_allCellCalled_recovery", "VALIDATED_EXISTING", index, "of", length(orig_levels), orig_ident)
  list(
    parameter = data.frame(
      orig.ident = orig_ident, dataset = unique(tab$dataset), n_cells_preQC = preqc_n,
      n_cells = current_n, expected_rate = expected_rate, pN = config$doubletfinder_pN,
      chosen_pK = chosen_pK, nExp_poi = nExp_poi,
      homotypic.prop = homotypic_reconstructed, nExp_adj = nExp_adj,
      homotypic_source = "reconstructed_from_rounded_nExp_adj_after_interrupted_parent_run",
      preliminary_HVG = NA_integer_, preliminary_PCs = "1:30", preliminary_resolution = config$cluster_resolution,
      JoinLayers_applied = TRUE, warning_count = 0L, elapsed_seconds = NA_real_, stringsAsFactors = FALSE
    ),
    summary = data.frame(orig.ident = orig_ident,
      doublet_number = sum(tab$df_classification == "Doublet"),
      singlet_number = sum(tab$df_classification == "Singlet"), stringsAsFactors = FALSE),
    cell_path = existing
  )
}

for (i in seq_along(orig_levels)) {
  result <- recover_existing(orig_levels[[i]], i)
  if (is.null(result)) result <- process_sample(orig_levels[[i]], i)
  parameter_rows[[i]] <- result$parameter
  summary_rows[[i]] <- result$summary
  classification_paths[[i]] <- result$cell_path
}

parameters <- do.call(rbind, parameter_rows)
summaries <- do.call(rbind, summary_rows)
summaries <- merge(
  parameters[, c("orig.ident", "dataset", "n_cells_preQC", "n_cells", "expected_rate",
                 "chosen_pK", "nExp_poi", "homotypic.prop", "nExp_adj")],
  summaries, by = "orig.ident", all = TRUE, sort = FALSE
)
write_tsv_new(parameters, file.path(df_output_dir, "doubletfinder_allCellCalled_parameters_recovered.tsv"))
write_tsv_new(summaries, file.path(df_output_dir, "doubletfinder_allCellCalled_summary_recovered.tsv"))

all_classifications <- do.call(rbind, lapply(classification_paths, function(path) {
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}))
if (nrow(all_classifications) != ncol(qc_object)) stop("Combined classification row count does not equal QC cells", call. = FALSE)
if (anyDuplicated(all_classifications$cell_id)) stop("Duplicate cell IDs in combined DoubletFinder results", call. = FALSE)
qc_object$DF_allCellCalled_classification <- all_classifications$df_classification[match(colnames(qc_object), all_classifications$cell_id)]
qc_object$DF_allCellCalled_pANN <- all_classifications$df_pANN[match(colnames(qc_object), all_classifications$cell_id)]
if (anyNA(qc_object$DF_allCellCalled_classification) || anyNA(qc_object$DF_allCellCalled_pANN)) {
  stop("Failed to map DoubletFinder results into all-cell-called object", call. = FALSE)
}
if (ncol(qc_object) != 410466L) stop("All-cell-called object cell count changed", call. = FALSE)

checkpoint <- timestamped_path(paths$checkpoints, "01c_DoubletFinder_allCellCalled", "qs2")
save_qs2_new(qc_object, checkpoint)
write_lines_new(checkpoint, file.path(paths$logs, paste0("01c_allCellCalled_checkpoint_path_", run_stamp, ".txt")))
write_lines_new(checkpoint, file.path(df_output_dir, "doubletfinder_allCellCalled_checkpoint_path_recovered.txt"))
log_message("05_doubletfinder_allCellCalled_recovery", "END", "all_cells=", ncol(qc_object),
            "singlets=", sum(all_classifications$df_classification == "Singlet"),
            "doublets=", sum(all_classifications$df_classification == "Doublet"),
            "checkpoint=", checkpoint,
            "elapsed_seconds=", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 3))
