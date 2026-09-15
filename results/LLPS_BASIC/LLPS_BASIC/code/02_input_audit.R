source(file.path(dirname(normalizePath(sys.frame(1)$ofile)), "00_config.R"))

start_time <- Sys.time()
log_message("02_input_audit", "START", "n_files=", length(config$input_files))

missing_files <- config$input_files[!file.exists(config$input_files)]
if (length(missing_files) > 0L) {
  stop("Input files not found: ", paste(missing_files, collapse = "; "), call. = FALSE)
}

dataset_from_file <- function(path) {
  x <- basename(path)
  x <- sub("\\.(rdata|qs2?|RData)$", "", x, ignore.case = TRUE)
  sub("_scRNA(_original)?$", "", x, ignore.case = TRUE)
}

read_objects <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rdata") {
    holder <- new.env(parent = emptyenv())
    object_names <- load(path, envir = holder, verbose = FALSE)
    return(setNames(lapply(object_names, function(nm) holder[[nm]]), object_names))
  }
  if (ext == "qs") {
    return(list(qs_object = qs::qread(path, nthreads = config$max_threads)))
  }
  if (ext == "qs2") {
    return(list(qs2_object = qs2::qs_read(path, nthreads = config$max_threads)))
  }
  stop("Unsupported input extension: ", path, call. = FALSE)
}

audit_one_file <- function(path) {
  log_message("02_input_audit", "READ_START", path)
  objects <- read_objects(path)
  rows <- list()
  metadata_rows <- list()
  orig_rows <- list()

  for (object_name in names(objects)) {
    object <- objects[[object_name]]
    is_seurat <- inherits(object, "Seurat")
    object_class <- paste(class(object), collapse = ";")
    if (!is_seurat) {
      rows[[length(rows) + 1L]] <- data.frame(
        file_path = path, file_name = basename(path), file_type = tolower(tools::file_ext(path)),
        dataset_name = dataset_from_file(path), object_name = object_name,
        object_class = object_class, is_seurat = FALSE, seurat_object_version = NA_character_,
        assays = NA_character_, has_RNA_assay = FALSE, RNA_layers = NA_character_,
        rownames_type = NA_character_, n_cells = NA_integer_, n_features = NA_integer_,
        duplicate_barcodes_within_object = NA_integer_, has_nFeature_RNA = NA,
        has_nCount_RNA = NA, has_orig_ident = NA, n_orig_ident = NA_integer_,
        existing_dataset_field = NA, audit_status = "NON_SEURAT_OBJECT",
        stringsAsFactors = FALSE
      )
      next
    }

    assays <- names(object@assays)
    has_rna <- "RNA" %in% assays
    layers <- if (has_rna) {
      paste(SeuratObject::Layers(object[["RNA"]]), collapse = ";")
    } else {
      NA_character_
    }
    features <- if (has_rna) rownames(object[["RNA"]]) else rownames(object)
    meta_names <- colnames(object[[]])
    orig_values <- if ("orig.ident" %in% meta_names) as.character(object$orig.ident) else character()

    rows[[length(rows) + 1L]] <- data.frame(
      file_path = path, file_name = basename(path), file_type = tolower(tools::file_ext(path)),
      dataset_name = dataset_from_file(path), object_name = object_name,
      object_class = object_class, is_seurat = TRUE,
      seurat_object_version = as.character(object@version), assays = paste(assays, collapse = ";"),
      has_RNA_assay = has_rna, RNA_layers = layers,
      rownames_type = typeof(features), n_cells = ncol(object), n_features = length(features),
      duplicate_barcodes_within_object = sum(duplicated(colnames(object))),
      has_nFeature_RNA = "nFeature_RNA" %in% meta_names,
      has_nCount_RNA = "nCount_RNA" %in% meta_names,
      has_orig_ident = "orig.ident" %in% meta_names,
      n_orig_ident = length(unique(orig_values)),
      existing_dataset_field = "dataset" %in% meta_names,
      audit_status = if (has_rna) "OK" else "BLOCKING_NO_RNA_ASSAY",
      stringsAsFactors = FALSE
    )

    metadata_rows[[length(metadata_rows) + 1L]] <- data.frame(
      dataset_name = dataset_from_file(path), object_name = object_name,
      metadata_field = meta_names, stringsAsFactors = FALSE
    )
    if (length(orig_values) > 0L) {
      tab <- table(orig_values, useNA = "ifany")
      orig_rows[[length(orig_rows) + 1L]] <- data.frame(
        dataset_name = dataset_from_file(path), object_name = object_name,
        orig_ident = names(tab), n_cells = as.integer(tab), stringsAsFactors = FALSE
      )
    }
  }
  log_message("02_input_audit", "READ_END", path, "objects=", length(objects))
  list(
    audit = do.call(rbind, rows),
    metadata = if (length(metadata_rows)) do.call(rbind, metadata_rows) else data.frame(),
    orig_ident = if (length(orig_rows)) do.call(rbind, orig_rows) else data.frame()
  )
}

results <- lapply(config$input_files, audit_one_file)
audit <- do.call(rbind, lapply(results, `[[`, "audit"))
metadata_fields <- do.call(rbind, lapply(results, `[[`, "metadata"))
orig_ident_counts <- do.call(rbind, lapply(results, `[[`, "orig_ident"))

seurat_manifest <- audit[audit$is_seurat, c(
  "file_path", "file_name", "file_type", "dataset_name", "object_name", "n_cells", "n_features"
), drop = FALSE]

if (nrow(seurat_manifest) != length(config$input_files)) {
  stop(
    "Expected exactly one Seurat object per selected file; found ", nrow(seurat_manifest),
    " across ", length(config$input_files), " files. Review audit before proceeding.",
    call. = FALSE
  )
}
if (anyDuplicated(seurat_manifest$dataset_name)) {
  stop("Dataset names derived from file names are not unique.", call. = FALSE)
}

write_tsv_new(seurat_manifest, file.path(paths$audit, "input_manifest.tsv"))
write_tsv_new(audit, file.path(paths$audit, "input_object_audit.tsv"))
write_tsv_new(metadata_fields, file.path(paths$audit, "metadata_fields_by_dataset.tsv"))
write_tsv_new(orig_ident_counts, file.path(paths$audit, "orig_ident_cell_counts.tsv"))

blocking <- audit$audit_status != "OK"
if (any(blocking)) {
  stop("Blocking input audit findings detected; review input_object_audit.tsv", call. = FALSE)
}
log_message("02_input_audit", "END", "seurat_objects=", nrow(seurat_manifest),
            "total_cells=", sum(seurat_manifest$n_cells),
            "elapsed_seconds=", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 3))
