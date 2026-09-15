source(file.path(dirname(normalizePath(sys.frame(1)$ofile)), "00_config.R"))

manifest <- read.delim(
  file.path(paths$audit, "input_manifest.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

read_selected_object <- function(path, object_name) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rdata") {
    holder <- new.env(parent = emptyenv())
    loaded <- load(path, envir = holder, verbose = FALSE)
    if (!object_name %in% loaded) stop("Manifest object missing from RData: ", object_name, call. = FALSE)
    return(holder[[object_name]])
  }
  if (ext == "qs") return(qs::qread(path, nthreads = config$max_threads))
  if (ext == "qs2") return(qs2::qs_read(path, nthreads = config$max_threads))
  stop("Unsupported extension: ", ext, call. = FALSE)
}

barcode_sets <- lapply(seq_len(nrow(manifest)), function(i) {
  log_message("02b_barcode_overlap_audit", "READ_START", manifest$file_path[[i]])
  object <- read_selected_object(manifest$file_path[[i]], manifest$object_name[[i]])
  cells <- colnames(object)
  log_message("02b_barcode_overlap_audit", "READ_END", manifest$file_path[[i]], "cells=", length(cells))
  cells
})
names(barcode_sets) <- manifest$dataset_name

pairs <- utils::combn(names(barcode_sets), 2L, simplify = FALSE)
overlap <- do.call(rbind, lapply(pairs, function(pair) {
  shared <- intersect(barcode_sets[[pair[[1L]]]], barcode_sets[[pair[[2L]]]])
  data.frame(
    dataset_1 = pair[[1L]], dataset_2 = pair[[2L]],
    n_shared_barcodes = length(shared),
    prefix_required = length(shared) > 0L,
    stringsAsFactors = FALSE
  )
}))

write_tsv_new(
  overlap,
  file.path(paths$audit, "barcode_overlap_by_dataset_20260914_151500.tsv")
)
log_message(
  "02b_barcode_overlap_audit", "END",
  "pairs_with_overlap=", sum(overlap$n_shared_barcodes > 0L),
  "max_shared=", max(overlap$n_shared_barcodes)
)
