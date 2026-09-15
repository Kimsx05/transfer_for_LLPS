# LLPS reusable scRNA-seq pipeline v1: centralized configuration

config <- list(
  run_id = "20260914_151017",
  seed = 123L,
  max_threads = 8L,
  isolated_library = "/home/jinshengxi/Rlibrary/LLPS_pipeline_v1_20260914_150459",
  shared_library = "/home/jinshengxi/Rlibrary/packages",
  run_dir = "/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017",
  input_files = c(
    "/home/jinshengxi/LLPS/1.1 单细胞数据/HRA003620_scRNA.qs",
    "/home/jinshengxi/LLPS/1.1 单细胞数据/GSE135337_scRNA_original.rdata",
    "/home/jinshengxi/LLPS/1.1 单细胞数据/GSE222315_scRNA_original.rdata",
    "/home/jinshengxi/LLPS/1.1 单细胞数据/CNP0000460_scRNA_original.rdata",
    "/home/jinshengxi/LLPS/1.1 单细胞数据/InhouseData_scRNA_original.rdata"
  ),
  hvg_n = 3000L,
  pcs = 1:30,
  cluster_resolution = 0.8,
  doubletfinder_pN = 0.25,
  qc_mad_multiplier = 3,
  mt_hard_ceiling = 25,
  hb_cutoff = 5,
  marker_only_pos = TRUE,
  marker_min_pct = 0.1,
  marker_logfc_threshold = 0.25
)

.libPaths(unique(c(config$isolated_library, config$shared_library, .libPaths())))
Sys.setenv(
  OMP_NUM_THREADS = config$max_threads,
  OPENBLAS_NUM_THREADS = config$max_threads,
  MKL_NUM_THREADS = config$max_threads,
  VECLIB_MAXIMUM_THREADS = config$max_threads,
  NUMEXPR_NUM_THREADS = config$max_threads,
  RCPP_PARALLEL_NUM_THREADS = config$max_threads
)
options(future.globals.maxSize = 100 * 1024^3)
set.seed(config$seed)

paths <- list(
  logs = file.path(config$run_dir, "00_logs"),
  audit = file.path(config$run_dir, "01_audit"),
  qc = file.path(config$run_dir, "02_qc"),
  doubletfinder = file.path(config$run_dir, "03_doubletfinder"),
  pca = file.path(config$run_dir, "04_pca"),
  harmony = file.path(config$run_dir, "05_harmony"),
  markers = file.path(config$run_dir, "06_markers"),
  annotation = file.path(config$run_dir, "07_annotation"),
  checkpoints = file.path(config$run_dir, "checkpoints"),
  code = file.path(config$run_dir, "code")
)

assert_new_file <- function(path) {
  if (file.exists(path) || dir.exists(path)) {
    stop("Refusing to overwrite existing path: ", path, call. = FALSE)
  }
  invisible(path)
}

write_tsv_new <- function(x, path) {
  assert_new_file(path)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

write_lines_new <- function(x, path) {
  assert_new_file(path)
  writeLines(x, path, useBytes = TRUE)
}

log_message <- function(step, ...) {
  path <- file.path(paths$logs, paste0(step, ".log"))
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), paste(..., collapse = " "), sep = "\t")
  cat(line, "\n", file = path, append = TRUE)
  message(line)
}
