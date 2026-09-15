source(file.path(dirname(normalizePath(sys.frame(1)$ofile)), "00_config.R"))

required_packages <- c(
  "Seurat", "SeuratObject", "harmony", "DoubletFinder", "SCP", "qs", "qs2",
  "Matrix", "dplyr", "ggplot2", "patchwork", "RhpcBLASctl", "future"
)

start_time <- Sys.time()
log_message("01_environment_check", "START")

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(config$max_threads)
  RhpcBLASctl::omp_set_num_threads(config$max_threads)
}
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = config$max_threads)
}
if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
}

pkg_rows <- lapply(required_packages, function(pkg) {
  location <- suppressWarnings(find.package(pkg, quiet = TRUE))
  installed <- length(location) == 1L && nzchar(location)
  data.frame(
    package = pkg,
    installed = installed,
    version = if (installed) as.character(utils::packageVersion(pkg, lib.loc = dirname(location))) else NA_character_,
    library_path = if (installed) location else NA_character_,
    stringsAsFactors = FALSE
  )
})
pkg_table <- do.call(rbind, pkg_rows)
write_tsv_new(pkg_table, file.path(paths$logs, "package_versions.tsv"))

missing <- pkg_table$package[!pkg_table$installed]
if (length(missing) > 0L) {
  stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

api_checks <- c(
  "DoubletFinder::paramSweep", "DoubletFinder::summarizeSweep",
  "DoubletFinder::find.pK", "DoubletFinder::doubletFinder",
  "DoubletFinder::modelHomotypic", "qs2::qs_save", "qs2::qs_read",
  "SCP::CellDimPlot", "harmony::RunHarmony"
)
api_ok <- vapply(strsplit(api_checks, "::", fixed = TRUE), function(parts) {
  parts[[2L]] %in% getNamespaceExports(parts[[1L]])
}, logical(1))
if (!all(api_ok)) {
  stop("Required APIs missing: ", paste(api_checks[!api_ok], collapse = ", "), call. = FALSE)
}

host <- Sys.info()[["nodename"]]
disk <- system2("df", c("-h", shQuote(config$run_dir)), stdout = TRUE, stderr = TRUE)
env_lines <- c(
  paste("start_time", format(start_time, "%Y-%m-%d %H:%M:%S %z"), sep = "\t"),
  paste("hostname", host, sep = "\t"),
  paste("user", Sys.info()[["user"]], sep = "\t"),
  paste("working_directory", getwd(), sep = "\t"),
  paste("R_version", R.version.string, sep = "\t"),
  paste("libPaths", paste(.libPaths(), collapse = " | "), sep = "\t"),
  paste("max_threads", config$max_threads, sep = "\t"),
  "disk_space",
  disk,
  "api_checks",
  paste(api_checks, api_ok, sep = "\t")
)
write_lines_new(env_lines, file.path(paths$logs, "environment_check.txt"))

session_path <- file.path(paths$logs, "sessionInfo_initial.txt")
assert_new_file(session_path)
session_lines <- capture.output(sessionInfo())
writeLines(session_lines, session_path, useBytes = TRUE)
log_message("01_environment_check", "END", "elapsed_seconds=", round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 3))
