source("/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/code/standardMAD_20260914_191846/00_config.R")
start_time <- Sys.time()
stamp <- format(start_time, "%Y%m%d_%H%M%S")
stage_dir <- file.path(config$run_dir, paste0("QC_DoubletFinder_overlap_", stamp))
if (file.exists(stage_dir) || dir.exists(stage_dir)) stop("Refusing existing stage directory: ", stage_dir)
dir.create(stage_dir, recursive = FALSE, showWarnings = FALSE)
if (!dir.exists(stage_dir)) stop("Could not create stage directory")

input_path <- "/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/checkpoints/01c_DoubletFinder_allCellCalled_20260914_231016.qs2"
message(format(Sys.time()), " START reading ", input_path)
obj <- qs2::qs_read(input_path, nthreads = config$max_threads)
if (ncol(obj) != 410466L) stop("Unexpected all-cell-called count: ", ncol(obj))
md <- obj[[]]
flags <- c("fail_nCount_low", "fail_nCount_high", "fail_nFeature_low", "fail_nFeature_high", "fail_mt", "fail_hb")
missing <- setdiff(c("dataset", "orig.ident", flags, "DF_allCellCalled_classification"), colnames(md))
if (length(missing)) stop("Missing metadata: ", paste(missing, collapse = ", "))
for (z in flags) if (anyNA(md[[z]])) stop("NA in ", z)
if (!all(md$DF_allCellCalled_classification %in% c("Singlet", "Doublet"))) stop("Invalid DF classifications")
md$qc_fail_any <- Reduce(`|`, md[flags])
md$df_doublet <- md$DF_allCellCalled_classification == "Doublet"
md$both_upper_high <- md$fail_nCount_high & md$fail_nFeature_high
md$upper_high_any <- md$fail_nCount_high | md$fail_nFeature_high
md$final_keep <- !md$qc_fail_any & !md$df_doublet

flag_summary <- do.call(rbind, lapply(c(flags, "both_upper_high", "upper_high_any", "qc_fail_any"), function(z) {
  hit <- md[[z]]
  data.frame(category=z, cells=sum(hit), doublets=sum(hit & md$df_doublet), singlets=sum(hit & !md$df_doublet),
             doublet_percentage=if (sum(hit)) 100*sum(hit & md$df_doublet)/sum(hit) else NA_real_)
}))
write_tsv_new(flag_summary, file.path(stage_dir, "QC_DoubletFinder_overlap.tsv"))

summarize_group <- function(d, level, value) data.frame(
  level=level, group=value, original_cells=nrow(d), QC_only_fail_cells=sum(d$qc_fail_any & !d$df_doublet),
  Doublet_only_cells=sum(!d$qc_fail_any & d$df_doublet), QC_Doublet_overlap_cells=sum(d$qc_fail_any & d$df_doublet),
  all_QC_fail_cells=sum(d$qc_fail_any), all_Doublets=sum(d$df_doublet), final_retained_cells=sum(d$final_keep),
  final_retention_percentage=100*mean(d$final_keep), stringsAsFactors=FALSE)
global <- summarize_group(md, "global", "all")
by_dataset <- do.call(rbind, lapply(split(md, md$dataset), function(d) summarize_group(d, "dataset", unique(as.character(d$dataset)))))
by_orig <- do.call(rbind, lapply(split(md, md$orig.ident), function(d) summarize_group(d, "orig.ident", unique(as.character(d$orig.ident)))))
write_tsv_new(global, file.path(stage_dir, "final_filter_summary_global.tsv"))
write_tsv_new(by_dataset, file.path(stage_dir, "final_filter_summary_by_dataset.tsv"))
write_tsv_new(by_orig, file.path(stage_dir, "final_filter_summary_by_orig_ident.tsv"))

report <- c(
  "# QC–DoubletFinder overlap summary", "",
  paste("- Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste("- Input:", input_path), paste("- Original cells:", nrow(md)),
  paste("- QC failures:", sum(md$qc_fail_any)), paste("- Doublets:", sum(md$df_doublet)),
  paste("- QC-only failures:", sum(md$qc_fail_any & !md$df_doublet)),
  paste("- Doublet-only cells:", sum(!md$qc_fail_any & md$df_doublet)),
  paste("- QC + Doublet overlap:", sum(md$qc_fail_any & md$df_doublet)),
  paste("- Final retained:", sum(md$final_keep)), "",
  "InhouseData has mt_qc_available=FALSE and fail_mt=FALSE; this is a documented assay limitation, not evidence of biological MT=0."
)
write_lines_new(report, file.path(stage_dir, "QC_DoubletFinder_overlap_summary.md"))

obj$qc_fail_any <- md$qc_fail_any
obj$DF_allCellCalled_doublet <- md$df_doublet
obj$final_keep <- md$final_keep
final_obj <- subset(obj, cells = rownames(md)[md$final_keep])
if (ncol(final_obj) != sum(md$final_keep)) stop("Final subset count mismatch")
checkpoint <- timestamped_path(paths$checkpoints, "01d_final_QC_DF_singlets", "qs2")
save_qs2_new(final_obj, checkpoint)
write_lines_new(checkpoint, file.path(stage_dir, "01d_checkpoint_path.txt"))
write_lines_new(c(report, "", paste("- Output checkpoint:", checkpoint), paste("- Elapsed seconds:", round(as.numeric(difftime(Sys.time(), start_time, units="secs")), 3))),
                file.path(stage_dir, "phase_completion_summary.md"))
message(format(Sys.time()), " END final cells=", ncol(final_obj), " checkpoint=", checkpoint)
