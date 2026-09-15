# Read-only standard-MAD sensitivity audit. Does not subset or save any object.
stamp <- "20260914_191846"
root <- "/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017"
checkpoint <- file.path(root, "checkpoints", "01_merged_full_20260914_153425.qs2")
outdir <- file.path(root, paste0("QC_sensitivity_standardMAD_", stamp))
.libPaths(unique(c("/home/jinshengxi/Rlibrary/LLPS_pipeline_v1_20260914_150459",
                   "/home/jinshengxi/Rlibrary/packages", .libPaths())))
threads <- 8L
Sys.setenv(OMP_NUM_THREADS=threads, OPENBLAS_NUM_THREADS=threads, MKL_NUM_THREADS=threads,
           RCPP_PARALLEL_NUM_THREADS=threads)
if (requireNamespace("RhpcBLASctl", quietly=TRUE)) {
  RhpcBLASctl::blas_set_num_threads(threads); RhpcBLASctl::omp_set_num_threads(threads)
}
if (!file.exists(checkpoint)) stop("Missing merge checkpoint", call.=FALSE)
if (file.exists(outdir) || dir.exists(outdir)) stop("Refusing existing output directory: ", outdir, call.=FALSE)
dir.create(outdir, recursive=FALSE)
newfile <- function(name) { p <- file.path(outdir,name); if(file.exists(p)) stop("Refusing overwrite: ",p); p }
wtsv <- function(x,name) write.table(x,newfile(name),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
qv <- function(x,p) as.numeric(quantile(x,p,na.rm=TRUE,names=FALSE))

message("Reading immutable checkpoint")
obj <- qs2::qs_read(checkpoint,nthreads=threads)
md <- obj[[]]
md$cell_id <- rownames(md); md$orig.ident <- as.character(md$orig.ident); md$dataset <- as.character(md$dataset)
md$log_count <- log1p(md$nCount_RNA); md$log_feature <- log1p(md$nFeature_RNA)

# Recompute MT/HB from counts layers and determine feature availability per group.
hb_list <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
mt_n <- setNames(numeric(nrow(md)),md$cell_id); hb_n <- setNames(numeric(nrow(md)),md$cell_id)
orig_mt <- setNames(vector("list",length(unique(md$orig.ident))),sort(unique(md$orig.ident)))
for(ly in grep("^counts",SeuratObject::Layers(obj[["RNA"]]),value=TRUE)) {
  mat <- SeuratObject::LayerData(obj[["RNA"]],layer=ly,fast=FALSE)
  ci <- match(colnames(mat),md$cell_id); if(anyNA(ci)) stop("Counts layer cell mapping failed")
  mt <- grep("^MT-",rownames(mat),value=TRUE); hb <- intersect(hb_list,rownames(mat))
  if(length(mt)) mt_n[colnames(mat)] <- mt_n[colnames(mat)] + Matrix::colSums(mat[mt,,drop=FALSE])
  if(length(hb)) hb_n[colnames(mat)] <- hb_n[colnames(mat)] + Matrix::colSums(mat[hb,,drop=FALSE])
  for(g in unique(md$orig.ident[ci])) orig_mt[[g]] <- union(orig_mt[[g]],mt)
}
md$percent.mt <- ifelse(md$nCount_RNA>0,100*mt_n[md$cell_id]/md$nCount_RNA,NA_real_)
md$percent.hb <- ifelse(md$nCount_RNA>0,100*hb_n[md$cell_id]/md$nCount_RNA,NA_real_)

groups <- sort(unique(md$orig.ident))
thr <- do.call(rbind,lapply(groups,function(g){
  x <- md[md$orig.ident==g,]; ds <- unique(x$dataset); if(length(ds)!=1) stop("orig.ident maps to multiple datasets: ",g)
  cm <- median(x$log_count,na.rm=TRUE); fm <- median(x$log_feature,na.rm=TRUE); mm <- median(x$percent.mt,na.rm=TRUE)
  c1 <- mad(x$log_count,center=cm,constant=1,na.rm=TRUE); f1 <- mad(x$log_feature,center=fm,constant=1,na.rm=TRUE)
  cnew <- mad(x$log_count,center=cm,constant=1.4826,na.rm=TRUE); fnew <- mad(x$log_feature,center=fm,constant=1.4826,na.rm=TRUE)
  mt1 <- mad(x$percent.mt,center=mm,constant=1,na.rm=TRUE)
  mtnew <- mad(x$percent.mt,center=mm,constant=1.4826,na.rm=TRUE)
  available <- length(orig_mt[[g]])>0 && ds!="InhouseData"
  data.frame(orig.ident=g,dataset=ds,n_cells=nrow(x),mt_qc_available=available,n_MT_features=length(orig_mt[[g]]),
    count_median_log=cm,count_MAD_old=c1,count_MAD_standard=cnew,
    A_count_lower=cm-3*c1,A_count_upper=cm+3*c1,B_count_lower=cm-3*c1,B_count_upper=cm+5*c1,
    D_count_lower=cm-3*cnew,D_count_upper=cm+5*cnew,E_count_lower=cm-3*cnew,E_count_upper=cm+3*cnew,
    feature_median_log=fm,feature_MAD_old=f1,feature_MAD_standard=fnew,
    A_feature_lower=fm-3*f1,A_feature_upper=fm+3*f1,B_feature_lower=fm-3*f1,B_feature_upper=fm+5*f1,
    D_feature_lower=fm-3*fnew,D_feature_upper=fm+5*fnew,E_feature_lower=fm-3*fnew,E_feature_upper=fm+3*fnew,
    mt_median=mm,mt_MAD_old=mt1,mt_MAD_standard=mtnew,
    mt_cutoff_old=if(available) min(mm+3*mt1,25) else NA_real_,
    mt_cutoff_standard=if(available) min(mm+3*mtnew,25) else NA_real_,hb_cutoff=5)
}))
wtsv(thr,"04_standardMAD_thresholds.tsv")

i <- match(md$orig.ident,thr$orig.ident)
hb_fail <- is.na(md$percent.hb)|md$percent.hb>5
mt_old <- ifelse(thr$mt_qc_available[i],is.na(md$percent.mt)|md$percent.mt>thr$mt_cutoff_old[i],FALSE)
mt_new <- ifelse(thr$mt_qc_available[i],is.na(md$percent.mt)|md$percent.mt>thr$mt_cutoff_standard[i],FALSE)
rule <- function(cl,cu,fl,fu,mtf) !(is.na(md$log_count)|md$log_count<cl[i]|md$log_count>cu[i]|
  is.na(md$log_feature)|md$log_feature<fl[i]|md$log_feature>fu[i]|mtf|hb_fail)
md$A <- rule(thr$A_count_lower,thr$A_count_upper,thr$A_feature_lower,thr$A_feature_upper,mt_old)
md$B <- rule(thr$B_count_lower,thr$B_count_upper,thr$B_feature_lower,thr$B_feature_upper,mt_old)
md$D <- rule(thr$D_count_lower,thr$D_count_upper,thr$D_feature_lower,thr$D_feature_upper,mt_new)
md$E <- rule(thr$E_count_lower,thr$E_count_upper,thr$E_feature_lower,thr$E_feature_upper,mt_new)

summary_row <- function(x,level,name){
  n<-nrow(x); vals<-vapply(c("A","B","D","E"),function(z)sum(x[[z]]),numeric(1))
  data.frame(level=level,group=name,initial_cells=n,
    Rule_A_old_retained=vals[1],Rule_A_old_retention_pct=100*vals[1]/n,
    Rule_B_old_retained=vals[2],Rule_B_old_retention_pct=100*vals[2]/n,
    Rule_D_new_retained=vals[3],Rule_D_new_retention_pct=100*vals[3]/n,
    Rule_E_reference_retained=vals[4],Rule_E_reference_retention_pct=100*vals[4]/n,
    D_minus_A=vals[3]-vals[1],D_minus_B=vals[3]-vals[2],stringsAsFactors=FALSE)
}
global <- summary_row(md,"global","ALL")
byds <- do.call(rbind,lapply(sort(unique(md$dataset)),function(g)summary_row(md[md$dataset==g,],"dataset",g)))
byorig <- do.call(rbind,lapply(groups,function(g)summary_row(md[md$orig.ident==g,],"orig.ident",g)))
byorig$dataset <- thr$dataset[match(byorig$group,thr$orig.ident)]
wtsv(global,"01_standardMAD_rule_comparison_global.tsv")
wtsv(byds,"02_standardMAD_rule_comparison_by_dataset.tsv")
wtsv(byorig,"03_standardMAD_rule_comparison_by_orig_ident.tsv")

# Tail flags for Rule D and integrity/pathology checks.
md$D_count_low <- md$log_count<thr$D_count_lower[i]; md$D_count_high <- md$log_count>thr$D_count_upper[i]
md$D_feature_low <- md$log_feature<thr$D_feature_lower[i]; md$D_feature_high <- md$log_feature>thr$D_feature_upper[i]
tail_global <- c(count_low=sum(md$D_count_low),count_high=sum(md$D_count_high),
                 feature_low=sum(md$D_feature_low),feature_high=sum(md$D_feature_high))
essential_cols <- c("count_MAD_standard","D_count_lower","D_count_upper",
                    "feature_MAD_standard","D_feature_lower","D_feature_upper")
bad_numeric <- !is.finite(as.matrix(thr[,essential_cols,drop=FALSE]))
bad_mt_numeric <- thr$mt_qc_available & !is.finite(thr$mt_cutoff_standard)
bad_mad <- thr$count_MAD_standard==0 | thr$feature_MAD_standard==0
near100 <- byorig$Rule_D_new_retention_pct>=99.5
dataset_near100 <- byds$Rule_D_new_retention_pct>=99.5
serious <- any(bad_numeric) || any(bad_mt_numeric) || any(bad_mad) || any(dataset_near100)
check <- data.frame(
  check=c("NA_or_Inf_numeric_thresholds","zero_standard_MAD","dataset_retention_ge_99.5pct","orig_ident_retention_ge_99.5pct","Inhouse_MT_unavailable"),
  count=c(sum(bad_numeric)+sum(bad_mt_numeric),sum(bad_mad),sum(dataset_near100),sum(near100),sum(!thr$mt_qc_available & thr$dataset=="InhouseData")),
  blocking=c(any(bad_numeric)||any(bad_mt_numeric),any(bad_mad),any(dataset_near100),FALSE,FALSE),stringsAsFactors=FALSE)
wtsv(check,"05_phase2_pathology_checks.tsv")

topchange <- head(byorig[order(-byorig$D_minus_A),c("dataset","group","initial_cells","Rule_A_old_retention_pct","Rule_B_old_retention_pct","Rule_D_new_retention_pct","D_minus_A","D_minus_B")],10)
report <- c("# Standard MAD sensitivity summary","",
  paste0("- Immutable checkpoint: `",checkpoint,"`"),
  "- Standard MAD is explicitly `constant = 1.4826`.",
  paste0("- Rule A_old retained: ",global$Rule_A_old_retained," (",round(global$Rule_A_old_retention_pct,3),"%)."),
  paste0("- Rule B_old retained: ",global$Rule_B_old_retained," (",round(global$Rule_B_old_retention_pct,3),"%)."),
  paste0("- Rule D_new retained: ",global$Rule_D_new_retained," (",round(global$Rule_D_new_retention_pct,3),"%)."),
  paste0("- Rule E_reference retained: ",global$Rule_E_reference_retained," (",round(global$Rule_E_reference_retention_pct,3),"%)."),
  paste0("- D minus A: ",global$D_minus_A," cells."),paste0("- D minus B: ",global$D_minus_B," cells."),
  paste0("- Rule D tail flags (overlap possible): count low=",tail_global[1],", count high=",tail_global[2],", feature low=",tail_global[3],", feature high=",tail_global[4],"."),
  "- InhouseData: mt_qc_available=FALSE; MT is not used for removal and zero percent.mt is not interpreted biologically.",
  paste0("- Phase 2 blocking pathology: ",serious,"."),"","## Largest changes versus old Rule A","",
  paste0("- ",topchange$dataset," / ",topchange$group,": +",topchange$D_minus_A," cells (D vs A)"))
writeLines(report,newfile("QC_standardMAD_summary.md"),useBytes=TRUE)
writeLines(c(capture.output(sessionInfo()),paste("blocking_pathology",serious,sep="\t")),newfile("sessionInfo.txt"))
message("AUDIT_COMPLETE blocking_pathology=",serious," outdir=",outdir)
