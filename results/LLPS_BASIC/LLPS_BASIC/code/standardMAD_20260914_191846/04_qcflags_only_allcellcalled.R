# Formal Rule D QC rerun from immutable merge checkpoint.
stamp <- "20260914_191846"
root <- "/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017"
checkpoint <- file.path(root,"checkpoints","01_merged_full_20260914_153425.qs2")
outdir <- file.path(root,"02_qc","QCflags_only_allCellCalled_20260914_195247")
.libPaths(unique(c("/home/jinshengxi/Rlibrary/LLPS_pipeline_v1_20260914_150459","/home/jinshengxi/Rlibrary/packages",.libPaths())))
threads<-8L; set.seed(123L)
Sys.setenv(OMP_NUM_THREADS=threads,OPENBLAS_NUM_THREADS=threads,MKL_NUM_THREADS=threads,RCPP_PARALLEL_NUM_THREADS=threads)
if(requireNamespace("RhpcBLASctl",quietly=TRUE)){RhpcBLASctl::blas_set_num_threads(threads);RhpcBLASctl::omp_set_num_threads(threads)}
if(file.exists(outdir)||dir.exists(outdir))stop("Refusing existing QC output directory: ",outdir)
dir.create(outdir,recursive=FALSE); dir.create(file.path(outdir,"plots"),recursive=FALSE)
newp<-function(dir,name){p<-file.path(dir,name);if(file.exists(p)||dir.exists(p))stop("Refusing overwrite: ",p);p}
wtsv<-function(x,name)write.table(x,newp(outdir,name),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
qv<-function(x,p)as.numeric(quantile(x,p,na.rm=TRUE,names=FALSE))

obj<-qs2::qs_read(checkpoint,nthreads=threads); md<-obj[[]]; md$cell_id<-rownames(md)
md$orig.ident<-as.character(md$orig.ident);md$dataset<-as.character(md$dataset)
hb_list<-c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
mt_n<-setNames(numeric(nrow(md)),md$cell_id);hb_n<-setNames(numeric(nrow(md)),md$cell_id)
orig_mt<-setNames(vector("list",length(unique(md$orig.ident))),sort(unique(md$orig.ident)))
for(ly in grep("^counts",SeuratObject::Layers(obj[["RNA"]]),value=TRUE)){
  mat<-SeuratObject::LayerData(obj[["RNA"]],layer=ly,fast=FALSE);ci<-match(colnames(mat),md$cell_id)
  if(anyNA(ci))stop("Counts layer cell mapping failed")
  mt<-grep("^MT-",rownames(mat),value=TRUE);hb<-intersect(hb_list,rownames(mat))
  if(length(mt))mt_n[colnames(mat)]<-mt_n[colnames(mat)]+Matrix::colSums(mat[mt,,drop=FALSE])
  if(length(hb))hb_n[colnames(mat)]<-hb_n[colnames(mat)]+Matrix::colSums(mat[hb,,drop=FALSE])
  for(g in unique(md$orig.ident[ci]))orig_mt[[g]]<-union(orig_mt[[g]],mt)
}
md$percent.mt<-ifelse(md$nCount_RNA>0,100*mt_n[md$cell_id]/md$nCount_RNA,NA_real_)
md$percent.hb<-ifelse(md$nCount_RNA>0,100*hb_n[md$cell_id]/md$nCount_RNA,NA_real_)
md$lc<-log1p(md$nCount_RNA);md$lf<-log1p(md$nFeature_RNA)
groups<-sort(unique(md$orig.ident))
thr<-do.call(rbind,lapply(groups,function(g){x<-md[md$orig.ident==g,];ds<-unique(x$dataset)
  cm<-median(x$lc,na.rm=TRUE);fm<-median(x$lf,na.rm=TRUE);mm<-median(x$percent.mt,na.rm=TRUE)
  cmad<-mad(x$lc,center=cm,constant=1.4826,na.rm=TRUE);fmad<-mad(x$lf,center=fm,constant=1.4826,na.rm=TRUE)
  mmad<-mad(x$percent.mt,center=mm,constant=1.4826,na.rm=TRUE);avail<-length(orig_mt[[g]])>0&&ds!="InhouseData"
  data.frame(orig.ident=g,dataset=ds,initial_cells=nrow(x),MAD_constant=1.4826,
   log_nCount_median=cm,log_nCount_MAD=cmad,log_nCount_lower=cm-3*cmad,log_nCount_upper=cm+5*cmad,
   nCount_lower=pmax(0,expm1(cm-3*cmad)),nCount_upper=expm1(cm+5*cmad),
   log_nFeature_median=fm,log_nFeature_MAD=fmad,log_nFeature_lower=fm-3*fmad,log_nFeature_upper=fm+5*fmad,
   nFeature_lower=pmax(0,expm1(fm-3*fmad)),nFeature_upper=expm1(fm+5*fmad),
   mt_qc_available=avail,n_MT_features=length(orig_mt[[g]]),percent_mt_median=mm,percent_mt_MAD=mmad,
   percent_mt_MAD_upper=if(avail)mm+3*mmad else NA_real_,mt_final_cutoff=if(avail)min(mm+3*mmad,25)else NA_real_,hb_cutoff=5)
}))
if(any(thr$log_nCount_MAD==0|thr$log_nFeature_MAD==0))stop("Zero MAD encountered")
essential<-c("log_nCount_lower","log_nCount_upper","log_nFeature_lower","log_nFeature_upper")
if(any(!is.finite(as.matrix(thr[,essential]))))stop("Non-finite QC threshold")
i<-match(md$orig.ident,thr$orig.ident)
md$fail_nCount_low<-is.na(md$lc)|md$lc<thr$log_nCount_lower[i]
md$fail_nCount_high<-!is.na(md$lc)&md$lc>thr$log_nCount_upper[i]
md$fail_nFeature_low<-is.na(md$lf)|md$lf<thr$log_nFeature_lower[i]
md$fail_nFeature_high<-!is.na(md$lf)&md$lf>thr$log_nFeature_upper[i]
md$mt_qc_available<-thr$mt_qc_available[i]
md$fail_mt<-ifelse(md$mt_qc_available,is.na(md$percent.mt)|md$percent.mt>thr$mt_final_cutoff[i],FALSE)
md$fail_hb<-is.na(md$percent.hb)|md$percent.hb>5
flags<-c("fail_nCount_low","fail_nCount_high","fail_nFeature_low","fail_nFeature_high","fail_mt","fail_hb")
md$QC_pass<-rowSums(md[,flags,drop=FALSE])==0

for(nm in c("percent.mt","percent.hb","fail_nCount_low","fail_nCount_high","fail_nFeature_low","fail_nFeature_high","mt_qc_available","fail_mt","fail_hb","QC_pass"))obj[[nm]]<-md[[nm]]
wtsv(thr,"qc_thresholds_standardMAD_by_orig_ident.tsv")
summ<-function(x,level,group)data.frame(level=level,group=group,dataset=if(level=="orig.ident")unique(x$dataset)else if(level=="dataset")group else "ALL",
 initial_cells=nrow(x),fail_nCount_low=sum(x$fail_nCount_low),fail_nCount_high=sum(x$fail_nCount_high),
 fail_nFeature_low=sum(x$fail_nFeature_low),fail_nFeature_high=sum(x$fail_nFeature_high),fail_mt=sum(x$fail_mt),fail_hb=sum(x$fail_hb),
 failure_overlap=sum(rowSums(x[,flags])>1),QC_pass=sum(x$QC_pass),removed=sum(!x$QC_pass),retention_pct=100*mean(x$QC_pass),
 mt_qc_available=all(x$mt_qc_available))
so<-do.call(rbind,lapply(groups,function(g)summ(md[md$orig.ident==g,],"orig.ident",g)))
sd<-do.call(rbind,lapply(sort(unique(md$dataset)),function(g)summ(md[md$dataset==g,],"dataset",g)))
wtsv(so,"qc_summary_standardMAD_by_orig_ident.tsv");wtsv(sd,"qc_summary_standardMAD_by_dataset.tsv")
combo<-as.data.frame(table(apply(md[,flags],1,function(z){n<-flags[z];if(!length(n))"PASS"else paste(n,collapse="+")})),stringsAsFactors=FALSE)
names(combo)<-c("failure_combination","n_cells");combo$percentage<-100*combo$n_cells/nrow(md);wtsv(combo,"qc_failure_overlap.tsv")

# Compact audit plots only.
long<-rbind(data.frame(dataset=md$dataset,metric="log1p(nCount_RNA)",value=md$lc),
 data.frame(dataset=md$dataset,metric="log1p(nFeature_RNA)",value=md$lf),
 data.frame(dataset=md$dataset,metric="percent.mt",value=md$percent.mt),
 data.frame(dataset=md$dataset,metric="percent.hb",value=md$percent.hb))
p<-ggplot2::ggplot(long,ggplot2::aes(dataset,value,fill=dataset))+ggplot2::geom_violin(scale="width",trim=TRUE,show.legend=FALSE)+
 ggplot2::facet_wrap(~metric,scales="free_y",ncol=1)+ggplot2::theme_bw(base_size=9)+ggplot2::theme(panel.grid=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=30,hjust=1))
for(ext in c("pdf","png")){pp<-newp(file.path(outdir,"plots"),paste0("qc_standardMAD_distributions.",ext));ggplot2::ggsave(pp,p,width=11,height=14,dpi=if(ext=="png")180 else 300,limitsize=FALSE)}
bar<-rbind(data.frame(dataset=sd$group,status="retained",cells=sd$QC_pass),data.frame(dataset=sd$group,status="removed",cells=sd$removed))
pb<-ggplot2::ggplot(bar,ggplot2::aes(dataset,cells,fill=status))+ggplot2::geom_col()+ggplot2::theme_bw()+ggplot2::theme(panel.grid=ggplot2::element_blank())
for(ext in c("pdf","png")){pp<-newp(file.path(outdir,"plots"),paste0("qc_standardMAD_retained_removed.",ext));ggplot2::ggsave(pp,pb,width=9,height=6,dpi=if(ext=="png")180 else 300)}

if(ncol(obj)!=410466L)stop("All-cell-called object cell count changed unexpectedly")
cp<-file.path(root,"checkpoints",paste0("01b_QCflags_only_noFiltering_",format(Sys.time(),"%Y%m%d_%H%M%S"),".qs2"))
if(file.exists(cp))stop("Refusing checkpoint overwrite")
qs2::qs_save(obj,cp,compress_level=4,shuffle=TRUE,nthreads=threads)
if(!file.exists(cp)||file.info(cp)$size<=0)stop("QC checkpoint write failed")
writeLines(cp,newp(outdir,"qc_standardMAD_checkpoint_path.txt"))
writeLines(c("Rule D QC flags only; no filtering","all cell-called cells retained = 410466","MAD constant = 1.4826","nCount/nFeature = lower 3 MAD, upper 5 MAD","MT = upper 3 standard MAD capped at 25%","InhouseData mt_qc_available = FALSE; fail_mt = FALSE","HB >5%",paste("eventual_QC_pass",sum(md$QC_pass))),newp(outdir,"qc_standardMAD_parameters.txt"))
message("QCFLAGS_ONLY_COMPLETE cells=",ncol(obj)," eventual_QC_pass=",sum(md$QC_pass)," checkpoint=",cp)
