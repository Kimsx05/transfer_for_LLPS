source("/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/code/standardMAD_20260914_191846/00_config.R")
start_time<-Sys.time(); stamp<-format(start_time,"%Y%m%d_%H%M%S")
input_path<-"/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/checkpoints/02_harmony_clustered_res0.8_standardMAD_20260915_092711.qs2"
out_dir<-file.path(paths$harmony,paste0("UMAP_blank_compatFallback_",stamp))
if(file.exists(out_dir)||dir.exists(out_dir)) stop("Existing output directory")
dir.create(out_dir,recursive=FALSE,showWarnings=FALSE)
obj<-qs2::qs_read(input_path,nthreads=config$max_threads)
cluster_col<-grep("RNA_snn_res\\.0\\.8$",colnames(obj[[]]),value=TRUE)
if(length(cluster_col)!=1L) stop("Cluster column not unique")
d<-as.data.frame(Seurat::Embeddings(obj,"umap")); colnames(d)[1:2]<-c("UMAP_1","UMAP_2")
d$cluster<-factor(obj[[cluster_col,drop=TRUE]],levels=sort(unique(obj[[cluster_col,drop=TRUE]])))
d$dataset<-factor(obj$dataset); d$orig.ident<-factor(obj$orig.ident)
theme_blank<-ggplot2::theme_void()+ggplot2::theme(plot.background=ggplot2::element_rect(fill="white",colour=NA),panel.background=ggplot2::element_rect(fill="white",colour=NA),legend.background=ggplot2::element_rect(fill="white",colour=NA))
base_plot<-function(group,title) ggplot2::ggplot(d,ggplot2::aes(x=UMAP_1,y=UMAP_2,colour=.data[[group]]))+ggplot2::geom_point(shape=16,size=0.05,alpha=0.65)+ggplot2::labs(colour=group,title=title)+theme_blank
plots<-list(UMAP_cluster=base_plot("cluster","Clusters"),UMAP_dataset=base_plot("dataset","Dataset"),UMAP_orig_ident=base_plot("orig.ident","orig.ident"),UMAP_split_by_dataset=base_plot("cluster","Clusters by dataset")+ggplot2::facet_wrap(~dataset,ncol=3)+ggplot2::theme(legend.position="none"))
for(nm in names(plots)){w<-if(nm=="UMAP_split_by_dataset")18 else 10;h<-if(nm=="UMAP_split_by_dataset")16 else 8;ggplot2::ggsave(file.path(out_dir,paste0(nm,".pdf")),plots[[nm]],width=w,height=h,bg="white",limitsize=FALSE);ggplot2::ggsave(file.path(out_dir,paste0(nm,".png")),plots[[nm]],width=w,height=h,dpi=200,bg="white",limitsize=FALSE)}
write_lines_new(c(paste("input",input_path,sep="\t"),"requested_backend\tSCP::CellDimPlot","fallback_backend\tggplot2 from Seurat UMAP embeddings","theme\tblank white","reason\tSCP 0.5.6 incompatible with installed ggplot2: numeric shape 19 translation error and missing internal gg_par","failed_SCP_outputs_preserved\tTRUE",paste("elapsed_seconds",round(as.numeric(difftime(Sys.time(),start_time,units="secs")),3),sep="\t")),file.path(out_dir,"UMAP_fallback_audit.txt"))
message(format(Sys.time())," END output=",out_dir)
