source('results/fig2_like_C0/phase1B_review/scripts/common.R')
suppressPackageStartupMessages({library(Seurat);library(qs)})
P2<-file.path(OUT,'phase2_external');RV<-file.path(P2,'reference_validation')
cfg<-readRDS(file.path(BASE,'cellchat_nichechat_fibro7/01_metadata/frozen_analysis_inputs.rds'))
m<-fread(file.path(CORE,'sample_metadata_reviewed.tsv'))
d<-as.data.table(cfg$metadata);d[,cell_type:=sub('^Fibro_','',cellchat_group)]
d<-merge(d,m[,.(sample,patient_key,primary_replication_allowed,tissue_origin)],by='sample',all.x=TRUE,sort=FALSE)
stopifnot(!anyNA(d$patient_key),!anyDuplicated(d$cell))
fwrite(d[,.(n_cells=.N),by=.(dataset,inclusion_status)],file.path(RV,'inherited_cell_inclusion.tsv'),sep='\t')
d<-d[!is.na(cell_type)&primary_replication_allowed==TRUE]
cat('Reading frozen full atlas counts:',cfg$major_object,'\n');flush.console()
o<-qs2::qs_read(cfg$major_object,nthreads=2);x<-LayerData(o,assay='RNA',layer='counts')
stopifnot(all(d$cell%in%colnames(x)))
# Aggregate directly from sparse raw counts; never persist another full single-cell object.
key<-paste(d$sample,d$cell_type,sep='::');lev<-sort(unique(key))
A<-sparseMatrix(i=match(d$cell,colnames(x)),j=match(key,lev),x=1,dims=c(ncol(x),length(lev)))
pb<-x%*%A;colnames(pb)<-lev
meta<-d[,.(n_cells=.N),by=.(sample,patient_key,dataset,tissue_origin,cell_type)]
meta[,library:=paste(sample,cell_type,sep='::')];meta<-meta[match(lev,library)]
stopifnot(identical(meta$library,colnames(pb)),all(Matrix::colSums(pb)>0),!anyDuplicated(rownames(pb)))
meta[,umi:=as.numeric(Matrix::colSums(pb))]
saveRDS(list(counts=pb,metadata=meta,types=c(CLUSTERS,'Epithelial','Endothelial','Myeloid','T_NK','B','Plasma','Mast')),
 file.path(RV,'patient_celltype_pseudobulk.rds'),compress='gzip')
fwrite(meta,file.path(RV,'sample_celltype_counts.tsv'),sep='\t')
cat('Complete:',nrow(pb),'genes',ncol(pb),'sample-type libraries;',sum(meta$n_cells),'cells\n')
