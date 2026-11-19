options(mc.cores=1L)
source('cellchat_nichechat_fibro7/scripts/nichenet_common.R')
setDTthreads(1)
library(edgeR);library(limma)
O<-'results/C5_Endothelial_NicheNet_spatial_v1'
w<-function(x,p)fwrite(x,file.path(O,p),sep='\t',na='NA')
manifest<-fread(file.path(NN_OUT,'verified_resource_manifest.tsv'))
manifest[,actual_sha256:=vapply(path,sha256,character(1))];stopifnot(all(manifest$actual_sha256==manifest$sha256));w(manifest,'00_input_audit/resource_verification.tsv')
mp<-fread('results/llps_C5_distribution_audit_20261110_v1/v2/patient_sample_mapping_final.csv')
s<-fread('results/llps_C5_distribution_audit_20261110_v1/v2/sample_celltype_counts.csv')
s[,Fib_total:=rowSums(.SD),.SDcols=paste0('C',0:6)];s[,fraction:=C5/Fib_total];w(s,'00_input_audit/specimen_sampling.tsv')
b<-readRDS('results/fig2_like_C0/phase2_external/reference_validation/patient_celltype_pseudobulk.rds');m<-as.data.table(b$metadata)
m[,patient_id:=mp$patient_id[match(sample,mp$specimen_id)]];stopifnot(!anyNA(m$patient_id),identical(m$library,colnames(b$counts)))
chk<-merge(m,s[,c('specimen_id','Endothelial',paste0('C',0:6)),with=FALSE],by.x='sample',by.y='specimen_id')
chk<-chk[cell_type%in%c('Endothelial',paste0('C',0:6))];chk[,expected:=mapply(function(t,i)s[[t]][i],cell_type,match(sample,s$specimen_id))];stopifnot(all(chk$n_cells==chk$expected));w(chk[,.(sample,cell_type,n_cells,expected)],'00_input_audit/pseudobulk_cache_count_check.tsv')
p<-s[tissue_origin=='Tumor',.(C5=sum(C5),Fib_total=sum(Fib_total),Endothelial=sum(Endothelial),specimens=paste(specimen_id,collapse=';'),stage=paste(sort(unique(clinical_group)),collapse=';')),by=.(patient_id,dataset)]
p[,`:=`(fraction=C5/Fib_total,eligible=Fib_total>=50 & Endothelial>=20)]
w(p,'01_receiver_pseudobulk/patient_sampling.tsv');print(p[,.(n=.N,eligible=sum(eligible),fraction_sd=sd(fraction[eligible])),by=dataset])
e<-m[cell_type=='Endothelial' & tissue_origin=='Tumor'];ids<-sort(unique(e$patient_id));a<-sparseMatrix(i=seq_len(nrow(e)),j=match(e$patient_id,ids),x=1,dims=c(nrow(e),length(ids)))
x<-b$counts[,match(e$library,colnames(b$counts)),drop=FALSE]%*%a;colnames(x)<-ids
p<-p[match(ids,patient_id)];stopifnot(!anyNA(p$patient_id));p[,library_UMI:=as.numeric(colSums(x))]
fitfun<-function(sel,label,extra=FALSE){
 d<-copy(p[sel]);if(nrow(d)<6 || sd(d$fraction)==0)return(NULL)
 d[,z:=as.numeric(scale(fraction))];d[,dataset:=factor(dataset)]
 rhs<-c(if(nlevels(d$dataset)>1)'dataset','z',if(extra)c('log1p(Fib_total)','log1p(Endothelial)'))
 design<-model.matrix(as.formula(paste('~',paste(rhs,collapse='+'))),d)
 if(qr(design)$rank<ncol(design)||nrow(design)-ncol(design)<3)return(NULL)
 y<-DGEList(as.matrix(x[,sel,drop=FALSE]));keep<-filterByExpr(y,design=design);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);v<-voom(y,design,plot=FALSE);f<-eBayes(lmFit(v,design),robust=TRUE)
 j<-which(colnames(design)=='z');t<-as.data.table(topTable(f,coef=j,n=Inf,sort.by='none'),keep.rownames='gene');t[,SE:=f$stdev.unscaled[,j]*sqrt(f$s2.post)];t[,`:=`(n_patients=nrow(d),residual_df=nrow(d)-ncol(design),fraction_sd=sd(d$fraction),model=label)]
 w(t,paste0('01_receiver_pseudobulk/DE_',label,'.tsv'));w(d,paste0('01_receiver_pseudobulk/design_',label,'.tsv'));t
}
main<-fitfun(which(p$eligible),'primary');fits<-list();for(ds in unique(p$dataset)){fits[[ds]]<-fitfun(which(p$eligible&p$dataset==ds),ds)}
fitfun(which(p$eligible),'sampling_adjusted',TRUE);fitfun(which(p$Fib_total>=100&p$Endothelial>=30),'stricter_sampling')
for(ds in unique(p$dataset))fitfun(which(p$eligible&p$dataset!=ds),paste0('without_',ds))
lt<-readRDS(file.path(NN_OUT,NN_RESOURCES[['ligand_target']]))
genes<-character();if(!is.null(main)){
 genes<-main[logFC>0&adj.P.Val<.10,gene]
 for(ds in c('InhouseData','GSE222315')){q<-fits[[ds]];genes<-intersect(genes,if(is.null(q))character() else q[logFC>0,gene])}
 genes<-intersect(genes,rownames(lt))
}
loo<-list();if(length(genes)){for(i in which(p$eligible)){q<-fitfun(setdiff(which(p$eligible),i),paste0('LOPO_',i));if(!is.null(q))loo[[length(loo)+1]]<-q[gene%in%genes,.(gene,logFC,omitted=p$patient_id[i])]};ll<-rbindlist(loo);w(ll,'08_patient_dataset_sensitivity/response_LOPO.tsv');stable<-ll[,.(positive_fraction=mean(logFC>0)),by=gene];genes<-intersect(genes,stable[positive_fraction>=.8,gene])}
counts<-p[eligible==TRUE,.N,by=dataset];gate<-nrow(p[eligible==TRUE])>=12&&sum(counts$N>=4)>=2&&length(genes)>=15
writeLines(if(gate)'PATIENT_ASSOCIATED' else 'CURATED_EXPLORATORY',file.path(O,'02_receiver_geneset/RESPONSE_MODE.txt'));w(data.table(gene=genes),'02_receiver_geneset/patient_response_genes.tsv')
w(data.table(n_eligible=sum(p$eligible),datasets_ge4=sum(counts$N>=4),positive_response_genes=length(genes),gate=gate),'02_receiver_geneset/gate.tsv')
# Small receiver pseudobulk only, not a duplicate single-cell matrix.
saveRDS(list(counts=x,metadata=p),file.path(O,'01_receiver_pseudobulk/Endothelial_patient_pseudobulk.rds'))
capture.output(sessionInfo(),file=file.path(O,'00_input_audit/sessionInfo.txt'))
cat('RESPONSE_COMPLETE',gate,length(genes),'\n')
