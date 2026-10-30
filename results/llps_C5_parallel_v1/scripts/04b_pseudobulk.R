# In-process raw-count aggregation; same formula and eligibility as Fig2 C0.
PB<-file.path(M1,'pseudobulk');dir.create(PB,showWarnings=FALSE)
write_tsv<-function(x,name)wt(x,file.path(PB,name))
meta_review<-fread(file.path(BASE,'results/fig2_like_C0/phase1B_review/sample_metadata_reviewed.tsv'))
counts<-SeuratObject::LayerData(o,assay='RNA',layer='counts')
groups<-paste(d$sample,ifelse(d$cluster=='C5','C5','Others'),sep='::')
lev<-sort(unique(groups));agg<-sparseMatrix(i=seq_along(groups),j=match(groups,lev),x=1,dims=c(nrow(d),length(lev)))
pc<-counts%*%agg;colnames(pc)<-lev
pm<-d[,.(n_cells=.N),by=.(sample,group=ifelse(cluster=='C5','C5','Others'))]
pm[,library:=paste(sample,group,sep='::')]
pm<-merge(pm,meta_review[,.(sample,dataset,tissue_origin,patient_uid,primary_replication_allowed)],by='sample')
pm<-pm[match(lev,library)];pb<-list(counts=pc,metadata=pm)
# Extract and reuse the actual C0 eligibility and fit_de code; no source-file mutation.
lines<-readLines(file.path(BASE,'results/fig2_like_C0/scripts/02_pseudobulk_DEG.R'))
src<-paste(lines,collapse='\n');a<-regexpr('eligible <- pm',src,fixed=TRUE)[1];b<-regexpr('writeLines(c(',src,fixed=TRUE)[1]
elig_code<-substr(src,a,b-1);elig_code<-gsub('C0','C5',elig_code,fixed=TRUE)
eval(parse(text=elig_code))
a<-regexpr('fit_de <- function',src,fixed=TRUE)[1];b<-regexpr("primary <- fit_de",src,fixed=TRUE)[1]
fit_code<-gsub('C0','C5',substr(src,a,b-1),fixed=TRUE)
suppressPackageStartupMessages(library(edgeR));eval(parse(text=fit_code))
audits<-list()
for(label in c('primary','min50','tumor_only')){
 ids<-if(label=='primary')sel else if(label=='min50')eligible[included&min_cells>=50,sample] else eligible[included&tissue_origin=='Tumor',sample]
 if(length(ids)<4){audits[[label]]<-data.table(model=label,n_samples=length(ids),status='NOT_FIT',reason='Fewer than four independent samples with required cells per arm');next}
 fit<-fit_de(ids,label);audits[[label]]<-cbind(fit$audit,status='FIT',reason='')
 write_tsv(fit$result,paste0('C5_vs_Rest_',label,'.tsv'))
}
write_tsv(rbindlist(audits,fill=TRUE),'model_summary.tsv')
