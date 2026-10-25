source('results/fig2_like_C0/music_optimization/scripts/common.R');CORE<-file.path(OPT,'06_signature_fallback')
stopifnot(readLines(file.path(CORE,'SIGNATURE_STATUS.txt'))=='PASS_EXPLORATORY_EXPRESSION_SCORE_ONLY')
uni<-readLines(file.path(CORE,'signature_gene_universe.txt'));sig<-readLines(file.path(CORE,'C0_identity_signature_genes.txt'))
suppressPackageStartupMessages(library(Seurat));obj<-qs::qread(OBJECT,nthreads=2);md<-as.data.table(obj[[]],keep.rownames='cell')
sm<-fread(file.path(OUT,'phase1B_review/sample_metadata_reviewed.tsv'));md[,patient_key:=sm$patient_key[match(orig.ident,sm$sample)]];md<-md[patient_key%in%TESTIDS]
x<-LayerData(obj,assay='RNA',layer='counts');u<-intersect(uni,rownames(x));x<-x[u,md$cell,drop=FALSE];stopifnot(all(sig%in%u));sigi<-match(sig,u);N<-length(uni);k<-length(sig)
score_sparse<-function(j){ix<-seq.int(x@p[j]+1,x@p[j+1]);if(x@p[j]==x@p[j+1])return(-.5+(N+1-k-1)/2/(N-k))
 vals<-x@x[ix];ind<-x@i[ix]+1L;positive<-vals>0;vals<-vals[positive];ind<-ind[positive];nz<-length(vals);nzero<-N-nz
 ranks<-rank(vals,ties.method='average')+nzero;target<-rep((nzero+1)/2,k);m<-match(sigi,ind);target[!is.na(m)]<-ranks[m[!is.na(m)]]
 (mean(target)-(k+1)/2)/(N-k)-.5}
md[,identity_signature_score:=vapply(seq_len(ncol(x)),score_sparse,numeric(1))]
z<-md[,.(n_cells=.N,median_identity_score=median(identity_signature_score),mean_identity_score=mean(identity_signature_score),median_existing_LLPS=median(LLPS77_AUC)),by=.(patient_key,dataset,cluster=subcluster)]
WT(z,'singlecell_signature_patient_cluster_summary.tsv')
p<-ggplot(z[n_cells>=5],aes(factor(cluster,levels=CLUSTERS),median_identity_score,fill=cluster))+geom_boxplot(width=.6,outlier.shape=NA)+geom_point(position=position_jitter(width=.12,seed=123),size=1.5)+scale_fill_manual(values=COLORS)+theme_pub()+theme(legend.position='none')+labs(title='C0 identity score at single-cell resolution',x=NULL,y='Patient median of single-cell signature score',caption='One point per reserved patient/cluster, at least 5 cells.\nOriginal LLPS AUCell values and cluster assignments were not modified.')
export_plot(p,'FigS_signature_singlecell_patient_medians',7,5)
cat('Single-cell identity-score check complete;',nrow(md),'cells summarized by patient.\n')
