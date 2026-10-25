source('results/fig2_like_C0/music_optimization/scripts/common.R');CORE<-file.path(OPT,'06_signature_fallback')
stopifnot(readLines(file.path(CORE,'SIGNATURE_STATUS.txt'))=='PASS_EXPLORATORY_EXPRESSION_SCORE_ONLY')
sig<-readLines(file.path(CORE,'C0_identity_signature_genes.txt'));uni<-readLines(file.path(CORE,'signature_gene_universe.txt'))
score<-function(v){N<-length(v);k<-length(sig);r<-rank(v,ties.method='average');(mean(r[match(sig,names(v))])-(k+1)/2)/(N-k)-.5}
audit<-list()
for(cohort in c('TCGA','GSE31684','GSE48075','GSE13507')){
 cat('Scoring',cohort,'\n');flush.console()
 clin<-fread(file.path(P2,'tables',paste0(cohort,'_clinical_harmonized.tsv')));clin<-clin[primary_tumor_eligible==TRUE]
 stopifnot(!anyDuplicated(clin$patient_unit),!anyDuplicated(clin$sample))
 if(cohort=='TCGA'){
  e<-fread(file.path(P2,'resources/TCGA_BLCA_HiSeqV2.gz'));setnames(e,1,'symbol');e<-e[symbol%in%uni,c('symbol',clin$sample),with=FALSE]
 }else{
  inv<-fread(file.path(P2,'tables/cohort_inventory.tsv'));nf<-inv$n_features[match(cohort,inv$cohort)]
  e<-fread(file.path(P2,'resources',paste0(cohort,'_series_matrix.txt.gz')),skip='ID_REF',nrows=nf);setnames(e,1,'probe')
  platform<-unique(clin$platform);map<-fread(file.path(P2,'tables',paste0(platform,'_probe_symbol.tsv')))
  e<-merge(map[symbol%in%uni],e,by='probe');e<-e[,c('symbol',clin$sample),with=FALSE]
 }
 # Probe-to-gene medians on supplied expression scale, fixed without outcome information.
 if(anyDuplicated(e$symbol))e<-e[,lapply(.SD,median,na.rm=TRUE),by=symbol]
 stopifnot(all(uni%in%e$symbol));e<-e[match(uni,symbol)];mat<-as.matrix(e[,-1]);rownames(mat)<-e$symbol
 stopifnot(identical(rownames(mat),uni),!anyNA(mat),all(is.finite(mat)),all(sig%in%rownames(mat)))
 sc<-apply(mat,2,score);clin[,C0_like_signature_score:=unname(sc[match(sample,names(sc))])]
 stopifnot(all(is.finite(clin$C0_like_signature_score)))
 WT(clin,paste0(cohort,'_signature_clinical.tsv'))
 WT(as.data.table(mat[sig,,drop=FALSE],keep.rownames='gene'),paste0(cohort,'_signature_gene_expression.tsv'))
 audit[[cohort]]<-data.table(cohort=cohort,n_primary=nrow(clin),n_unique_patients=uniqueN(clin$patient_unit),signature_genes_present=length(sig),universe_genes=length(uni),missing_signature_genes=0,score_min=min(sc),score_max=max(sc),method='fixed-universe normalized mean rank; expression score, not abundance')
 rm(e,mat);gc()
}
WT(rbindlist(audit),'external_signature_coverage.tsv');cat('External scoring complete\n')
