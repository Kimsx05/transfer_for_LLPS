# Zero and missingness handling is shared across network, pathway and LR summaries.
edge_summary<-function(eligible,weight,inference_status) {
  use<-which(eligible)
  w<-weight[use]
  n<-length(use);missing<-sum(is.na(w));observed<-w[!is.na(w)]
  detected<-sum(observed>0)
  complete<-n>0L && missing==0L
  nz<-observed[observed>0]
  list(n_eligible=n,n_successful_eligible=length(observed),n_missing_eligible=missing,n_detected=detected,
       prevalence=if(complete)detected/n else NA_real_,
       prevalence_lower_bound=if(n)detected/n else NA_real_,prevalence_upper_bound=if(n)(detected+missing)/n else NA_real_,
       consensus_strength=if(complete)median(observed) else NA_real_,
       mean_weight=if(complete)mean(observed) else NA_real_,IQR=if(complete)IQR(observed) else NA_real_,
       median_nonzero_weight=if(length(nz))median(nz) else NA_real_,max_weight=if(length(observed))max(observed) else NA_real_,
       observed_median_weight=if(length(observed))median(observed) else NA_real_,aggregation_complete=complete,
       retained=complete && n>=3L && detected/n>=.30)
}
summarize_edges<-function(d,by=c('source','target')) {
  r<-d[,edge_summary(eligible,weight,inference_status),by=by]
  r[,display_strength:=ifelse(retained,consensus_strength,0)]
  r
}
matrix_from_statistics<-function(stat,column) {
  m<-blank_matrix()
  m[cbind(match(stat$source,GROUPS),match(stat$target,GROUPS))]<-stat[[column]]
  m
}
write_matrix<-function(m,path)write_tsv(data.table(source=rownames(m),as.data.table(m)),path)
