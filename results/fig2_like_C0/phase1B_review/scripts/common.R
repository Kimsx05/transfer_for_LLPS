source('results/fig2_like_C0/scripts/common.R')
PHASE1 <- CORE
CORE <- file.path(OUT,'phase1B_review')
set.seed(123)
WT <- function(x,name) fwrite(x,file.path(CORE,name),sep='\t',na='NA')
export_plot <- function(p,stem,w=7,h=5) {
 ggsave(file.path(CORE,paste0(stem,'.pdf')),p,width=w,height=h,device=cairo_pdf,bg='white')
 ggsave(file.path(CORE,paste0(stem,'.png')),p,width=w,height=h,dpi=350,device=ragg::agg_png,bg='white')
}
focus <- c('EPITHELIAL_MESENCHYMAL_TRANSITION','TGF_BETA_SIGNALING','INFLAMMATORY_RESPONSE','HYPOXIA','TNFA_SIGNALING_VIA_NFKB')
short <- c('EMT','TGF-beta','Inflammatory response','Hypoxia','TNF-alpha / NF-kB')
names(short) <- paste0('HALLMARK_',focus)
# Patient-cluster sandwich covariance, finite-sample correction; no cell-level inference.
cluster_lm <- function(data,formula,term='tumor',label='model') {
 mf<-model.frame(formula,data=data,na.action=na.pass); ok<-complete.cases(mf)
 z<-data[ok];fit<-lm(formula,z);X<-model.matrix(fit);n<-nrow(X);k<-ncol(X)
 G<-uniqueN(z$patient_key)
 if(qr(X)$rank<k || G<=k+1 || n<=k || !term%in%names(coef(fit)))
  return(data.table(model=label,n=n,n_patients=G,status='NOT_ESTIMABLE',beta=NA_real_,SE=NA_real_,lower=NA_real_,upper=NA_real_,p_value=NA_real_))
 meat<-crossprod(rowsum(X*residuals(fit),z$patient_key));bread<-solve(crossprod(X))
 V<-G/(G-1)*(n-1)/(n-k)*bread%*%meat%*%bread
 stopifnot(min(diag(V))>-1e-10)
 se<-sqrt(pmax(diag(V),0))[match(term,colnames(X))];b<-coef(fit)[term];df<-G-k
 data.table(model=label,n=n,n_patients=G,status='ESTIMABLE',beta=unname(b),SE=se,
 lower=unname(b)-qt(.975,df)*se,upper=unname(b)+qt(.975,df)*se,p_value=2*pt(-abs(b/se),df),df=df)
}
paired <- function(d,value) {
 z<-d[!is.na(patient_uid)&is.finite(get(value))]
 ids<-z[,.(n=.N,nt=uniqueN(tissue_origin)),by=patient_uid][n==2&nt==2,patient_uid]
 dcast(z[patient_uid%in%ids],patient_uid+dataset~tissue_origin,value.var=value)
}
