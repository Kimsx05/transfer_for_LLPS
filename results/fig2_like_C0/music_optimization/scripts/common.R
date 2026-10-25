source('results/fig2_like_C0/phase1B_review/scripts/common.R')
OPT<-file.path(OUT,'music_optimization');P2<-file.path(OUT,'phase2_external');RV<-file.path(P2,'reference_validation')
.libPaths(c(file.path(P2,'Rlib'),.libPaths()));suppressPackageStartupMessages(library(nnls))
source(file.path(P2,'resources/MuSiC_analysis.R'))
DSCOL<-c(CNP0000460='#0072B2',GSE135337='#E69F00',GSE222315='#009E73',HRA003620='#CC79A7',InhouseData='#D55E00')
WT<-function(x,name)fwrite(x,file.path(CORE,name),sep='\t',na='NA')
cor_safe<-function(x,y,method='spearman'){ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok];if(length(x)<3||sd(x)==0||sd(y)==0)NA_real_ else cor(x,y,method=method)}
metric<-function(t,p){ok<-is.finite(t)&is.finite(p);t<-t[ok];p<-p[ok];e<-p-t;c(n=length(t),rho=cor_safe(t,p),MAE=mean(abs(e)),RMSE=sqrt(mean(e^2)),bias=mean(e))}
b<-readRDS(file.path(RV,'patient_celltype_pseudobulk.rds'));X<-b$counts;M<-as.data.table(b$metadata);TYPES<-b$types
old<-fread(file.path(RV,'natural_mixture_predictions.tsv'));TEST<-old[scenario=='Full_RNA'&primary_patient_specimen==TRUE&n_stromal>=50];TESTIDS<-TEST$patient_key
# Patient aggregate statistics only; no cells or samples from reserved patients enter new references.
basis<-function(x,m,mincells=5,balanced=FALSE){
 k<-paste(m$patient_key,m$cell_type,sep='::');lv<-unique(k)
 A<-sparseMatrix(i=seq_along(k),j=match(k,lv),x=1,dims=c(length(k),length(lv)))
 z<-as.matrix(x%*%A);md<-m[,.(n_cells=sum(n_cells),dataset=dataset[1]),by=.(key=paste(patient_key,cell_type,sep='::'),patient_key,cell_type)];md<-md[match(lv,key)]
 ok<-md$n_cells>=mincells&colSums(z)>0;md<-md[ok];z<-z[,ok,drop=FALSE]
 th<-sweep(z,2,colSums(z),'/');types<-unique(m$cell_type)
 meanmat<-varmat<-matrix(NA_real_,nrow(z),length(types),dimnames=list(rownames(z),types));ss<-setNames(rep(NA_real_,length(types)),types)
 for(ct in types){j<-which(md$cell_type==ct);if(!length(j))next;w<-rep(1/length(j),length(j))
  if(balanced){ds<-md$dataset[j];w<-1/as.numeric(table(ds)[ds]);w<-w/sum(w)}
  meanmat[,ct]<-as.vector(th[,j,drop=FALSE]%*%w);ss[ct]<-sum(colSums(z[,j,drop=FALSE])/md$n_cells[j]*w)
  if(length(j)>1)varmat[,ct]<-rowSums(sweep((th[,j,drop=FALSE]-meanmat[,ct])^2,2,w,'*'))/(1-sum(w^2))
 }
 cov<-md[,.(n_patients=.N,n_cells=sum(n_cells),n_datasets=uniqueN(dataset)),by=cell_type]
 list(D=sweep(meanmat,2,ss,'*'),M=meanmat,Sigma=varmat,S=ss,coverage=cov,theta=th,metadata=md)
}
fit_music<-function(y,bas,genes){
 g<-intersect(genes,names(y));y<-y[g];nz<-is.finite(y)&y>0;g<-g[nz];y<-y[nz]/sum(y[nz])
 fit<-music.iter(y,bas$D[g,,drop=FALSE],bas$S,bas$Sigma[g,,drop=FALSE],iter.max=1000,nu=.0001,eps=.01)
 fit$p<-setNames(fit$p.weight,colnames(bas$D));fit$genes<-g;fit$Y<-y;fit
}
