options(mc.cores=1L)
.libPaths(c('results/C5_spatial_colocalization_esda_v1/.Rlib',.libPaths()))
suppressPackageStartupMessages({library(spdep);library(data.table);library(Matrix)})
setDTthreads(1);o='results/C5_spatial_colocalization_esda_v1';partners=c('Endothelial','Epithelial','Myeloid','T_NK');args=commandArgs(TRUE)
meta=fread(file.path(o,'00_input_audit/spatial_tissue_patient_audit.tsv'));pil=fread(file.path(o,'00_input_audit/pilot_prespecified.tsv'))$section_id
ss=if(length(args)&&args[1]=='pilot')pil else sort(meta$section_id)
for(sid in ss){dest=file.path(o,'04_lee_L_sensitivity',paste0(sid,'.tsv'));if(file.exists(dest))next
 d=fread(file.path(o,'figures',paste0(sid,'_plot_data.tsv.gz')));d=d[analysis_included==TRUE];ed=fread(file.path(o,'01_spatial_weights',paste0(sid,'_edges.tsv.gz')));n=nrow(d)
 A=sparseMatrix(i=c(ed$source,ed$target)+1L,j=c(ed$target,ed$source)+1L,x=1,dims=c(n,n));lw=mat2listw(A,style='W',zero.policy=FALSE);stopifnot(length(lw$neighbours)==n,all(abs(vapply(lw$weights,sum,numeric(1))-1)<1e-12));z=list()
 for(p in partners){for(method in c('raw','adjusted')){x=if(method=='raw')d$five_UCell else d[[paste0(p,'__adjusted_five')]];y=if(method=='raw')d[[paste0('major_mean__',p)]] else d[[paste0(p,'__adjusted_abundance')]];set.seed(123);a=lee.mc(x,y,lw,nsim=999,alternative='two.sided',zero.policy=FALSE);z[[length(z)+1]]=data.table(section_id=sid,partner=p,method=method,n=n,L=unname(a$statistic),p=a$p.value,permutations=999,seed=123)}}
 fwrite(rbindlist(z),dest,sep='\t');message('Lee DONE ',sid)
}
f=rbindlist(lapply(list.files(file.path(o,'04_lee_L_sensitivity'),pattern='\\.tsv$',full.names=TRUE),fread));f[,FDR:=p.adjust(p,method='BH'),by=method];fwrite(f,file.path(o,'tables/lee_results.tsv'),sep='\t')
