# Reuse the actual C0 common functions, then redirect outputs only.
source('results/llps_C0_functional_signature_v2/scripts/common.R')
OUT<-file.path(BASE,'results/llps_C5_parallel_v1')
M1<-file.path(OUT,'04_C5_vs_Rest_DEG');M2<-file.path(OUT,'05_GO_KEGG');M3<-file.path(OUT,'06_Hallmark_GSEA')
load_frozen_C0<-load_frozen
load_frozen<-function(){z<-load_frozen_C0();z$metadata[,group:=ifelse(cluster=='C5','C5','Rest')];z}
future::plan(future::sequential)
