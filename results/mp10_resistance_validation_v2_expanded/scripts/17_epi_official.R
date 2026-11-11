.libPaths(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/refdir/Rlib',.libPaths()))
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
suppressPackageStartupMessages({library(Matrix);library(data.table);library(AUCell);library(UCell);library(Seurat)})
setDTthreads(1);B='/home/data/t070721/codex_workspace/Bladder Metabolism/results/mp10_resistance_validation_v2_expanded';args=commandArgs(TRUE);input=if(length(args))args[1] else file.path(B,'02_epi/pilot');set.seed(20261111)
readbin=function(p,type,size){con=file(p,'rb');on.exit(close(con));readBin(con,type,n=file.info(p)$size/size,size=size,endian='little')}
genes=readLines(file.path(input,'genes.txt'));cells=readLines(file.path(input,'cells.txt'));x=new('dgCMatrix',i=readbin(file.path(input,'indices.bin'),'integer',4),p=readbin(file.path(input,'indptr.bin'),'integer',4),x=readbin(file.path(input,'data.bin'),'numeric',8),Dim=as.integer(c(length(genes),length(cells))),Dimnames=list(genes,cells));map=fread(file.path(B,'00_plan_audit/epi_signature_feature_mapping.tsv'));cv=fread(file.path(B,'00_plan_audit/epi_coverage.tsv'));sets=split(map[!is.na(input_gene)&signature%in%cv[status=='READY',signature],input_gene],map[!is.na(input_gene)&signature%in%cv[status=='READY',signature],signature]);started=proc.time()[3];out=list()
if(!'--cycle-only'%in%args){
for(a in seq(1,ncol(x),256)){
 ii=a:min(a+255,ncol(x));set.seed(20261111+a);rank=AUCell_buildRankings(x[,ii,drop=FALSE],plotStats=FALSE,splitByBlocks=TRUE,BPPARAM=BiocParallel::SerialParam(),verbose=FALSE);auc=AUCell_calcAUC(sets,rank,aucMaxRank=4501,normAUC=TRUE,nCores=1,verbose=FALSE);rm(rank)
 uc=ScoreSignatures_UCell(x[,ii,drop=FALSE],features=sets,maxRank=1500,ncores=1,BPPARAM=BiocParallel::SerialParam(),chunk.size=256);aa=as.data.table(t(getAUC(auc)),keep.rownames='cell_id');aa=melt(aa,id.vars='cell_id',variable.name='signature',value.name='score');aa[,method:='official_AUCell'];uu=as.data.table(uc,keep.rownames='cell_id');setnames(uu,names(uu),sub('_UCell$','',names(uu)));uu=melt(uu,id.vars='cell_id',variable.name='signature',value.name='score');uu[,method:='UCell'];out[[length(out)+1]]=rbind(aa,uu);rm(auc,uc);gc(FALSE)
}
fwrite(rbindlist(out),file.path(input,'official_scores.tsv.gz'),sep='\t');cat('ranking_seconds',proc.time()[3]-started,'n',ncol(x),'\n')
}
# Standard Seurat RNA cycle scoring; per-specimen control bins, full original background.
cc=readRDS(file.path(B,'01_signatures/cell_cycle_Seurat2019.rds'));obj=CreateSeuratObject(x,min.cells=0,min.features=0);obj=NormalizeData(obj,verbose=FALSE);set.seed(20261111);obj=CellCycleScoring(obj,s.features=cc$s.genes,g2m.features=cc$g2m.genes,set.ident=FALSE,seed=20261111);cycle=obj@meta.data[,c('S.Score','G2M.Score','Phase')];fwrite(as.data.table(cycle,keep.rownames='cell_id'),file.path(input,'cycle_scores.tsv'),sep='\t');cat('total_seconds',proc.time()[3]-started,'\n');capture.output(sessionInfo(),file=file.path(input,'sessionInfo.txt'))
