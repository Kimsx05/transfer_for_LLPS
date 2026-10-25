source('results/fig2_like_C0/phase1B_review/scripts/common.R')
suppressPackageStartupMessages({library(edgeR);library(fgsea);library(clusterProfiler);library(org.Hs.eg.db)})
p<-readRDS(file.path(CORE,'stromal_cluster_pseudobulk.rds'));meta<-p$metadata
original_samples<-fread(file.path(PHASE1,'pseudobulk_sample_eligibility.tsv'))[included==TRUE,sample]
sets<-gmtPathways(file.path(OUT,'resources/h.all.v2025.1.Hs.symbols.gmt'))
k<-fread(file.path(OUT,'resources/KEGG_hsa_gene_pathway.tsv'),header=FALSE,col.names=c('gene','term'))
k[,`:=`(gene=sub('hsa:','',gene,fixed=TRUE),term=sub('path:','',term,fixed=TRUE))]
kn<-fread(file.path(OUT,'resources/KEGG_hsa_pathway_names.tsv'),header=FALSE,col.names=c('term','name'))
kn[,name:=sub(' - Homo sapiens \\(human\\)$','',name)]
run_model<-function(comparators,name,restrict_samples=NULL) {
 mm<-copy(meta[cluster%in%c('C0',comparators)&sample%in%original_samples]);mm[,group:=ifelse(cluster=='C0','C0','Others')]
 n<-mm[,.(n_cells=sum(n_cells)),by=.(sample,group)]
 elig<-n[,.(both=.N==2,minimum=min(n_cells)),by=sample][both&minimum>=20,sample]
 if(!is.null(restrict_samples))elig<-intersect(elig,restrict_samples)
 stopifnot(length(elig)>=4)
 mm<-mm[sample%in%elig];key<-paste(mm$sample,mm$group,sep='::');lev<-sort(unique(key))
 ag<-sparseMatrix(i=seq_along(key),j=match(key,lev),x=1,dims=c(length(key),length(lev)))
 pb<-as.matrix(p$counts[,mm$library]%*%ag);colnames(pb)<-lev
 md<-unique(mm[,.(sample,dataset,tissue_origin,patient_key,group)]);md[,library:=paste(sample,group,sep='::')];md<-md[match(lev,library)]
 md[,`:=`(sample=factor(sample),group=factor(group,levels=c('Others','C0')))]
 design<-model.matrix(~sample+group,md);stopifnot(qr(design)$rank==ncol(design),!anyDuplicated(md[group=='C0',patient_key]))
 y<-DGEList(pb);y<-y[filterByExpr(y,group=md$group),,keep.lib.sizes=FALSE];y<-calcNormFactors(y)
 y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE);q<-glmQLFTest(fit,coef='groupC0')
 a<-as.data.table(topTags(q,n=Inf,sort.by='none')$table,keep.rownames='gene');a[,`:=`(C0_up=FDR<.05&logFC>.25,rank_statistic=sign(logFC)*sqrt(pmax(F,0)))]
 WT(a,paste0(name,'_DEG.tsv'));WT(md,paste0(name,'_samples.tsv'))
 r<-a[order(-rank_statistic,gene)];ranks<-setNames(r$rank_statistic,r$gene);set.seed(123)
 gs<-fgseaMultilevel(pathways=sets,stats=ranks,minSize=10,maxSize=500,eps=0,BPPARAM=BiocParallel::SerialParam())
 gs[,leadingEdge:=vapply(leadingEdge,paste,character(1),collapse=';')];WT(gs,paste0(name,'_GSEA.tsv'))
 valid<-intersect(a$gene,AnnotationDbi::keys(org.Hs.eg.db,keytype='SYMBOL'))
 mp<-as.data.table(AnnotationDbi::select(org.Hs.eg.db,keys=valid,keytype='SYMBOL',columns='ENTREZID'))
 mp<-unique(mp[!is.na(ENTREZID)]);mp[,n_ids:=uniqueN(ENTREZID),by=SYMBOL];mp<-mp[n_ids==1]
 ma<-merge(a,mp,by.x='gene',by.y='SYMBOL');setorder(ma,-logCPM,gene);ma<-ma[!duplicated(ENTREZID)]
 up<-ma[C0_up==TRUE,ENTREZID];univ<-ma$ENTREZID
 go<-as.data.table(as.data.frame(enrichGO(up,universe=univ,OrgDb=org.Hs.eg.db,keyType='ENTREZID',ont='BP',
  pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,minGSSize=10,maxGSSize=500,readable=TRUE)))
 ke<-as.data.table(as.data.frame(enricher(up,universe=univ,TERM2GENE=k[,.(term,gene)],TERM2NAME=kn,
  pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,minGSSize=10,maxGSSize=500)))
 WT(go,paste0(name,'_GO.tsv'));WT(ke,paste0(name,'_KEGG.tsv'))
 WT(ma[,.(gene,ENTREZID,C0_up)],paste0(name,'_foreground_background.tsv'))
 list(DEG=a,GSEA=gs,GO=go,KEGG=ke,samples=elig,audit=data.table(model=name,comparator=paste(comparators,collapse='+'),
  n_samples=length(elig),n_Tumor=sum(md[group=='C0',tissue_origin]=='Tumor'),n_Normal=sum(md[group=='C0',tissue_origin]=='Normal'),
  tested=nrow(a),up=sum(a$C0_up)))
}
fib<-run_model(c('C4','C5'),'C0_vs_fibro45')
cat('Fibro45 complete\n')
inc2<-run_model(c('C2','C4','C5'),'C0_vs_fibro245')
cat('Fibro245 complete\n')
matched<-run_model(paste0('C',1:6),'C0_vs_all_matched_fibro45',fib$samples)
WT(rbindlist(list(fib$audit,inc2$audit,matched$audit)),'functional_sensitivity_models.tsv')
original<-list(DEG=fread(file.path(PHASE1,'C0_vs_others_DEG.tsv')),GSEA=fread(file.path(PHASE1,'C0_GSEA_results.tsv')),
 GO=fread(file.path(PHASE1,'C0_GO_enrichment.tsv')),KEGG=fread(file.path(PHASE1,'C0_KEGG_enrichment.tsv')))
models<-list(original_all=original,all_same11=matched,fibro_C4_C5=fib,fibro_C2_C4_C5=inc2)
for(db in c('GSEA','GO','KEGG')){
 combined<-rbindlist(lapply(names(models),function(n){d<-copy(models[[n]][[db]]);d[,model:=n];d}),fill=TRUE)
 WT(combined,paste0(db,'_all_comparator_results.tsv'))
}
gs<-fread(file.path(CORE,'GSEA_all_comparator_results.tsv'))
focusgs<-gs[pathway%in%names(short)];focusgs[,pathway_label:=factor(short[pathway],levels=rev(unname(short)))]
focusgs[,model:=factor(model,levels=names(models))]
pplot<-ggplot(focusgs,aes(model,pathway_label,fill=NES))+geom_tile(color='white',linewidth=.7)+
 geom_text(aes(label=sprintf('%.2f\nq=%s',NES,vapply(padj,ptext,character(1)))),size=3)+
 scale_fill_gradient2(low='#3B6A95',mid='white',high='#BF4C52',midpoint=0,name='NES')+
 scale_x_discrete(labels=c('Original\nall others','Same 11 samples\nall others','Fibro-like\nC4 + C5','Include unresolved\nC2 + C4 + C5'))+
 theme_pub()+theme(axis.line=element_blank(),axis.ticks=element_blank())+
 labs(title='C0 functional signals: comparator sensitivity',x=NULL,y=NULL,caption='All pathways shown irrespective of significance. q = Hallmark-wide BH FDR.')
export_plot(pplot,'FigS_functional_comparator_GSEA',9,5)
WT(focusgs,'core_GSEA_comparator_comparison.tsv')
for(db in c('GO','KEGG')){
 showids<-fread(file.path(PHASE1,paste0('C0_',db,'_displayed_terms.tsv')))$ID
 comb<-fread(file.path(CORE,paste0(db,'_all_comparator_results.tsv')))
 WT(comb[ID%in%showids],paste0(db,'_original_top_terms_comparison.tsv'))
}

# Leading-edge membership is inherited from the completed primary GSEA.
le<-rbindlist(lapply(names(short),function(path){row<-original$GSEA[pathway==path];data.table(pathway=path,gene=strsplit(row$leadingEdge,';',fixed=TRUE)[[1]])}))
valid<-intersect(unique(le$gene),AnnotationDbi::keys(org.Hs.eg.db,keytype='SYMBOL'))
ann<-as.data.table(AnnotationDbi::select(org.Hs.eg.db,keys=valid,keytype='SYMBOL',columns=c('GOALL','ONTOLOGYALL')))
ecm<-unique(ann[GOALL%in%c('GO:0031012','GO:0062023'),SYMBOL])
le[,`:=`(GO_ECM=gene%in%ecm,COL_or_FN1=grepl('^COL[0-9]',gene)|gene=='FN1')]
le<-merge(le,original$DEG[,.(gene,primary_logFC=logFC,primary_FDR=FDR,rank_statistic)],by='gene')
le<-merge(le,fib$DEG[,.(gene,fibro45_logFC=logFC,fibro45_FDR=FDR)],by='gene',all.x=TRUE)
freq<-le[,.(n_pathways=uniqueN(pathway)),by=gene];le<-merge(le,freq,by='gene')
WT(le,'core_GSEA_leading_edge_gene_pathway.tsv')
summ<-le[,.(n_leading=.N,n_GO_ECM=sum(GO_ECM),n_COL_FN1=sum(COL_or_FN1),fraction_GO_ECM=mean(GO_ECM),
 n_shared_pathways=sum(n_pathways>1),n_fibro45_tested=sum(!is.na(fibro45_FDR)),
 n_fibro45_up=sum(fibro45_logFC>.25&fibro45_FDR<.05,na.rm=TRUE)),by=pathway]
WT(summ,'leading_edge_specificity_summary.tsv')
overlap<-CJ(pathway1=names(short),pathway2=names(short));overlap[,`:=`(n_overlap=0L,Jaccard=0)]
for(i in seq_len(nrow(overlap))){a<-le[pathway==overlap$pathway1[i],gene];b<-le[pathway==overlap$pathway2[i],gene];overlap[i,`:=`(n_overlap=length(intersect(a,b)),Jaccard=length(intersect(a,b))/length(union(a,b)))]}
WT(overlap,'leading_edge_pathway_overlap.tsv')
chosen<-unique(unlist(lapply(names(short),function(path)head(le[pathway==path][order(-rank_statistic,gene),gene],6))))
topshared<-freq[order(-n_pathways,gene)][n_pathways>=3,gene];chosen<-unique(c(chosen,head(topshared,8)))
plotd<-le[gene%in%chosen];plotd[,gene:=factor(gene,levels=rev(sort(chosen)))];plotd[,label:=factor(short[pathway],levels=unname(short))]
pplot<-ggplot(plotd,aes(label,gene))+geom_point(aes(color=primary_logFC,shape=GO_ECM),size=3.3)+
 scale_color_gradient(low='#E9DAB1',high='#BE4654',name='Primary log2FC')+
 scale_shape_manual(values=c('FALSE'=16,'TRUE'=17),name='GO extracellular\nmatrix annotation')+
 theme_pub()+theme(axis.text.x=element_text(angle=25,hjust=1),axis.text.y=element_text(size=9))+
 labs(title='C0 leading-edge genes across core pathways',x=NULL,y=NULL,
 caption='Six highest-ranked leading-edge genes per pathway plus shared genes (>=3 pathways).\nComplete membership and comparator-specific DE are available in the tables.')
export_plot(pplot,'FigS_leading_edge_gene_pathway',8.4,max(6,.21*length(chosen)+2))
WT(data.table(gene=chosen),'leading_edge_displayed_genes.tsv')
capture.output(sessionInfo(),file=file.path(CORE,'sessionInfo.txt'))
cat('Functional specificity review complete\n')
