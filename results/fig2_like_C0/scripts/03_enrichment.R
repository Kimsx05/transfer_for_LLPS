source('results/fig2_like_C0/scripts/common.R')
suppressPackageStartupMessages({library(clusterProfiler);library(org.Hs.eg.db);library(fgsea)})
a <- fread(file.path(CORE,'C0_vs_others_DEG.tsv'))
stopifnot(nrow(a)>0,any(a$C0_up))
status('G_GO','Mapping tested genes and using the tested-gene background')
# Exact SYMBOL mapping only. Ambiguous symbols and duplicate Entrez IDs cannot
# inflate the foreground or background; choose highest logCPM for duplicate IDs.
valid <- intersect(a$gene,AnnotationDbi::keys(org.Hs.eg.db,keytype='SYMBOL'))
mp <- as.data.table(AnnotationDbi::select(org.Hs.eg.db,keys=valid,keytype='SYMBOL',columns='ENTREZID'))
mp <- unique(mp[!is.na(ENTREZID)])
mp[,n_entrez:=uniqueN(ENTREZID),by=SYMBOL]
mp <- mp[n_entrez==1]
mapped <- merge(a,mp,by.x='gene',by.y='SYMBOL')
setorder(mapped,-logCPM,gene);mapped[,id_representative:=!duplicated(ENTREZID)]
audit <- merge(a[,.(gene,C0_up)],mapped[,.(gene,ENTREZID,id_representative)],by='gene',all.x=TRUE)
audit[,mapping_status:=fcase(is.na(ENTREZID),'UNMAPPED_OR_AMBIGUOUS',id_representative,'USED',default='DUPLICATE_ENTREZ_NONREPRESENTATIVE')]
write_tsv(audit,'enrichment_gene_mapping_audit.tsv')
mapped <- mapped[id_representative==TRUE]
universe <- mapped$ENTREZID;up <- mapped[C0_up==TRUE,ENTREZID]
write_tsv(mapped[,.(gene,ENTREZID,C0_up)],'enrichment_foreground_background.tsv')
go <- enrichGO(gene=up,universe=universe,OrgDb=org.Hs.eg.db,keyType='ENTREZID',ont='BP',
 pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,minGSSize=10,maxGSSize=500,readable=TRUE)
got <- as.data.table(as.data.frame(go));stopifnot(nrow(got)>0)
write_tsv(got,'C0_GO_enrichment.tsv')

# Greedy overlap pruning affects display only. Complete results remain unchanged.
representative_terms <- function(tab,n=12) {
 tab <- copy(tab[p.adjust<.05]);setorder(tab,p.adjust,pvalue,ID)
 selected <- integer();sets <- strsplit(tab$geneID,'/',fixed=TRUE)
 for(i in seq_len(nrow(tab))) {
  redundant <- any(vapply(selected,function(j)length(intersect(sets[[i]],sets[[j]]))/length(union(sets[[i]],sets[[j]]))>=.6,logical(1)))
  if(!redundant)selected<-c(selected,i)
  if(length(selected)>=n)break
 }
 tab[selected]
}
bar_enrichment <- function(tab,kind) {
 top <- representative_terms(tab)
 write_tsv(top,paste0('C0_',kind,'_displayed_terms.tsv'))
 if(!nrow(top)) {
  p <- ggplot()+annotate('text',x=0,y=0,label='No terms pass BH FDR < 0.05',size=4)+theme_void()+labs(title=paste('C0 upregulated genes:',kind))
 } else {
  top[,label:=vapply(Description,function(x)paste(strwrap(x,48),collapse='\n'),character(1))]
  top[,label:=factor(label,levels=rev(unique(label)))]
  p <- ggplot(top,aes(Count,label,fill=-log10(p.adjust)))+geom_col(width=.72)+
   scale_fill_gradient(low='#F4DA93',high='#BE4B5A',name=expression(-log[10]('FDR')))+
   theme_pub()+theme(axis.text.y=element_text(size=9),legend.position='right')+
   labs(title=paste('C0 upregulated genes:',if(kind=='GO')'GO Biological Process' else 'KEGG'),
    x='Number of upregulated genes',y=NULL,caption='BH FDR < 0.05; redundant terms omitted from display only.')
 }
 export_plot(p,paste0('Fig2',if(kind=='GO')'G' else 'H','_C0_',kind,'_barplot'),8.2,6.4)
}
bar_enrichment(got,'GO');status('G_GO','COMPLETE')

status('H_KEGG','Using the same foreground and tested-gene background as GO')
if(file.exists(file.path(OUT,'resources/KEGG_hsa_gene_pathway.tsv.part')))stop('KEGG download is incomplete; do not use partial annotation')
k <- fread(file.path(OUT,'resources/KEGG_hsa_gene_pathway.tsv'),header=FALSE,col.names=c('gene','term'))
kn <- fread(file.path(OUT,'resources/KEGG_hsa_pathway_names.tsv'),header=FALSE,col.names=c('term','name'))
stopifnot(nrow(k)>30000,all(grepl('^hsa:[0-9]+$',k$gene)),all(grepl('^path:hsa[0-9]+$',k$term)))
k[,`:=`(gene=sub('^hsa:','',gene),term=sub('^path:','',term))]
kn[,`:=`(term=sub('^path:','',term),name=sub(' - Homo sapiens \\(human\\)$','',name))]
stopifnot(all(k$term%in%kn$term))
ke <- enricher(gene=up,universe=universe,TERM2GENE=k[,.(term,gene)],TERM2NAME=kn,
 pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,minGSSize=10,maxGSSize=500)
kt <- as.data.table(as.data.frame(ke));stopifnot(nrow(kt)>0)
write_tsv(kt,'C0_KEGG_enrichment.tsv');bar_enrichment(kt,'KEGG')
write_tsv(data.table(database=c('GO_BP','KEGG'),n_mapped_tested=length(universe),n_mapped_up=length(up),
 n_annotated_background=c(length(go@universe),length(ke@universe)),n_result_terms=c(nrow(got),nrow(kt)),
 n_FDR05=c(sum(got$p.adjust<.05),sum(kt$p.adjust<.05))), 'enrichment_analysis_summary.tsv')
status('H_KEGG','COMPLETE')

status('I_GSEA','Hallmark v2025.1.Hs; signed sqrt(QL F); all tested genes')
gmtfile <- file.path(OLD,'00_preflight/h.all.v2025.1.Hs.symbols.gmt')
file.copy(gmtfile,file.path(OUT,'resources/h.all.v2025.1.Hs.symbols.gmt'),overwrite=TRUE)
pathways <- gmtPathways(gmtfile)
r <- a[order(-rank_statistic,gene)];ranks <- setNames(r$rank_statistic,r$gene)
stopifnot(!anyDuplicated(names(ranks)),all(is.finite(ranks)),length(pathways)==50)
set.seed(20261025)
gs <- fgseaMultilevel(pathways=pathways,stats=ranks,minSize=10,maxSize=500,eps=0,
 BPPARAM=BiocParallel::SerialParam())
setorder(gs,padj)
gs_out <- copy(gs);gs_out[,leadingEdge:=vapply(leadingEdge,paste,character(1),collapse=';')]
write_tsv(gs_out,'C0_GSEA_results.tsv')
preferred <- paste0('HALLMARK_',c('EPITHELIAL_MESENCHYMAL_TRANSITION','ANGIOGENESIS','TGF_BETA_SIGNALING',
 'INFLAMMATORY_RESPONSE','HYPOXIA','CHOLESTEROL_HOMEOSTASIS','FATTY_ACID_METABOLISM'))
selected <- preferred[preferred%in%gs[padj<.05,pathway]]
# Prefer significant requested pathways; fill to five with positive, then negative,
# FDR-ranked alternatives. Never force a non-significant pathway into a figure.
if(length(selected)<5)selected <- unique(c(selected,head(gs[padj<.05&NES>0&!pathway%in%selected,pathway],5-length(selected))))
if(length(selected)<5)selected <- unique(c(selected,head(gs[padj<.05&!pathway%in%selected,pathway],5-length(selected))))
aliases <- c(HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION='EMT',HALLMARK_ANGIOGENESIS='ANGIOGENESIS',HALLMARK_TGF_BETA_SIGNALING='TGFB')
selection <- merge(data.table(pathway=unique(c(preferred,selected))),gs_out[,.(pathway,NES,pval,padj)],by='pathway',all.x=TRUE)
selection[,`:=`(exported=pathway%in%selected,significant_FDR05=padj<.05,
 direction=fifelse(NES>0,'positive_NES_C0','negative_NES_Others'))]
write_tsv(selection,'GSEA_pathway_plot_selection.tsv')
if(!length(selected))stop('No Hallmark pathway passes FDR < 0.05; no forced GSEA plots')
for(path in selected) {
 row <- gs[pathway==path]
 short <- if(path%in%names(aliases))unname(aliases[path])else sub('^HALLMARK_','',path)
 label <- gsub('_',' ',sub('^HALLMARK_','',path))
 hits <- names(ranks)%in%pathways[[path]]
 w <- abs(ranks);running <- cumsum(ifelse(hits,w/sum(w[hits]),-1/sum(!hits)))
 curve <- data.table(rank=seq_along(ranks),gene=names(ranks),statistic=as.numeric(ranks),hit=hits,running_ES=as.numeric(running))
 # A single standard enrichment plot: red running curve and a separate hit strip
 # on the same rank axis. No multi-pathway panels or combined figure sheets.
 span <- diff(range(running));hit_bottom <- min(running)-.18*span
 p <- ggplot(curve,aes(rank,running_ES))+
  geom_hline(yintercept=0,color='grey65',linewidth=.35,linetype='dashed')+
  geom_line(color='#C43B3B',linewidth=.85)+
  geom_segment(data=curve[hit==TRUE],aes(x=rank,xend=rank),y=hit_bottom,yend=hit_bottom+.075*span,
   inherit.aes=FALSE,linewidth=.25,color='#222222')+
  expand_limits(y=hit_bottom-.025*span)+
  scale_y_continuous(breaks=pretty(range(running),5))+
  scale_x_continuous(expand=expansion(mult=c(.005,.005)))+theme_pub()+
  labs(title=label,subtitle=sprintf('NES = %.2f | P = %s | BH FDR = %s',row$NES,ptext(row$pval),ptext(row$padj)),
   x='Rank in ordered gene list (C0 to Others)',y='Running enrichment score',
   caption='Rank: sign(log2FC) x sqrt(QL F). Black ticks: pathway genes.')
 export_plot(p,paste0('Fig2I_GSEA_',short),7,4.7)
 write_tsv(curve,paste0('GSEA_',short,'_curve_data.tsv.gz'))
}
capture.output(sessionInfo(),file=file.path(OUT,'summary/sessionInfo_enrichment.txt'))
capture.output(AnnotationDbi::metadata(org.Hs.eg.db),file=file.path(OUT,'summary/org.Hs.eg.db_metadata.txt'))
status('I_GSEA',sprintf('COMPLETE: %d significant Hallmarks; %d individual plots',sum(gs$padj<.05),length(selected)))
