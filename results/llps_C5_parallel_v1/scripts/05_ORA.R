source('results/llps_C5_parallel_v1/scripts/common.R')
suppressPackageStartupMessages({library(clusterProfiler);library(org.Hs.eg.db)})
direction<-Sys.getenv('ORA_DIRECTION','up'); if(direction=='down'){M2<-file.path(M2,'down');dir.create(M2,showWarnings=FALSE)}
el<-fread(file.path(M1,'ORA_detection_eligibility.tsv'));up<-fread(file.path(M1,paste0('markers_C5_',direction,'.tsv')))
# Detectability eligibility uses the specified 25% threshold in either group,
# independently of direction or fold change. Conditioning the universe on positive
# logFC would reuse the foreground selection and is intentionally avoided.
sym<-unique(el[detectable_universe==TRUE,symbol]);valid<-intersect(sym,AnnotationDbi::keys(org.Hs.eg.db,keytype='SYMBOL'))
mp<-as.data.table(AnnotationDbi::select(org.Hs.eg.db,keys=valid,keytype='SYMBOL',columns='ENTREZID'));mp<-unique(mp[!is.na(ENTREZID)])
mp[,n_ids:=uniqueN(ENTREZID),by=SYMBOL];amb<-mp[n_ids>1];mp<-mp[n_ids==1]
upsym<-unique(el[gene%in%up$gene,symbol]);ids<-unique(mp$ENTREZID);upid<-unique(mp[SYMBOL%in%upsym,ENTREZID])
wt(merge(el[detectable_universe==TRUE],mp,by.x='symbol',by.y='SYMBOL',all.x=TRUE),file.path(M2,'background_mapping.tsv'))
wt(data.table(symbol=setdiff(sym,mp$SYMBOL)),file.path(M2,'unmapped_or_ambiguous_symbols.tsv'))
ann<-as.data.table(AnnotationDbi::select(org.Hs.eg.db,keys=ids,keytype='ENTREZID',columns=c('GOALL','ONTOLOGYALL')))
gouniv<-unique(ann[ONTOLOGYALL=='BP'&!is.na(GOALL),ENTREZID]);gofg<-intersect(upid,gouniv)
go<-enrichGO(gene=gofg,universe=gouniv,OrgDb=org.Hs.eg.db,keyType='ENTREZID',ont='BP',pAdjustMethod='BH',pvalueCutoff=.05,qvalueCutoff=.05,minGSSize=10,maxGSSize=500,readable=TRUE)
stopifnot(!is.null(go));wt(as.data.table(go@result),file.path(M2,'GO_BP_all.tsv'))
gosig<-as.data.table(as.data.frame(go));wt(gosig,file.path(M2,'GO_BP_significant.tsv'))
# Run the requested enrichKEGG(hsa) against an explicitly restored project KEGG
# REST snapshot, preserving reproducibility and avoiding a fresh annotation version.
k<-fread(file.path(BASE,'results/fig2_like_C0/resources/KEGG_hsa_gene_pathway.tsv'),header=FALSE,col.names=c('gene','term'));k[,`:=`(gene=sub('hsa:','',gene,fixed=TRUE),term=sub('path:','',term,fixed=TRUE))]
kn<-fread(file.path(BASE,'results/fig2_like_C0/resources/KEGG_hsa_pathway_names.tsv'),header=FALSE,col.names=c('term','name'));kn[,name:=sub(' - Homo sapiens \\(human\\)$','',name)]
ke<-clusterProfiler:::get_KEGG_Env();assign('organism','hsa',envir=ke);assign('_type_','KEGG',envir=ke)
assign('KEGGPATHID2EXTID',as.data.frame(k[,.(term,gene)]),envir=ke);assign('KEGGPATHID2NAME',as.data.frame(kn),envir=ke)
kuniv<-intersect(ids,k$gene);kfg<-intersect(upid,kuniv)
kk<-enrichKEGG(gene=kfg,organism='hsa',keyType='kegg',universe=kuniv,pAdjustMethod='BH',pvalueCutoff=.05,qvalueCutoff=.05,minGSSize=10,maxGSSize=500)
stopifnot(!is.null(kk));kk<-setReadable(kk,OrgDb=org.Hs.eg.db,keyType='ENTREZID')
wt(as.data.table(kk@result),file.path(M2,'KEGG_all.tsv'));ksig<-as.data.table(as.data.frame(kk));wt(ksig,file.path(M2,'KEGG_significant.tsv'))
wt(data.table(database=c('GO_BP','KEGG'),detectable_original_features=sum(el$detectable_universe),detectable_unique_symbols=length(sym),symbol_to_Entrez_success=length(unique(mp$SYMBOL)),unmapped_or_ambiguous_symbols=length(setdiff(sym,mp$SYMBOL)),mapped_unique_Entrez=length(ids),effective_background=c(length(gouniv),length(kuniv)),input_marker_features=nrow(up),input_unique_marker_symbols=length(upsym),mapped_marker_Entrez=length(upid),effective_foreground=c(length(gofg),length(kfg)),reported_terms=c(nrow(go@result),nrow(kk@result)),significant_terms=c(nrow(gosig),nrow(ksig))),file.path(M2,'ORA_summary.tsv'))
wt(data.table(database=c('GO_BP','KEGG'),source=c(paste0('org.Hs.eg.db ',packageVersion('org.Hs.eg.db')),'Project KEGG REST human snapshot'),source_file=c(find.package('org.Hs.eg.db'),file.path(BASE,'results/fig2_like_C0/resources/KEGG_hsa_gene_pathway.tsv')),snapshot_mtime=c(NA_character_,as.character(file.info(file.path(BASE,'results/fig2_like_C0/resources/KEGG_hsa_gene_pathway.tsv'))$mtime))),file.path(M2,'ORA_annotation_sources.tsv'))
wt(data.table(ENTREZID=gouniv,in_foreground=gouniv%in%gofg),file.path(M2,'GO_background_foreground.tsv'))
wt(data.table(ENTREZID=kuniv,in_foreground=kuniv%in%kfg),file.path(M2,'KEGG_background_foreground.tsv'))
# Simplify significant GO terms only; preserve the untouched complete output above.
go_for_simplify<-go;go_for_simplify@result<-as.data.frame(gosig)
if(nrow(gosig)>1){cat('GO simplify on',nrow(gosig),'significant terms\n');flush.console();gos<-as.data.table(as.data.frame(simplify(go_for_simplify,cutoff=.7,by='p.adjust',select_fun=min,measure='Wang')))}else gos<-gosig
wt(gos,file.path(M2,'GO_BP_simplified.tsv'))
cat('ORA complete\n')
