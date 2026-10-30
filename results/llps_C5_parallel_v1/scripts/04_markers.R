source('results/llps_C5_parallel_v1/scripts/common.R')
z<-load_frozen();o<-z$object;d<-z$metadata;x<-z$expression
wt(d[,.(cell,sample,patient_key,dataset,cluster,group,LLPS77_AUC)],file.path(M1,'cell_metadata.tsv.gz'))
wt(d[,.N,by=.(dataset,sample,patient_key,cluster)],file.path(M1,'sample_cluster_coverage.tsv'))
wt(data.table(input=OBJECT,n_genes=nrow(x),n_cells=ncol(x),n_C5=sum(d$group=='C5'),n_Rest=sum(d$group=='Rest'),expression='RNA/data',normalization=z$norm$normalization.method,scale_factor=z$norm$scale.factor,max_norm_check_error=z$max_norm_error),file.path(M1,'input_summary.tsv'))
if(file.exists(file.path(M1,'gene_ranking_all.tsv')) && file.exists(file.path(M1,'markers_C5_up.tsv'))){
r<-fread(file.path(M1,'gene_ranking_all.tsv'));a<-fread(file.path(M1,'markers_all.tsv'));up<-fread(file.path(M1,'markers_C5_up.tsv'))
}else{
cat('FindMarkers: complete unfiltered signed comparison\n');flush.console()
r<-FindMarkers(o,ident.1='C5',ident.2=setdiff(CLUSTERS,'C5'),assay='RNA',slot='data',test.use='wilcox',min.pct=0,logfc.threshold=0,only.pos=FALSE,random.seed=123,verbose=FALSE)
r<-as.data.table(r,keep.rownames='gene');r[,`:=`(expression_layer='RNA/data',min_pct_threshold=0,logFC_prefilter=0)]
wt(r,file.path(M1,'gene_ranking_all.tsv'))
cat('FindMarkers: prespecified positive candidates\n');flush.console()
a<-FindMarkers(o,ident.1='C5',ident.2=setdiff(CLUSTERS,'C5'),assay='RNA',slot='data',test.use='wilcox',min.pct=.25,logfc.threshold=.75,only.pos=TRUE,random.seed=123,verbose=FALSE)
a<-as.data.table(a,keep.rownames='gene');a[,`:=`(expression_layer='RNA/data',criteria='max(pct.1,pct.2)>=0.25; avg_log2FC>0.75; Bonferroni p_val_adj<0.05')]
setorder(a,-avg_log2FC,p_val_adj,gene);wt(a,file.path(M1,'markers_all.tsv'))
up<-a[p_val_adj<.05&avg_log2FC>.75&pmax(pct.1,pct.2)>=.25];wt(up,file.path(M1,'markers_C5_up.tsv'))
wt(head(up,30),file.path(M1,'top30_markers.tsv'))
}
# Mapping is annotation only: preserve all original normalized expression rows.
gm<-fread(file.path(BASE,'fig2_C0_multimethod_deconvolution/03_methods/CIBERSORTx_web/resources/gencode_v36_gene_map.tsv'));gm[,id:=sub('\\..*$','',id)];gm<-unique(gm[,.(id,gene)]);gm[,n_symbols:=uniqueN(gene),by=id];gm<-gm[n_symbols==1]
mp<-data.table(gene=rownames(x),mean_normalized=Matrix::rowMeans(x),n_detected=Matrix::rowSums(x>0))
mp[,symbol:=gene];ei<-grep('^ENSG',mp$gene);mapped<-gm$gene[match(sub('\\..*$','',mp$gene[ei]),gm$id)];mp[ei[!is.na(mapped)],symbol:=mapped[!is.na(mapped)]]
# If identifiers converge, use the most detected row, then highest overall mean,
# then original name; no selection by p-value or directional fold change.
setorder(mp,symbol,-n_detected,-mean_normalized,gene);mp[,representative:=!duplicated(symbol)]
wt(mp,file.path(M1,'gene_identifier_mapping.tsv'))
rr<-merge(r,mp,by='gene');rr[,detectable_universe:=pmax(pct.1,pct.2)>=.25&n_detected>0&is.finite(avg_log2FC)&is.finite(p_val)]
wt(rr[,.(gene,symbol,representative,n_detected,detectable_universe,pct.1,pct.2)],file.path(M1,'ORA_detection_eligibility.tsv'))
# Full finite, expressed ranking; all-zero and non-testable rows are explicitly flagged.
rr[,rank_eligible:=representative&n_detected>0&is.finite(avg_log2FC)&is.finite(p_val)]
rank<-rr[rank_eligible==TRUE];setorder(rank,-avg_log2FC,symbol);rank[,rank_position:=.I]
wt(rank[,.(gene,symbol,avg_log2FC,p_val,p_val_adj,rank_position)],file.path(M1,'GSEA_preranked_genes.tsv'))
wt(data.table(total_expression_rows=nrow(x),returned_full_comparison=nrow(r),ranked_symbols=nrow(rank),duplicate_symbol_rows=sum(!mp$representative),nonexpressed_rows=sum(mp$n_detected==0),tied_rank_rows=sum(duplicated(rank$avg_log2FC)|duplicated(rank$avg_log2FC,fromLast=TRUE)),ranking='avg_log2FC descending; exact ties ordered by symbol; no numerical jitter',markers_candidate=nrow(a),significant_up=nrow(up)),file.path(M1,'marker_and_ranking_QC.tsv'))
# Top markers + prespecified existing C5 genes. No sample-count exclusion.
focus<-unique(c(head(up$gene,30),c('POSTN','COL11A1','COL10A1')));focus<-intersect(focus,rownames(x))
xx<-x[focus,,drop=FALSE]
review<-list()
for(level in c('sample','patient_key','dataset')){
 for(id in unique(d[[level]]))for(g in c('C5','Rest')){
  ii<-which(d[[level]]==id&d$group==g);if(!length(ii))next
  review[[length(review)+1]]<-data.table(level=level,unit=id,group=g,gene=focus,n_cells=length(ii),pct_detected=Matrix::rowMeans(xx[,ii,drop=FALSE]>0),mean_log_normalized=Matrix::rowMeans(xx[,ii,drop=FALSE]),mean_linear_normalized=Matrix::rowMeans(expm1(xx[,ii,drop=FALSE])))
 }
}
rv<-rbindlist(review);wt(rv,file.path(M1,'marker_sample_review.tsv'))
wide<-dcast(rv,level+unit+gene~group,value.var=c('n_cells','pct_detected','mean_log_normalized','mean_linear_normalized'))
wide[,`:=`(both_groups=is.finite(n_cells_C5)&is.finite(n_cells_Rest),mean_log_difference=mean_log_normalized_C5-mean_log_normalized_Rest,descriptive_log2_ratio=log2((mean_linear_normalized_C5+1e-9)/(mean_linear_normalized_Rest+1e-9)))]
wt(wide,file.path(M1,'marker_unit_contrasts.tsv'))
consistent<-wide[both_groups==TRUE,.(n_units=.N,n_positive=sum(mean_log_difference>0),fraction_positive=mean(mean_log_difference>0),median_log_difference=median(mean_log_difference),min_C5_cells=min(n_cells_C5),min_Rest_cells=min(n_cells_Rest)),by=.(level,gene)]
wt(consistent,file.path(M1,'marker_unit_consistency.tsv'))
wt(r[gene%in%c('POSTN','COL11A1','COL10A1')],file.path(M1,'prespecified_C5_gene_review.tsv'))
flags<-copy(up);flags[,`:=`(mitochondrial=grepl('^MT-',gene),ribosomal=grepl('^RP[SL][0-9]',gene),stress_candidate=grepl('^(HSP|DNAJ)',gene)|gene%in%c('FOS','FOSB','JUN','JUNB','JUND','ATF3','DDIT3','DDIT4','EGR1','HSPA5','XBP1'))]
wt(flags[,.(gene,avg_log2FC,p_val_adj,mitochondrial,ribosomal,stress_candidate)],file.path(M1,'marker_QC_flags.tsv'))
# Dot plot displays frozen expression, not additional inference.
showgenes<-head(up$gene,20);dots<-rbindlist(lapply(CLUSTERS,function(cl){ii<-which(d$cluster==cl);data.table(gene=showgenes,cluster=cl,pct=100*Matrix::rowMeans(x[showgenes,ii,drop=FALSE]>0),mean=Matrix::rowMeans(x[showgenes,ii,drop=FALSE]))}))
dots[,mean_Z:=as.numeric(scale(mean)),by=gene];dots[,`:=`(cluster=factor(cluster,levels=CLUSTERS),gene=factor(gene,levels=rev(showgenes)))]
wt(dots,file.path(M1,'top20_marker_dotplot_data.tsv'))
pp<-ggplot(dots,aes(cluster,gene))+geom_point(aes(size=pct,colour=mean_Z))+scale_size_continuous(range=c(.3,5),name='Detected (%)')+scale_color_gradient2(low=COOL,mid='#F4F0EB',high=WARM,name='Mean expression\nrow Z-score')+theme_pub()+labs(title='C5 vs Rest All: top 20 markers',subtitle='Ranked by avg_log2FC among prespecified significant C5-up markers',x=NULL,y=NULL,caption='Frozen RNA LogNormalize data; all 20,146 cells. Cell-level inference is exploratory.')
export(pp,file.path(M1,'C5_top20_marker_dotplot'),8,8)
capture.output(sessionInfo(),file=file.path(OUT,'sessionInfo.txt'))
cat('Markers complete:',nrow(up),'significant up; ranked',nrow(rank),'symbols\n')

# Symmetric down selection from the complete unfiltered comparison; preserve Bonferroni rule.
down<-r[p_val_adj<.05&avg_log2FC<(-.75)&pmax(pct.1,pct.2)>=.25]
wt(down,file.path(M1,'markers_C5_down.tsv'))
all_de<-copy(r);all_de[,`:=`(P=p_val,adjusted_P=p_val_adj,FDR_BH=p.adjust(p_val,'BH'),adjustment_used='Bonferroni; Seurat total RNA feature count')]
wt(all_de,file.path(M1,'DEG_all.tsv'))
source('results/llps_C5_parallel_v1/scripts/04b_pseudobulk.R')
