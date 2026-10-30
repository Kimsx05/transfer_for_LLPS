source('results/llps_C5_parallel_v1/scripts/common.R')
V2<-file.path(BASE,'results/llps_C0_functional_signature_v2');M2<-M3
suppressPackageStartupMessages({library(fgsea);library(patchwork)})
sets<-readRDS(file.path(V2,'resources/gene_sets.rds'))$H
full<-fread(file.path(OUT,'04_C5_vs_Rest_DEG/gene_ranking_all.tsv'))
r<-fread(file.path(OUT,'04_C5_vs_Rest_DEG/GSEA_preranked_genes.tsv'))
r<-merge(r,full[,.(gene,pct.1,pct.2)],by='gene',all.x=TRUE)
stopifnot(!anyNA(r$pct.1),!anyDuplicated(r$symbol),all(is.finite(r$avg_log2FC)))
r[,max_detection:=pmax(pct.1,pct.2)];setorder(r,-avg_log2FC,symbol)
results<-list();diag<-list();leads<-list()
for(threshold in c(.10,.05,.25)){
 tag<-paste0(round(threshold*100),'pct');rr<-r[max_detection>=threshold];rr[,rank_position:=seq_len(.N)]
 wt(rr,file.path(M2,paste0('preranked_',tag,'.tsv')));ranks<-setNames(rr$avg_log2FC,rr$symbol)
 result_path<-file.path(M2,paste0(if(threshold==.10)'GSEA_primary_' else 'GSEA_sensitivity_',tag,'.tsv'))
 if(file.exists(result_path)){
  g<-fread(result_path)[is.finite(pval)]
  g[,leadingEdge:=strsplit(leadingEdge,';',fixed=TRUE)]
  for(extra in intersect(c('matched_genes','analysis_status'),names(g)))g[,(extra):=NULL]
 }else{
  set.seed(123)
  g<-fgseaMultilevel(pathways=sets,stats=ranks,minSize=10,maxSize=500,eps=0,nPermSimple=10000,BPPARAM=BiocParallel::SerialParam(RNGseed=123))
 }
 eligible_n<-vapply(sets,function(gs)sum(names(ranks)%in%gs),integer(1))
 stopifnot(setequal(g$pathway,names(sets)[eligible_n>=10&eligible_n<=500]),all(is.finite(g$pval)))
 leads[[tag]]<-rbindlist(lapply(seq_len(nrow(g)),function(i)data.table(threshold=threshold,pathway=g$pathway[i],gene=g$leadingEdge[[i]],NES=g$NES[i],FDR=g$padj[i])))
 g[,leadingEdge:=vapply(leadingEdge,paste,character(1),collapse=';')]
 g<-merge(data.table(pathway=names(sets),matched_genes=eligible_n),g,by='pathway',all.x=TRUE)
 g[,analysis_status:=ifelse(is.finite(pval),'TESTED','NOT_TESTED_GENESET_SIZE_OUTSIDE_10_500')]
 g[,detection_threshold:=threshold];setorder(g,-NES)
 wt(g,file.path(M2,paste0(if(threshold==.10)'GSEA_primary_' else 'GSEA_sensitivity_',tag,'.tsv')));results[[tag]]<-g
 diag[[tag]]<-data.table(threshold=threshold,n_ranks=length(ranks),n_pathways_tested=sum(is.finite(g$pval)),nonfinite=0,duplicate_symbols=0,n_rows_in_ties=sum(duplicated(ranks)|duplicated(ranks,fromLast=TRUE)),min_logFC=min(ranks),max_logFC=max(ranks),head100_median_detection=median(head(rr$max_detection,100)),n_positive=sum(g$padj<.05&g$NES>0,na.rm=TRUE),n_negative=sum(g$padj<.05&g$NES<0,na.rm=TRUE))
 cat(tag,':',nrow(rr),'genes;',sum(g$padj<.05&g$NES>0,na.rm=TRUE),'positive;',sum(g$padj<.05&g$NES<0,na.rm=TRUE),'negative\n');flush.console()
}
wt(rbindlist(diag),file.path(M2,'filter_and_result_summary.tsv'));wt(rbindlist(leads),file.path(M2,'leading_edge.tsv'))
long<-rbindlist(results,use.names=TRUE,fill=TRUE);wt(long,file.path(M2,'GSEA_all_thresholds.tsv'))
wide<-dcast(long,pathway~detection_threshold,value.var=c('ES','NES','padj'))
wide[,`:=`(same_direction_all_thresholds=sign(ES_0.05)==sign(ES_0.1)&sign(ES_0.1)==sign(ES_0.25),significant_all_thresholds=padj_0.05<.05&padj_0.1<.05&padj_0.25<.05)]
wt(wide,file.path(M2,'threshold_direction_comparison.tsv'))
g<-copy(results[['10pct']]);rr<-r[max_detection>=.10];ranks<-setNames(rr$avg_log2FC,rr$symbol)
# Plotting below reuses the verified v2 running-ES construction, without rerunning DE.
# Summary retains all tested Hallmarks, separating statistically significant bars.
g[,significance:=ifelse(padj<.05,'FDR < 0.05','Not significant')];g[,label:=factor(pathway,levels=rev(pathway))]
pp<-ggplot(g,aes(NES,label,fill=NES>0,alpha=significance))+geom_col(width=.75)+
 scale_fill_manual(values=c('TRUE'=WARM,'FALSE'=COOL),guide='none')+scale_alpha_manual(values=c('FDR < 0.05'=1,'Not significant'=.3),name=NULL)+
 scale_y_discrete(labels=setNames(label_term(g$pathway,52),g$pathway))+theme_pub()+
 labs(title='C5 vs Rest All: Hallmark GSEA',subtitle='Primary: detection ≥10% in either group; avg_log2FC ranking',x='Normalized enrichment score (NES)',y=NULL,
 caption='Gene-set enrichment of a cell-level contrast; not independent-patient validation.')
export(pp,file.path(M2,'figures/Hallmark_GSEA_summary'),10.5,15)
# Prespecified objective selection: four largest positive NES and two smallest
# negative NES, conditional on BH FDR <0.05; never choose by mechanism preference.
chosen<-rbind(head(g[padj<.05&NES>0][order(-NES)],4),head(g[padj<.05&NES<0][order(NES)],2))
wt(chosen,file.path(M2,'displayed_GSEA_pathways.tsv'))
for(i in seq_len(nrow(chosen))){
 row<-chosen[i];hits<-which(names(ranks)%in%sets[[row$pathway]])
 hit<-seq_along(ranks)%in%hits;weights<-abs(unname(ranks));increment<-ifelse(hit,weights/sum(weights[hit]),-1/sum(!hit))
 curve<-data.table(rank=0:length(ranks),running_ES=c(0,cumsum(increment)))
 observed_ES<-if(max(curve$running_ES)>-min(curve$running_ES))max(curve$running_ES) else min(curve$running_ES)
 stopifnot(abs(observed_ES-row$ES)<1e-8)
 wt(data.table(rank=seq_along(ranks),gene=names(ranks),ranking=unname(ranks),hit=hit,running_ES=curve$running_ES[-1]),file.path(M2,paste0(sub('HALLMARK_','',row$pathway),'_curve_data.tsv.gz')))
 p1<-ggplot(curve,aes(rank,running_ES))+geom_line(colour=if(row$NES>0)WARM else COOL,linewidth=.8)+scale_x_continuous(limits=c(0,length(ranks)),expand=expansion(mult=0))+theme_pub()+labs(title=gsub('_',' ',sub('HALLMARK_','',row$pathway)),subtitle=paste0('NES = ',sprintf('%.2f',row$NES),' | ES = ',sprintf('%.2f',row$ES),' | BH FDR = ',ptxt(row$padj)),x=NULL,y='Enrichment score')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
 p2<-ggplot(data.table(rank=hits),aes(x=rank))+geom_segment(aes(xend=rank,y=0,yend=1),linewidth=.25)+scale_x_continuous(limits=c(0,length(ranks)),expand=expansion(mult=0))+theme_void()+theme(plot.background=element_rect(fill='transparent',colour=NA),plot.margin=margin(0,18,0,10))+labs(y='Hits')
 metric<-data.table(rank=seq_along(ranks),metric=unname(ranks))
 p3<-ggplot(metric,aes(rank,metric))+geom_area(fill='#A2A2A2')+scale_x_continuous(limits=c(0,length(ranks)),expand=expansion(mult=0))+theme_pub()+labs(x='Rank in expression-filtered gene list',y='avg_log2FC',caption='C5-enriched genes at left; Rest-enriched genes at right.')
 pp<-(p1/p2/p3)+plot_layout(heights=c(3,.6,1.5));pp<-pp&theme(plot.background=element_rect(fill='transparent',colour=NA))
 export(pp,file.path(M2,'figures',paste0('GSEA_',sub('HALLMARK_','',row$pathway))),8,6.3)
}
cat('GSEA significant:',sum(g$padj<.05&g$NES>0,na.rm=TRUE),'positive;',sum(g$padj<.05&g$NES<0,na.rm=TRUE),'negative\n')
long[,threshold_label:=factor(paste0(round(detection_threshold*100),'%'),levels=c('5%','10%','25%'))]
long[,label:=factor(pathway,levels=rev(g$pathway))];long[,star:=ifelse(!is.finite(padj),'NA',ifelse(padj<.05,'*',''))]
p<-ggplot(long,aes(threshold_label,label,fill=NES))+geom_tile(colour='white',linewidth=.2)+geom_text(aes(label=star),size=4)+scale_y_discrete(labels=setNames(label_term(g$pathway,52),g$pathway))+scale_fill_gradient2(low=COOL,mid='#FAF7F2',high=WARM,midpoint=0,name='NES',na.value='#D8D8D8')+theme_pub()+theme(axis.line=element_blank(),axis.ticks=element_blank())+labs(title='Hallmark GSEA: detection-threshold sensitivity',subtitle='10% primary; * BH FDR < 0.05; NA: fewer than 10 matched genes',x='Minimum detection in either C5 or Rest',y=NULL)
export(p,file.path(M2,'figures/NES_threshold_comparison'),10,15)
