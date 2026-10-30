source('results/llps_C5_parallel_v1/scripts/common.R')
tab<-fread(file.path(M1,'DEG_all.tsv'))
tab[,status:=fcase(p_val_adj<.05&avg_log2FC>.75&pmax(pct.1,pct.2)>=.25,'C5-up',p_val_adj<.05&avg_log2FC<(-.75)&pmax(pct.1,pct.2)>=.25,'C5-down',default='Not significant')]
setorder(tab,p_val_adj,-avg_log2FC);vv<-copy(tab);vv[,neglog10FDR:=-log10(pmax(p_val_adj,.Machine$double.xmin))]
label<-unique(rbind(tab[gene%in%c('POSTN','COL11A1','COL10A1')],head(tab[status=='C5-up'][order(-avg_log2FC)],5),head(tab[status=='C5-down'][order(avg_log2FC)],5)),by='gene')
label[,neglog10FDR:=-log10(pmax(p_val_adj,.Machine$double.xmin))]
# Original C0 v1 volcano geometry; substitute v2 marker thresholds and labels.
src<-paste(readLines('results/llps_C0_functional_state_v1/scripts/04_DEG_LLPS.R'),collapse='\n')
a<-regexpr('pp<-ggplot(vv',src,fixed=TRUE)[1];b<-regexpr('export_plot(pp',src,fixed=TRUE)[1];code<-substr(src,a,b-1)
code<-gsub('Tumor-up','C5-up',code,fixed=TRUE);code<-gsub('Tumor-down','C5-down',code,fixed=TRUE);code<-gsub('logFC','avg_log2FC',code,fixed=TRUE)
code<-sub('c(-.25,.25)','c(-.75,.75)',code,fixed=TRUE)
code<-substr(code,1,regexpr(' labs(title=',code,fixed=TRUE)[1]-1)
eval(parse(text=paste0(code,"labs(title='C5 vs Rest All',subtitle=sprintf('%s up; %s down | Seurat RNA/data Wilcoxon',sum(tab$status=='C5-up'),sum(tab$status=='C5-down')),x='avg_log2FC (C5 / Rest)',y=expression(-log[10]('Bonferroni-adjusted P')),color=NULL,caption='max detection >=25%; |avg_log2FC| >0.75; Bonferroni-adjusted P <0.05.\\nCell-level inference is exploratory; patient-level pseudobulk validation reported separately.')")))
export(pp,file.path(M1,'C5_vs_Rest_volcano'),8.8,6.2)
wt(vv,file.path(M1,'volcano_data.tsv'));wt(label,file.path(M1,'volcano_labels.tsv'))
pb<-fread(file.path(M1,'pseudobulk/C5_vs_Rest_primary.tsv'));cmp<-merge(tab,pb,by='gene',suffixes=c('_cell','_pseudobulk'))
cmp[,same_direction:=sign(avg_log2FC)==sign(logFC)];wt(cmp,file.path(M1,'cell_vs_pseudobulk_comparison.tsv'))
wt(cmp[status=='C5-up',.(n_cell_up_tested=.N,n_positive_pb=sum(logFC>0),n_positive_significant_pb=sum(logFC>0&FDR<.05),fraction_positive_pb=mean(logFC>0))],file.path(M1,'pseudobulk_marker_validation_summary.tsv'))
wt(cmp[gene%in%c('POSTN','COL11A1','COL10A1')],file.path(M1,'prespecified_gene_pseudobulk_validation.tsv'))
# Keep prespecified genes even when the original pseudobulk expression filter excludes them.
focus<-merge(tab[gene%in%c('POSTN','COL11A1','COL10A1')],pb,by='gene',all.x=TRUE)
focus[,pseudobulk_status:=ifelse(is.finite(FDR),'TESTED','NOT_TESTED: excluded by original edgeR filterByExpr')]
wt(focus,file.path(M1,'prespecified_gene_pseudobulk_validation.tsv'))
