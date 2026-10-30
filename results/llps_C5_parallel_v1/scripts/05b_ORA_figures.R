source('results/llps_C5_parallel_v1/scripts/common.R')
GO_IDS<-c('GO:0030198','GO:0032963','GO:0007155','GO:0001944','GO:0001525','GO:0030574','GO:0032964','GO:0061448','GO:0030036','GO:0030111')
# Reuse v3 plotting expression; replace the fixed C0 IDs with verified significant C5 IDs.
src<-paste(readLines('results/llps_C0_functional_signature_v3/scripts/01_curated_ORA.R'),collapse='\n')
start<-regexpr(' a[,label:=factor',src,fixed=TRUE)[1];end<-regexpr(" export(pp,",src,fixed=TRUE)[1]
plot_code<-substr(src,start,end-1)
plot_code<-gsub('C0','C5',plot_code,fixed=TRUE)
plot_code<-sub('Biologically curated selection from significant marker-based ORA','Representative significant marker-based ORA',plot_code,fixed=TRUE)
plot_code<-sub('C5 vs Rest All (C1–C6). Statistics reused unchanged from v2; BH FDR < 0.05.','C5 vs C0+C1+C2+C3+C4+C6; new ORA using the original v2 algorithm; BH FDR < 0.05.',plot_code,fixed=TRUE)
for(db in c('GO_BP','KEGG')){
 orig<-fread(file.path(M2,paste0(db,'_all.tsv')));sig<-fread(file.path(M2,paste0(db,'_significant.tsv')))
 a<-if(db=='GO_BP')orig[match(GO_IDS,ID)]else head(sig[order(p.adjust,pvalue,ID)],10)
 stopifnot(!anyNA(a$ID),all(a$ID%in%sig$ID),all(a$p.adjust<.05))
 a[,`:=`(display_order=seq_len(.N),FDR=p.adjust,selection=if(db=='GO_BP')'Significant ECM/collagen, adhesion, vascular, connective tissue, cytoskeletal and Wnt representatives'else'Top significant pathways; fewer than ten available')]
 stopifnot(all(a$Count==orig$Count[match(a$ID,orig$ID)]),identical(a$GeneRatio,orig$GeneRatio[match(a$ID,orig$ID)]))
 wt(a,file.path(M2,paste0('selected_',db,'.tsv')))
 eval(parse(text=plot_code));export(pp,file.path(M2,paste0('C5_',db,'_barplot')),10,7.1)
}
# Requested themes checked in the full result, including nonsignificant terms.
go<-fread(file.path(M2,'GO_BP_all.tsv'))
focus<-go[grepl('migration|angiogen|vascula|extracellular matrix|cell adhesion|inflamm|immune|lipid|cholesterol',Description,ignore.case=TRUE)]
wt(focus,file.path(M2,'requested_theme_review_all.tsv'))
