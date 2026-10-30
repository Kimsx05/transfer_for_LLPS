source('results/llps_C5_parallel_v1/scripts/common.R')
r<-fread(file.path(M1,'gene_ranking_all.tsv'));up<-fread(file.path(M1,'markers_C5_up.tsv'));down<-fread(file.path(M1,'markers_C5_down.tsv'))
stopifnot(!anyDuplicated(r$gene),setequal(up$gene,r[p_val_adj<.05&avg_log2FC>.75&pmax(pct.1,pct.2)>=.25,gene]),setequal(down$gene,r[p_val_adj<.05&avg_log2FC<(-.75)&pmax(pct.1,pct.2)>=.25,gene]))
m<-fread(file.path(OUT,'01_Normal_Tumor/C5_fraction_per_sample_reviewed.tsv'))
stopifnot(sum(m$n_C5)==716,sum(m$n_all_Fibro_Pericyte)==20146,max(abs(m[n_all_Fibro_Pericyte>0,C5_proportion-n_C5/n_all_Fibro_Pericyte]))<1e-12,all(is.na(m[n_all_Fibro_Pericyte==0,C5_proportion])))
for(db in c('GO_BP','KEGG')){
 all<-fread(file.path(M2,paste0(db,'_all.tsv')));s<-fread(file.path(M2,paste0('selected_',db,'.tsv')));ix<-match(s$ID,all$ID)
 stopifnot(!anyNA(ix),all(s$p.adjust==all$p.adjust[ix]),all(s$Count==all$Count[ix]),all(s$GeneRatio==all$GeneRatio[ix]),all(s$p.adjust<.05))
}
for(tag in c('5','10','25')){
 rr<-fread(file.path(M3,paste0('preranked_',tag,'pct.tsv')))
 g<-fread(file.path(M3,paste0(if(tag=='10')'GSEA_primary_'else'GSEA_sensitivity_',tag,'pct.tsv')))
 stopifnot(!anyDuplicated(rr$symbol),all(diff(rr$avg_log2FC)<=0),all(pmax(rr$pct.1,rr$pct.2)>=as.numeric(tag)/100),nrow(g)==50)
 test<-g[is.finite(pval)];stopifnot(max(abs(test$padj-p.adjust(test$pval,'BH')))<1e-12)
}
cx<-fread(file.path(OUT,'03_Bulk_Survival/C0_C5_Cox_comparison.tsv'))
old<-fread(file.path(BASE,'results/cibersortx_phase2_clinical_v1/all_cox_results.tsv'))[endpoint=='OS']
z<-merge(cx[cluster=='C0'],old,by=c('dataset','model'),suffixes=c('_new','_old'))
stopifnot(nrow(z)==8,max(abs(z$HR_new-z$HR_old))<1e-12,all(z$N_new==z$N_old))
counts<-m[analysis_included==TRUE,.(N=.N,n_zero=sum(n_C5==0),median_C5=median(C5_proportion)),by=tissue_origin]
wt(counts,file.path(OUT,'01_Normal_Tumor/all_reviewed_tissue_summary.tsv'))
writeLines(c('PASS: C5=716; frozen Fibro/Pericyte=20146; zero samples preserved.',paste('PASS: full DEG rows',nrow(r),'up',nrow(up),'down',nrow(down)),
 'PASS: displayed ORA IDs/FDR/Count/GeneRatio exactly match complete tables.',
 'PASS: GSEA new C5 rankings, unique symbols, correct thresholds and within-threshold BH.',
 'PASS: all eight reused C0 Full14 Cox results match previous results (HR and N).'),file.path(OUT,'validation.txt'))
