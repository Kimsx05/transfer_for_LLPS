source('results/fig2_like_C0/scripts/common.R')
suppressPackageStartupMessages(library(edgeR))
status('G_DEG','C0 versus pooled C1-C6; raw-count pseudobulk')
pb <- readRDS(file.path(CORE,'C0_Others_pseudobulk_counts.rds'))
pm <- pb$metadata
eligible <- pm[,.(n_groups=.N,min_cells=min(n_cells)),by=sample]
eligible <- merge(eligible,unique(pm[,.(sample,dataset,tissue_origin,patient_uid,primary_replication_allowed)]),by='sample')
# Rules are declared before testing expression. A minimum of 20 cells in EACH arm
# keeps sample aggregates usable; 50-cell and tumor-only fits are sensitivity analyses.
eligible[,candidate:=n_groups==2 & min_cells>=20 & primary_replication_allowed==TRUE]
eligible[,patient_key:=fifelse(is.na(patient_uid),sample,patient_uid)]
sel <- eligible[candidate==TRUE][order(-min_cells,sample)][!duplicated(patient_key),sample]
eligible[,included:=sample%in%sel]
eligible[,reason:=fcase(included,'INCLUDED',!candidate,'FEWER_THAN_20_CELLS_PER_ARM_OR_ALTERNATE_SPECIMEN',
 default='OTHER_SPECIMEN_FROM_SAME_KNOWN_PATIENT')]
write_tsv(eligible,'pseudobulk_sample_eligibility.tsv')
writeLines(c('Primary: C0 versus pooled C1-C6, within the same sample.',
 'At least 20 cells in each arm; retain one specimen per known patient by largest minimum arm size, tie by sample name.',
 'Honor existing primary_replication_allowed flag for the user-confirmed T4463209a/b relationship.',
 'edgeR TMM, filterByExpr(group), robust estimateDisp and robust glmQLFit; design ~sample + group, Others reference.',
 'BH FDR < 0.05 AND log2FC > 0.25 defines C0 upregulated genes for BOTH GO and KEGG.',
 'Sample effects absorb sample-level dataset/tissue intercepts; pooled effects need not be identical in tumors and normal tissue.',
 'Sensitivity: >=50 cells per arm and tumor-only; same sample-independent gene filtering. These do not replace primary p values.',
 'Ranks = sign(logFC)*sqrt(QL F), all expressed/tested genes; no arbitrary top-N gene cutoff.'),
 file.path(OUT,'summary/DEG_rules.txt'))
fit_de <- function(samples,label) {
 meta <- copy(pm[sample%in%samples]);setorder(meta,sample,group)
 meta[,sample:=factor(sample)];meta[,group:=factor(group,levels=c('Others','C0'))]
 design <- model.matrix(~sample+group,meta)
 stopifnot(length(samples)>=4,qr(design)$rank==ncol(design),nrow(design)-ncol(design)>=3)
 y <- DGEList(as.matrix(pb$counts[,meta$library,drop=FALSE]))
 keep <- filterByExpr(y,group=meta$group)
 y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM')
 y <- estimateDisp(y,design,robust=TRUE)
 fit <- glmQLFit(y,design,robust=TRUE)
 q <- glmQLFTest(fit,coef='groupC0')
 res <- as.data.table(topTags(q,n=Inf,sort.by='none')$table,keep.rownames='gene')
 res[,rank_statistic:=sign(logFC)*sqrt(pmax(F,0))]
 res[,C0_up:=FDR<.05 & logFC>.25]
 write_tsv(data.table(library=meta$library,design),paste0('pseudobulk_design_',label,'.tsv'))
 write_tsv(data.table(library=meta$library,sample=meta$sample,group=meta$group,n_cells=meta$n_cells,
 library_size=y$samples$lib.size,TMM_factor=y$samples$norm.factors),paste0('pseudobulk_normalization_',label,'.tsv'))
 list(result=res,audit=data.table(model=label,n_samples=length(samples),n_libraries=nrow(meta),
 n_Normal=uniqueN(meta[tissue_origin=='Normal',sample]),n_Tumor=uniqueN(meta[tissue_origin=='Tumor',sample]),
 n_datasets=uniqueN(meta$dataset),n_genes=nrow(res),n_UP=sum(res$C0_up),rank=qr(design)$rank,
 residual_df=nrow(design)-ncol(design)))
}
primary <- fit_de(sel,'primary')
res <- primary$result;setorder(res,FDR)
write_tsv(res,'C0_vs_others_DEG.tsv');write_tsv(res[C0_up==TRUE],'C0_upregulated_genes.tsv')
write_tsv(res[order(-rank_statistic,gene),.(gene,rank_statistic,logFC,F,PValue,FDR)],'C0_preranked_gene_list.tsv')
audits <- list(primary$audit)
for(label in c('min50','tumor_only')) {
 ids <- if(label=='min50') eligible[included&min_cells>=50,sample] else eligible[included&tissue_origin=='Tumor',sample]
 if(length(ids)<4) stop('Insufficient samples for planned sensitivity: ',label)
 fit <- fit_de(ids,label);audits[[label]]<-fit$audit
 write_tsv(fit$result,paste0('C0_vs_others_DEG_',label,'_sensitivity.tsv'))
 idx <- match(res$gene,fit$result$gene)
 res[,(paste0(label,'_logFC')):=fit$result$logFC[idx]]
 res[,(paste0(label,'_FDR')):=fit$result$FDR[idx]]
}
write_tsv(res,'C0_DEG_sensitivity_comparison.tsv')
write_tsv(rbindlist(audits),'pseudobulk_model_summary.tsv')
capture.output(sessionInfo(),file=file.path(OUT,'summary/sessionInfo_DEG.txt'))
status('G_DEG',sprintf('COMPLETE: %d tested; %d upregulated',nrow(res),sum(res$C0_up)))
