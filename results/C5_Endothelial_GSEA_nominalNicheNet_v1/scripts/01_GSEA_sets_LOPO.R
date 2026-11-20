source('results/C5_Endothelial_GSEA_nominalNicheNet_v1/scripts/common.R')
library(fgsea);library(edgeR);library(limma)
manifest<-rd('00_input_audit/resource_verification.tsv');manifest[,current_sha256:=vapply(path,sha256,character(1))];stopifnot(all(manifest$current_sha256==manifest$sha256));w(manifest,'00_input_audit/resource_verification.tsv')
files<-c(file.path(OLD,'01_receiver_pseudobulk',paste0('DE_',models,'.tsv')),file.path(OLD,'01_receiver_pseudobulk/Endothelial_patient_pseudobulk.rds'),file.path(OLD,'02_receiver_geneset/expressed_background.tsv'),file.path(OLD,'04_NicheNet_ligand_activity/all_program_activities.tsv'),'results/llps_C0_functional_signature_v2/resources/gene_sets.rds','results/llps_C0_functional_signature_v2/resources/gene_set_sources.tsv')
current_hashes<-data.table(path=files,sha256=vapply(files,sha256,character(1)));hash_path<-file.path(O,'00_input_audit/input_hashes.tsv');if(file.exists(hash_path)){prior_hashes<-fread(hash_path);stopifnot(identical(current_hashes$path,prior_hashes$path),identical(current_hashes$sha256,prior_hashes$sha256))}else w(current_hashes,'00_input_audit/input_hashes.tsv')
w(rbindlist(lapply(names(sets),function(n)data.table(pathway=n,gene=sets[[n]]))),'00_input_audit/Hallmark_2026_frozen.tsv')
gall<-list();qa<-list()
for(m in models){
 d<-copy(DE[[m]]);qa[[m]]<-data.table(model=m,input_genes=nrow(d),nonfinite=sum(!is.finite(d$t)),duplicate_symbols=sum(duplicated(d$gene)),ENSG=sum(grepl('^ENSG',d$gene)),min_DE_FDR=min(d$adj.P.Val,na.rm=TRUE))
 stopifnot(!anyDuplicated(d$gene));d<-d[is.finite(t)];setorder(d,-t,gene);w(d[,.(gene,t,logFC,P.Value,adj.P.Val)],paste0('01_Endothelial_Hallmark_GSEA/ranking_',m,'.tsv'))
 path<-file.path(O,paste0('01_Endothelial_Hallmark_GSEA/GSEA_',m,'.tsv'))
 if(file.exists(path)){g<-fread(path)}else{
 set.seed(123);g<-fgseaMultilevel(sets,setNames(d$t,d$gene),minSize=10,maxSize=500,eps=0,nPermSimple=10000,BPPARAM=BiocParallel::SerialParam(RNGseed=123))
 g[,leadingEdge:=vapply(leadingEdge,paste,character(1),collapse=';')];g[,`:=`(model=m,dataset=if(m%in%c('InhouseData','GSE222315'))m else 'pooled',ranked_gene_count=nrow(d),resource_version='MSigDB 2026.1.Hs')];g[,padj:=p.adjust(pval,'BH')];w(g,paste0('01_Endothelial_Hallmark_GSEA/GSEA_',m,'.tsv'))
 };gall[[m]]<-g;cat('GSEA',m,nrow(g),'FDR10',sum(g$padj<.1),'\n');flush.console()
}
gall<-rbindlist(gall);w(gall,'01_Endothelial_Hallmark_GSEA/all_models.tsv');w(rbindlist(qa),'00_input_audit/DE_rank_QC.tsv')
wide<-dcast(gall,pathway~model,value.var=c('NES','padj'));wide[,both_major_positive:=NES_GSE222315>0 & NES_InhouseData>0];wide[,both_major_same_as_primary:=sign(NES_GSE222315)==sign(NES_primary)&sign(NES_InhouseData)==sign(NES_primary)];w(wide,'01_Endothelial_Hallmark_GSEA/model_comparison.tsv')
gate<-wide[pathway=='HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',NES_primary>0 & padj_primary<.1 & both_major_positive];writeLines(if(gate)'COHORT_SUPPORTED_EXPLORATORY' else 'REFERENCE_ONLY',file.path(O,'05_EMT_program_NicheNet/EMT_STATUS.txt'))
a<-copy(DE$primary);for(m in models[-1]){z<-DE[[m]][,.(gene,logFC,P.Value,adj.P.Val)];setnames(z,names(z)[-1],paste0(names(z)[-1],'_',m));a<-merge(a,z,by='gene',all.x=TRUE)}
a[,`:=`(in_background=gene%in%bg$gene,ECM_marker=gene%in%cfg$ECM_marker_exclusion)];a[,both_cohort_positive:=!is.na(logFC_GSE222315)&!is.na(logFC_InhouseData)&logFC_GSE222315>0&logFC_InhouseData>0]
a[,cohort_availability:=ifelse(is.na(logFC_GSE222315)|is.na(logFC_InhouseData),'EXPRESSION_FILTERED_NOT_TESTED','BOTH_TESTED')]
defs<-list();counts<-list()
for(n in names(cfg$nominal_p)){thr<-cfg$nominal_p[[n]];g<-selectgenes(DE$primary,DE$InhouseData,DE$GSE222315,thr);a[,(paste0('Set_',n)):=gene%in%g];defs[[n]]<-data.table(program=paste0('Set_',n),gene=g);defs[[paste0(n,'_noECM')]]<-data.table(program=paste0('Set_',n,'_noECM'),gene=setdiff(g,cfg$ECM_marker_exclusion));counts[[n]]<-data.table(set=n,p_threshold=thr,positive_primary=sum(a$logFC>0&a$P.Value<thr),positive_both_cohorts=sum(a$logFC>0&a$P.Value<thr&a$both_cohort_positive),matched=length(g),ECM=sum(g%in%cfg$ECM_marker_exclusion));w(a[gene%in%g],paste0('02_nominal_response_genesets/Set_',n,'_metrics.tsv'))}
defs[['EMT']]<-data.table(program='EMT',gene=intersect(sets[['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']],bg$gene));w(rbindlist(defs),'02_nominal_response_genesets/frozen_genesets.tsv');w(rbindlist(counts),'02_nominal_response_genesets/set_counts.tsv');w(a,'02_nominal_response_genesets/all_gene_audit.tsv');w(data.table(gene=cfg$ECM_marker_exclusion),'02_nominal_response_genesets/prespecified_ECM_markers.tsv')
# Refit only the small existing endothelial patient pseudobulk, not single-cell objects.
b<-readRDS(file.path(OLD,'01_receiver_pseudobulk/Endothelial_patient_pseudobulk.rds'));p<-as.data.table(b$metadata);x<-b$counts;stopifnot(identical(p$patient_id,colnames(x)))
fit<-function(sel){d<-copy(p[sel]);if(nrow(d)<6||sd(d$fraction)==0)return(NULL);d[,z:=as.numeric(scale(fraction))];d[,dataset:=factor(dataset)];design<-model.matrix(if(nlevels(d$dataset)>1)~dataset+z else ~z,d);if(qr(design)$rank<ncol(design)||nrow(design)-ncol(design)<3)return(NULL);y<-DGEList(as.matrix(x[,sel,drop=FALSE]));keep<-filterByExpr(y,design);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);v<-voom(y,design,plot=FALSE);f<-eBayes(lmFit(v,design),robust=TRUE);as.data.table(topTable(f,coef='z',n=Inf,sort.by='none'),keep.rownames='gene')}
# Reproduction check before LOPO.
f0<-fit(which(p$eligible));j<-merge(f0,DE$primary,by='gene');stopifnot(nrow(j)==nrow(DE$primary),max(abs(j$t.x-j$t.y))<1e-7);w(data.table(n=nrow(j),max_t_difference=max(abs(j$t.x-j$t.y))),'00_input_audit/model_reproduction.tsv')
loo<-list();lmetrics<-list();lst<-list();baseg<-defs[['B']]$gene
for(i in which(p$eligible)){
 id<-p$patient_id[i];fn<-paste0('04_rank_sensitivity/LOPO/',id,'_primary.tsv.gz');if(file.exists(file.path(O,fn))){d<-fread(file.path(O,fn));ih<-fread(file.path(O,paste0('04_rank_sensitivity/LOPO/',id,'_InhouseData.tsv.gz')));gs<-fread(file.path(O,paste0('04_rank_sensitivity/LOPO/',id,'_GSE222315.tsv.gz')))}else{
 sel<-setdiff(which(p$eligible),i);d<-fit(sel);ih<-fit(sel[p$dataset[sel]=='InhouseData']);gs<-fit(sel[p$dataset[sel]=='GSE222315']);stopifnot(!is.null(d),!is.null(ih),!is.null(gs));w(d,fn);w(ih,paste0('04_rank_sensitivity/LOPO/',id,'_InhouseData.tsv.gz'));w(gs,paste0('04_rank_sensitivity/LOPO/',id,'_GSE222315.tsv.gz'))}
 g<-selectgenes(d,ih,gs);loo[[id]]<-data.table(program=paste0('LOPO_',id),gene=g);lst[[id]]<-data.table(patient=id,dataset=p$dataset[i],n_genes=length(g),jaccard=length(intersect(g,baseg))/length(union(g,baseg)),status=if(length(g)>=15)'EVALUABLE' else 'NOT_EVALUABLE_TOO_FEW_GENES');lmetrics[[id]]<-d[gene%in%baseg,.(gene,omitted_patient=id,logFC,t,P.Value,adj.P.Val,selected=gene%in%g)];cat('LOPO',id,length(g),'\n');flush.console()
}
w(rbindlist(loo),'04_rank_sensitivity/LOPO_genesets.tsv');w(rbindlist(lst),'04_rank_sensitivity/LOPO_status.tsv');w(rbindlist(lmetrics),'04_rank_sensitivity/SetB_gene_LOPO_metrics.tsv');capture.output(sessionInfo(),file=file.path(O,'00_input_audit/sessionInfo_GSEA.txt'));cat('STAGE1_COMPLETE\n')
