source('results/fig2_like_C0/phase1B_review/scripts/common.R')
suppressPackageStartupMessages({library(Seurat);library(qs)})
m<-fread(file.path(PHASE1,'C0_fraction_per_sample.tsv'))
# Pair mappings are explicit rows in Chen et al. Table 1, NOT inferred names.
review<-data.table(code=c('01','04','05','11','24','27','30'),
 tumor_accession=paste0('HRS',499831:499837),normal_accession=paste0('HRS',499838:499844))
review[,`:=`(tumor_sample=paste0('HRA003620_BC-',code),normal_sample=paste0('HRA003620_NC-',code),
 patient_uid=paste0('HRA003620::Chen_Table1_patient_',code),verified=TRUE,
 source='https://pmc.ncbi.nlm.nih.gov/articles/PMC9893503/#T1',
 evidence='Table 1 explicitly groups BC and NC specimens by patient; study methods confirm seven corresponding adjacent-normal specimens')]
WT(review,'HRA003620_pairing_review.tsv')
for(i in seq_len(nrow(review)))m[sample%in%c(review$tumor_sample[i],review$normal_sample[i]),patient_uid:=review$patient_uid[i]]
# Keep prior eligibility as explicitly historical audit fields.
setnames(m,c('known_TN_pair','main_independent_included','exclusion_reason'),
 c('phase1_known_TN_pair','phase1_main_independent_included','phase1_exclusion_reason'))
m[dataset=='HRA003620',patient_mapping_status:='VERIFIED_CHEN_TABLE1_PATIENT_PAIR']
m[,patient_key:=fifelse(is.na(patient_uid),sample,patient_uid)]
m[,tumor:=as.numeric(tissue_origin=='Tumor')]
m[,tissue_origin:=factor(tissue_origin,levels=c('Normal','Tumor'))]
m[,analysis_included:=primary_replication_allowed&n_all_Fibro_Pericyte>0]
WT(m,'sample_metadata_reviewed.tsv')
writeLines(c('# HRA003620 pairing review',
 'Chen et al., Frontiers in Oncology 2023, Table 1 explicitly places BC-01/NC-01, BC-04/NC-04, BC-05/NC-05, BC-11/NC-11, BC-24/NC-24, BC-27/NC-27 and BC-30/NC-30 under seven patients.',
 'Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC9893503/#T1',
 'Independent NGDC support: BC-01 run HRR854156 and NC-01 experiment HRX607517 both identify HRI318417; BC-30 HRS499837 and NC-30 HRS499844 both identify HRI318423.',
 'All 7 biological pairs verified. Two tumor specimens have no retained stromal cells, so only 5 HRA pairs have defined C0 proportions. Zero denominator remains NA.',
 'Earlier Phase 1 explicitly treated these identities as unknown. Its independent-sample P is superseded for inference by this review; historical figures are retained unchanged, with the refinement documented here.',
 'Study-patient keys group related specimens for all later validation splits. Names alone were not used as pairing evidence.'),file.path(CORE,'HRA003620_pairing_review.md'))
obj<-qread(OBJECT,nthreads=2)
d<-as.data.table(obj[[]],keep.rownames='cell')
frozen<-fread(file.path(PHASE1,'cell_metadata_umap_LLPSscore.tsv.gz'))
stopifnot(identical(d$cell,frozen$cell),identical(as.character(d$subcluster),frozen$cluster),
 max(abs(d$LLPS77_AUC-frozen$LLPS77_AUC))<1e-14)
q<-d[subcluster=='C0',.(n_C0=.N,median_LLPS=median(LLPS77_AUC),median_nFeature=median(nFeature_RNA),
 median_nCount=median(nCount_RNA),median_percent_mt=median(percent.mt)),by=.(sample=orig.ident)]
q<-merge(m,q,by='sample',all.x=TRUE,suffixes=c('','_check'))
stopifnot(all(q$n_C0[!is.na(q$n_C0_check)]==q$n_C0_check[!is.na(q$n_C0_check)]))
q[,`:=`(LLPS_primary=primary_replication_allowed&n_C0>=10,LLPS_anycapture=primary_replication_allowed&n_C0>=1,
 log10_features=log10(median_nFeature),log10_counts=log10(median_nCount))]
WT(q,'C0_LLPS_per_specimen.tsv')
composition<-d[,.N,by=.(sample=orig.ident,cluster=subcluster)]
composition<-merge(composition,m[,.(sample,dataset,tissue_origin,patient_uid,n_all_Fibro_Pericyte)],by='sample')
composition[,fraction:=N/n_all_Fibro_Pericyte]
WT(composition,'all_cluster_composition_per_specimen.tsv')
counts<-LayerData(obj,assay='RNA',layer='counts')
key<-paste(d$orig.ident,d$subcluster,sep='::');lev<-sort(unique(key))
agg<-sparseMatrix(i=seq_along(key),j=match(key,lev),x=1,dims=c(nrow(d),length(lev)))
pb<-counts%*%agg;colnames(pb)<-lev
meta<-d[,.(n_cells=.N),by=.(sample=orig.ident,cluster=subcluster)];meta[,library:=paste(sample,cluster,sep='::')]
meta<-merge(meta,m[,.(sample,dataset,tissue_origin,patient_uid,patient_key,primary_replication_allowed)],by='sample')
meta<-meta[match(lev,library)];stopifnot(identical(meta$library,colnames(pb)))
saveRDS(list(counts=pb,metadata=meta),file.path(CORE,'stromal_cluster_pseudobulk.rds'),compress='gzip')
identity<-fread(file.path(FIBROOT,'marker_discovery_review_v3/cluster_identity_and_state_review.csv'))
identity[,comparator_role:=fcase(cluster%in%c('C4','C5'),'PRIMARY_FIBRO_LIKE_COMPARATOR',cluster=='C0','TARGET',
 cluster=='C2','UNRESOLVED_EXCLUDED_PRIMARY_INCLUDED_SECONDARY',default='MURAL_EXCLUDED')]
WT(identity,'fibro_like_comparator_definition.tsv')
cat('Inputs, pairing, C0 sample medians and cluster aggregates complete\n')
