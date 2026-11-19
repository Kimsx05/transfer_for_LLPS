options(mc.cores=1L);source('cellchat_nichechat_fibro7/scripts/nichenet_common.R');setDTthreads(1)
O<-'results/C5_Endothelial_NicheNet_spatial_v1';w<-function(x,p)fwrite(x,file.path(O,p),sep='\t',na='NA')
mp<-fread('results/llps_C5_distribution_audit_20261110_v1/v2/patient_sample_mapping_final.csv');man<-fread(file.path(ROOT,'01_metadata/per_sample_input_manifest.tsv'));lr<-readRDS(file.path(NN_OUT,'lr_network_human_21122021.rds'));print(names(lr));genes<-unique(c(lr$from,lr$to));rows<-list();bg<-list()
for(i in seq_len(nrow(man))){
 s<-man$sample[i];map<-mp[specimen_id==s];if(map$tissue_origin!='Tumor')next
 cat('EXPRESSION',s,'\n');flush.console();z<-qs2::qs_read(man$input_path[i],nthreads=1);x<-z$expression;m<-as.data.table(z$metadata);stopifnot(identical(colnames(x),m$cell))
 for(g in c(paste0('Fibro_C',0:6),'Endothelial')){
  ix<-which(m$cellchat_group==g);if(!length(ix))next
  gs<-intersect(genes,rownames(x));a<-x[gs,ix,drop=FALSE]
  rows[[length(rows)+1]]<-data.table(gene=gs,group=g,sample=s,patient_id=map$patient_id,dataset=map$dataset,n_cells=length(ix),n_detected=as.numeric(rowSums(a>0)),sum_log_normalized=as.numeric(rowSums(a)))
  if(g=='Endothelial'){a<-x[,ix,drop=FALSE];bg[[length(bg)+1]]<-data.table(gene=rownames(a),sample=s,patient_id=map$patient_id,dataset=map$dataset,n_cells=length(ix),n_detected=as.numeric(rowSums(a>0)))}
 };rm(z,x);gc(FALSE)
}
r<-rbindlist(rows);w(r,'03_candidate_ligands/expression_specimen.tsv.gz');r<-r[,.(n_cells=sum(n_cells),n_detected=sum(n_detected),sum_log_normalized=sum(sum_log_normalized)),by=.(gene,group,patient_id,dataset)];r[,`:=`(pct=n_detected/n_cells,mean_log_normalized=sum_log_normalized/n_cells)];w(r,'03_candidate_ligands/expression_patient.tsv.gz')
b<-rbindlist(bg);b<-b[,.(n_cells=sum(n_cells),n_detected=sum(n_detected)),by=.(gene,patient_id,dataset)];b[,pct:=n_detected/n_cells];w(b,'03_candidate_ligands/Endothelial_background_patient.tsv.gz');cat('EXPRESSION_COMPLETE\n')
