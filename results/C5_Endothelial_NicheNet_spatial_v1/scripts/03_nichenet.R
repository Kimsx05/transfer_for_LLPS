options(mc.cores=1L);source('cellchat_nichechat_fibro7/scripts/nichenet_common.R');setDTthreads(1);library(nichenetr)
O<-'results/C5_Endothelial_NicheNet_spatial_v1';w<-function(x,p)fwrite(x,file.path(O,p),sep='\t',na='NA')
r<-fread(file.path(O,'03_candidate_ligands/expression_patient.tsv.gz'));b<-fread(file.path(O,'03_candidate_ligands/Endothelial_background_patient.tsv.gz'))
sumexpr<-function(r,threshold=.1){r[,.(pooled_pct=sum(n_detected)/sum(n_cells),mean_log_expression=sum(sum_log_normalized)/sum(n_cells),patients_any=sum(n_detected>0),support_patients=sum(pct>=threshold & n_cells>=ifelse(group[1]=='Endothelial',20,5)),support_datasets=uniqueN(dataset[pct>=threshold & n_cells>=ifelse(group[1]=='Endothelial',20,5)]),n_patients=.N),by=.(gene,group)]}
su<-sumexpr(r);w(su,'03_candidate_ligands/expression_summary.tsv');lr<-as.data.table(readRDS(file.path(NN_OUT,'lr_network_human_21122021.rds')));lr<-unique(lr[,.(ligand=from,receptor=to)])
ls<-su[group=='Fibro_C5'];rs<-su[group=='Endothelial'];pairs<-merge(lr,ls,by.x='ligand',by.y='gene');pairs<-merge(pairs,rs,by.x='receptor',by.y='gene',suffixes=c('_C5','_Endo'));pairs[,eligible:=pooled_pct_C5>=.1 & pooled_pct_Endo>=.1 & support_patients_C5>=2 & support_patients_Endo>=2];w(pairs,'03_candidate_ligands/broad_LR_all.tsv')
lt<-readRDS(file.path(NN_OUT,'ligand_target_matrix_nsga2r_final.rds'));ligs<-intersect(unique(pairs[eligible==TRUE,ligand]),colnames(lt));bg<-b[n_cells>=20,.(pooled_pct=sum(n_detected)/sum(n_cells),support=sum(pct>=.1)),by=gene];bg<-bg[pooled_pct>=.1 & support>=2 & gene%in%rownames(lt)];w(bg,'02_receiver_geneset/expressed_background.tsv')
# Common-patient expression contrasts; equal patient weight, no cell-level inferential p-values.
sp<-list();for(other in c('Fibro_C0','Fibro_C1','Fibro_C4','REST')){
 a<-r[group=='Fibro_C5' & n_cells>=5];q<-if(other=='REST')r[group%in%paste0('Fibro_C',c(0:4,6)),.(n_cells=sum(n_cells),n_detected=sum(n_detected),sum_log_normalized=sum(sum_log_normalized)),by=.(gene,patient_id,dataset)] else r[group==other];q<-q[n_cells>=5];q[,mean_other:=sum_log_normalized/n_cells];j<-merge(a,q[,.(gene,patient_id,mean_other)],by=c('gene','patient_id'));j[,difference:=mean_log_normalized-mean_other];j[,comparison:=other];sp[[other]]<-j[,.(gene,patient_id,dataset,comparison,difference,mean_C5=mean_log_normalized,mean_other)]
};sp<-rbindlist(sp);w(sp,'03_candidate_ligands/common_patient_specificity.tsv.gz');w(sp[,.(n_patients=.N,median_difference=median(difference),positive_fraction=mean(difference>0)),by=.(gene,comparison)],'03_candidate_ligands/specificity_summary.tsv')
mode<-readLines(file.path(O,'02_receiver_geneset/RESPONSE_MODE.txt'));programs<-list()
if(mode=='PATIENT_ASSOCIATED')programs[['PATIENT_ASSOCIATED']]<-fread(file.path(O,'02_receiver_geneset/patient_response_genes.tsv'))$gene else {
 gmt<-'results/fig2_like_C0/resources/h.all.v2025.1.Hs.symbols.gmt';lines<-strsplit(readLines(gmt),'\t');for(z in lines)if(z[1]%in%c('HALLMARK_ANGIOGENESIS','HALLMARK_TGF_BETA_SIGNALING'))programs[[z[1]]]<-z[-c(1,2)]
 w(data.table(path=gmt,sha256=sha256(gmt),version='MSigDB human 2025.1',interpretation='curated functional programs, not observed C5-induced DEGs'),'02_receiver_geneset/program_provenance.tsv')
}
cc<-rbindlist(lapply(c(10,20,30),function(t)fread(paste0('results/llps_C5_cellchat_v2/04_patient_level_summary/C5_LR_summary_t',t,'.csv.gz'))));cc<-cc[stratification=='tumor_only' & sender=='Fibro_C5' & receiver=='Endothelial'];w(cc,'06_CellChat_integration/tumor_LR_all_thresholds.tsv');ccmain<-cc[min_cells==20 & positive_patients>0]
all<-list();links<-list();defs<-list();random<-list()
for(nm in names(programs)){
 g<-intersect(programs[[nm]],bg$gene);defs[[nm]]<-data.table(program=nm,gene=programs[[nm]],in_background=programs[[nm]]%in%bg$gene)
 if(length(g)<15){cat('PROGRAM_NOT_EVALUABLE',nm,length(g),'\n');next}
 cat('ACTIVITY',nm,length(g),length(ligs),'\n');flush.console()
 a<-as.data.table(predict_ligand_activities(g,bg$gene,lt,ligs));setorder(a,-aupr_corrected,-pearson);a[,`:=`(rank=seq_len(.N),program=nm,geneset_n=length(g),background_n=nrow(bg),CellChat_supported=test_ligand%in%ccmain$ligand)];all[[nm]]<-a
 for(l in head(a$test_ligand,10)){z<-as.data.table(get_weighted_ligand_target_links(l,g,lt,n=250));z[,program:=nm];links[[length(links)+1]]<-z}
 # 20 size/detection-bin matched random programs: calibration diagnostic, not ligand p-values.
 cuts<-unique(quantile(bg$pooled_pct,seq(0,1,.1)));bins<-cut(bg$pooled_pct,cuts,include.lowest=TRUE);names(bins)<-bg$gene
 for(k in 1:20){gg<-unlist(lapply(split(g,bins[g]),function(gs){pool<-names(bins)[bins==bins[gs[1]]];sample(pool,length(gs))}));ra<-as.data.table(predict_ligand_activities(gg,bg$gene,lt,ligs));ra[,`:=`(program=nm,replicate=k)];random[[length(random)+1]]<-ra}
}
a<-rbindlist(all);w(a,'04_NicheNet_ligand_activity/all_program_activities.tsv');w(rbindlist(defs),'02_receiver_geneset/frozen_program_genes.tsv');w(rbindlist(links,fill=TRUE),'05_ligand_target_links/top10_links.tsv');w(rbindlist(random),'04_NicheNet_ligand_activity/matched_random_programs.tsv.gz')
a[,normalized_rank:=rank/.N,by=program];ranked<-a[,.(mean_normalized_rank=mean(normalized_rank),best_rank=min(rank),programs=paste(program,collapse=';')),by=test_ligand];setorder(ranked,mean_normalized_rank,best_rank);w(ranked,'04_NicheNet_ligand_activity/aggregate_priority.tsv')
selected<-head(ranked[test_ligand%in%ccmain$ligand,test_ligand],5);sel<-ccmain[ligand%in%c(selected,'FN1','COL1A1','COL1A2')];sel[,role:=ifelse(ligand%in%selected,'candidate','ECM_control')];setorder(sel,ligand,-patient_condition_median);sel<-sel[,head(.SD,1),by=ligand];w(sel,'07_spatial_expression_validation/frozen_axes.tsv')
integ<-merge(a,ccmain,by.x='test_ligand',by.y='ligand',all=TRUE,allow.cartesian=TRUE);integ[,evidence_class:=fifelse(is.na(rank),'CellChat_only',fifelse(is.na(feature),'NicheNet_only','both'))];w(integ,'06_CellChat_integration/full_evidence.tsv')
w(a[test_ligand%in%c('ADM','MDK','SEMA3C','WNT2','FN1','COL1A1','COL1A2')],'04_NicheNet_ligand_activity/requested_ligand_ranks.tsv')
cat('NICHENET_COMPLETE',mode,length(ligs),nrow(bg),'\n')
