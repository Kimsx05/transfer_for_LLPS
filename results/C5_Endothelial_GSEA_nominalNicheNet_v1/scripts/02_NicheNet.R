source('results/C5_Endothelial_GSEA_nominalNicheNet_v1/scripts/common.R');library(nichenetr)
lt<-readRDS(file.path(NN_OUT,'ligand_target_matrix_nsga2r_final.rds'));stopifnot(all(bg$gene%in%rownames(lt)),all(ligs%in%colnames(lt)))
defs<-fread(file.path(O,'02_nominal_response_genesets/frozen_genesets.tsv'));loo<-fread(file.path(O,'04_rank_sensitivity/LOPO_genesets.tsv'));defs<-rbind(defs,loo)
cc<-rd('06_CellChat_integration/tumor_LR_all_thresholds.tsv');support<-unique(cc[min_cells==20 & positive_patients>0,ligand]);all<-list();links<-list();status<-list()
run<-function(g,nm){
 fn<-file.path(O,'03_nominal_NicheNet/checkpoints',paste0(nm,'.tsv'))
 if(file.exists(fn))return(fread(fn))
 set.seed(123);z<-as.data.table(predict_ligand_activities(g,bg$gene,lt,ligs));setorder(z,-aupr_corrected,-pearson,test_ligand);z[,`:=`(rank=seq_len(.N),program=nm,geneset_n=length(g),background_n=nrow(bg),CellChat_supported=test_ligand%in%support)];stopifnot(nrow(z)==246,all(is.finite(z$aupr_corrected)));fwrite(z,fn,sep='\t');cat('ACTIVITY_COMPLETE',nm,length(g),'\n');flush.console();z
}
for(nm in unique(defs$program)){
 g<-defs[program==nm,gene];status[[nm]]<-data.table(program=nm,n_genes=length(g),status=if(length(g)>=15)'EVALUABLE' else 'NOT_EVALUABLE_TOO_FEW_GENES');if(length(g)<15)next
 a<-run(g,nm);all[[nm]]<-a
 if(!grepl('^LOPO_',nm)){for(l in head(a$test_ligand,10)){z<-as.data.table(get_weighted_ligand_target_links(l,g,lt,n=250));if(nrow(z)){z[,program:=nm];links[[length(links)+1]]<-z}}}
 w(rbindlist(all),'03_nominal_NicheNet/all_activities.tsv');w(rbindlist(links,fill=TRUE),'03_nominal_NicheNet/top10_ligand_target_links.tsv')
}
w(rbindlist(status),'03_nominal_NicheNet/program_status.tsv')
# Fixed same-size/detection-bin matched random background programs, 50 draws.
g<-defs[program=='Set_B',gene];cuts<-unique(quantile(bg$pooled_pct,seq(0,1,.1)));bins<-cut(bg$pooled_pct,cuts,include.lowest=TRUE);names(bins)<-bg$gene
set.seed(123);randomsets<-lapply(1:50,function(k){unlist(lapply(split(g,bins[g],drop=TRUE),function(gs){pool<-names(bins)[bins==bins[gs[1]]];sample(pool,length(gs),replace=FALSE)}),use.names=FALSE)})
w(rbindlist(lapply(seq_along(randomsets),function(k){data.table(program=sprintf('RANDOM_%02d',k),gene=randomsets[[k]],detection_bin=as.character(bins[randomsets[[k]]]))})), '04_rank_sensitivity/random_genesets.tsv')
ra<-list();for(k in seq_along(randomsets)){ra[[k]]<-run(randomsets[[k]],sprintf('RANDOM_%02d',k));w(rbindlist(ra),'04_rank_sensitivity/random_activities.tsv.gz')}
capture.output(sessionInfo(),file=file.path(O,'00_input_audit/sessionInfo_NicheNet.txt'));cat('NICHENET_COMPLETE\n')
