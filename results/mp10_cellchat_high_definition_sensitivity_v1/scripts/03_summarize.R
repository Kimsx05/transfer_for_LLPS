source('mp10_cellchat_high_definition_sensitivity_v1/scripts/common.R')
p<-if('--partial'%in%commandArgs(TRUE))rbindlist(lapply(list.files(tb('membership'),pattern='__status.tsv$',full.names=TRUE),fread),fill=TRUE)else fread(tb('network_status.tsv'));stopifnot(all(p$status%in%c('completed','completed_reused','NE')));dict<-fread(tb('final_label_dictionary.tsv'));tmes<-dict[parent!='Epithelial',celltype];domain<-fread(tb('fixed_LR_domain.tsv'));nc<-fread(tb('node_counts.tsv'));edges<-list();paths<-list();qc<-list();lk<-0L
for(i in seq_len(nrow(p))){z<-p[i];if(z$status=='NE')next;o<-qread(file.path(root,'objects/networks',paste0(z$network_id,'.qs')),nthreads=1);stopifnot(o$cache_key==z$cache_key);nd<-names(o$counts);lrs<-rownames(o$LR);stopifnot(identical(lrs,domain$interaction_id));sparse<-list();k<-0L
 for(di in c('output','input'))for(t in tmes){computed<-t%in%nd;nn<-if(computed)as.integer(o$counts[t])else 0L
  getvec<-function(st,slot){if(!computed)return(rep(NA_real_,length(lrs)));if(di=='output')o[[slot]][st,t,]else o[[slot]][t,st,]}
  hp<-getvec('Epi High','prob');lp<-getvec('Epi Low','prob');hP<-getvec('Epi High','pval');lP<-getvec('Epi Low','pval');hw<-hp;lw<-lp;hw[which(!(hp>0 & hP<.05))]<-0;lw[which(!(lp>0 & lP<.05))]<-0
  if(computed)for(st in c('High','Low')){pr<-if(st=='High')hp else lp;pv<-if(st=='High')hP else lP;sel<-which(pr!=0 | pv!=1);k<-k+1L;sparse[[k]]<-data.table(filter_stage='unfiltered',min_cells=NA_integer_,lr_domain='fixed_full',background_id='not_applicable_raw',network_id=z$network_id,cohort=z$cohort,scope=z$scope,high=z$high,direction=di,target=t,epi_state=st,sender=if(di=='output')paste('Epi',st)else t,receiver=if(di=='input')paste('Epi',st)else t,LR_order=sel,interaction_id=lrs[sel],pathway=o$LR$pathway_name[sel],ligand=o$LR$ligand[sel],receptor=o$LR$receptor[sel],probability=pr[sel],internal_P=pv[sel],supported=pr[sel]>0 & pv[sel]<.05,node_n=nn)}
  for(dom in c('fixed_full','legacy_1875_intersection')){ix<-if(dom=='fixed_full')seq_along(lrs)else which(domain$in_legacy_1875);for(th in c(10L,30L)){
   ok<-computed & nn>th;reason<-if(!computed)'NE_not_modeled_N_le_10'else if(!ok)'NE_postfilter_N_le_threshold'else'computed'
   lk<-lk+1L;edges[[lk]]<-data.table(network_id=z$network_id,cohort=z$cohort,scope=z$scope,high=z$high,min_cells=th,lr_domain=dom,direction=di,target=t,node_n=nn,computed=ok,status=reason,H_W=if(ok)sum(hw[ix])else NA_real_,L_W=if(ok)sum(lw[ix])else NA_real_,H_R=if(ok)sum(hp[ix])else NA_real_,L_R=if(ok)sum(lp[ix])else NA_real_,H_C=if(ok)sum(hp[ix]>0 & hP[ix]<.05)else NA_integer_,L_C=if(ok)sum(lp[ix]>0 & lP[ix]<.05)else NA_integer_)
   v<-data.table(pathway=o$LR$pathway_name[ix],H_W=hw[ix],L_W=lw[ix],H_R=hp[ix],L_R=lp[ix],H_C=as.integer(hp[ix]>0 & hP[ix]<.05),L_C=as.integer(lp[ix]>0 & lP[ix]<.05));v<-v[,lapply(.SD,sum),by=pathway];if(!ok)for(n in setdiff(names(v),'pathway'))set(v,j=n,value=NA_real_);v[,`:=`(network_id=z$network_id,cohort=z$cohort,scope=z$scope,high=z$high,min_cells=th,lr_domain=dom,direction=di,target=t,computed=ok,status=reason)];paths[[lk]]<-v
  }}
 }
 put(rbindlist(sparse),paste0('lr/',z$network_id,'__raw_sparse.tsv.gz'));qc[[i]]<-data.table(network_id=z$network_id,n_at_P_005=sum(o$pval==.05 & o$prob>0),n_prob_nonzero=sum(o$prob>0),all_prob_zero_P_one=all(o$pval[o$prob==0]==1));cat('Summarized',z$network_id,'\n')
}
ed<-rbindlist(edges);pa<-rbindlist(paths);put(ed,'epi_tme_edges.tsv.gz');put(pa,'epi_tme_pathways.tsv.gz');put(rbindlist(qc),'support_boundary_QC.tsv')
# Build explicit shared background universes. All three High definitions use identical node counts.
node_set<-function(co,sc,th)sort(nc[cohort==co & scope==sc & high=='H50' & N>th & node%in%tmes,node])
configs<-list();k<-0L
for(co in unique(p$cohort))for(sc in unique(p[cohort==co,scope]))for(th in c(10L,30L)){
 native<-node_set(co,sc,th);sets<-list(native=native,min10_min30_common=node_set(co,sc,30))
 if(sc=='ALL')sets$C10_C30_common<-intersect(node_set('C10','ALL',th),node_set('C30','ALL',th))
 if(co=='C10'){dsets<-setdiff(unique(p[cohort=='C10',scope]),'ALL');if(sc=='ALL')for(ds in dsets)sets[[paste0('ALL_dataset_common:',ds)]]<-intersect(native,node_set(co,ds,th)) else sets[[paste0('ALL_dataset_common:',sc)]]<-intersect(native,node_set(co,'ALL',th))}
 for(ct in names(sets)){ss<-sort(intersect(native,sets[[ct]]));k<-k+1L;configs[[k]]<-data.table(cohort=co,scope=sc,min_cells=th,background_id=paste0('bg_',substr(seqhash(ss),1,16)),background_context=ct,universe_members=paste(ss,collapse=';'),n_universe=length(ss))}
}
bg<-rbindlist(configs);put(bg,'background_comparison_map.tsv');bgu<-bg[,.(background_contexts=paste(background_context,collapse='|'),universe_members=universe_members[1],n_universe=n_universe[1]),by=.(cohort,scope,min_cells,background_id)]
effects<-list();k<-0L
for(ii in seq_len(nrow(bgu))){b<-bgu[ii];members<-strsplit(b$universe_members,';',fixed=TRUE)[[1]];v<-copy(ed[cohort==b$cohort & scope==b$scope & min_cells==b$min_cells]);v[,`:=`(background_id=b$background_id,background_contexts=b$background_contexts,universe_members=b$universe_members)];
 # Work group by group to exclude each target from the same universe.
 for(hi in c('H50','H10','H05'))for(dom in unique(v$lr_domain))for(di in c('output','input'))for(metric in c('W','R')){
  u<-copy(v[high==hi & lr_domain==dom & direction==di]);u[,`:=`(metric=metric,H=get(paste0('H_',metric)),L=get(paste0('L_',metric)),target_in_background_universe=target%in%members)]
  bh<-sum(u[target%in%members,H]);bl<-sum(u[target%in%members,L]);u[,`:=`(B_H=bh-ifelse(target_in_background_universe,H,0),B_L=bl-ifelse(target_in_background_universe,L,0),background_members=vapply(target,function(t)paste(setdiff(members,t),collapse=';'),character(1)),n_background=length(members)-as.integer(target_in_background_universe))]
  u[,`:=`(raw_difference=H-L,raw_ratio=ifelse(L>0,H/L,NA_real_),r_B=ifelse(B_L>0 & n_background>0,B_H/B_L,NA_real_))]
  u[,`:=`(RR=ifelse(L>0 & r_B>0,H/L/r_B,NA_real_),E=H-r_B*L,new_supported_contact=metric=='W' & computed & H>0 & L==0)]
  u[target_in_background_universe==FALSE,`:=`(RR=NA_real_,E=NA_real_,r_B=NA_real_)]
  u[,RR_status:=fcase(!computed,status,!target_in_background_universe,'NE_outside_comparison_universe',n_background==0,'undefined_no_background',B_L==0,'undefined_background_L_zero',B_H==0,'undefined_background_H_zero',L==0,'undefined_target_L_zero',default='defined')]
  u[,log2_RR:=ifelse(RR>0,log2(RR),NA_real_)];u[,relative_enhanced:=!is.na(RR)&RR>1];u[,rank_W_or_R:=frank(-RR,ties.method='min',na.last='keep')];k<-k+1L;effects[[k]]<-u
 }
}
ef<-rbindlist(effects);put(ef,'background_adjusted_effects.tsv.gz')
keys<-c('network_id','min_cells','lr_domain','direction','target')
# Pathway E uses overall target-specific r_B, not a separately normalized pathway background.
pe<-merge(pa,ef[,c(keys,'background_id','background_contexts','metric','r_B','RR','E','RR_status'),with=FALSE],by=keys,allow.cartesian=TRUE)
pe[,`:=`(H_path=ifelse(metric=='W',H_W,H_R),L_path=ifelse(metric=='W',L_W,L_R))];pe[,E_path:=H_path-r_B*L_path];put(pe,'pathway_background_contributions.tsv.gz')
chk<-pe[is.finite(E),.(E=sum(E_path)-E[1]),by=c(keys,'background_id','metric')];stopifnot(all(abs(chk$E)<1e-10));put(data.table(check='pathway_contribution_additivity',max_abs_error=max(abs(chk$E)),pass=TRUE),'pathway_additivity_QC.tsv')
# Main rankings and predeclared + ALL top3 union.
native<-ef[grepl('(^|\\|)native(\\||$)',background_contexts)];put(native,'native_effects_rankings.tsv.gz');top<-native[scope=='ALL' & min_cells==10 & lr_domain=='fixed_full' & relative_enhanced==TRUE & metric=='W'][order(rank_W_or_R),head(.SD,3),by=.(cohort,high,direction)]
put(top,'ALL_top3_positive.tsv');old6<-c('Naive/Memory B','Treg','cDC1','Plasma','Capillary EC','Arterial EC');focus<-unique(c(old6,top$target));put(data.table(target=focus,predeclared=focus%in%old6),'focus_candidate_union.tsv')
sel<-pe[cohort=='C10' & scope=='ALL' & high=='H05' & min_cells==10 & lr_domain=='fixed_full' & metric=='W' & grepl('(^|\\|)native(\\||$)',background_contexts) & target%in%focus & E_path>0][order(-E_path),head(.SD,5),by=.(target,direction)]
put(sel[,.(target,direction,pathway,selection_E_path=E_path,overall_RR=RR)],'focus_pathway_selection.tsv')
# Comparable background tables: min10 vs min30 must be numerically identical on common nodes.
cm<-merge(ef,bg[,.(cohort,scope,min_cells,background_id,comparison_context=background_context)],by=c('cohort','scope','min_cells','background_id'),allow.cartesian=TRUE)
put(cm,'comparable_background_effects.tsv.gz')
v<-cm[comparison_context=='min10_min30_common' & target_in_background_universe==TRUE];pair<-dcast(v,cohort+scope+high+lr_domain+direction+target+metric+background_id~min_cells,value.var=c('H','L','RR','E'));errs<-vapply(c('H','L','RR','E'),function(n){a<-pair[[paste0(n,'_10')]];b<-pair[[paste0(n,'_30')]];stopifnot(identical(is.na(a),is.na(b)));max(c(0,abs(a-b)),na.rm=TRUE)},numeric(1));stopifnot(all(errs<1e-12));put(data.table(check=names(errs),max_abs_error=errs,pass=TRUE),'min_filter_common_background_QC.tsv')
cat('SUMMARIZATION DONE\n')
