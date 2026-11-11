# Source common.R first. Restores only computed cells in the frozen LR domain.
recover_lr<-function(id,direction,target,epi_state='High',min_cells=10L,lr_domain='fixed_full'){
 stopifnot(direction%in%c('input','output'),epi_state%in%c('High','Low'),min_cells%in%c(10L,30L),lr_domain%in%c('fixed_full','legacy_1875_intersection'))
 dom<-fread(tb('fixed_LR_domain.tsv'));if(lr_domain=='legacy_1875_intersection')dom<-dom[in_legacy_1875==TRUE]
 nc<-fread(tb('node_counts.tsv'));nn<-nc[network_id==id & node==target,N];stopifnot(length(nn)==1)
 dom[,`:=`(network_id=id,min_cells=min_cells,lr_domain=lr_domain,direction=direction,target=target,epi_state=epi_state,computed=nn>min_cells,probability=NA_real_,internal_P=NA_real_)]
 if(nn>min_cells){sp<-fread(tb(paste0('lr/',id,'__raw_sparse.tsv.gz')));di<-direction;ta<-target;st<-epi_state;sp<-sp[direction==di & target==ta & epi_state==st];dom[,`:=`(probability=0,internal_P=1)];ix<-match(dom$interaction_id,sp$interaction_id);sel<-which(!is.na(ix));dom[sel,`:=`(probability=sp$probability[ix[sel]],internal_P=sp$internal_P[ix[sel]])]}
 dom[,supported:=ifelse(computed,probability>0 & internal_P<.05,NA)];dom
}
