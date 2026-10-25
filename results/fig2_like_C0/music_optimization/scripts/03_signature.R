source('results/fig2_like_C0/music_optimization/scripts/common.R');CORE<-file.path(OPT,'06_signature_fallback')
uni<-readLines(file.path(OPT,'02_model_comparison/fixed_shared_gene_universe.txt'))
# Gene IDs alone are read; no external stage or outcome is consulted.
con<-gzfile(file.path(P2,'resources/TCGA_BLCA_HiSeqV2.gz'),'rt');tc<-character();readLines(con,1)
repeat{ll<-readLines(con,2000);if(!length(ll))break;tc<-c(tc,sub('\t.*','',ll))};close(con)
uni<-intersect(uni,tc);writeLines(sort(uni),file.path(CORE,'signature_gene_universe.txt'))
tr<-which(!M$patient_key%in%TESTIDS);ba<-basis(X[,tr,drop=FALSE],M[tr],5,TRUE)
stopifnot(!any(ba$metadata$patient_key%in%TESTIDS))
L<-log2(ba$M[uni,,drop=FALSE]*1e6+1);other<-setdiff(colnames(L),'C0');fc<-L[,'C0']-apply(L[,other,drop=FALSE],1,max)
pat<-ba$metadata;th<-ba$theta[uni,,drop=FALSE]*1e6
# Paired donor profiles: a donor contributes once regardless of the number of specimens.
contr<-list();dscon<-list()
for(id in unique(pat$patient_key)){j0<-which(pat$patient_key==id&pat$cell_type=='C0');jo<-which(pat$patient_key==id&pat$cell_type%in%setdiff(CLUSTERS,'C0'))
 if(length(j0)&&length(jo))contr[[id]]<-log2(th[,j0]+1)-rowMeans(log2(th[,jo,drop=FALSE]+1))}
for(ds in unique(pat$dataset)){j0<-which(pat$dataset==ds&pat$cell_type=='C0');jo<-which(pat$dataset==ds&pat$cell_type%in%setdiff(CLUSTERS,'C0'))
 if(length(j0)&&length(jo))dscon[[ds]]<-rowMeans(log2(th[,j0,drop=FALSE]+1))-rowMeans(log2(th[,jo,drop=FALSE]+1))}
gmt<-strsplit(readLines(file.path(OUT,'resources/h.all.v2025.1.Hs.symbols.gmt')),'\t')
state<-unique(unlist(lapply(gmt,function(a)if(a[1]%in%paste0('HALLMARK_',c('INFLAMMATORY_RESPONSE','HYPOXIA','TNFA_SIGNALING_VIA_NFKB')))a[-c(1,2)]else NULL)))
cand<-data.table(gene=uni,C0_CPM=ba$M[uni,'C0']*1e6,worst_other_log2FC=fc,paired_patient_positive=rowMeans(do.call(cbind,contr)>0),n_paired_patients=length(contr),positive_datasets=rowSums(do.call(cbind,dscon)>0),n_evaluable_datasets=length(dscon),state_or_generic_ECM=uni%in%state|grepl('^COL',uni)|uni=='FN1')
# Historical DE only annotates the frozen training-selected list; it cannot affect inclusion/order.
cand[,eligible:=C0_CPM>=5&worst_other_log2FC>=1&paired_patient_positive>=.7&positive_datasets>=2&!state_or_generic_ECM]
setorder(cand,-worst_other_log2FC);genes<-head(cand[eligible==TRUE,gene],20)
deg<-fread(file.path(OUT,'phase1_core/C0_vs_others_DEG.tsv'));cand<-merge(cand,deg[,.(gene,historical_logFC=logFC,historical_FDR=FDR)],by='gene',all.x=TRUE,sort=FALSE)
cand[,selected:=gene%in%genes];WT(cand,'identity_signature_candidate_audit.tsv');writeLines(genes,file.path(CORE,'C0_identity_signature_genes.txt'));WT(cand[state_or_generic_ECM==TRUE&worst_other_log2FC>0],'state_ECM_candidates_not_in_identity.tsv')
cat('Training-selected genes:',paste(genes,collapse=', '),'\n');flush.console()
if(length(genes)<5){writeLines('FAIL_TOO_FEW_TRAINING_SPECIFIC_GENES',file.path(CORE,'SIGNATURE_STATUS.txt'));quit(status=0)}
# Unidirectional normalized average-rank score, fixed universe and average ties.
score<-function(v){v<-v[uni];rr<-rank(v,ties.method='average');N<-length(v);k<-length(genes);(mean(rr[match(genes,uni)])-(k+1)/2)/(N-k)-.5}
stopifnot(abs(score(setNames(seq_along(uni),uni)))<=.5)
# Reserved-patient type profiles, same aggregation unit as reference statistics.
kt<-paste(M$patient_key,M$cell_type,sep='::');lv<-unique(kt);AA<-sparseMatrix(i=seq_along(kt),j=match(kt,lv),x=1,dims=c(length(kt),length(lv)))
xx<-X%*%AA;md<-M[,.(n_cells=sum(n_cells),dataset=dataset[1]),by=.(key=paste(patient_key,cell_type,sep='::'),patient_key,cell_type)];md<-md[match(lv,key)]
j<-which(md$patient_key%in%TESTIDS&md$n_cells>=5);ctrl<-copy(md[j]);ctrl[,score:=vapply(j,function(k)score(setNames(as.numeric(xx[,k]),rownames(xx))),numeric(1))]
WT(ctrl,'held_out_patient_type_scores.tsv')
auc<-function(pos,neg){mean(outer(pos,neg,`>`))+.5*mean(outer(pos,neg,`==`))}
a<-auc(ctrl[cell_type=='C0',score],ctrl[cell_type%in%setdiff(CLUSTERS,'C0'),score]);aall<-auc(ctrl[cell_type=='C0',score],ctrl[cell_type!='C0',score])
paired<-dcast(ctrl[cell_type%in%c('C0','C5')],patient_key~cell_type,value.var='score');paired<-paired[complete.cases(paired)];delta<-median(paired$C0-paired$C5)
WT(data.table(AUC_C0_vs_other_stroma=a,AUC_C0_vs_all_other_types=aall,n_C0_profiles=sum(ctrl$cell_type=='C0'),n_other_stroma_profiles=sum(ctrl$cell_type%in%setdiff(CLUSTERS,'C0')),n_C0_C5_pairs=nrow(paired),median_paired_C0_minus_C5=delta),'signature_identity_validation.tsv')
mix<-copy(TEST);mix[,signature_score:=vapply(sample,function(s)score(Matrix::rowSums(X[,M$sample==s,drop=FALSE])),numeric(1))];mix[,true_total_stroma:=n_stromal/total_cells]
WT(mix[,.(sample,patient_key,dataset,true_C0,true_total_stroma,signature_score)],'held_out_mixture_scores.tsv')
rmix<-cor_safe(mix$true_C0,mix$signature_score);rstro<-cor_safe(mix$true_total_stroma,mix$signature_score)
rpartial<-cor(resid(lm(rank(signature_score)~rank(true_total_stroma)+factor(dataset),mix)),resid(lm(rank(true_C0)~rank(true_total_stroma)+factor(dataset),mix)))
byds<-mix[,.(n=.N,rho_C0=cor_safe(true_C0,signature_score),rho_total_stroma=cor_safe(true_total_stroma,signature_score)),by=dataset];WT(byds,'signature_mixture_dataset_correlations.tsv')
WT(data.table(rho_C0_fraction=rmix,rho_total_stroma=rstro,partial_rank_rho_C0_adjusted_stroma_dataset=rpartial),'signature_mixture_correlations.tsv')
pass<-a>=.75&delta>0&rmix>0&all(byds[n>=5,rho_C0]>0)
writeLines(if(pass)'PASS_EXPLORATORY_EXPRESSION_SCORE_ONLY'else'FAIL_IDENTITY_OR_MIXTURE_VALIDATION',file.path(CORE,'SIGNATURE_STATUS.txt'))
p<-ggplot(ctrl,aes(factor(cell_type,levels=TYPES),score,fill=cell_type))+geom_boxplot(width=.6,outlier.shape=NA)+geom_point(position=position_jitter(width=.1,seed=123),size=1.3)+scale_fill_manual(values=c(COLORS,setNames(rep('grey80',7),setdiff(TYPES,CLUSTERS))))+theme_pub()+theme(legend.position='none',axis.text.x=element_text(angle=35,hjust=1))+labs(title='C0 identity signature in held-out patient profiles',subtitle=sprintf('%d genes | C0 vs other stroma AUC = %.3f',length(genes),a),x=NULL,y='C0-like signature score',caption='Patient-type aggregates; reserved patients excluded from signature selection.')
export_plot(p,'FigS_signature_celltype_specificity',9,5.3)
p<-ggplot(mix,aes(true_C0,signature_score,color=dataset))+geom_point(size=2.7)+scale_color_manual(values=DSCOL)+theme_pub()+theme(legend.position='bottom')+labs(title='C0-like expression signal in held-out mixtures',subtitle=sprintf('rho with C0 fraction = %.3f | rho with total stroma = %.3f',rmix,rstro),x='Observed C0 / all Fibro + Pericyte',y='C0-like signature score',color=NULL,caption='Expression score is not a cell-fraction estimate; total stromal content may contribute.')
export_plot(p,'FigS_signature_mixture_validation',7.5,5.5)
cat('Identity AUC',a,'paired delta',delta,'mixture rho',rmix,'stroma rho',rstro,'pass',pass,'\n')
