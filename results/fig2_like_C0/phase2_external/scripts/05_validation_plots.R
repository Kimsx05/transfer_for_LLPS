source('results/fig2_like_C0/phase1B_review/scripts/common.R')
P2<-file.path(OUT,'phase2_external');CORE<-file.path(P2,'reference_validation')
pred<-fread(file.path(CORE,'natural_mixture_predictions.tsv'))
metrics<-fread(file.path(CORE,'C0_validation_metrics.tsv'))
d<-pred[scenario=='Full_RNA'&primary_patient_specimen==TRUE&n_stromal>=50]
a<-metrics[scenario=='Full_RNA'&dataset=='OVERALL']
p<-ggplot(d,aes(true_C0,pred_C0,shape=dataset))+geom_abline(slope=1,intercept=0,linetype='dashed',color='grey50')+
 geom_point(size=3,color=COLORS['C0'],stroke=1)+scale_shape_manual(values=c(16,17,15,18,8))+coord_equal(xlim=c(0,1),ylim=c(0,1))+
 scale_x_continuous(labels=scales::percent)+scale_y_continuous(labels=scales::percent)+theme_pub()+
 theme(legend.position='bottom')+labs(title='C0 recovery in held-out patients',subtitle=sprintf('Spearman rho = %.2f | MAE = %.2f | RMSE = %.2f | n = %d',a$spearman,a$MAE,a$RMSE,a$n_patients),
 x='Observed C0 / all Fibro + Pericyte',y='Estimated C0 / all Fibro + Pericyte',shape='Held-out dataset',
 caption='One specimen per patient; at least 50 stromal cells. Dashed line: perfect agreement.\nGSE222315 fold is not estimable: insufficient training patients for C6.')
export_plot(p,'FigS_C0_deconvolution_validation',7.5,6.6)
ctrl<-fread(file.path(CORE,'held_out_pure_type_predictions.tsv'))
c0c5<-ctrl[cell_type%in%c('C0','C5')]
# Patient-average diagnostic profiles: no patient contributes repeated weight through paired samples.
long<-melt(c0c5,id.vars=c('patient_key','cell_type','pred_stromal'),measure.vars=paste0('estimate_',CLUSTERS),variable.name='predicted',value.name='coefficient')
long[,within_stroma:=coefficient/pred_stromal];long[,predicted:=sub('estimate_','',predicted)]
long<-long[,.(within_stroma=mean(within_stroma,na.rm=TRUE)),by=.(patient_key,cell_type,predicted)]
heat<-long[,.(median=median(within_stroma,na.rm=TRUE),n_patients=.N),by=.(cell_type,predicted)]
fwrite(heat,file.path(CORE,'C0_C5_confusion_patient_summary.tsv'),sep='\t')
heat[,predicted:=factor(predicted,levels=CLUSTERS)]
p<-ggplot(heat,aes(predicted,cell_type,fill=median))+geom_tile(color='white',linewidth=.7)+geom_text(aes(label=sprintf('%.0f%%',100*median)),size=3.6)+
 scale_fill_gradient(low='white',high=COLORS['C1'],limits=c(0,1),labels=scales::percent,name='Median\nassignment')+
 theme_pub()+theme(axis.line=element_blank(),axis.ticks=element_blank())+
 labs(title='C0–C5 identity recovery',x='Estimated stromal cluster',y='Held-out pure cell type',caption='Patient-averaged profiles; assignment normalized over C0–C6.\nComponentwise medians need not sum to 100%.\nPure-type controls are diagnostic and do not replace mixed-tissue validation.')
export_plot(p,'FigS_C0_C5_deconvolution_confusion',7.5,3.9)
z<-metrics[dataset=='OVERALL'];z[,scenario:=factor(scenario,levels=c('Full_RNA','GPL570','GPL6947','GPL6102','GPL6947_gene_efficiency_SD0.5'))]
p<-ggplot(z,aes(scenario,MAE))+geom_hline(yintercept=.15,linetype='dashed',color='grey50')+geom_col(width=.65,fill=COLORS['C0'],color=COLORS['C1'],linewidth=.4)+
 geom_text(aes(label=sprintf('%.3f',MAE)),vjust=-.5,size=3.5)+scale_x_discrete(labels=c('All RNA genes','GPL570\ncoverage','GPL6947\ncoverage','GPL6102\ncoverage','GPL6947 +\ngene-efficiency shift'))+
 scale_y_continuous(expand=expansion(mult=c(0,.15)))+theme_pub()+labs(title='Platform sensitivity of C0 recovery',x=NULL,y='Mean absolute error in C0 fraction',
 caption='Same held-out patients; coverage uses actual GEO platform annotations.\nEfficiency shift: fixed gene-wise log2 SD = 0.5; a stress test, not matched microarray validation.\nDashed line: prespecified MAE gate (0.15).')
export_plot(p,'FigS_C0_deconvolution_platform_sensitivity',8.3,5.4)
cv<-fread(file.path(CORE,'training_patient_coverage.tsv'));cv[,cell_type:=factor(cell_type,levels=rev(c(CLUSTERS,'Epithelial','Endothelial','Myeloid','T_NK','B','Plasma','Mast')))]
p<-ggplot(cv,aes(held_out_dataset,cell_type,fill=pmin(n_patients,10)))+geom_tile(color='white')+geom_text(aes(label=ifelse(n_patients<3,paste0(n_patients,' !'),n_patients)),size=3.1)+
 scale_fill_gradient(low='white',high=COLORS['C1'],limits=c(0,10),name='Training\npatients',breaks=c(0,3,6,10),labels=c('0','3','6','10+'))+
 theme_pub()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=20,hjust=1))+
 labs(title='Patient coverage of the 14-type reference',x='Held-out dataset',y=NULL,caption='At least 20 cells per training patient/type. ! marks fewer than three patients.\nA deficient type invalidates the fold; it is not removed or merged.')
export_plot(p,'FigS_C0_reference_patient_coverage',8,6.7)
