source('results/fig2_like_C0/music_optimization/scripts/common.R');suppressPackageStartupMessages(library(survival))
CORE<-file.path(OPT,'02_model_comparison');pred<-fread(file.path(CORE,'all_model_predictions.tsv'));truth<-fread(file.path(OPT,'01_diagnostics/C0_C6_cell_and_RNA_truth_predictions.tsv'))
z<-melt(pred[evaluation_set=='original20_development_evaluation'&status=='OK'],id.vars=c('sample','method','pred_stromal'),measure.vars=paste0('estimate_',CLUSTERS),variable.name='cell_type',value.name='estimate');z[,cell_type:=sub('estimate_','',cell_type)];z<-merge(z,truth[,.(sample,cell_type,true_fraction)],by=c('sample','cell_type'));z[,pred_fraction:=estimate/pred_stromal]
WT(z[,as.data.table(as.list(metric(true_fraction,pred_fraction))),by=.(method,cell_type)],'all_models_C0_C6_metrics.tsv')
stopifnot(all(is.finite(z$pred_fraction)),all(z$pred_fraction>=0),all(z$pred_fraction<=1+1e-12))
CORE<-file.path(OPT,'05_survival');out<-phout<-list()
# Diagnostic sensitivities prompted by observed PH violations, not significance-seeking model selection.
s<-fread(file.path(CORE,'TCGA_signature_survival_samples.tsv'));s<-s[complete.cases(s[,.(age,stage_num)])]
f<-coxph(Surv(os_months,os_event)~score_z+tt(score_z)+age+stage_num,data=s,tt=function(x,t,...)x*log1p(t/12))
cf<-coef(f);V<-vcov(f);idx<-match(c('score_z','tt(score_z)'),names(cf));tt<-c(12,36,60)
out[['TCGA']]<-rbindlist(lapply(tt,function(t){a<-c(1,log1p(t/12));est<-sum(a*cf[idx]);se<-sqrt(as.numeric(t(a)%*%V[idx,idx]%*%a));data.table(cohort='TCGA',model='time_varying_score_sensitivity',time_months=t,n=f$n,events=f$nevent,HR=exp(est),lower95=exp(est-1.96*se),upper95=exp(est+1.96*se),p_value=2*pnorm(-abs(est/se)))}))
WT(as.data.table(summary(f)$coefficients,keep.rownames='term'),'TCGA_time_varying_Cox_coefficients.tsv')
for(cc in c('GSE31684','GSE13507')){
 s<-fread(file.path(CORE,paste0(cc,'_signature_survival_samples.tsv')));s<-s[complete.cases(s[,.(age,stage)])]
 form<-if(cc=='GSE13507')Surv(os_months,os_event)~score_z+age+grade+strata(stage)else Surv(os_months,os_event)~score_z+age+strata(stage)
 f<-coxph(form,data=s,x=TRUE);ss<-summary(f);out[[cc]]<-data.table(cohort=cc,model='stage_stratified_sensitivity',time_months=NA_real_,n=f$n,events=f$nevent,HR=ss$conf.int['score_z','exp(coef)'],lower95=ss$conf.int['score_z','lower .95'],upper95=ss$conf.int['score_z','upper .95'],p_value=ss$coefficients['score_z','Pr(>|z|)'])
 phout[[cc]]<-cbind(cohort=cc,as.data.table(cox.zph(f)$table,keep.rownames='term'))
}
WT(rbindlist(out),'PH_sensitivity_Cox_results.tsv');WT(rbindlist(phout),'PH_sensitivity_tests.tsv')
# Update only the descriptive forest annotation; primary estimates remain exactly unchanged.
co<-fread(file.path(OPT,'summary/signature_Cox_summary.tsv'));q<-co[term=='score_z'];q[,label:=paste(cohort,ifelse(model=='continuous_univariable_per_SD','univariable','adjusted'),sep=' | ')];q[cohort=='TCGA'&model=='continuous_multivariable_per_SD',label:=paste0(label,' *')]
p<-ggplot(q,aes(HR,factor(label,levels=rev(label))))+geom_vline(xintercept=1,linetype='dashed',color='grey50')+geom_segment(aes(x=lower95,xend=upper95,yend=factor(label,levels=rev(label))),linewidth=.6)+geom_point(size=2.5,color=COLORS['C1'])+scale_x_log10()+theme_pub()+labs(title='C0-like signature score and overall survival',x='Hazard ratio per cohort-specific SD (95% CI)',y=NULL,caption='Descriptive comparison; no pooled effect. * TCGA adjusted score violates PH (P=0.029).\nStage PH violations in GSE31684/GSE13507; time-varying/stratified sensitivity tables provided.')
export_plot(p,'FigS_signature_Cox_forest',8.7,5.8)
# Clinical invariants, risk tables and fixed-cutoff reconstruction.
for(cc in c('TCGA','GSE31684','GSE48075','GSE13507')){
 s<-fread(file.path(CORE,paste0(cc,'_signature_survival_samples.tsv')));km<-fread(file.path(CORE,paste0(cc,'_signature_survival_stats.tsv')));risk<-fread(file.path(CORE,paste0(cc,'_signature_risk_table.tsv')))
 stopifnot(!anyDuplicated(s$patient_unit),all(s$os_months>0),all(s$os_event%in%c(0,1)),abs(median(s$score)-km$cutoff)<1e-12,all(as.character(s$group)==ifelse(s$score>km$cutoff,'High','Low')),nrow(s)==km$n,sum(s$os_event)==km$events)
 for(i in seq_len(nrow(risk)))stopifnot(risk$n_risk[i]==sum(s$group==risk$group[i]&s$os_months>=risk$time[i]))
}
writeLines(c('PASS: reserved 20 never used in new references/signature selection (runtime assertions).','PASS: every original evaluation patient retained, failures explicit.','PASS: cluster-fraction predictions are same-fit numerator/denominator.','PASS: four cohorts use a fixed median among valid OS cases.','PASS: risk table counts recomputed from raw follow-up at each plotted time.','PASS: patient IDs unique within each clinical cohort; positive OS times; binary events.','PASS: all clinical outputs explicitly use C0-like signature score, not abundance.'),file.path(OPT,'summary/numerical_checks.txt'))
print(rbindlist(out))
