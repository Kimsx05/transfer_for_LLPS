source('results/fig2_like_C0/music_optimization/scripts/common.R')
suppressPackageStartupMessages(library(survival))
SIG<-file.path(OPT,'06_signature_fallback');stopifnot(readLines(file.path(SIG,'SIGNATURE_STATUS.txt'))=='PASS_EXPLORATORY_EXPRESSION_SCORE_ONLY')
GRCOL<-c(Low=GROUP_COLORS['Normal'],High=GROUP_COLORS['Tumor']);names(GRCOL)<-c('Low','High')
stage_cols<-c('I'='#D9EAF2','II'='#A6CEE3','III'='#4F98C5','IV'='#1F78B4','Ta'='#D9EAF2','Tis'='#C9DAD0','T1'='#A6CEE3','T2'='#71AFCE','T3'='#4F98C5','T4'='#1F78B4')
allstage<-allcox<-allkm<-allph<-allmodels<-inclusions<-list()
for(cohort in c('TCGA','GSE31684','GSE48075','GSE13507')){
 cat('Clinical analysis',cohort,'\n');flush.console()
 d<-fread(file.path(SIG,paste0(cohort,'_signature_clinical.tsv')));d[,score:=C0_like_signature_score]
 if(cohort=='TCGA'){
  d[,stage:=sub('Stage ','',stage_group)];lev<-c('I','II','III','IV');label<-'Pathological AJCC stage';d[,stage_num:=match(stage,lev)]
 }else{
  d[,stage:=T_component];lev<-c('Ta','Tis','T1','T2','T3','T4');label<-switch(cohort,GSE31684='Pathological T at cystectomy',GSE48075='Clinical T',GSE13507='Reported TNM T component')
  if(cohort=='GSE48075')d[grepl('pathological',stage_definition),stage:=NA_character_]
  d[,stage_num:=as.numeric(c(Ta=0,T1=1,T2=2,T3=3,T4=4)[stage])]
 }
 d[,stage:=factor(stage,levels=lev)];ds<-d[!is.na(stage)&is.finite(score)];ds[,stage:=droplevels(stage)]
 inc<-d[,.(sample,patient_unit,stage_raw,stage,stage_num,os_months,os_event,score)];inc[,`:=`(stage_included=!is.na(stage),survival_included=is.finite(os_months)&os_months>0&os_event%in%c(0,1))];inclusions[[cohort]]<-cbind(cohort=cohort,inc)
 CORE<-file.path(OPT,'04_stage')
 kw<-kruskal.test(score~stage,ds);trend<-ds[is.finite(stage_num)];ct<-cor.test(trend$stage_num,trend$score,method='spearman',exact=FALSE)
 st<-data.table(cohort=cohort,estimand='C0-like signature score',stage_definition=label,n=nrow(ds),n_groups=nlevels(ds$stage),KW_chisq=unname(kw$statistic),KW_p=kw$p.value,trend_n=nrow(trend),trend_rho=unname(ct$estimate),trend_p=ct$p.value,trend_note='Tis excluded from ordinal trend; T is not overall stage')
 WT(st,paste0(cohort,'_signature_stage_stats.tsv'));WT(ds[,.(n=.N,median=median(score),q25=quantile(score,.25),q75=quantile(score,.75)),by=stage],paste0(cohort,'_signature_stage_groups.tsv'));allstage[[cohort]]<-st
 counts<-ds[,.N,by=stage];xl<-setNames(paste0(counts$stage,'\n(n=',counts$N,')'),counts$stage)
 p<-ggplot(ds,aes(stage,score,fill=stage))+geom_boxplot(width=.55,outlier.shape=NA,linewidth=.4)+geom_point(position=position_jitter(width=.13,seed=123),size=1.1,alpha=.6)+scale_fill_manual(values=stage_cols)+scale_x_discrete(labels=xl)+theme_pub()+theme(legend.position='none')+
  labs(title=paste(cohort,'| C0-like signature score'),subtitle=sprintf('Kruskal–Wallis P = %s | n = %d',ptext(kw$p.value),nrow(ds)),x=label,y='C0-like signature score',caption='Exploratory expression-based association; not a cell-fraction estimate.')
 export_plot(p,paste0('Fig2E_signature_',cohort,'_stage_C0'),6.6,5.3)
 # Fixed median split among valid OS cases; retained without outcome-optimized cutpoints.
 s<-d[is.finite(os_months)&os_months>0&os_event%in%c(0,1)&is.finite(score)];stopifnot(!anyDuplicated(s$patient_unit))
 cutoff<-median(s$score);s[,group:=factor(ifelse(score>cutoff,'High','Low'),levels=c('Low','High'))];sdscore<-sd(s$score);s[,score_z:=(score-mean(score))/sdscore]
 fit<-survfit(Surv(os_months,os_event)~group,data=s,conf.type='log-log');lr<-survdiff(Surv(os_months,os_event)~group,data=s);lrp<-pchisq(lr$chisq,df=1,lower.tail=FALSE)
 binary<-coxph(Surv(os_months,os_event)~group,data=s,x=TRUE);univ<-coxph(Surv(os_months,os_event)~score_z,data=s,x=TRUE)
 extract<-function(fit,model,terms){ss<-summary(fit);j<-match(terms,rownames(ss$coefficients));data.table(cohort=cohort,model=model,term=terms,n=fit$n,events=fit$nevent,HR=ss$conf.int[j,'exp(coef)'],lower95=ss$conf.int[j,'lower .95'],upper95=ss$conf.int[j,'upper .95'],p_value=ss$coefficients[j,'Pr(>|z|)'],logHR=ss$coefficients[j,'coef'],SE=ss$coefficients[j,'se(coef)'])}
 co<-rbind(extract(binary,'binary_high_vs_low','groupHigh'),extract(univ,'continuous_univariable_per_SD','score_z'))
 # Covariates are chosen only from availability and prespecified sample/event criteria.
 s[,age:=as.numeric(age_years)];s[,grade:=tolower(grade_raw)];s[grepl('high',grade),grade:='High'];s[grepl('low',grade),grade:='Low'];s[!grade%in%c('Low','High'),grade:=NA_character_];s[,grade:=factor(grade,levels=c('Low','High'))]
 vars<-c('age','stage_num');cc<-s[complete.cases(s[,..vars])]
 if(nrow(cc)<.7*nrow(s)||sum(cc$os_event)<30){vars<-'age';cc<-s[is.finite(age)]}
 if(nrow(cc)>=.7*nrow(s)&&sum(cc$os_event)>=10*(1+length(vars))){
  # Grade is optional only with >=20 per level, >=5 events/level and correlation <0.8 with stage.
  cg<-cc[!is.na(grade)];gt<-cg[,.(n=.N,events=sum(os_event)),by=grade]
  if(nrow(gt)==2&&all(gt$n>=20)&all(gt$events>=5)&&nrow(cg)>=.9*nrow(cc)&&sum(cg$os_event)>=10*(2+length(vars))&&is.finite(cor_safe(as.numeric(cg$grade),cg$stage_num))&&abs(cor_safe(as.numeric(cg$grade),cg$stage_num))<.8){vars<-c(vars,'grade');cc<-cg}
  # Drop constant covariates by identifiability only, never by P value.
  vars<-vars[vapply(vars,function(v)uniqueN(cc[[v]])>1,logical(1))]
  fm<-as.formula(paste('Surv(os_months,os_event)~score_z+',paste(vars,collapse='+')))
  multi<-coxph(fm,data=cc,x=TRUE);co<-rbind(co,extract(multi,'continuous_multivariable_per_SD',names(coef(multi))))
  models<-list(binary=binary,univariable=univ,multivariable=multi)
  allmodels[[cohort]]<-data.table(cohort=cohort,n_OS=nrow(s),events_OS=sum(s$os_event),n_multivariable=nrow(cc),events_multivariable=sum(cc$os_event),covariates=paste(vars,collapse='+'),purity='unavailable; not imputed',grade_decision=if('grade'%in%vars)'included'else'excluded: sparse/missing/correlated',score_SD=sdscore)
 }else{models<-list(binary=binary,univariable=univ);allmodels[[cohort]]<-data.table(cohort=cohort,n_OS=nrow(s),events_OS=sum(s$os_event),n_multivariable=0,covariates='NOT_ESTIMABLE')}
 ph<-rbindlist(lapply(names(models),function(nm){z<-cox.zph(models[[nm]],transform='km');cbind(cohort=cohort,model=nm,as.data.table(z$table,keep.rownames='term'))}))
 CORE<-file.path(OPT,'05_survival');WT(ph,paste0(cohort,'_signature_PH_tests.tsv'));WT(co,paste0(cohort,'_signature_Cox_stats.tsv'))
 bco<-co[model=='binary_high_vs_low'];km<-data.table(cohort=cohort,n=nrow(s),events=sum(s$os_event),n_low=sum(s$group=='Low'),n_high=sum(s$group=='High'),cutoff=cutoff,cutoff_rule='median valid OS cases; Low<=median, High>median',time_unit='months',logrank_p=lrp,HR_high_vs_low=bco$HR,lower95=bco$lower95,upper95=bco$upper95)
 WT(km,paste0(cohort,'_signature_survival_stats.tsv'));WT(s[,.(sample,patient_unit,os_months,os_event,score,score_z,group,stage,stage_num,age,grade)],paste0(cohort,'_signature_survival_samples.tsv'))
 # One graphic: survival curve and aligned at-risk rows share the same time axis.
 ss<-summary(fit,censored=TRUE);curve<-data.table(time=ss$time,survival=ss$surv,lower=ss$lower,upper=ss$upper,group=sub('group=','',ss$strata),n_censor=ss$n.censor)
 curve<-rbind(data.table(time=0,survival=1,lower=1,upper=1,group=c('Low','High'),n_censor=0),curve);setorder(curve,group,time)
 xmax<-max(s$os_months);ticks<-pretty(c(0,xmax),n=5);ticks<-ticks[ticks>=0&ticks<=xmax]
 rs<-summary(fit,times=ticks,extend=TRUE);risk<-data.table(time=rs$time,n_risk=rs$n.risk,group=sub('group=','',rs$strata));risk[,y:=ifelse(group=='Low',-.15,-.24)]
 WT(risk,paste0(cohort,'_signature_risk_table.tsv'));WT(curve,paste0(cohort,'_signature_KM_curve.tsv'))
 labs_grp<-c(Low=sprintf('Low (n=%d)',km$n_low),High=sprintf('High (n=%d)',km$n_high))
 p<-ggplot(curve,aes(time,survival,color=group,fill=group))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.13,colour=NA,show.legend=FALSE)+geom_step(linewidth=.7)+
  geom_point(data=curve[n_censor>0],shape=3,size=1.4,show.legend=FALSE)+geom_text(data=risk,aes(x=time,y=y,label=n_risk,color=group),inherit.aes=FALSE,size=3.1,show.legend=FALSE)+
  annotate('text',x=-.07*xmax,y=c(-.15,-.24),label=c('Low','High'),hjust=1,size=3.2,color=GRCOL)+
  annotate('text',x=0,y=-.07,label='Number at risk',hjust=0,size=3.1,color='black')+
  scale_color_manual(values=GRCOL,labels=labs_grp,name=NULL)+scale_fill_manual(values=GRCOL)+scale_x_continuous(breaks=ticks,limits=c(-.13*xmax,xmax),expand=expansion(mult=c(0,.02)))+
  scale_y_continuous(breaks=c(0,.25,.5,.75,1),labels=scales::percent,limits=c(-.29,1.03),expand=c(0,0))+theme_pub()+theme(legend.position='top')+
  labs(title=paste(cohort,'| C0-like signature score'),subtitle=sprintf('OS | Log-rank P = %s | HR High/Low = %.2f (%.2f–%.2f)',ptext(lrp),bco$HR,bco$lower95,bco$upper95),x='Time (months)',y='Overall survival',caption=sprintf('Fixed median cutoff = %.4f; shaded bands: 95%% CI; +: censoring.\nExploratory expression score, not abundance. HR is an unadjusted Cox estimate.',cutoff))
 export_plot(p,paste0('Fig2F_signature_',cohort,'_survival_C0'),7.6,6.2)
 allcox[[cohort]]<-co;allkm[[cohort]]<-km;allph[[cohort]]<-ph
}
CORE<-file.path(OPT,'summary');st<-rbindlist(allstage);st[,KW_FDR_across_cohorts:=p.adjust(KW_p,'BH')];st[,trend_FDR_across_cohorts:=p.adjust(trend_p,'BH')];WT(st,'signature_stage_summary.tsv')
co<-rbindlist(allcox);co[term=='score_z',FDR_across_cohorts:=p.adjust(p_value,'BH'),by=model];WT(co,'signature_Cox_summary.tsv')
km<-rbindlist(allkm);km[,logrank_FDR_across_cohorts:=p.adjust(logrank_p,'BH')];WT(km,'signature_KM_summary.tsv');WT(rbindlist(allph),'signature_PH_summary.tsv');WT(rbindlist(allmodels,fill=TRUE),'clinical_model_specifications.tsv');WT(rbindlist(inclusions,fill=TRUE),'clinical_inclusion_manifest.tsv')
# Descriptive forest: comparable per-SD units, but different stages/platforms and selected cohorts.
f<-co[term=='score_z'];f[,label:=paste(cohort,ifelse(model=='continuous_univariable_per_SD','univariable','adjusted'),sep=' | ')]
p<-ggplot(f,aes(HR,factor(label,levels=rev(label))))+geom_vline(xintercept=1,linetype='dashed',color='grey50')+geom_segment(aes(x=lower95,xend=upper95,yend=factor(label,levels=rev(label))),linewidth=.6)+geom_point(size=2.5,color=COLORS['C1'])+scale_x_log10()+theme_pub()+labs(title='C0-like signature score and overall survival',x='Hazard ratio per cohort-specific SD (95% CI)',y=NULL,caption='Descriptive comparison; no pooled effect. Cohorts differ in stage composition and platform.\nProportional-hazards tests and exact adjustment sets are reported separately.')
CORE<-file.path(OPT,'05_survival');export_plot(p,'FigS_signature_Cox_forest',8.1,5.7)
cat('All four signature clinical analyses complete\n')
