from analyze import *
def interval(g,value):
 v=g[value].dropna().to_numpy();rng=np.random.default_rng(123);groups=[x[value].dropna().to_numpy() for _,x in g.groupby('dataset')];groups=[x for x in groups if len(x)];boot=np.median(np.concatenate([rng.choice(a,(2000,len(a)),replace=True) for a in groups],axis=1),axis=1) if groups else np.array([np.nan]);return dict(n_patients=len(v),median=np.median(v) if len(v) else np.nan,ci_low=np.quantile(boot,.025),ci_high=np.quantile(boot,.975),positive_patients=int((v>0).sum()),negative_patients=int((v<0).sum()),n_datasets=len(groups))
def patient_summary(d,keys,value,prefix):
 known=d[d.patient_uid.notna()].copy();pat=known.groupby(['patient_uid','dataset']+keys,dropna=False)[value].mean().reset_index();save(pat,O/'06_patient_dataset_summary'/f'{prefix}_patient_effects.tsv');summ=[];lodo=[]
 for key,g in pat.groupby(keys,dropna=False):
  key=key if isinstance(key,tuple) else (key,);tags=dict(zip(keys,key));summ.append(dict(**tags,**interval(g,value)))
  for ds in sorted(d.dataset.unique()):
   h=g[g.dataset!=ds];lodo.append(dict(**tags,omitted_dataset=ds,removed_confirmed_patients=len(g)-len(h),**interval(h,value)))
 save(summ,O/'06_patient_dataset_summary'/f'{prefix}_patient_CI.tsv');save(lodo,O/'06_patient_dataset_summary'/f'{prefix}_LODO.tsv')
 ds=d.groupby(['dataset']+keys,dropna=False).agg(median=(value,'median'),minimum=(value,'min'),maximum=(value,'max'),n_sections=(value,'size'),positive_sections=(value,lambda x:int((x>0).sum())),negative_sections=(value,lambda x:int((x<0).sum()))).reset_index();save(ds,O/'06_patient_dataset_summary'/f'{prefix}_dataset.tsv')
 # All-source LODO is descriptive: unknown identities cannot be bootstrapped as patients.
 ld=[]
 for key,g in d.groupby(keys,dropna=False):
  key=key if isinstance(key,tuple) else (key,)
  for dsx in sorted(d.dataset.unique()):ld.append(dict(zip(keys,key),omitted_dataset=dsx,n_sections=int((g.dataset!=dsx).sum()),median=g.loc[g.dataset!=dsx,value].median()))
 save(ld,O/'06_patient_dataset_summary'/f'{prefix}_LODO_sections_descriptive.tsv');return pat,pd.DataFrame(summ),ds
if __name__=='__main__':
 g=pd.read_csv(O/'tables/global_results.tsv',sep='\t');assert len(g)==2520
 pat,ci,ds=patient_summary(g,['partner','method'],'I','global')
 ring=pd.read_csv(O/'tables/ring_results.tsv',sep='\t');patient_summary(ring,['partner','method','ring'],'I','rings')
 lee=pd.read_csv(O/'tables/lee_results.tsv',sep='\t').merge(META.reset_index()[['section_id','dataset','patient_uid']],on='section_id');patient_summary(lee,['partner','method'],'L','lee')
 local=pd.read_csv(O/'tables/local_summary.tsv',sep='\t').merge(META.reset_index()[['section_id','dataset','patient_uid','tissue']],on='section_id');local['HH_fraction']=local.HH_FDR/local.n_tests;local['HH_at_least_10']=local.HH_FDR>=10
 components=[]
 for sid in META.index:
  d=pd.read_parquet(O/'03_local_bivariate_moran'/f'{sid}.parquet');d=d[d.analysis_included].reset_index(drop=True);w=weights(graph(d))
  for p in PARTNERS:
   hh=(d[p+'__class']=='HH').to_numpy();nh=int(hh.sum());nc,lab=connected_components(w.sparse[hh][:,hh],directed=False) if nh else (0,np.array([],int));sizes=np.bincount(lab) if nh else np.array([],int);components.append(dict(section_id=sid,partner=p,HH_components=nc,HH_largest_component=int(sizes.max()) if nh else 0,HH_components_ge10=int((sizes>=10).sum())))
 local=local.merge(pd.DataFrame(components),on=['section_id','partner'])
 save(local,O/'tables/local_summary_annotated.tsv');patient_summary(local,['partner'],'HH_fraction','local_HH')
 save(local.groupby(['partner','dataset']).agg(n_sections=('section_id','size'),sections_any_HH=('HH_FDR',lambda a:int((a>0).sum())),sections_ge10_HH=('HH_FDR',lambda a:int((a>=10).sum())),median_HH_fraction=('HH_fraction','median'),total_HH_spots=('HH_FDR','sum'),sections_contiguous_HH_ge10=('HH_largest_component',lambda a:int((a>=10).sum()))).reset_index(),O/'06_patient_dataset_summary/local_HH_dataset_support.tsv')
 # Ring 0 has no spatial lag: report Pearson as scale-matched counterpart and Spearman.
 same0=[]
 for sid in META.index:
  d=pd.read_parquet(O/'03_local_bivariate_moran'/f'{sid}.parquet');d=d[d.analysis_included]
  for partner in PARTNERS:
   for method in ['raw','adjusted']:
    x=d.five_UCell if method=='raw' else d[partner+'__adjusted_five'];y=d['major_mean__'+partner] if method=='raw' else d[partner+'__adjusted_abundance'];same0.append(dict(section_id=sid,partner=partner,method=method,pearson=stats.pearsonr(x,y).statistic,spearman=stats.spearmanr(x,y).statistic))
 same0=pd.DataFrame(same0);save(same0,O/'tables/ring0_same_spot.tsv')
 rings0=ring.merge(same0,on=['section_id','partner','method']);rings0['opposite_to_same_spot_Pearson']=np.sign(rings0.I)!=np.sign(rings0.pearson);save(rings0,O/'tables/ring0_to_ring123_comparison.tsv')
 primary=g[g.method=='primary'];wide=g.pivot(index=['section_id','partner'],columns='method',values='I').reset_index();delta=[]
 for p in PARTNERS:
  w=wide[wide.partner==p]
  for method in ['reverse','AUCell','log1p_abundance','major_fraction','q05','adjusted_base','adjusted_full','A26','Core15']:
   delta.append(dict(partner=p,method=method,spearman_vs_primary=stats.spearmanr(w.primary,w[method]).statistic,sign_agreement=np.mean(np.sign(w.primary)==np.sign(w[method])),median_difference=np.median(w[method]-w.primary)))
 save(delta,O/'tables/method_agreement.tsv');cl=primary.merge(lee[lee.method=='raw'],on=['section_id','partner'],suffixes=('_moran','_lee'));save(cl,O/'tables/moran_lee_comparison.tsv');save([dict(partner=p,spearman=stats.spearmanr(z.I,z.L).statistic,sign_agreement=np.mean(np.sign(z.I)==np.sign(z.L))) for p,z in cl.groupby('partner')],O/'tables/moran_lee_agreement.tsv')
 # Patient paired raw-adjusted difference has a patient-level interval.
 ww=wide.merge(META.reset_index()[['section_id','dataset','patient_uid']],on='section_id');ww['adjusted_minus_raw']=ww.adjusted_full-ww.primary;patient_summary(ww,['partner'],'adjusted_minus_raw','adjustment_change')
 same=primary.copy();patient_summary(same,['partner'],'same_spot_spearman','same_spot')
 # Same-spot negative / adjacent positive is explicitly counted, not inferred from HH maps.
 peripheral=pd.read_csv(O/'tables/peripheral_patterns.tsv',sep='\t').merge(primary[['section_id','partner','I','same_spot_spearman']],on=['section_id','partner']).merge(META.reset_index()[['section_id','dataset','patient_uid']],on='section_id');peripheral=peripheral.merge(same0[same0.method=='raw'][['section_id','partner','pearson']],on=['section_id','partner']);peripheral['section_Pearson_negative_neighbor_positive']=(peripheral.I>0)&(peripheral.pearson<0);peripheral['section_same_negative_neighbor_positive']=(peripheral.I>0)&(peripheral.same_spot_spearman<0);save(peripheral,O/'tables/peripheral_patterns_annotated.tsv');patient_summary(peripheral,['partner'],'fraction_high_same_low_neighbor_high','peripheral')
 # Tissue is descriptive only (Normal patients unconfirmed); no fixed-high fraction test.
 save(g[g.method.isin(['primary','adjusted_full'])].groupby(['dataset','tissue','partner','method']).agg(n_sections=('I','size'),median_I=('I','median')).reset_index(),O/'06_patient_dataset_summary/tissue_descriptive.tsv')
 # CellChat existing tumor-only strengths and eligible denominators, no refitting.
 ccroot=Path('results/llps_C5_cellchat_v2');cc=pd.read_csv(ccroot/'05_threshold_comparison/C5_partner_network_all_thresholds.csv');cc=cc[cc.stratification.eq('tumor_only')].copy();cc['partner']=np.where(cc.sender=='Fibro_C5',cc.receiver,cc.sender);cc['direction']=np.where(cc.sender=='Fibro_C5','outgoing','incoming');cc=cc[cc.partner.isin(PARTNERS)];save(cc,O/'07_cellchat_integration/tumor_partner_strength_all_thresholds.tsv')
 pairs=[]
 for threshold,folder in [(10,'01_threshold10'),(20,'02_threshold20'),(30,'03_threshold30')]:
  c=pd.read_csv(ccroot/folder/'C5_vs_fibro_common_specimen_strength.csv');c=c[c.other.eq('Fibro_C0')&c.partner.isin(PARTNERS)&c.clinical_group.isin(['NMIBC','MIBC'])];c=c.groupby(['patient_id','dataset','clinical_group','partner','direction']).difference.median().reset_index();c['min_cells']=threshold;pairs.append(c)
 pairs=pd.concat(pairs,ignore_index=True);save(pairs,O/'07_cellchat_integration/C5_C0_tumor_paired_patient_condition.tsv');paired_summary=pairs.groupby(['min_cells','partner','direction']).agg(paired_median_difference=('difference','median'),patient_condition_units=('difference','size'),independent_patients=('patient_id','nunique'),datasets=('dataset','nunique'),positive_units=('difference',lambda x:int((x>0).sum()))).reset_index();save(paired_summary,O/'07_cellchat_integration/C5_C0_tumor_paired_summary.tsv')
 integration=[]
 for p in PARTNERS:
  z=primary[primary.partner==p];l=local[local.partner==p];dsp=ds[(ds.partner==p)&(ds.method=='primary')];q=dict(partner=p,section_median_Moran=z.I.median(),positive_sections=int((z.I>0).sum()),positive_global_FDR_sections=int(((z.I>0)&(z.FDR<.05)).sum()),negative_global_FDR_sections=int(((z.I<0)&(z.FDR<.05)).sum()),HH_any_sections=int((l.HH_FDR>0).sum()),HH_ge10_sections=int((l.HH_FDR>=10).sum()),HH_contiguous_ge10_sections=int((l.HH_largest_component>=10).sum()),HH_any_datasets=l.loc[l.HH_FDR>0,'dataset'].nunique(),positive_dataset_medians=int((dsp['median']>0).sum()))
  for method in ['primary','adjusted_full']:
   zci=ci[(ci.partner==p)&(ci.method==method)].iloc[0]
   for col in ['median','ci_low','ci_high','n_patients']:q[method+'_'+col]=zci[col]
  for direction in ['outgoing','incoming']:
   c=cc[(cc.min_cells==20)&(cc.partner==p)&(cc.direction==direction)].iloc[0];pc=paired_summary[(paired_summary.min_cells==20)&(paired_summary.partner==p)&(paired_summary.direction==direction)].iloc[0]
   q[direction+'_CellChat_strength']=c.patient_condition_median;q[direction+'_eligible_patients']=c.eligible_patients;q[direction+'_C5_minus_C0']=pc.paired_median_difference;q[direction+'_paired_patients']=pc.independent_patients
  integration.append(q)
 save(integration,O/'07_cellchat_integration/CellChat_spatial_integration.tsv')
 print(pd.DataFrame(integration).to_string(index=False))
