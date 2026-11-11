from env import *
# Full raw scoring summaries, including separately labelled PARP treatment arms.
a=read('03_cell_models/sample_scores.tsv.gz');keys=['dataset','study','background','state','baseline_include','signature']+(['assay_treatment'] if 'assay_treatment' in a else []);tab=a.groupby(keys,dropna=False).agg(n_RNA=('score','count'),mean=('score','mean'),median=('score','median'),SD=('score','std'),min=('score','min'),max=('score','max')).reset_index();save(tab,'03_cell_models/line_state_score_summary.tsv')
# The one ratio missing value has a biological detection reason, not a missing platform gene.
pb=read('02_epi/pseudobulk_counts.tsv.gz').set_index('gene');den=pb.loc[['LRRC8A','LRRC8B','LRRC8C','LRRC8D','LRRC8E']].sum();rows=[]
for sample in den[den==0].index:rows.append(dict(sample=sample,NPEPPS_count=pb.loc['NPEPPS',sample],VRAC_denominator_total_count=0,score=np.nan,reason='all5 denominator genes mapped but zero count in this sample; native ratio NA; no pseudocount/subset substitution'))
save(pd.DataFrame(rows),'02_epi/pseudobulk_ratio_zero_audit.tsv')
# Distinguish detected member counts from fractions in existing and future Epi effect tables.
p=B/'02_epi/sample_effects.tsv';ef=pd.read_csv(p,sep='\t');c=read('00_plan_audit/epi_coverage.tsv').set_index('signature').n_present
if 'mean_detection' in ef:ef=ef.rename(columns={'mean_detection':'mean_detected_genes'})
ef['mean_detection_rate']=ef.mean_detected_genes/ef.signature.map(c);save(ef,'02_epi/sample_effects.tsv')
p=B/'scripts/23_epi_statistics.py';s=p.read_text().replace("mean_detection=d.get(sn+'_detected',pd.Series(dtype=float)).mean()", "mean_detected_genes=d.get(sn+'_detected',pd.Series(dtype=float)).mean(),mean_detection_rate=d.get(sn+'_detected',pd.Series(dtype=float)).mean()/float(cov.loc[cov.signature==sn,'n_present'].iloc[0])");p.write_text(s)
# Explicit v2 implementation amendment score comparison.
old=read('00_plan_audit/amendments/before_DMSO_reference_GSE285648_scores.tsv');new=read('00_plan_audit/scores/GSE285648_scores.tsv');cmp=old[['matrix_id','signature','score']].merge(new[['matrix_id','signature','score']],on=['matrix_id','signature'],suffixes=('_before','_after'));cmp['difference']=cmp.score_after-cmp.score_before;save(cmp,'00_plan_audit/amendments/GSE285648_reference_before_after.tsv');assert np.max(np.abs(cmp.loc[cmp.signature=='MP10','difference']))<1e-8
# Seal actual data-dependent implementation corrections, preserving original freeze files.
paths=['00_plan_audit/amendments/GSE285648_DMSO_reference.json','00_plan_audit/amendments/literal_pTN_fix.json','00_plan_audit/amendments/immediate_full_scoring.json'];seal=[]
for x in paths:
 p=B/x
 if p.exists():seal.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
save(pd.DataFrame(seal),'00_plan_audit/amendment_manifest.tsv')
# Actual dataset-by-signature coverage completeness, including all blocked definitions.
manifest=read('01_signatures/signature_manifest.tsv');coverage=read('00_plan_audit/coverage.tsv');assert coverage.groupby('dataset').signature.nunique().eq(len(manifest)).all();save(pd.DataFrame([dict(file='GSE285648_Figure1.jpg',actual_original_asset='MOL2-20-779-g001.jpg',actual_figure='Figure4',used_for_IC50=False),dict(file='GSE285648_IC50_Figure1.jpg',actual_original_asset='MOL2-20-779-g007.jpg',actual_figure='Figure1A',used_for_IC50=True)]),'00_plan_audit/source_figure_identity.tsv');print('COMPLETION_TABLES_READY')
