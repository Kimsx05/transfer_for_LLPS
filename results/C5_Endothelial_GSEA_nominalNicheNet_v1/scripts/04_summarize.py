from pathlib import Path
import pandas as pd,numpy as np,json,itertools
from scipy.stats import spearmanr
O=Path('results/C5_Endothelial_GSEA_nominalNicheNet_v1');OLD=Path('results/C5_Endothelial_NicheNet_spatial_v1')
def read(p,old=False):return pd.read_csv((OLD if old else O)/p,sep='\t')
def write(d,p):d.to_csv(O/p,sep='\t',index=False)
a=read('03_nominal_NicheNet/all_activities.tsv');defs=read('02_nominal_response_genesets/frozen_genesets.tsv');ranks=a.pivot(index='test_ligand',columns='program',values='rank')
rows=[]
for x,y in itertools.combinations(['Set_A','Set_B','Set_C','Set_A_noECM','Set_B_noECM','Set_C_noECM','EMT'],2):
 if x not in ranks or y not in ranks:continue
 rows.append(dict(program1=x,program2=y,rank_spearman=spearmanr(ranks[x],ranks[y]).statistic,top10_intersection=int(((ranks[x]<=10)&(ranks[y]<=10)).sum()),top20_intersection=int(((ranks[x]<=20)&(ranks[y]<=20)).sum())))
write(pd.DataFrame(rows),'04_rank_sensitivity/rank_comparisons.tsv');write(ranks.reset_index(),'04_rank_sensitivity/all_ligand_ranks.tsv')
base=a[a.program.eq('Set_B')].copy();loo=a[a.program.str.startswith('LOPO_')];ls=[]
for prog,z in loo.groupby('program'):
 j=base.merge(z,on='test_ligand',suffixes=('_B','_LOPO'));ls.append(dict(program=prog,n_genes=z.geneset_n.iloc[0],rank_spearman=spearmanr(j.rank_B,j.rank_LOPO).statistic,top10_overlap=((j.rank_B<=10)&(j.rank_LOPO<=10)).sum(),top20_overlap=((j.rank_B<=20)&(j.rank_LOPO<=20)).sum()))
write(pd.DataFrame(ls),'04_rank_sensitivity/LOPO_rank_summary.tsv')
lr=loo.groupby('test_ligand').agg(LOPO_rank_median=('rank','median'),LOPO_rank_min=('rank','min'),LOPO_rank_max=('rank','max'),LOPO_top10_fraction=('rank',lambda x:(x<=10).mean()),LOPO_top20_fraction=('rank',lambda x:(x<=20).mean()),LOPO_n=('rank','size')).reset_index();write(lr,'04_rank_sensitivity/LOPO_ligand_stability.tsv')
null=read('04_rank_sensitivity/random_activities.tsv.gz');nc=[]
for _,z in base.iterrows():
 q=null[null.test_ligand.eq(z.test_ligand)];nc.append(dict(test_ligand=z.test_ligand,rank_SetB=z['rank'],aupr_corrected=z.aupr_corrected,random_n=len(q),null_auc_median=q.aupr_corrected.median(),null_auc_q025=q.aupr_corrected.quantile(.025),null_auc_q975=q.aupr_corrected.quantile(.975),empirical_exceedance=(1+(q.aupr_corrected>=z.aupr_corrected).sum())/(1+len(q)),random_rank_median=q['rank'].median(),random_top10_fraction=(q['rank']<=10).mean()))
nc=pd.DataFrame(nc)
from scipy.stats import false_discovery_control
nc['empirical_BH_across246']=false_discovery_control(nc.empirical_exceedance,method='bh');write(nc,'04_rank_sensitivity/random_program_calibration.tsv')
metrics=read('02_nominal_response_genesets/all_gene_audit.tsv');summ=[]
for s in ['A','B','C']:
 d=read(f'02_nominal_response_genesets/Set_{s}_metrics.tsv');summ.append(dict(set=s,n=len(d),min_FDR=d['adj.P.Val'].min(),median_primary_P=d['P.Value'].median(),median_logFC=d.logFC.median(),cohort_logFC_spearman=spearmanr(d.logFC_GSE222315,d.logFC_InhouseData).statistic,stage_positive_fraction=(d.logFC_stage_adjusted>0).mean(),sampling_positive_fraction=(d.logFC_sampling_adjusted>0).mean(),ECM_n=d.ECM_marker.sum(),ECM_genes=';'.join(d.loc[d.ECM_marker,'gene'])))
write(pd.DataFrame(summ),'02_nominal_response_genesets/set_diagnostics.tsv')
gloo=read('04_rank_sensitivity/SetB_gene_LOPO_metrics.tsv');st=gloo.groupby('gene').agg(n_tested=('logFC','size'),positive_fraction=('logFC',lambda x:(x>0).mean()),selected_fraction=('selected','mean'),min_logFC=('logFC','min'),max_logFC=('logFC','max'),max_nominal_P=('P.Value','max')).reset_index();write(st,'04_rank_sensitivity/SetB_gene_influence_summary.tsv')
co=[]
for ds in ['InhouseData','GSE222315']:
 d=read(f'01_receiver_pseudobulk/DE_{ds}.tsv',True);other='GSE222315' if ds=='InhouseData' else 'InhouseData';q=read(f'01_receiver_pseudobulk/DE_{other}.tsv',True).set_index('gene');d=d[(d.logFC>0)&(d['P.Value']<.05)].copy();d['other_cohort_logFC']=d.gene.map(q.logFC);d['other_cohort_P']=d.gene.map(q['P.Value']);d['discovery_cohort']=ds;co.append(d)
write(pd.concat(co),'02_nominal_response_genesets/cohort_nominal_discoveries.tsv')
# Identity: patient-level aggregation retains sample counts and cannot diagnose doublets/ambient RNA definitively.
idn=read('02_nominal_response_genesets/identity_co_detection_specimen.tsv');cols=[c for c in idn if c.startswith('n_')];pat=idn.groupby(['patient_id','dataset'])[cols].sum().reset_index()
for c in cols[1:]:pat[c.replace('n_','fraction_',1)]=pat[c]/pat.n_endothelial
write(pat,'02_nominal_response_genesets/identity_co_detection_patient.tsv')
# Fixed old resources; all ranked ligands retained, not only intersection positives.
expr=read('03_candidate_ligands/expression_summary.tsv',True);c5=expr[expr.group.eq('Fibro_C5')].drop(columns='group').rename(columns={'gene':'test_ligand',**{c:'C5_'+c for c in expr if c not in ['gene','group']}})
cc=read('06_CellChat_integration/tumor_LR_all_thresholds.tsv',True);cm=cc[(cc.min_cells==20)&(cc.positive_patients>0)].copy();cm=cm.rename(columns={'ligand':'test_ligand'})
e=base.merge(c5,on='test_ligand',how='left').merge(cm,on='test_ligand',how='outer',suffixes=('','_CC'));e['evidence_class']=np.where(e['rank'].isna(),'CellChat_only',np.where(e.feature.isna(),'NicheNet_only','both'))
spec=read('03_candidate_ligands/specificity_summary.tsv',True)
for comparison in ['Fibro_C0','Fibro_C1','Fibro_C4','REST']:
 s=spec[spec.comparison.eq(comparison)].drop(columns='comparison').rename(columns={'gene':'test_ligand',**{c:comparison+'_'+c for c in spec if c not in ['gene','comparison']}});e=e.merge(s,on='test_ligand',how='left')
# Reuse previously audited receptor complex definitions, with expression support for all existing CC axes.
old=read('06_CellChat_integration/complete_expression_activity_LR_target_evidence.tsv',True);rec=old[['test_ligand','receptor','required_subunits','Endo_min_subunit_pct','Endo_min_subunit_support_patients','NicheNet_expressed_receptors']].drop_duplicates();e=e.merge(rec,on=['test_ligand','receptor'],how='left');e=e.merge(lr,on='test_ligand',how='left').merge(nc[['test_ligand','empirical_exceedance','empirical_BH_across246']],on='test_ligand',how='left')
links=read('03_nominal_NicheNet/top10_ligand_target_links.tsv');tar=links[links.program.eq('Set_B')].groupby('ligand').target.agg(lambda x:';'.join(sorted(set(x.dropna()))));e['nominal_program_predicted_targets']=e.test_ligand.map(tar);e['nominal_target_link_status']=np.where(e['rank']>10,'NOT_EXTRACTED_BEYOND_TOP10',np.where(e.nominal_program_predicted_targets.fillna('').str.len()>0,'TOP250_PROGRAM_INTERSECTION','NO_TOP250_PROGRAM_INTERSECTION'))
pairs=read('00_input_audit/official_human_v2_LR_pairs.tsv');edges=set(zip(pairs.ligand,pairs.receptor));e['NicheNet_direct_subunit_edges']=e.apply(lambda z:';'.join(g for g in str(z.required_subunits).split(';') if (z.test_ligand,g) in edges) if pd.notna(z.required_subunits) else '',axis=1);e['all_CC_required_subunits_have_NicheNet_direct_edges']=e.apply(lambda z:all((z.test_ligand,g) in edges for g in str(z.required_subunits).split(';')) if pd.notna(z.required_subunits) else np.nan,axis=1)
sp=read('08_patient_dataset_sensitivity/patient_bootstrap_LODO.tsv',True);sp=sp[(sp.top_fraction==.1)&(sp.ring==1)&(sp.metric=='adjusted_rho')&(sp.omitted_dataset=='NONE')&(sp.role!='curated_program')];sp=sp.rename(columns={'ligand':'test_ligand','median_effect':'existing_spatial_adjusted_rho','ci_low':'spatial_ci_low','ci_high':'spatial_ci_high','BH_FDR':'spatial_BH_FDR','n_patients':'spatial_n_patients'})
e=e.merge(sp[['test_ligand','receptor','existing_spatial_adjusted_rho','spatial_ci_low','spatial_ci_high','spatial_BH_FDR','spatial_n_patients']],on=['test_ligand','receptor'],how='left');e['spatial_status']=np.where(e.spatial_n_patients.notna(),'EXISTING_AXIS_EVALUATED','NOT_PREVIOUSLY_EVALUATED_AS_THIS_AXIS');e=e.sort_values(['rank','test_ligand']);write(e,'06_CellChat_spatial_integration/SetB_complete_evidence.tsv');write(cc[cc.ligand.isin(base.loc[base['rank']<=20,'test_ligand'])],'06_CellChat_spatial_integration/top20_CellChat_all_thresholds.tsv')
# Standard program comparisons: rank only. Explicit distinct MSigDB versions.
prior=read('04_NicheNet_ligand_activity/all_program_activities.tsv',True);prior['resource_version']='2025.1.Hs';new=a[a.program.isin(['Set_A','Set_B','Set_C','EMT'])].copy();new['resource_version']=np.where(new.program.eq('EMT'),'2026.1.Hs','observed nominal-DEG');comp=pd.concat([prior,new],ignore_index=True);write(comp,'05_EMT_program_NicheNet/program_rank_comparison.tsv');write(comp.pivot(index='test_ligand',columns='program',values='rank').reset_index(),'05_EMT_program_NicheNet/program_rank_wide.tsv')
focus=['ADM','MDK','SEMA3C','WNT2','FN1','COL1A1','COL1A2','TGFB1','VEGFA'];write(a[a.test_ligand.isin(focus)&~a.program.str.startswith('LOPO_')],'03_nominal_NicheNet/requested_ligand_ranks.tsv')
# Precommitted optional lightweight audit selection, before reading new spatial expression.
oldaxes=read('07_spatial_expression_validation/frozen_axes_subunits.tsv',True);tested=set(oldaxes.ligand);candidates=e[(e['rank']<=20)&(e.evidence_class=='both')&~e.test_ligand.isin(tested)].sort_values(['rank','patient_condition_median'],ascending=[True,False]).drop_duplicates('test_ligand').head(3);frozen=read('06_CellChat_spatial_integration/frozen_new_spatial_axes.tsv');assert list(zip(candidates.test_ligand,candidates.receptor))==list(zip(frozen.test_ligand,frozen.receptor));write(candidates,'06_CellChat_spatial_integration/new_axes_full_evidence.tsv')
print('SUMMARY_COMPLETE',base[['test_ligand','rank']].head(10).to_dict('records'))
# Additional transparent influence and target membership tables.
ti=read('04_rank_sensitivity/target_influence_activities.tsv');ts=[]
for gene,z in ti.groupby('removed_gene'):
 j=base.merge(z,on='test_ligand',suffixes=('_B','_drop'));ts.append(dict(removed_gene=gene,rank_spearman=spearmanr(j.rank_B,j.rank_drop).statistic,top10_overlap=((j.rank_B<=10)&(j.rank_drop<=10)).sum(),max_top10_rank_change=(j.loc[j.rank_B<=10,'rank_drop']-j.loc[j.rank_B<=10,'rank_B']).abs().max()))
write(pd.DataFrame(ts),'04_rank_sensitivity/target_influence_summary.tsv');write(ti.merge(base[['test_ligand','rank']],on='test_ligand',suffixes=('_drop','_original')),'04_rank_sensitivity/target_influence_rank_changes.tsv')
h=read('00_input_audit/Hallmark_2026_frozen.tsv');write(links.merge(h,on='target' if 'target' in h else 'gene') if 'target' in h else links.merge(h.rename(columns={'gene':'target'}),on='target',how='left'),'03_nominal_NicheNet/target_Hallmark_membership.tsv')
write(a[a.program.eq('EMT')],'05_EMT_program_NicheNet/EMT_activities_REFERENCE_ONLY.tsv');write(links[links.program.eq('EMT')],'05_EMT_program_NicheNet/EMT_top10_target_links.tsv');write(pd.DataFrame([dict(program='EMT',original_genes=h[h.pathway.eq('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')].gene.nunique(),matched_genes=defs[defs.program.eq('EMT')].gene.nunique(),status=(O/'05_EMT_program_NicheNet/EMT_STATUS.txt').read_text().strip(),version='2026.1.Hs')]),'05_EMT_program_NicheNet/EMT_coverage.tsv')
# Patient aggregation of the lightweight expression audit: descriptive only, not inferential replication.
light=read('06_CellChat_spatial_integration/new_axes_expression_audit.tsv');known=light[light.patient_id_status.isin(['VERIFIED_WITHIN_STUDY','AUTHOR_CROSSWALK_ACCEPTED_BY_USER_SOURCE_PRIORITY'])&light.patient_uid.notna()].copy()
if len(known):write(known.groupby(['ligand','receptor','dataset','patient_uid']).agg(n_sections=('section_id','size'),ligand_detection=('ligand_detection_sender','mean'),receptor_joint_detection=('receptor_joint_detection_receiver','mean'),edge_detection=('edge_joint_detection','mean')).reset_index(),'06_CellChat_spatial_integration/new_axes_known_patient_description.tsv')
# Selection must agree with the frozen pre-spatial file (no post-spatial reselection).
worst=[]
for ligand in base.loc[base['rank']<=10,'test_ligand']:
 z=loo[loo.test_ligand.eq(ligand)].sort_values(['rank','program'],ascending=[False,True]);worst.append(z.iloc[0].to_dict())
write(pd.DataFrame(worst),'04_rank_sensitivity/top10_worst_LOPO_patient.tsv')
program_stats=[]
for p,z in comp.groupby('program'):
 j=base[['test_ligand','rank']].merge(z[['test_ligand','rank']],on='test_ligand',suffixes=('_B','_other'))
 program_stats.append(dict(program=p,rank_spearman_with_SetB=spearmanr(j.rank_B,j.rank_other).statistic,top10_overlap=((j.rank_B<=10)&(j.rank_other<=10)).sum(),top20_overlap=((j.rank_B<=20)&(j.rank_other<=20)).sum()))
write(pd.DataFrame(program_stats),'05_EMT_program_NicheNet/program_agreement_with_SetB.tsv')
# Exact status codes; unresolved or user-assumed identifiers never become patients.
write(light[['section_id','dataset','patient_uid','patient_id_status','tissue']].drop_duplicates(),'00_input_audit/reused_spatial_identity.tsv')
staticcols=['test_ligand','receptor','feature','positive_patients','eligible_patients','positive_datasets','eligible_datasets','patient_condition_median','required_subunits','Endo_min_subunit_pct','Endo_min_subunit_support_patients','NicheNet_expressed_receptors','NicheNet_direct_subunit_edges','all_CC_required_subunits_have_NicheNet_direct_edges','existing_spatial_adjusted_rho','spatial_ci_low','spatial_ci_high','spatial_n_patients']+[c for c in e if c.startswith(('C5_','Fibro_C','REST_'))]
static=e[staticcols].drop_duplicates();write(a[a.program.isin(['Set_A','Set_B','Set_C','EMT'])].merge(static,on='test_ligand',how='left'),'06_CellChat_spatial_integration/all_primary_sets_expression_CC_evidence.tsv.gz')
