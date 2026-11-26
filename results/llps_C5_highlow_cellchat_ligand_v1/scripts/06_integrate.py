import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import pandas as pd,numpy as np,hashlib
B=Path('/home/jinshengxi/bladder_metabolism/results');O=B/'llps_C5_highlow_cellchat_ligand_v1';D=O/'06_NicheNet_spatial_integration'
focus=pd.read_csv(O/'00_input_audit/candidate_overlap.csv').gene.tolist()
source=[]
def read(rel):
 p=B/rel;source.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size));return pd.read_csv(p,sep='\t')
old=read('C5_Endothelial_NicheNet_spatial_v1/04_NicheNet_ligand_activity/all_program_activities.tsv')
new=read('C5_Endothelial_GSEA_nominalNicheNet_v1/03_nominal_NicheNet/all_activities.tsv')
lopo=read('C5_Endothelial_GSEA_nominalNicheNet_v1/04_rank_sensitivity/LOPO_ligand_stability.tsv')
spec=read('C5_Endothelial_NicheNet_spatial_v1/03_candidate_ligands/specificity_summary.tsv')
spatial=read('C5_Endothelial_NicheNet_spatial_v1/tables/final_axis_summary.tsv')
axes=read('C5_Endothelial_GSEA_nominalNicheNet_v1/06_CellChat_spatial_integration/new_axes_full_evidence.tsv') if (B/'C5_Endothelial_GSEA_nominalNicheNet_v1/06_CellChat_spatial_integration/new_axes_full_evidence.tsv').exists() else pd.read_csv(next((B/'C5_Endothelial_GSEA_nominalNicheNet_v1').glob('06*/new_axes_full_evidence.tsv')),sep='\t')
cor=pd.read_csv(O/'05_LLPS_ligand_correlation/pseudobulk_continuous_associations.csv.gz');cor=cor.query("min_cells==10 and stratum=='ALL' and score_type=='score_mean'")
within=pd.read_csv(O/'05_LLPS_ligand_correlation/patient_effect_tests.csv.gz');within=within.query("sensitivity=='original' and stratum=='ALL'")
cc=pd.read_csv(O/'04_pathway_LR_analysis/within_median_t10_statistics.csv.gz');cc=cc.query("kind=='LR' and stratum=='Tumor_ALL' and direction=='outgoing' and partner=='Endothelial'")
cc=cc.merge(pd.read_csv(O/'00_input_audit/DB_interaction_gene_components.csv')[['feature','ligand_genes','receptor_genes']],on='feature',how='left',validate='many_to_one')
depth=pd.read_csv(O/'07_sensitivity_analysis/pseudobulk_dataset_depth_adjusted.csv').query('min_cells==10')
rows=[]
for gene in focus:
 row={'gene':gene,'interpretation':'Observational association only; no causal or C5-specific spatial attribution'}
 z=spec[(spec.gene==gene)&(spec.comparison=='Fibro_C0')]
 if len(z):row.update({f'C5_vs_C0_{k}':z.iloc[0][k] for k in ['n_patients','median_difference','positive_fraction']})
 z=depth[depth.gene==gene]
 if len(z):row.update({f'depth_adjusted_pseudobulk_{k}':z.iloc[0][k] for k in ['adjusted_beta','FDR_focus','residual_df']})
 z=cor[cor.gene==gene]
 if len(z):row.update({f'pseudobulk_{k}':z.iloc[0][k] for k in ['n_patients','pearson','spearman','adjusted_beta','adjusted_p','adjusted_p_FDR_focus','residual_df']})
 for metric in ['adjusted_r','spearman','expression_difference','detection_difference']:
  z=within[(within.gene==gene)&(within.metric==metric)]
  if len(z):row.update({f'within_{metric}_{k}':z.iloc[0][k] for k in ['n','median','ci_low','ci_high','FDR_focus']})
 for label,tab in [('first_NicheNet',old),('nominal_NicheNet',new)]:
  z=tab[tab.test_ligand==gene]
  for _,r in z.iterrows():
   for k in ['rank','aupr_corrected','pearson']:row[f'{label}_{r.program}_{k}']=r[k]
 z=lopo[lopo.test_ligand==gene]
 if len(z):row.update(z.drop(columns='test_ligand').iloc[0].to_dict())
 z=spatial[spatial.ligand==gene]
 if len(z):row.update({f'spatial_{k}':v for k,v in z.drop(columns='ligand').iloc[0].to_dict().items()})
 z=cc[cc.ligand==gene].sort_values('difference_median',key=abs,ascending=False)
 row['CellChat_Endothelial_all_axes']=';'.join(z.feature.astype(str));row['CellChat_Endothelial_axis_count']=len(z)
 # Preserve every receptor-specific value in a linked long table; do not cherry-pick one axis for evidence tier.
 row['CellChat_Endothelial_FDR05_axes']=';'.join(z.loc[z.FDR<.05,'feature'].astype(str))
 row['CellChat_Endothelial_positive_effect_axes']=';'.join(z.loc[z.difference_median>0,'feature'].astype(str))
 row['CellChat_Endothelial_negative_effect_axes']=';'.join(z.loc[z.difference_median<0,'feature'].astype(str))
 row['CellChat_Endothelial_n_patients']=z.n_patients.max() if len(z) else np.nan
 row['primary_association_depth_sensitive']=bool(row.get('pseudobulk_adjusted_p_FDR_focus',1)<.05 and not row.get('depth_adjusted_pseudobulk_FDR_focus',1)<.05)
 row['limitations']='NicheNet target-program prior, not independent ligand-target validation; spatial mixtures and C5 localization NO-GO; no causal interpretation'
 if row.get('pseudobulk_adjusted_p_FDR_focus',1)<.05 and row.get('within_adjusted_r_FDR_focus',1)<.05:row['association_evidence']='Both within-patient and between-patient association; inspect signs and LODO'
 elif row.get('pseudobulk_adjusted_p_FDR_focus',1)<.05 or row.get('within_adjusted_r_FDR_focus',1)<.05:row['association_evidence']='One association level supported; other level not established'
 else:row['association_evidence']='No FDR-supported focused LLPS association'
 rows.append(row)
pd.DataFrame(rows).to_csv(D/'ligand_evidence_wide.csv',index=False)
cc[cc.ligand.isin(focus)].to_csv(D/'all_focused_Endothelial_LR_changes.csv',index=False)
for name,tab in [('historical_NicheNet_first',old),('historical_NicheNet_nominal',new),('historical_NicheNet_LOPO',lopo),('historical_C5_specificity',spec),('historical_spatial_summary',spatial),('historical_new_axes_evidence',axes)]:
 col='test_ligand' if 'test_ligand' in tab else 'gene' if 'gene' in tab else 'ligand';tab[tab[col].isin(focus)].to_csv(D/(name+'.csv'),index=False)
pd.DataFrame(source).to_csv(D/'source_provenance.csv',index=False)
print('INTEGRATION COMPLETE')
