from pathlib import Path
import pandas as pd,numpy as np
from scipy.stats import rankdata,spearmanr
O=Path('results/C5_fivegene_LLPS77_spatial_v1');P=Path('results/C5_spatial_exploratory_all63_v1');labels=pd.read_csv(O/'00_input_audit/spatial_tissue_patient_audit.tsv',sep='\t');out=[]
def rho(x,y):return float(spearmanr(x,y).statistic) if len(x)>3 and np.std(x)>0 and np.std(y)>0 else np.nan
def residual(y,cs):
 X=np.column_stack([np.ones(len(y))]+[rankdata(c)/len(y) for c in cs]);v=rankdata(y)/len(y);return v-X@np.linalg.lstsq(X,v,rcond=None)[0]
for row in labels.itertuples():
 sid=row.section_id;s=pd.read_parquet(P/'01_spatial_scores'/f'{sid}.parquet');d=pd.read_parquet(O/'01_LLPS77_spatial_scores'/f'{sid}.parquet');assert list(s.spot_id)==list(d.spot_id);info=dict(section_id=sid,dataset=row.dataset,patient_uid=row.patient_uid,tissue=row.tissue);n=len(s);five=s.five_UCell.to_numpy();ll=d.LLPS77_UCell.to_numpy();fib=s.major_mean__Fibroblast.to_numpy();dep=s.nCount.to_numpy();nf=s.nFeature.to_numpy();epi=s.major_mean__Epithelial.to_numpy();cs=[fib,dep,nf,epi]
 for method in ['UCell','AUCell']:
  out.append(dict(**info,analysis='LLPS_original_common_score_agreement',method=method,effect=rho(d['LLPS77_'+method],d['LLPS77_common_'+method]),coverage=d.LLPS77_covered_genes.iloc[0]/77,n_spots=n))
 mask=fib>=np.quantile(fib,.75)
 for name,variables,mask0 in [('Fib_high_adjusted',cs,mask),('C2C4_program_adjusted',cs+[s.C2C4_program_UCell.to_numpy()],np.ones(n,bool)),('all_major_adjusted',[s[c].to_numpy() for c in s if c.startswith('major_mean__')]+[dep,nf],np.ones(n,bool)),('q05_Fib_Epi_adjusted',[s.major_q05__Fibroblast.to_numpy(),s.major_q05__Epithelial.to_numpy(),dep,nf],np.ones(n,bool)),('detection_counts_adjusted',cs+[s.five_detected_genes.to_numpy(),d.LLPS77_detected_genes.to_numpy()],np.ones(n,bool))]:
  out.append(dict(**info,analysis=name,method='UCell',effect=rho(residual(five[mask0],[v[mask0] for v in variables]),residual(ll[mask0],[v[mask0] for v in variables])),coverage=d.LLPS77_covered_genes.iloc[0]/77,n_spots=mask0.sum()))
 for fraction in [.05,.1,.2]:
  h=five>np.quantile(five,1-fraction);fh=fib>=np.quantile(fib,.75);fl=fib<=np.quantile(fib,.25);out.append(dict(**info,analysis='Five_high_fib_distribution',method='UCell',top_fraction=fraction,high_in_Fib_top_quartile=(h&fh).sum()/max(h.sum(),1),high_in_Fib_bottom_quartile=(h&fl).sum()/max(h.sum(),1),n_spots=h.sum()))
df=pd.DataFrame(out);df.to_csv(O/'tables/coverage_composition_sensitivity.tsv',sep='\t',index=False);rng=np.random.default_rng(123);rows=[]
for (analysis,method),g in df.dropna(subset=['effect']).groupby(['analysis','method']):
 p=g.dropna(subset=['patient_uid']).groupby(['dataset','patient_uid']).effect.mean().reset_index();v=p.effect.to_numpy();n=len(v);boot=np.median(rng.choice(v,(1000,n)),axis=1) if n>=3 else np.array([np.nan]);rows.append(dict(analysis=analysis,method=method,median_section=g.effect.median(),n_patients=n,median_patient=np.median(v) if n else np.nan,ci_low=np.quantile(boot,.025),ci_high=np.quantile(boot,.975)))
pd.DataFrame(rows).to_csv(O/'07_dataset_heterogeneity/additional_sensitivity_summary.tsv',sep='\t',index=False)
# LOO associations summarized at patient level; no independent feature optimization.
g=pd.read_csv(O/'tables/gene_contribution.tsv',sep='\t');g=g[g.program=='five'];rows=[]
for gene,d in g.groupby('gene'):
 for metric in ['LOO_LLPS_rho','LOO_LLPS_adjusted']:
  v=d.dropna(subset=['patient_uid']).groupby(['dataset','patient_uid'])[metric].mean().dropna().to_numpy();n=len(v);boot=np.median(rng.choice(v,(1000,n)),axis=1) if n>=3 else np.array([np.nan]);rows.append(dict(gene_removed=gene,metric=metric,n_patients=n,median_patient=np.median(v) if n else np.nan,ci_low=np.quantile(boot,.025),ci_high=np.quantile(boot,.975)))
pd.DataFrame(rows).to_csv(O/'07_dataset_heterogeneity/Five_LOO_patient_associations.tsv',sep='\t',index=False)
print('Supplementary composition, coverage and LOO summaries complete')
