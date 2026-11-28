from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
O=Path(__file__).resolve().parents[1];D=O/'08_clinical_association';p=pd.read_csv(D/'section_niche_proportions_complete.tsv',sep='\t');meta=pd.read_csv(O/'00_input_audit/section_metadata.tsv',sep='\t');meta.groupby(['dataset','tissue','source_sample_type'],dropna=False).agg(n_sections=('section_id','nunique'),n_patients=('patient_uid','nunique')).reset_index().to_csv(D/'clinical_sample_sizes.tsv',sep='\t',index=False)
rows=[]
for variable,groups in [('tissue',['Normal','Tumor']),('source_sample_type',['NMIBC','MIBC'])]:
 for dataset,g in p.groupby('dataset'):
  for niche,h in g.groupby('niche'):
   a=h[h[variable]==groups[0]];b=h[h[variable]==groups[1]];n1=len(a);n2=len(b)
   # No across-dataset case/control comparisons; eligible formal test only with >=5 verified patients each.
   aa=a[a.patient_uid.notna()].groupby('patient_uid').proportion.mean();bb=b[b.patient_uid.notna()].groupby('patient_uid').proportion.mean();paired=aa.index.intersection(bb.index);valid=len(aa)>=5 and len(bb)>=5 and not len(paired)
   pv=stats.mannwhitneyu(aa,bb,alternative='two-sided').pvalue if valid else np.nan
   rows.append(dict(variable=variable,dataset=dataset,niche=niche,group_a=groups[0],group_b=groups[1],n_sections_a=n1,n_sections_b=n2,n_patients_a=len(aa),n_patients_b=len(bb),n_paired_patients=len(paired),mean_a=a.proportion.mean(),mean_b=b.proportion.mean(),difference=b.proportion.mean()-a.proportion.mean(),p=pv,status='exploratory_unpaired_patient_test' if valid else 'descriptive_only_insufficient_verified_independent_groups'))
r=pd.DataFrame(rows);r['FDR']=np.nan
for _,g in r.groupby(['variable','dataset']):
 idx=g.index[g.p.notna()]
 if len(idx):r.loc[idx,'FDR']=multipletests(r.loc[idx,'p'],method='fdr_bh')[1]
r.to_csv(D/'clinical_associations.tsv',sep='\t',index=False)
p[p.patient_uid.notna()].groupby(['dataset','patient_uid','niche']).agg(proportion=('proportion','mean'),present=('n_spots',lambda x:bool((x>0).any()))).reset_index().to_csv(D/'patient_niche_proportions.tsv',sep='\t',index=False)
p.assign(present=p.n_spots>0).groupby(['dataset','tissue','niche'],dropna=False).agg(mean_proportion=('proportion','mean'),presence_fraction=('present','mean'),n_sections=('section_id','size')).reset_index().to_csv(D/'clinical_descriptive_summary.tsv',sep='\t',index=False)
(D/'CLINICAL_LIMITATIONS.md').write_text('Only verified original tissue and exact source_sample_type labels used. No inferred T stage or patient names. Comparisons stratified by dataset to avoid source-confounded pooled disease tests. Formal exploratory unpaired patient tests require >=5 verified independent patients in each group and no paired identities; otherwise descriptive NA. Proportions sum to1; single-niche differences are exploratory, not independent outcomes. Insufficient verified groups preclude robust compositional regression and T-stage modeling. Absence is explicitly0 for proportions; absent niche scores remain NA.\n')
