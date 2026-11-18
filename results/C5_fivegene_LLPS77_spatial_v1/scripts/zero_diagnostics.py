from pathlib import Path
import pandas as pd
O=Path('results/C5_fivegene_LLPS77_spatial_v1');P=Path('results/C5_spatial_exploratory_all63_v1');rows=[]
for row in pd.read_csv(O/'00_input_audit/section_mapping.tsv',sep='\t').itertuples():
 s=pd.read_parquet(P/'01_spatial_scores'/f'{row.section_id}.parquet')
 for sig in ['five','A26','core15']:
  for method in ['UCell','AUCell']:
   zero=s[sig+'_'+method]==0;detected=s[sig+'_detected_genes']>0;rows.append(dict(section_id=row.section_id,dataset=row.dataset,signature=sig,method=method,score_zero_fraction=zero.mean(),count_zero_fraction=(~detected).mean(),score_zero_but_counts_positive_fraction=(zero&detected).mean(),detected_gene_median=s[sig+'_detected_genes'].median()))
df=pd.DataFrame(rows);df.to_csv(O/'tables/zero_score_vs_zero_counts.tsv',sep='\t',index=False);print(df.query("method=='UCell'").groupby('signature')[['score_zero_fraction','count_zero_fraction','score_zero_but_counts_positive_fraction','detected_gene_median']].median().to_string())
