from pathlib import Path
import pandas as pd,json,gzip
O=Path('results/C5_fivegene_LLPS77_spatial_v1');T=O/'00_input_audit/pilot_statistics';C=O/'00_input_audit/statistic_checkpoints';C.mkdir(exist_ok=True);tables={p.stem:pd.read_csv(p,sep='\t') for p in T.glob('*.tsv')};tables['pairs']=pd.read_parquet(T/'matched_spot_pairs.parquet')
for sid in pd.read_csv(O/'00_input_audit/pilot_prespecified.tsv',sep='\t').section_id:
 payload={k:d[d.section_id==sid].to_dict('records') for k,d in tables.items()}
 with gzip.open(C/(sid+'.json.gz'),'wt') as f:json.dump(payload,f,default=lambda v:v.item() if hasattr(v,'item') else str(v))
print('Pilot checkpoints retained')
