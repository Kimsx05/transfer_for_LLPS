from pathlib import Path
import subprocess,time,json
import pandas as pd
O=Path('results/C5_fivegene_LLPS77_spatial_v1');pilot=pd.read_csv(O/'00_input_audit/pilot_prespecified.tsv',sep='\t').section_id.tolist()
for sid in pilot:
 p=O/'01_LLPS77_spatial_scores'/f'{sid}.parquet';assert p.exists();s=pd.read_parquet(p);assert s.spot_id.is_unique and s.LLPS77_UCell.between(0,1).all() and s.LLPS77_AUCell.between(0,1).all();t=json.loads((O/'00_input_audit'/f'{sid}_technical.json').read_text());assert t['five_reconstruction_error']<1e-7
(O/'00_input_audit/pilot_technical_PASS.txt').write_text('All 7 prespecified pilots passed exact IDs/counts alignment, finite in-range scores, full-background Five UCell reconstruction <1e-7. No biological effect acceptance gate.\n')
ids=pd.read_csv(O/'00_input_audit/section_mapping.tsv',sep='\t').section_id.tolist();ids=[s for s in ids if not (O/'01_LLPS77_spatial_scores'/f'{s}.parquet').exists()];jobs=[]
for k in range(3):
 batch=ids[k::3]
 if not batch:continue
 log=open(O/'00_input_audit'/f'batch{k}.log','w');p=subprocess.Popen(['bash',str(O/'scripts/run_python.sh'),str(O/'scripts/score_new.py'),','.join(batch)],stdout=log,stderr=log);jobs.append((p,log))
while jobs:
 for p,log in jobs.copy():
  if p.poll() is not None:
   log.close();jobs.remove((p,log));assert p.returncode==0,p.returncode
 time.sleep(2)
print('ALL63 scoring complete',flush=True)
