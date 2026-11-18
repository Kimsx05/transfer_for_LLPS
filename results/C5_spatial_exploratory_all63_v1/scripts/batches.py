"""Resume non-overlapping section batches; at most three serial workers within <=4 CPUs."""
import subprocess,sys,time
from pathlib import Path
import pandas as pd
O=Path('results/C5_spatial_exploratory_all63_v1')
ids=pd.read_csv(O/'00_inputs_and_QC/section_mapping.tsv',sep='\t').section_id
todo=[s for s in ids if not (O/'01_spatial_scores'/f'{s}.parquet').exists()]
queue=[todo[i:i+4] for i in range(0,len(todo),4)];running=[]
try:
 while queue or running:
  while queue and len(running)<3:
   chunk=queue.pop(0);log=open(O/'00_inputs_and_QC'/f'batch_{chunk[0]}.log','w');p=subprocess.Popen([sys.executable,str(O/'scripts/score_all.py'),'explicit',','.join(chunk)],stdout=log,stderr=log);running.append((p,log,chunk))
  for p,log,chunk in list(running):
   code=p.poll()
   if code is not None:
    log.close();running.remove((p,log,chunk))
    if code:raise RuntimeError(f'Batch failed {chunk}: {code}; completed outputs remain reusable')
    print('DONE BATCH',chunk,flush=True)
  if running:time.sleep(1)
finally:
 for p,log,chunk in running:
  if p.poll() is None:p.terminate()
  log.close()
