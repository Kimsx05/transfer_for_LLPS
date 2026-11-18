from pathlib import Path
import subprocess,time,pandas as pd
O=Path('results/C5_fivegene_LLPS77_spatial_v1');ids=pd.read_csv(O/'00_input_audit/section_mapping.tsv',sep='\t').section_id.tolist()
while not all((O/'01_LLPS77_spatial_scores'/f'{sid}.parquet').exists() for sid in ids):time.sleep(5)
subprocess.run(['bash',str(O/'scripts/run_python.sh'),str(O/'scripts/prepare_plot_data.py')],check=True)
jobs=[]
for k in range(2):
 log=open(O/'00_input_audit'/f'maps{k}.log','w');p=subprocess.Popen(['Rscript',str(O/'scripts/maps.R'),','.join(ids[k::2])],stdout=log,stderr=log);jobs.append((p,log))
while jobs:
 for p,log in jobs.copy():
  if p.poll() is not None:
   log.close();jobs.remove((p,log));assert p.returncode==0,p.returncode
 time.sleep(2)
print('All63 maps exported',flush=True)
