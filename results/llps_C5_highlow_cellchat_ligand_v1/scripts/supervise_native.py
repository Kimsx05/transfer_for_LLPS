import os,subprocess,time,json,signal
from pathlib import Path
O=Path('/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1')
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS','OMP_THREAD_LIMIT']:os.environ[k]='1'
os.sched_setaffinity(0,{3});jobs=[]
for w in range(3):
 f=open(O/f'02_CellChat_highlow/worker_{w}.log','a');p=subprocess.Popen(['taskset','-c',str(w),'Rscript','-e',"source('"+str(O/'scripts/02_cellchat.R')+"')",str(w)],stdout=f,stderr=subprocess.STDOUT);jobs.append(p)
try:
 while any(p.poll() is None for p in jobs):
  rows=[]
  for p in jobs:
   r={'pid':p.pid,'exit_code':p.poll()}
   if p.poll() is None:
    lines=Path(f'/proc/{p.pid}/status').read_text().splitlines();r.update({s.split(':')[0]:s.split(':')[1].strip() for s in lines if s.startswith(('VmRSS:','Threads:','Cpus_allowed_list:','State:'))});r['cpu_ticks']=sum(map(int,Path(f'/proc/{p.pid}/stat').read_text().split()[13:15]))
   rows.append(r)
  status={'time':time.strftime('%Y-%m-%d %H:%M:%S'),'complete':len(list((O/'02_CellChat_highlow').glob('*/COMPLETE.csv'))),'errors':len(list((O/'02_CellChat_highlow').glob('*/ERROR.csv'))),'workers':rows}
  with open(O/'00_input_audit/native_runtime.jsonl','a') as f:f.write(json.dumps(status)+'\n')
  print(json.dumps(status),flush=True)
  rss=sum(float(r.get('VmRSS','0 kB').split()[0]) for r in rows)/1048576
  if rss>7.4:raise RuntimeError('Aggregate R RSS exceeded 7.4 GiB; checkpoints retained. Reduce concurrency before resuming.')
  if status['errors']:raise RuntimeError('A sample recorded an ERROR.csv; inspect before resuming.')
  time.sleep(30)
finally:
 for p in jobs:
  if p.poll() is None:p.terminate()
print('NATIVE WORKERS FINISHED', [p.returncode for p in jobs],flush=True)
if any(p.returncode for p in jobs):raise SystemExit(1)
