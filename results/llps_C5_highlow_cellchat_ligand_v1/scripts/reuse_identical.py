import os,time,csv,gzip,hashlib,shutil,sys
from pathlib import Path
O=Path('/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1');jobs=list(csv.DictReader(open(O/'02_CellChat_highlow/job_manifest.csv')));jobs=[r for r in jobs if r['runnable']=='TRUE'];labels={}
with gzip.open(O/'01_LLPS_stratification/cell_assignments.csv.gz','rt') as f:
 for r in csv.DictReader(f):labels.setdefault((r['specimen_id'],r['method']),[]).append((r['cell'],r['state']))
groups={}
for j in jobs:
 h=hashlib.sha256(repr(sorted(labels[(j['specimen_id'],j['method'])])).encode()).hexdigest();groups.setdefault((j['specimen_id'],h),[]).append(j['job_id'])
records=[]
for key,ids in groups.items():
 for i in ids:records.append(dict(job_id=i,assignment_sha256=key[1],identical_jobs=';'.join(ids)))
p=O/'00_input_audit/identical_assignment_audit.csv'
with p.open('w') as f:w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
if '--audit-only' in sys.argv:sys.exit(0)
# Reuse only exact barcode-state identity, within this new analysis. No old probabilities reused.
while True:
 for ids in groups.values():
  ready=[i for i in ids if (O/'02_CellChat_highlow'/i/'unfiltered.rds').exists()]
  if not ready:continue
  src=O/'02_CellChat_highlow'/ready[0]/'unfiltered.rds'
  for i in ids:
   dest=O/'02_CellChat_highlow'/i
   if dest.exists():continue
   dest.mkdir();tmp=dest/'reuse.tmp';shutil.copyfile(src,tmp);tmp.rename(dest/'unfiltered.rds');(dest/'identical_assignment_reuse.txt').write_text(str(src)+'\nExact barcode/state labels and expression, seed, parameters identical. Threshold exports remain separate.\n');print('REUSE',i,'FROM',ready[0],flush=True)
 if all((O/'02_CellChat_highlow'/j['job_id']/'COMPLETE.csv').exists() or (O/'02_CellChat_highlow'/j['job_id']/'ERROR.csv').exists() for j in jobs):break
 time.sleep(15)
