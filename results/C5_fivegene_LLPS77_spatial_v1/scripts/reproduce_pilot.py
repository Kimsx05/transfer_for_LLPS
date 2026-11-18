from pathlib import Path
import tempfile,subprocess,pandas as pd,numpy as np,json
O=Path('results/C5_fivegene_LLPS77_spatial_v1');sid='GSE171351_A1';p=O/'01_LLPS77_spatial_scores'/f'{sid}.parquet';old=pd.read_parquet(p)
with tempfile.TemporaryDirectory(dir=O/'00_input_audit') as td:
 backup=Path(td)/p.name;p.rename(backup)
 try:subprocess.run(['bash',str(O/'scripts/run_python.sh'),str(O/'scripts/score_new.py'),sid],check=True,stdout=open(O/'00_input_audit/reproduce_pilot.log','w'),stderr=subprocess.STDOUT)
 except BaseException:
  if not p.exists():backup.rename(p)
  raise
 new=pd.read_parquet(p);assert list(old.spot_id)==list(new.spot_id);rows=[]
 for col in old:
  if col.endswith(('_UCell','_AUCell')):
   err=float(np.max(abs(old[col]-new[col])));rows.append(dict(section_id=sid,score=col,max_absolute_difference=err));assert err<1e-12
 pd.DataFrame(rows).to_csv(O/'00_input_audit/independent_rerun_reproducibility.tsv',sep='\t',index=False)
print('Independent complete rerun identical')
