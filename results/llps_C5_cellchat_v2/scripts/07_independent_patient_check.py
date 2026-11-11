import csv,gzip,statistics,collections,math
from pathlib import Path
O=Path('/home/jinshengxi/bladder_metabolism/results/llps_C5_cellchat_v2')
manifest=list(csv.DictReader((O/'00_audit/manifest.csv').open()))
for t,folder in [(10,'01_threshold10'),(20,'02_threshold20'),(30,'03_threshold30')]:
 elig=collections.defaultdict(dict);positive=collections.defaultdict(dict);annot={}
 for r in manifest:
  d=O/folder/r['run_id']
  for z in csv.DictReader((d/'network_all_14_groups.csv').open()):
   if z['clinical_group']=='Normal' or z['sender']==z['receiver'] or 'Fibro_C5' not in (z['sender'],z['receiver']):continue
   if z['eligible']=='TRUE':elig[z['sender'],z['receiver']][z['specimen_id']]=z
  for z in csv.DictReader((d/'LR_detected_all_groups.csv').open()):
   if z['clinical_group']=='Normal' or z['sender']==z['receiver'] or 'Fibro_C5' not in (z['sender'],z['receiver']):continue
   key=z['sender'],z['receiver'],z['interaction_id'];positive[key][z['specimen_id']]=float(z['weight']);annot[key]=z
 rows=[]
 for key,obs in positive.items():
  ne=elig[key[:2]];units=collections.defaultdict(list)
  for sp,r in ne.items():units[r['patient_id'],r['clinical_group']].append(obs.get(sp,0.0))
  patients={r['patient_id'] for r in ne.values()};posp={ne[sp]['patient_id'] for sp in obs};pd=collections.defaultdict(set)
  for sp in obs:pd[ne[sp]['dataset']].add(ne[sp]['patient_id'])
  a=annot[key];rows.append(dict(min_cells=t,sender=key[0],receiver=key[1],feature=key[2],ligand=a['ligand'],receptor=a['receptor'],pathway=a['pathway'],annotation=a['annotation'],eligible_specimens=len(ne),eligible_patients=len(patients),positive_specimens=len(obs),positive_patients=len(posp),positive_datasets=len(pd),patient_prevalence=len(posp)/len(patients),patient_condition_median=statistics.median(statistics.median(v) for v in units.values()),dataset_support=';'.join(f'{d}:{len(s)}' for d,s in sorted(pd.items()))))
 rows.sort(key=lambda r:(-r['positive_patients'],-r['positive_datasets'],-r['patient_prevalence'],-r['patient_condition_median']))
 with (O/f'00_audit/independent_LR_patient_counts_t{t}.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 summary=O/f'04_patient_level_summary/C5_LR_summary_t{t}.csv.gz'
 if summary.exists():
  lookup={(r['sender'],r['receiver'],r['feature']):r for r in rows};n=0
  for z in csv.DictReader(gzip.open(summary,'rt')):
   if z['stratification']!='tumor_only' or int(z['positive_patients'])==0:continue
   a=lookup[z['sender'],z['receiver'],z['feature']]
   for c in ['eligible_specimens','eligible_patients','positive_specimens','positive_patients','positive_datasets']:assert int(z[c])==a[c],(t,c,z,a)
   for c in ['patient_prevalence','patient_condition_median']:assert math.isclose(float(z[c]),a[c],rel_tol=1e-10,abs_tol=1e-12),(t,c,z,a)
   n+=1
  assert n==len(rows),(n,len(rows));print('Independent patient counts/medians PASS',t,n)
 if t==20:
  print('TOP t20 tumor LR by independent patient support')
  for r in rows[:12]:print({k:r[k] for k in ['sender','receiver','ligand','receptor','pathway','positive_patients','eligible_patients','patient_condition_median','dataset_support']})
