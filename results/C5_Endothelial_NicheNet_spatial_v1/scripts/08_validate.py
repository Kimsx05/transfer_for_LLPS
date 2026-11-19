from pathlib import Path
import pandas as pd,numpy as np,hashlib,struct,subprocess
O=Path('results/C5_Endothelial_NicheNet_spatial_v1');rows=[]
def check(name,ok,detail=''):rows.append(dict(check=name,passed=bool(ok),detail=str(detail)))
files=list((O/'07_spatial_expression_validation').glob('*/expression.parquet'));check('63 spatial sections',len(files)==63,len(files));n=0
for f in files:
 d=pd.read_parquet(f);n+=len(d);check(f.parent.name+' unique spots',not d.spot_id.duplicated().any());check(f.parent.name+' done',(f.parent/'DONE').exists())
check('265002 spots',n==265002,n)
a=pd.read_csv(O/'04_NicheNet_ligand_activity/all_program_activities.tsv',sep='\t');check('3 evaluable curated programs',a.program.nunique()==3);check('246 ligands per program',(a.groupby('program').size()==246).all());check('finite activity',np.isfinite(a[['aupr_corrected','aupr','auroc','pearson']]).all().all())
for r in pd.read_csv(O/'00_input_audit/source_checksums.tsv',sep='\t').itertuples():
 h=hashlib.sha256()
 with open(r.path,'rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 check('source unchanged '+r.path,h.hexdigest()==r.sha256)
for f in (O/'figures').glob('*.png'):
 data=f.read_bytes();pos=data.find(b'pHYs');ppm=struct.unpack('>IIB',data[pos+4:pos+13]) if pos>=0 else None;check(f.name+' 300dpi',ppm is not None and abs(ppm[0]*.0254-300)<1,ppm)
for f in (O/'figures').glob('*.pdf'):
 q=subprocess.run(['pdfimages','-list',str(f)],capture_output=True,text=True);check(f.name+' vector',q.returncode==0 and len(q.stdout.strip().splitlines())==2)
check('report',(O/'COMPREHENSIVE_REPORT.md').exists());pd.DataFrame(rows).to_csv(O/'00_input_audit/final_validation.tsv',sep='\t',index=False);print('VALIDATION',sum(r['passed'] for r in rows),len(rows));assert all(r['passed'] for r in rows)
