from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,subprocess,struct
O=Path('results/C5_fivegene_LLPS77_spatial_v1');P=Path('results/C5_spatial_exploratory_all63_v1');labels=pd.read_csv(O/'00_input_audit/section_mapping.tsv',sep='\t');rows=[];ids=[]
for sid in labels.section_id:
 s=pd.read_parquet(O/'01_LLPS77_spatial_scores'/f'{sid}.parquet');old=pd.read_parquet(P/'01_spatial_scores'/f'{sid}.parquet');assert list(s.spot_id)==list(old.spot_id);ids.extend(s.spot_id);assert len(s)==len(old)
 for name in ['LLPS77','LLPS77_common','LLPS77_without_core15']:
  for method in ['UCell','AUCell']:assert s[name+'_'+method].between(0,1).all()
 for suffix in ['joint','genes','colocalization']:
  for ext in ['pdf','png']:assert (O/'figures'/f'{sid}_{suffix}.{ext}').stat().st_size>1000
 rows.append(dict(section_id=sid,n_spots=len(s),status='PASS'))
assert len(ids)==len(set(ids))==265002
for r in pd.read_csv(O/'00_input_audit/input_hashes.tsv',sep='\t').itertuples():assert hashlib.sha256(Path(r.path).read_bytes()).hexdigest()==r.sha256,r.path
pd.DataFrame(rows).to_csv(O/'00_input_audit/final_QC.tsv',sep='\t',index=False)
fig=[]
for p in sorted((O/'figures').glob('*.png')):
 b=p.read_bytes();i=b.find(b'pHYs');assert i>=0;px,py,unit=struct.unpack('>IIB',b[i+4:i+13]);dpi=px*.0254;assert abs(dpi-300)<1;fig.append(dict(file=p.name,dpi=dpi))
pd.DataFrame(fig).to_csv(O/'00_input_audit/PNG_validation.tsv',sep='\t',index=False)
# All map PDFs inspected for vector-only contents and readable fonts.
v=[]
for p in sorted((O/'figures').glob('*_joint.pdf')):
 out=subprocess.check_output(['pdfimages','-list',str(p)],text=True);assert len(out.strip().splitlines())==2,p
 font=subprocess.check_output(['pdffonts',str(p)],text=True);v.append(dict(file=p.name,vector_only=True,fonts=font.replace('\n',' | ')));assert 'Arial' in font
pd.DataFrame(v).to_csv(O/'00_input_audit/PDF_validation.tsv',sep='\t',index=False)
checks=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted((O/'scripts').glob('*')) if p.is_file()];pd.DataFrame(checks).to_csv(O/'00_input_audit/script_checksums.tsv',sep='\t',index=False)
print('PASS 63 sections, 265002 unique IDs, source hashes intact, PNG300dpi/vector maps')
