from pathlib import Path
import subprocess,struct,json,pandas as pd
O=Path('results/C5_spatial_exploratory_all63_v1');rows=[]
for directory in ['02_spatial_maps','figures']:
 for p in sorted((O/directory).glob('*.png')):
  with p.open('rb') as f:
   assert f.read(8)==b'\x89PNG\r\n\x1a\n';width=height=dpi=None
   while True:
    head=f.read(8)
    if not head:break
    length,kind=struct.unpack('>I4s',head)
    if kind in [b'IHDR',b'pHYs']:
     v=f.read(length)
     if kind==b'IHDR':width,height=struct.unpack('>II',v[:8])
     elif kind==b'pHYs':x,y,unit=struct.unpack('>IIB',v);dpi=x*.0254 if unit==1 else None
    else:f.seek(length,1)
    f.read(4)
    if kind==b'IEND':break
  assert dpi is not None and abs(dpi-300)<.1,(p,dpi)
  rows.append(dict(path=str(p.relative_to(O)),kind='PNG',status='PASS',width=width,height=height,dpi=dpi))
 for p in sorted((O/directory).glob('*.pdf')):
  info=subprocess.run(['pdfinfo',str(p)],capture_output=True,text=True,check=True).stdout
  assert 'Pages:' in info
  fonts=subprocess.run(['pdffonts',str(p)],capture_output=True,text=True,check=True).stdout.splitlines()[2:]
  names=[x.split()[0] for x in fonts if x.strip()];assert all('Arial' in x for x in names),(p,names)
  if p.name.endswith('_overview.pdf'):
   images=subprocess.run(['pdfimages','-list',str(p)],capture_output=True,text=True,check=True).stdout.splitlines();assert len(images)==2,(p,'unexpected raster in primary map')
  rows.append(dict(path=str(p.relative_to(O)),kind='PDF',status='PASS',fonts=';'.join(names)))
pd.DataFrame(rows).to_csv(O/'00_inputs_and_QC/figure_validation.tsv',sep='\t',index=False)
print('Validated',len(rows),'figures',flush=True)
