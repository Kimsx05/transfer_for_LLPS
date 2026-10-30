from pathlib import Path
from PIL import Image
import csv, subprocess
root=Path(__file__).resolve().parents[1]
rows=[]
for p in sorted(root.rglob('*.pdf')):
    png=p.with_suffix('.png')
    assert png.exists()
    im=Image.open(png)
    assert im.mode=='RGBA' and im.getchannel('A').getextrema()==(0,255)
    rasters=subprocess.check_output(['pdfimages','-list',str(p)],text=True).splitlines()[2:]
    assert not rasters, (p,rasters)
    assert 'Pages:           1' in subprocess.check_output(['pdfinfo',str(p)],text=True)
    rows.append(dict(path=str(p.relative_to(root)),bytes=p.stat().st_size,png_pixels=f'{im.width}x{im.height}',transparent=True,PDF_rasters=0,status='PASS'))
assert len(rows)>0
with (root/'figure_QC.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
print(f'{len(rows)} PDF/PNG figure pairs passed transparency/vector/size checks')
