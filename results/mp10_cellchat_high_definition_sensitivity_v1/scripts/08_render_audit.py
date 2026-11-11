#!/usr/bin/python3
from pathlib import Path
import csv, json, subprocess, hashlib, re
from PIL import Image, ImageChops, ImageStat, ImageDraw
root=Path(__file__).resolve().parents[1]
def read(p):
 with p.open() as h:return list(csv.DictReader(h,delimiter='\t'))
def write(p,rows):
 with p.open('w') as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
manifest=read(root/'config/figure_manifest.tsv');checks=[]
for row in manifest:
 p=Path(row['file']);assert p.is_file(),p
 rec=dict(figure_id=row['figure_id'],part=row['part'],file=str(p),exists=True,n_bytes=p.stat().st_size,vector_PDF='NA',dpi='NA',render_ok=True)
 if p.suffix=='.png':
  with Image.open(p) as im:
   dpi=im.info.get('dpi',(0,0));assert abs(dpi[0]-300)<1 and abs(dpi[1]-300)<1,(p,dpi)
   assert im.width==round(float(row['width_in'])*300) and im.height==round(float(row['height_in'])*300)
   diff=ImageChops.difference(im.convert('RGB'),Image.new('RGB',im.size,'white'));assert diff.getbbox() is not None
   rec['dpi']=dpi[0]
 else:
  raw=p.read_bytes();assert raw.startswith(b'%PDF');assert b'/Subtype /Image' not in raw,p
  rec['vector_PDF']=True
  result=subprocess.run(['/usr/bin/mutool','draw','-q','-F','txt','-o','-',str(p)],capture_output=True,text=True)
  assert result.returncode==0,(p,result.stderr);assert len(result.stdout.strip())>0,p
 checks.append(rec)
write(root/'qc/render_format_checks.tsv',checks)
# Contact sheets of full PNG outputs aid visual review; retain original figures unchanged.
full=[Path(x['file']) for x in manifest if x['part']=='png'];sheetpaths=[]
for offset in range(0,len(full),12):
 batch=full[offset:offset+12];canvas=Image.new('RGB',(1800,4*510),'white');draw=ImageDraw.Draw(canvas)
 for j,p in enumerate(batch):
  with Image.open(p) as im:
   im=im.convert('RGB');im.thumbnail((580,470));x=(j%3)*600+(600-im.width)//2;y=(j//3)*510+25;canvas.paste(im,(x,y));draw.text(((j%3)*600+8,(j//3)*510+5),p.stem,fill='black')
 dest=root/'qc'/f'contact_sheet_{offset//12+1}.png';canvas.save(dest);sheetpaths.append(str(dest))
# Render main result PDFs independently of saved PNG.
for stem in ['C10_ALL_input_effects','C10_ALL_output_effects','C30_ALL_input_effects','candidate_rank_input','dataset_direction_input']:
 p=root/'figures'/f'{stem}.pdf'
 if p.exists():
  dest=root/'qc'/f'{stem}_PDF_preview.png';result=subprocess.run(['/usr/bin/mutool','draw','-q','-r','100','-o',str(dest),str(p)],capture_output=True,text=True);assert result.returncode==0,(p,result.stderr)
with (root/'FIGURE_INDEX.md').open('w') as h:
 h.write('# Figures\n\nEach full PNG has a vector PDF. Figures with legends also have panel-only and legend-only PDFs. Direct plot data and scales are registered in `config/figure_manifest.tsv`.\n\n')
 for p in full:h.write(f'- [{p.stem}](figures/{p.name})\n')
rows=[]
for p in sorted(root.rglob('*')):
 if p.is_file() and p.suffix not in {'.so','.o'} and 'rcpp_cache' not in str(p) and 'logs' not in p.relative_to(root).parts and p.name not in {'bundle_checksums.tsv','run.lock'}:
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
  rows.append(dict(path=str(p.relative_to(root)),bytes=p.stat().st_size,sha256=h.hexdigest()))
write(root/'config/bundle_checksums.tsv',rows)
print(f'Validated {len(checks)} rendered files for {len(full)} figures. Contact sheets: {len(sheetpaths)}')
