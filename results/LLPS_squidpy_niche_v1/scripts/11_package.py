from pathlib import Path
import json,hashlib,zipfile,html
import pandas as pd
O=Path(__file__).resolve().parents[1];D=O/'delivery';D.mkdir(exist_ok=True)
def gallery(paths,title):
 items=[]
 for p in paths:
  rel=p.relative_to(O);items.append(f'<section><h3>{html.escape(str(rel))}</h3><img loading="lazy" src="{html.escape(str(rel))}" style="max-width:900px;width:100%"></section>')
 return '<!doctype html><meta charset="utf-8"><title>'+title+'</title><body style="font-family:Arial;background:white"><h1>'+title+'</h1>'+''.join(items)
allpng=sorted(p for p in (O/'09_figures').rglob('*.png') if not p.is_symlink());(O/'GALLERY.html').write_text(gallery(allpng,'LLPS spatial niche clustering v1: all figures'))
sids=set(pd.read_csv(O/'09_figures/representative_sections.tsv',sep='\t').section_id);selected=[]
for p in sorted(O.rglob('*')):
 if not p.is_file() or p.is_symlink():continue
 rel=p.relative_to(O)
 if any(x.startswith('.') or x=='__pycache__' for x in rel.parts) or rel.parts[0]=='delivery':continue
 if p.name in ['GALLERY.html','DELIVERY_GALLERY.html','storage_summary.json','11_package.log','resources_11_package.json']:continue
 if p.suffix in ['.pyc','.npy','.npz'] or p.name.endswith('_plot_data.tsv.gz') or p.name.endswith('_umap_data.tsv.gz'):continue
 if rel.parts[0]=='09_figures':
  if p.suffix not in ['.png','.pdf','.md','.json','.tsv']:continue
  if 'spatial' in rel.parts and not any(p.name.startswith(sid+'_') for sid in sids):continue
  if 'graph_QC' in rel.parts and not p.name.startswith(('R1_','R2_','R3_')):continue
  if 'sections' in rel.parts:continue
 if p.name=='audited_spots.parquet':continue
 selected.append(p)
(O/'DELIVERY_GALLERY.html').write_text(gallery([p for p in selected if p.suffix=='.png'],'LLPS spatial niches: delivery figure selection'));selected.append(O/'DELIVERY_GALLERY.html')
rows=[dict(path=str(p.relative_to(O)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in selected];pd.DataFrame(rows).to_csv(D/'output_manifest.tsv',sep='\t',index=False)
archive=D/'LLPS_squidpy_niche_v1_lightweight.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in selected:z.write(p,'results/LLPS_squidpy_niche_v1/'+str(p.relative_to(O)))
 z.write(D/'output_manifest.tsv','results/LLPS_squidpy_niche_v1/output_manifest.tsv')
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
(D/'package_sha256.txt').write_text(hashlib.sha256(archive.read_bytes()).hexdigest()+'  '+archive.name+'\n')
size=sum(p.stat().st_size for p in O.rglob('*') if p.is_file() and not p.is_symlink());(O/'10_reports/storage_summary.json').write_text(json.dumps(dict(new_disk_bytes=size,delivery_zip_bytes=archive.stat().st_size,input_expression_copied_bytes=0,HE_source_copied_bytes=0,environment_included_in_total=True,archive_files=len(selected)+1,all_local_pngs=len(allpng),note='Storage summary and this packaging process own logs excluded from archive to avoid self-referential changing checksums.'),indent=2))
print('files',len(selected),'zip_MB',archive.stat().st_size/1e6,'CRC_validation_PASS',flush=True)
