from analyze import *
import html,struct,subprocess
if __name__=='__main__':
 checks=[]
 def check(name,ok,detail=''):checks.append(dict(check=name,passed=bool(ok),detail=str(detail)));assert ok,(name,detail)
 for hfile in ['source_hashes.tsv','context_source_hashes.tsv']:
  h=pd.read_csv(O/'00_input_audit'/hfile,sep='\t')
  for row in h.itertuples():check('input_unchanged',hashlib.sha256(Path(row.path).read_bytes()).hexdigest()==row.sha256,row.path)
 g=pd.read_csv(O/'tables/global_results.tsv',sep='\t');check('global_complete',len(g)==2520);check('primary_252',len(g[g.method=='primary'])==252);check('all_global_OK',g.status.eq('OK').all());check('primary_BH',np.allclose(g[g.method=='primary'].FDR,bh(g[g.method=='primary'].p_two_sided)))
 total=0
 for sid in META.index:
  print('VERIFY',sid,flush=True)
  src=pd.read_parquet(S/'01_spatial_scores'/f'{sid}.parquet');d=pd.read_parquet(O/'03_local_bivariate_moran'/f'{sid}.parquet');check('spot_order',np.array_equal(src.spot_id,d.spot_id),sid);check('Five_UCell_unchanged',np.array_equal(src.five_UCell,d.five_UCell),sid);check('Five_AUCell_unchanged',np.array_equal(src.five_AUCell,d.five_AUCell),sid);total+=len(d)
  v=d[d.analysis_included].reset_index(drop=True);w=weights(graph(v));p=np.concatenate([v[z+'__p_two_sided'].to_numpy() for z in PARTNERS]);q=bh(p).reshape(4,len(v))
  for j,z in enumerate(PARTNERS):
   check('local_BH',np.allclose(q[j],v[z+'__FDR']),sid+' '+z);check('local_native_P_resolution',np.allclose(v[z+'__p_sim']*1000,np.round(v[z+'__p_sim']*1000)),sid+' '+z)
   m=esda.Moran_BV(v.five_UCell,v['major_mean__'+z],w,permutations=0);expect=g[(g.section_id==sid)&(g.partner==z)&(g.method=='primary')].I.iloc[0];check('Moran_package_recheck',np.isclose(m.I,expect,rtol=1e-12),sid+' '+z)
  for kind in ['local','descriptive_quadrants']:
   for suffix in ['.pdf','.png','_panel_only.pdf','_legend_only.pdf']:check('map_exists',(O/'figures'/f'{sid}_{kind}{suffix}').is_file(),sid+kind+suffix)
 check('spot_total',total==265002,total);check('lee_504',len(pd.read_csv(O/'tables/lee_results.tsv',sep='\t'))==504);check('rings_1512',len(pd.read_csv(O/'tables/ring_results.tsv',sep='\t'))==1512)
 t=pd.read_csv(O/'tables/translation_sensitivity.tsv',sep='\t');check('translation_504',len(t)==504)
 # Inspect actual PNG physical resolution chunks, not filename metadata.
 for f in (O/'figures').glob('*.png'):
  b=f.read_bytes();i=8;dpi=None
  while i+12<len(b):
   n=struct.unpack('>I',b[i:i+4])[0];tag=b[i+4:i+8]
   if tag==b'pHYs':xx,yy,unit=struct.unpack('>IIB',b[i+8:i+17]);dpi=xx*.0254 if unit==1 else None;break
   i+=12+n
  check('PNG_300dpi',dpi is not None and abs(dpi-300)<.1,f.name)
 for f in [O/'figures'/f'{PILOT[0]}_local.pdf',O/'figures/Figure_A_global.pdf']:
  txt=subprocess.check_output(['pdfimages','-list',str(f)],text=True);check('vector_PDF_no_raster',len(txt.strip().splitlines())==2,f.name);font=subprocess.check_output(['pdffonts',str(f)],text=True);check('Arial_font','Arial' in font,f.name)
 save(checks,O/'00_input_audit/final_validation.tsv')
 # Objective representative list uses original seven, not current spatial effect.
 pilot=pd.read_csv(O/'00_input_audit/pilot_prespecified.tsv',sep='\t');pilot['map_path']=pilot.section_id.map(lambda s:'figures/'+s+'_local.pdf');save(pilot,O/'figures/Figure_B_representative_index.tsv')
 body=['<!doctype html><meta charset="utf-8"><title>Five-gene spatial colocalization</title><style>body{font:16px Arial;margin:30px}img{max-width:620px}article{display:inline-block;vertical-align:top;width:650px;margin:10px}h2{clear:both}</style><h1>Five-gene-associated spatial colocalization</h1><p>Exploratory, not verified C5 abundance or communication. Full figures retain original coordinates. HH maps use per-section BH across four targets; descriptive quadrants are NOT significance maps.</p><p><a href="COMPREHENSIVE_REPORT.md">Report</a> | <a href="07_cellchat_integration/CellChat_spatial_integration.tsv">CellChat–spatial integration</a></p><h2>Prespecified representatives (Figure B)</h2>']
 for sid in PILOT:body.append(f'<a href="figures/{sid}_local.pdf">{html.escape(sid)}</a><br>')
 for dataset,ms in META.reset_index().sort_values(['dataset','tissue','section_id']).groupby('dataset'):
  body.append('<h2>'+html.escape(dataset)+'</h2>')
  for r in ms.itertuples():
   s=r.section_id;body.append(f'<article><h3>{html.escape(s)} ({r.tissue})</h3><a href="figures/{s}_local.pdf"><img loading="lazy" src="figures/{s}_local.png"></a><p><a href="figures/{s}_local.pdf">Full local FDR map</a> | <a href="figures/{s}_descriptive_quadrants.pdf">Descriptive quadrants</a> | <a href="03_local_bivariate_moran/{s}.parquet">Spot results</a></p></article>')
 (O/'GALLERY.html').write_text('\n'.join(body));save([dict(path=str(f.relative_to(O)),bytes=f.stat().st_size) for f in O.rglob('*') if f.is_file() and '.Rlib' not in str(f)],O/'00_input_audit/output_manifest.tsv');save([dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in (O/'scripts').glob('*') if f.is_file()],O/'00_input_audit/script_hashes.tsv');print('FINAL VALIDATION',len(checks),'PASS')
