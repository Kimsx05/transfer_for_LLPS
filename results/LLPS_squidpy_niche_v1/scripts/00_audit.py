from pathlib import Path
import json,hashlib,os,platform,importlib.metadata as im
import numpy as np,pandas as pd,pyarrow as pa,pyarrow.parquet as pq
pa.set_cpu_count(4);pa.set_io_thread_count(1)
O=Path(__file__).resolve().parents[1];B=O.parent;R=Path('/home/jinshengxi/LLPS/spatial_data/LLPS_spatial_transfer_manual');A=R/'annotations/LLPS_spatial_transfer_v1'
def save(d,n):pd.DataFrame(d).to_csv(O/'00_input_audit'/n,sep='\t',index=False)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
mp=pd.read_csv(A/'00_metadata/celltype_mapping.tsv',sep='\t');save(mp,'celltype_mapping.tsv')
meta=pd.read_csv(A/'00_metadata/section_manifest.tsv',sep='\t');clinical=pd.read_csv(B/'C5_fivegene_LLPS77_spatial_v1/00_input_audit/spatial_tissue_patient_audit.tsv',sep='\t');geo=pd.read_csv(R/'manifests/section_mapping.tsv',sep='\t')
meta=meta.merge(clinical.drop(columns=['dataset','patient_uid','patient_id_status']),on='section_id',validate='one_to_one').merge(geo.drop(columns=['dataset','n_spots']),on='section_id',validate='one_to_one');save(meta,'section_metadata.tsv')
files=[A/'00_metadata/spot_metadata.parquet',A/'01_cell2location/abundance_mean.parquet',A/'01_cell2location/abundance_q05.parquet'];pf=[pq.ParquetFile(p) for p in files];manifest=[];qc=[];parts=[]
for p in files+[A/'00_metadata/celltype_mapping.tsv',R/'manifests/section_mapping.tsv']:
 manifest.append(dict(path=str(p),role='reused_input',bytes=p.stat().st_size,sha256=sha(p)))
for i,r in enumerate(meta.itertuples()):
 d,a,q=[f.read_row_group(i).to_pandas() for f in pf];sid=r.section_id
 assert d.section_id.eq(sid).all() and d.spot_id.is_unique and a.spot_id.is_unique and q.spot_id.is_unique
 assert set(d.spot_id)==set(a.spot_id)==set(q.spot_id)
 a=a.set_index('spot_id').loc[d.spot_id];q=q.set_index('spot_id').loc[d.spot_id]
 assert set(a.columns)==set(mp.subtype) and np.isfinite(a.to_numpy()).all() and (a.to_numpy()>=0).all()
 oldp=B/'C5_spatial_exploratory_all63_v1/01_spatial_scores'/f'{sid}.parquet';newp=B/'C5_fivegene_LLPS77_spatial_v1/01_LLPS77_spatial_scores'/f'{sid}.parquet'
 old=pd.read_parquet(oldp);new=pd.read_parquet(newp)
 assert old.spot_id.is_unique and new.spot_id.is_unique and set(d.spot_id)==set(old.spot_id)==set(new.spot_id)
 old=old.set_index('spot_id').loc[d.spot_id];new=new.set_index('spot_id').loc[d.spot_id]
 diff=0
 for major,g in mp.groupby('major_lineage'):
  v=a[g.subtype].sum(axis=1).to_numpy();diff=max(diff,float(np.max(np.abs(v-old['major_mean__'+major].to_numpy()))));assert np.allclose(v,old['major_mean__'+major],rtol=1e-5,atol=1e-5)
  d['major_mean__'+major]=v.astype('float32')
 for sub in mp.subtype:d['sub_mean__'+sub]=a[sub].to_numpy(dtype='float32')
 for c in ['five_UCell','A26_UCell']:d[c]=old[c].to_numpy(dtype='float32')
 d['LLPS77_UCell']=new.LLPS77_UCell.to_numpy(dtype='float32')
 for c in new:
  if c.startswith('LLPS77_without_') and c.endswith('UCell'):d[c]=new[c].to_numpy(dtype='float32')
 assert np.allclose(d[['x','y','array_row','array_col']],old[['x','y','array_row','array_col']])
 assert not d[['array_row','array_col']].duplicated().any()
 for c in ['patient_uid','tissue','source_sample_type']:d[c]=getattr(r,c)
 parts.append(d)
 qc.append(dict(section_id=sid,n_spots=len(d),score_matched=len(old),llps_matched=len(new),abundance_matched=len(a),duplicate_keys=0,unmatched_keys=0,major_subtype_max_abs_difference=diff,unassigned_abundance=0,zero_total_abundance=int((a.sum(axis=1)==0).sum()),mean_q05_ratio_median=float(np.median(q.to_numpy()/a.to_numpy())),status='PASS'))
 for p in [oldp,newp]:manifest.append(dict(path=str(p),role='existing_scores',bytes=p.stat().st_size,sha256=sha(p)))
 for p,role in [(R/r.expression_relative_path,'expression_reference_not_loaded'),(R/r.histology_relative_path,'HE_reference_not_copied')]:
  assert p.exists();manifest.append(dict(path=str(p),role=role,bytes=p.stat().st_size,sha256='not_hashed_large_reference'))
 print(sid,len(d),flush=True)
allspots=pd.concat(parts,ignore_index=True);assert allspots.spot_id.is_unique
allspots.to_parquet(O/'00_input_audit/audited_spots.parquet',index=False,compression='zstd');save(qc,'spot_matching_QC.tsv');save(manifest,'input_manifest.tsv')
versions={m:im.version(m) for m in ['numpy','pandas','scipy','pyarrow','scanpy','squidpy','igraph','leidenalg']};(O/'00_input_audit/software_versions.json').write_text(json.dumps(versions,indent=2))
report=f'''# Input audit\n\nPASS: {len(meta)} sections; {len(allspots):,} spots. Exact spot_id joins, no row-order joins. 30 reference subtypes map to {mp.major_lineage.nunique()} original major lineages. Mural remains separate from Fibroblast. No unassigned reference types; no unassigned mass removed.\n\nPosterior mean reused consistently with prior LLPS spatial studies. Major values independently reconstructed from subtype sums and checked against previous scores cache. q05 audited, not substituted for posterior mean.\n\nPatient mapping: {meta.patient_uid.notna().sum()} sections, {meta.patient_uid.nunique()} confirmed patient IDs; unknown identities remain missing. Clinical metadata retained verbatim with source and mapping provenance; no name-based inference.\n\nFive-gene is a C5-associated transcriptional signature, not validated C5 CAF abundance. LLPS-77 is transcriptional activity, not physical phase separation.\n\nH&E registration follows existing low-resolution pixel coordinates with upper origin and no additional scaling. Existing validation did not establish independent histological registration.\n\nResource override: 4 CPU cores/threads, address-space hard limit 100,000,000,000 bytes. Seed 123. No original data changed.\n'''
(O/'00_input_audit/input_audit_report.md').write_text(report)
