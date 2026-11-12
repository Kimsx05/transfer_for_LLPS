#!/usr/bin/env python3
"""Read-only source audit, per-section export; --validate-only needs only portable package.
Dependencies: numpy pandas pyarrow h5py. No cell2location, torch, R or AnnData needed.
"""
import os
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
import argparse,hashlib,json,shutil,sys
from pathlib import Path
from urllib.parse import quote
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import h5py
pa.set_cpu_count(1);pa.set_io_thread_count(1)
def elem(g):
    if isinstance(g,h5py.Dataset):
        return g.asstr()[()] if h5py.check_string_dtype(g.dtype) else g[()]
    if g.attrs.get('encoding-type')=='categorical':
        c=elem(g['codes']);v=elem(g['categories']);return np.array([v[i] if i>=0 else None for i in c])
    if 'values' in g and 'mask' in g:
        v=elem(g['values']).astype(object);v[elem(g['mask'])]=None;return v
    raise ValueError(g.name)
def obs(f):
    g=f['obs']; idx=elem(g[g.attrs.get('_index','_index')]);return pd.DataFrame({k:elem(g[k]) for k in g if k!=g.attrs.get('_index','_index')},index=idx)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()
def tsv(df,p):df.to_csv(p,sep='\t',index=False)
def validate(out):
    sm=pd.read_csv(out/'00_metadata/section_manifest.tsv',sep='\t');mp=pd.read_csv(out/'00_metadata/celltype_mapping.tsv',sep='\t');labels=mp.subtype.tolist()
    files=[pq.ParquetFile(out/p) for p in ['00_metadata/spot_metadata.parquet','01_cell2location/abundance_mean.parquet','01_cell2location/abundance_q05.parquet']]
    assert len(labels)==len(set(labels))==30
    assert all(f.num_row_groups==len(sm) for f in files)
    seen=set();rows=[]
    for i,r in enumerate(sm.itertuples()):
        m,a,b=[f.read_row_group(i).to_pandas() for f in files];ids=m.spot_id.tolist()
        assert len(m)==len(a)==len(b)==r.n_spots
        assert ids==a.spot_id.tolist()==b.spot_id.tolist() and len(set(ids))==len(ids) and not seen.intersection(ids)
        seen.update(ids)
        assert m.section_id.eq(r.section_id).all() and m.dataset.eq(r.dataset).all()
        assert np.isfinite(m[['x','y','array_row','array_col']].to_numpy()).all()
        for x in (a,b):
            assert x.columns.tolist()==['spot_id']+labels
            assert np.isfinite(x[labels].to_numpy()).all() and (x[labels].to_numpy()>=0).all()
        rows.append(dict(section_id=r.section_id,n_spots=len(m),portable_validation='PASS'))
    return len(seen),rows
ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path);ap.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--validate-only',action='store_true');args=ap.parse_args();O=args.output.resolve()
if args.validate_only:
    n,rows=validate(O)
    inv=pd.read_csv(O/'03_validation/file_manifest.tsv',sep='\t')
    for r in inv.itertuples():assert sha(O/r.relative_path)==r.sha256, r.relative_path
    print(f'PASS: {len(rows)} sections, {n} spots; all package checksums match');sys.exit(0)
if args.project is None:ap.error('--project is required for export')
P=args.project.resolve();R=P/'cell2location_full_run/v2_30type';G=P/'cell2location_final_validation_gate_v1';E=P/'spatial_ecology_top20_v2/00_inputs';A=P/'neighborhood_v1/00_inputs_and_audit'
assert O!=P
if (O/'01_cell2location/abundance_mean.parquet').exists():raise RuntimeError('Refusing overwrite; choose a new --output')
free=shutil.disk_usage(O.parent).free
if free<2*1024**3:raise RuntimeError('Require >=2 GiB free for safe lightweight export')
for d in ['00_metadata/evidence','01_cell2location','02_spatial_inputs','03_validation','scripts']:(O/d).mkdir(parents=True,exist_ok=True)
if Path(__file__).resolve()!=O/'scripts/export_spatial_annotations.py':shutil.copyfile(__file__,O/'scripts/export_spatial_annotations.py')
manifest=pd.read_csv(R/'01_manifest/01_spatial_manifest_frozen.csv').fillna('');cfg=json.loads((R/'00_config/frozen_run_config.json').read_text());labels=pd.read_csv(R/'00_config/final_30_labels.csv').iloc[:,0].tolist()
assert manifest.section_id.is_unique and manifest.eligible_for_cell2location.all()
gate=pd.read_csv(G/'00_final_run_manifest.csv').set_index('section');assert set(gate.index)==set(manifest.section_id)
mp=pd.read_csv(G/'05_major_sanity/05_subtype_to_major_mapping.csv').set_index('subtype').loc[labels].reset_index();mp['annotation_version']='v2_30type';tsv(mp,O/'00_metadata/celltype_mapping.tsv')
lat=pd.read_csv(A/'lattice_recovery/03_recovered_coordinates_by_section.csv.gz');assert not lat.duplicated(['section','barcode']).any()
hist=pd.read_csv(E/'histology_manifest.tsv',sep='\t').set_index('section_id')
identity=pd.read_csv(P/'spatial_ecology_top20_v2/05_07_shared_K12/v1/section_identity_overlay.tsv',sep='\t').fillna('').set_index('section_id')
evidence=[R/'00_config/frozen_run_config.json',R/'00_config/FROZEN_RUN_CONFIG.md',R/'code/01_run_one_section.py',G/'CELL2LOCATION_FINAL_VALIDATION_GATE_REPORT.md',G/'09_final_status/DOWNSTREAM_USAGE_TABLE.csv',G/'00_final_run_manifest.csv',A/'qc/phase1_compute_provenance.csv',A/'code/01_recover_exact_major_q05.py',P/'cell2location_preflight_v1/code/preflight/04_prepare_visium_inputs.R',P/'cell2location_setup_v1/requirements_frozen_20260908.txt']
registry=[]
def register(p,role):
    registry.append(dict(path=str(p),role=role,bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,sha256=sha(p)))
for p in evidence:
    shutil.copyfile(p,O/'00_metadata/evidence'/p.name);register(p,'provenance_evidence')
for p in [R/'01_manifest/01_spatial_manifest_frozen.csv',A/'lattice_recovery/03_recovered_coordinates_by_section.csv.gz',E/'histology_manifest.tsv',P/'spatial_ecology_top20_v2/05_07_shared_K12/v1/section_identity_overlay.tsv',G/'05_major_sanity/05_subtype_to_major_mapping.csv']:register(p,'metadata_source')
sections=[];inputs=[];checks=[];writers={};source_bytes=0
for r in manifest.itertuples():
    sec=r.section_id;d=R/'02_sections'/sec;raw=Path(r.raw_input_path);post=d/'spatial_posterior.h5ad';print(sec,flush=True)
    assert (d/'COMPLETED_AND_VALID').is_file()
    mc=json.loads((d/'model_config.json').read_text());assert mc['labels']==labels
    for k in ['N_cells_per_location','detection_alpha','max_epochs','batch_size','train_size','seed','posterior_num_samples','posterior_batch_size_max']:assert mc[k]==cfg[k]
    assert mc['lr']==cfg['learning_rate']
    with h5py.File(raw,'r') as f,h5py.File(post,'r') as pf:
        ob=obs(f);po=obs(pf);bar=ob.index.astype(str).tolist();xy=f['obsm/spatial'][()]
        assert bar==po.index.astype(str).tolist() and len(set(bar))==len(bar)
        assert ob.dataset.eq(r.dataset).all() and ob.section_id.eq(sec).all()
        assert np.array_equal(xy,pf['obsm/spatial'][()]) and np.isfinite(xy).all()
        ids=[quote(r.dataset,safe='')+'::'+quote(sec,safe='')+'::'+quote(b,safe='') for b in bar]
        sublat=lat[lat.section.eq(sec)].set_index('barcode');assert set(sublat.index)==set(bar) and sublat.dataset.eq(r.dataset).all()
        im=hist.loc[sec];impath=Path(im.image_path);assert impath.is_file()
        prov=f['uns/input_provenance'];sf={k:float(elem(v)) for k,v in prov['scale_factors'].items()};coordinates=list(elem(prov['coordinate_columns']))
        md=pd.DataFrame(dict(spot_id=ids,dataset=r.dataset,section_id=sec,sample_id=r.specimen,barcode=bar,source_obs_id=bar,array_row=sublat.loc[bar,'array_row'].to_numpy(),array_col=sublat.loc[bar,'array_col'].to_numpy(),x=xy[:,0],y=xy[:,1],coordinate_system=str(im.coordinate_space),in_tissue=ob.in_tissue.to_numpy(),in_tissue_provenance='preflight assigned TRUE to hard-QC retained spots; not original Space Ranger flag',source_sample_id=identity.loc[sec,'source_sample_id'],patient_uid=identity.loc[sec,'patient_uid']))
        frames={'metadata':md};maxdiff={}
        for metric,key,prefix in [('mean','means_cell_abundance_w_sf','meanscell_abundance_w_sf_'),('q05','q05_cell_abundance_w_sf','q05cell_abundance_w_sf_')]:
            src=d/f'abundance_{metric}.csv.gz';tab=pd.read_csv(src,index_col=0);assert tab.index.tolist()==bar and tab.columns.tolist()==labels
            v=tab.to_numpy();assert np.isfinite(v).all() and (v>=0).all()
            g=pf['obsm'][key];pbar=elem(g[g.attrs.get('_index','_index')]);assert list(pbar)==bar
            orig=np.column_stack([elem(g[(prefix+l).replace('/','__SLASH__')]) for l in labels]);diff=float(np.max(np.abs(v-orig)));assert np.allclose(v,orig,rtol=1e-6,atol=1e-7)
            maxdiff[metric]=diff;tab=tab.reset_index(drop=True);tab.insert(0,'spot_id',ids);frames[metric]=tab;source_bytes+=src.stat().st_size;register(src,'abundance_'+metric)
        # Stream sparse values; do not load the expression matrix.
        x=f['X'];assert x.attrs['encoding-type'] in ('csr_matrix','csc_matrix');counts=x['data'];count_ok=True
        for start in range(0,len(counts),1000000):
            z=counts[start:start+1000000];count_ok &= bool(np.isfinite(z).all() and (z>=0).all() and np.equal(z,np.floor(z)).all())
        assert count_ok
        genes=elem(f['var'][f['var'].attrs.get('_index','_index')]);ens=sum(str(g).startswith('ENSG') for g in genes)
        assert len(set(genes))==len(genes)
        inputs.append(dict(section_id=sec,dataset=r.dataset,expression_path=str(raw),format='H5AD sparse X',expression_bytes=raw.stat().st_size,n_genes=len(genes),gene_identifier_type='original exact feature names; gene_symbol column; mixed if Ensembl_count>0',Ensembl_count=ens,raw_counts_available=True,counts_validation='finite nonnegative integer values',coordinate_path=str(raw)+'::obsm/spatial',coordinate_columns=json.dumps(coordinates),coordinate_system=im.coordinate_space,scale_factors_path=str(raw)+'::uns/input_provenance/scale_factors',scale_factors_json=json.dumps(sf),image_path=str(impath),image_bytes=impath.stat().st_size,image_scope='lowres HE PNG; optional for scoring/neighborhoods',source_seurat_path=str(P/'ST/all_ST_combined_list_SeuratObject_updated.qs'),posterior_path=str(post),posterior_bytes=post.stat().st_size,match='PASS exact obs names/order and coordinates',migration_required='counts for expression scores/DE/NicheNet; coordinates already packaged; image optional',minimal_transfer='copy per-section raw H5AD; omit posterior/model/reference/combined Seurat; optionally repack sparse counts+features+barcodes'))
    for k,frame in frames.items():
        rel={'metadata':'00_metadata/spot_metadata.parquet','mean':'01_cell2location/abundance_mean.parquet','q05':'01_cell2location/abundance_q05.parquet'}[k]
        table=pa.Table.from_pandas(frame,preserve_index=False)
        if k not in writers:writers[k]=pq.ParquetWriter(O/rel,table.schema,compression='zstd')
        writers[k].write_table(table)
        # independent readback after closing below, one row group per section
    assert len(bar)==gate.loc[sec,'n_spots']
    sections.append(dict(dataset=r.dataset,section_id=sec,sample_id=r.specimen,source_sample_id=identity.loc[sec,'source_sample_id'],patient_uid=identity.loc[sec,'patient_uid'],patient_id_status=identity.loc[sec,'patient_id_status'],n_spots=len(bar),raw_input_path=str(raw),annotation_path=str(d),source_qc_status=r.source_qc_status,QC_status='PASS_WITH_WARNING: original subtype identifiability',annotation_version='v2_30type'))
    checks.append(dict(section_id=sec,n_metadata=len(bar),n_mean=len(frames['mean']),n_q05=len(frames['q05']),n_types=len(labels),coordinates='PASS',counts='PASS',mean_vs_posterior_max_abs=maxdiff['mean'],q05_vs_posterior_max_abs=maxdiff['q05'],config='PASS',original_gate_spots='PASS'))
for w in writers.values():w.close()
tsv(pd.DataFrame(sections),O/'00_metadata/section_manifest.tsv');tsv(pd.DataFrame(inputs),O/'02_spatial_inputs/spatial_input_manifest.tsv')
n,portable=validate(O)
# Compare every Parquet value against source CSV after independent disk readback.
for metric in ['mean','q05']:
    pqf=pq.ParquetFile(O/f'01_cell2location/abundance_{metric}.parquet')
    for i,r in enumerate(manifest.itertuples()):
        src=pd.read_csv(R/'02_sections'/r.section_id/f'abundance_{metric}.csv.gz',index_col=0)
        assert np.array_equal(src.to_numpy(),pqf.read_row_group(i).to_pandas()[labels].to_numpy())
for c in checks:c['parquet_vs_csv_exact']='PASS'
tsv(pd.DataFrame(checks),O/'03_validation/section_validation.tsv');tsv(pd.DataFrame(registry),O/'03_validation/source_file_manifest.tsv')
summary=dict(n_sections=len(sections),n_spots=n,expected_sections_hint=63,expected_spots_hint=262520,spot_difference=n-262520,source_abundance_bytes=source_bytes,expression_bytes=sum(x['expression_bytes'] for x in inputs),optional_image_bytes=sum(x['image_bytes'] for x in inputs),posterior_bytes=sum(x['posterior_bytes'] for x in inputs),free_bytes_before=free,free_bytes_after=shutil.disk_usage(O).free,export_float_conversion_error=0,python=sys.version,dependencies={k:__import__(k).__version__ for k in ['numpy','pandas','pyarrow','h5py']},threads=1)
reference_checks={}
for key in ['reference_input_path','reference_model_path','signature_matrix_path']:
    rp=Path(cfg[key]);reference_checks[key]={'exists':rp.is_file(),'bytes':rp.stat().st_size if rp.is_file() else None}
    if key=='reference_input_path' and rp.is_file():
        with h5py.File(rp,'r') as rf:reference_checks[key]['shape']=[int(x) for x in rf['X'].attrs['shape']]
summary['reference_checks']=reference_checks
(O/'03_validation/audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
import subprocess
report_script=Path(__file__).with_name('finalize_reports.py')
if report_script.resolve() != (O/'scripts/finalize_reports.py').resolve():shutil.copyfile(report_script,O/'scripts/finalize_reports.py')
subprocess.run([sys.executable,str(O/'scripts/finalize_reports.py'),str(O)],check=True)
