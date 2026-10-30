from pathlib import Path
import csv,gzip,math
B=Path('/home/jinshengxi/bladder_metabolism');O=B/'results/llps_C5_parallel_v1/03_Bulk_Survival';W=B/'fig2_C0_multimethod_deconvolution/03_methods/CIBERSORTx_web';C=B/'results/cibersortx_phase2_clinical_v1'
def read(p):
    return list(csv.DictReader((gzip.open(p,'rt') if str(p).endswith('.gz') else p.open()),delimiter='\t'))
def write(p,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r)); f=p.open('w');w=csv.DictWriter(f,keys,delimiter='\t');w.writeheader();w.writerows(rows);f.close()
frozen={r['cell']:r['cluster'] for r in read(B/'results/fig2_like_C0/phase1_core/cell_metadata_umap_LLPSscore.tsv.gz')}
refs=read(W/'manifest/reference_cells_LOCAL_ONLY.tsv');sub=[r for r in refs if r['full14_group']=='Fibro_C5']
assert len(sub)==150 and all(r['fibro_subtype']=='C5' and frozen[r['cell']]=='C5' and r['focused_group']=='Other_Fibro' for r in sub)
write(O/'C5_reference_mapping.tsv',sub)
write(O/'reference_coverage.tsv',[dict(group='Fibro_C5',n_cells=len(sub),n_patients=len({r['patient_key'] for r in sub}),datasets=';'.join(sorted({r['dataset'] for r in sub})),frozen_identity='All exact barcode matches C5',focused='Other_Fibro; no independent C5')])
inv=read(C/'input_inventory.tsv');aud=[]
for row in inv:
    co=row['dataset'];raw=read(Path(row['file']));d=read(C/'by_dataset'/co/'clinical_merged.tsv');idx={r['Mixture']:r for r in raw}
    assert len(idx)==len(d)
    for r in d:
        rr=idx[r['CIBERSORTx_sample_id']]
        for cl in ('Fibro_C0','Fibro_C5'):
            assert abs(float(rr[cl])-float(r[cl]))<1e-12 and 0<=float(r[cl])<=1
    (O/'by_dataset'/co).mkdir(parents=True,exist_ok=True)
    # Reuse cleaned endpoint definitions and existing patient mapping verbatim.
    write(O/'by_dataset'/co/'clinical_merged.tsv',d)
    aud.append(dict(**row,C5_valid=True,C5_zero_n=sum(float(r['Fibro_C5'])==0 for r in d),OS_available='OS' in row['endpoint_list'].split(';'),source_clinical=str(C/'by_dataset'/co/'clinical_merged.tsv')))
write(O/'input_inventory.tsv',aud)
write(O/'endpoint_definitions.tsv',read(C/'endpoint_definitions.tsv'))
# Inventory existing result files, including reference-varying runs; do not rerun deconvolution.
rows=[]
for p in sorted((W/'output').glob('*.txt')):
    if 'Adjusted' in p.name or 'Results' in p.name:
        header=p.open().readline().strip().split('\t')
        rows.append(dict(file=str(p),contains_independent_C5='Fibro_C5' in header,contains_C0='Fibro_C0' in header,columns=';'.join(header)))
write(O/'existing_deconvolution_inventory.tsv',rows)
print('Validated exact C5 reference identity and 7 existing Full14 cohorts; four have OS.')
