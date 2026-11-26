from pathlib import Path
import csv,hashlib,zipfile,json
O=Path('/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1')
D=O/'delivery';D.mkdir(exist_ok=True)
(D/'UPLOAD_STATUS.md').write_text('''# Upload status

Local analysis and QC are complete. GitHub upload is NOT complete.
Target: Kimsx05/transfer_for_LLPS/results/llps_C5_highlow_cellchat_ligand_v1/.
Both Git Data blob and Contents create-file APIs returned HTTP 403 (Resource not accessible by integration). No remote files were created.
Automatic approval review rejected inspecting host cached Git credentials; no further credential access was attempted. Restore the connector repository write permission before uploading.

The ZIP contains reports, gallery, all figures and plotting data, scripts, audits, full statistics across all methods/thresholds, ligand and sensitivity statistics, and primary within_median_t10 specimen/patient paired data. Original expression matrices, Seurat/CellChat RDS caches, invalid provisional networks, and redundant sensitivity per-specimen/per-patient dense pair exports are excluded. Full local results remain intact. See upload_manifest.csv for included SHA256 hashes. Extract at repository root to retain the requested results directory.
''')
files=[]
for p in sorted(O.rglob('*')):
 if not p.is_file():continue
 r=p.relative_to(O)
 if r.parts[0] in ['delivery','02_CellChat_highlow','02_CellChat_highlow_provisional_merged_not_for_analysis']:continue
 if '__pycache__' in r.parts or p.suffix.lower() in ['.rds','.qs','.qs2','.pyc']:continue
 if r.parts[0]=='03_patient_level_comparison' and not (p.name.startswith('within_median_t10_') or p.name.endswith('DONE.csv')):continue
 if p.name=='output_manifest.csv':continue
 files.append(p)
files.append(D/'UPLOAD_STATUS.md')
rows=[{'relative_path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
with (D/'upload_manifest.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['relative_path','bytes','sha256']);w.writeheader();w.writerows(rows)
files.append(D/'upload_manifest.csv')
zpath=D/'llps_C5_highlow_light_results.zip'
with zipfile.ZipFile(zpath,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,'results/'+O.name+'/'+str(p.relative_to(O)))
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
(D/'archive_sha256.txt').write_text(hashlib.sha256(zpath.read_bytes()).hexdigest()+'  '+zpath.name+'\n')
rows=[{'path':str(p),'relative_path':str(p.relative_to(O)),'bytes':p.stat().st_size} for p in sorted(O.rglob('*')) if p.is_file() and p.name!='output_manifest.csv']
with (O/'output_manifest.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['path','relative_path','bytes']);w.writeheader();w.writerows(rows)
print(json.dumps({'archive':str(zpath),'files':len(files),'MiB':round(zpath.stat().st_size/2**20,2),'status':'LOCAL_COMPLETE_UPLOAD_BLOCKED'}))
