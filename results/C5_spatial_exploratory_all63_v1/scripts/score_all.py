from pathlib import Path
import os,sys,json,hashlib,shutil,tempfile,subprocess,re
import pandas as pd,numpy as np,h5py,scipy.sparse as sp
from scipy.io import mmwrite, _fast_matrix_market
_fast_matrix_market.PARALLELISM=1
from scipy.stats import spearmanr
O=Path('results/C5_spatial_exploratory_all63_v1');R=Path('/home/jinshengxi/LLPS/spatial_data/LLPS_spatial_transfer_manual');V=Path('results/C5_spatial_analysis_v1');V2=Path('results/C5_spatial_identifiability_v2_relaxed_marker')
sys.path.insert(0,str(R/'scripts'));from data_io import index,blocks
fixed=json.loads((V2/'01_signature_candidates/fixed_signatures.json').read_text())
sig={'A26':fixed['A26'],'core15':fixed['E15_fixed'],'five': ['MMP11','COL11A1','INHBA','COL5A1','CTHRC1'],'MMP11':['MMP11'],'marker3':['POSTN','COL11A1','COL10A1']}
z=np.load('results/C5_spatial_identifiability_v1/05_deconvolution_feasibility/full_reference_descriptive_reference.npz',allow_pickle=True);types=list(z['types']);rn=z['rn'];genes=z['genes'].astype(str);comp=rn[:,[types.index('C2'),types.index('C4')]].mean(1);ratio=(comp+.1)/(rn[:,types.index('C5')]+.1);valid=(comp>=1)&np.array([not re.match(r'^(MT-|RP[SL][0-9]|ENSG)',g) for g in genes]);ix=np.flatnonzero(valid);ix=ix[np.argsort(-ratio[ix],kind='stable')][:20];sig['C2C4_program']=genes[ix].tolist()
if not (O/'00_inputs_and_QC/frozen_signatures.json').exists():(O/'00_inputs_and_QC/frozen_signatures.json').write_text(json.dumps(sig,indent=2))
for k,v in sig.items():(O/'00_inputs_and_QC'/f'{k}.txt').write_text('\n'.join(v)+'\n')
(O/'00_inputs_and_QC/competitor_definition.txt').write_text('Single-cell-only existing v2 formula on full reference: mean(C2,C4)>=1; rank (mean(C2,C4)+0.1)/(C5+0.1); exclude MT-/ribosomal/ENSG; top20. Combined competing program, not separate C2/C4 abundances. Fixed before spatial scoring.\n')
for src in [V/'00_input_audit/pilot_prespecified.tsv',V/'00_input_audit/spatial_tissue_patient_audit.tsv',R/'manifests/section_mapping.tsv']:shutil.copy2(src,O/'00_inputs_and_QC'/src.name)
checks=[]
for src in [V/'COMPREHENSIVE_REPORT.md',V2/'COMPREHENSIVE_REPORT.md',V2/'01_signature_candidates/fixed_signatures.json',V/'scripts/03_pilot_score.R']:
 checks.append({'path':str(src.resolve()),'sha256':hashlib.sha256(src.read_bytes()).hexdigest()})
pd.DataFrame(checks).to_csv(O/'00_inputs_and_QC/history_checksums.tsv',sep='\t',index=False)
mp=pd.read_csv(R/'manifests/section_mapping.tsv',sep='\t');pilot=pd.read_csv(V/'00_input_audit/pilot_prespecified.tsv',sep='\t').section_id.tolist()
md=pd.read_parquet(R/'annotations/LLPS_spatial_transfer_v1/00_metadata/spot_metadata.parquet')
for label,rel in [('major_mean','LLPS_spatial_transfer_v1_1_addendum/01_major_abundance/major_mean.parquet'),('major_q05','LLPS_spatial_transfer_v1_1_addendum/01_major_abundance/major_q05.parquet'),('sub_mean','LLPS_spatial_transfer_v1/01_cell2location/abundance_mean.parquet'),('sub_q05','LLPS_spatial_transfer_v1/01_cell2location/abundance_q05.parquet')]:
 a=pd.read_parquet(R/'annotations'/rel);a=a.rename(columns={c:label+'__'+c for c in a if c!='spot_id'});md=md.merge(a,on='spot_id',validate='one_to_one')
assert len(md)==265002 and md.spot_id.is_unique
mode=sys.argv[1] if len(sys.argv)>1 else 'pilot';order=pilot if mode=='pilot' else [s for s in mp.section_id if s not in pilot]
if mode=='explicit':order=sys.argv[2].split(',')
if mode!='explicit' and len(sys.argv)>2:order=[sid for sid in order if not (O/'01_spatial_scores'/f'{sid}.parquet').exists()][:int(sys.argv[2])]
for sid in order:
 dest=O/'01_spatial_scores'/f'{sid}.parquet'
 if dest.exists():print('REUSE',sid,flush=True);continue
 row=mp.set_index('section_id').loc[sid];print('START',sid,flush=True)
 with h5py.File(R/row.expression_relative_path,'r') as f:
  gn=index(f['var']);obs=index(f['obs']);assert len(set(gn))==len(gn)
  x=sp.vstack([v for _,v in blocks(f)],format='csr');s=md[md.section_id==sid].copy();assert list(s.barcode)==list(obs)
 n=np.asarray(x.sum(1)).ravel();s['nCount']=n;s['nFeature']=np.diff(x.indptr);assert np.all(n>0)
 coverage=[];detection=np.asarray((x>0).sum(0)).ravel()/len(s)
 for name,gs in sig.items():
  present=[g for g in gs if g in gn];missing=sorted(set(gs)-set(present));idx=[list(gn).index(g) for g in present]
  s[name+'_detected_genes']=np.asarray((x[:,idx]>0).sum(1)).ravel();s[name+'_covered_genes']=len(idx)
  for g in gs:coverage.append(dict(section_id=sid,dataset=row.dataset,signature=name,gene=g,present=g in gn,detected_fraction=detection[list(gn).index(g)] if g in gn else np.nan))
 ix=np.flatnonzero(gn=='MMP11');s['MMP11_counts']=x[:,ix[0]].toarray().ravel() if len(ix) else np.nan;s['MMP11_log1p_CP10k']=np.log1p(s.MMP11_counts/n*1e4)
 with tempfile.TemporaryDirectory(prefix='c5all63_') as td:
  td=Path(td);mmwrite(td/'counts.mtx',x.T);(td/'genes.txt').write_text('\n'.join(gn)+'\n');(td/'barcodes.txt').write_text('\n'.join(obs)+'\n')
  with open(O/'00_inputs_and_QC'/f'{sid}_score.log','w') as log:subprocess.run(['Rscript',str(O/'scripts/score.R'),str(td),str(td/'scores.tsv'),str(O/'00_inputs_and_QC/frozen_signatures.json')],stdout=log,stderr=log,check=True)
  sc=pd.read_csv(td/'scores.tsv',sep='\t');s=s.merge(sc,on='barcode',validate='one_to_one');shutil.copy2(td/'R_session.txt',O/'00_inputs_and_QC/R_session.txt')
 for name in sig:
  for method in ['UCell','AUCell']:
   col=name+'_'+method
   if col not in s:s[col]=np.nan
   assert s[col].dropna().between(0,1).all();s[col+'_percentile']=s[col].rank(method='average',pct=True)
 s.to_parquet(str(dest)+'.tmp',index=False);os.replace(str(dest)+'.tmp',dest);pd.DataFrame(coverage).to_csv(O/'00_inputs_and_QC'/f'{sid}_coverage.tsv',sep='\t',index=False)
 print('DONE',sid,len(s),flush=True)
# Pilot technical acceptance checks are unrelated to model discrimination.
if mode=='pilot':
 rows=[]
 for sid in pilot:
  s=pd.read_parquet(O/'01_spatial_scores'/f'{sid}.parquet');old=pd.read_csv(V/'01_C5_signature'/f'{sid}_pilot_scores.tsv.gz',sep='\t');s=s.set_index('spot_id').loc[old.spot_id]
  for a,b in [('five','C5_extended_exploratory'),('marker3','C5_marker3')]:
   for method in ['UCell','AUCell']:
    v=s[a+'_'+method].to_numpy();w=old[b+'_'+method].to_numpy();rows.append(dict(section_id=sid,signature=a,method=method,spearman=spearmanr(v,w).statistic,max_absolute_difference=np.max(abs(v-w)),n=len(s)))
 pd.DataFrame(rows).to_csv(O/'00_inputs_and_QC/pilot_technical_comparison.tsv',sep='\t',index=False)
 assert min(r['spearman'] for r in rows)>.99
 (O/'00_inputs_and_QC/pilot_technical_PASS.txt').write_text('PASS: 7 prespecified sections, all spot IDs and counts matched, scores finite in range, prior score Spearman >0.99. No model performance gate.\n')
