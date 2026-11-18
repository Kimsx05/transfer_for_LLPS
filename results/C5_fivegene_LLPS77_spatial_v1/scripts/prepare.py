from pathlib import Path
import sys,json,hashlib,shutil,h5py,pandas as pd
O=Path('results/C5_fivegene_LLPS77_spatial_v1');P=Path('results/C5_spatial_exploratory_all63_v1');R=Path('/home/jinshengxi/LLPS/spatial_data/LLPS_spatial_transfer_manual')
for d in ['00_input_audit','01_LLPS77_spatial_scores','02_Fivegene_spatial_distribution','03_Fivegene_controls','04_Fivegene_LLPS77_association','05_spatial_colocalization','06_epithelial_and_TME_neighborhood','07_dataset_heterogeneity','figures','tables','scripts']:(O/d).mkdir(exist_ok=True)
sys.path.insert(0,str(R/'scripts'));from data_io import index
ll=Path('results/C5_spatial_analysis_v1/00_input_audit/LLPS77_exact.txt').read_text().splitlines();old=Path('/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/results/tables/LLPS77_matched_genes.txt').read_text().splitlines();assert len(ll)==len(set(ll))==77 and set(ll)==set(old)
sig=json.loads((P/'00_inputs_and_QC/frozen_signatures.json').read_text());mp=pd.read_csv(R/'manifests/section_mapping.tsv',sep='\t');common=set(ll);rows=[]
for r in mp.itertuples():
 with h5py.File(R/r.expression_relative_path) as f:gs=set(index(f['var']))
 common &= gs;rows.append(dict(section_id=r.section_id,dataset=r.dataset,covered=len(set(ll)&gs),missing=';'.join(sorted(set(ll)-gs))))
sets={'LLPS77':ll,'LLPS77_common':[g for g in ll if g in common]}
for name in ['five','A26','core15']:
 overlap=sorted(set(ll)&set(sig[name]));rows.append(dict(section_id='OVERLAP:'+name,covered=len(overlap),missing=';'.join(overlap)))
 if overlap:sets['LLPS77_without_'+name]=[g for g in ll if g not in sig[name]]
(O/'00_input_audit/scoring_sets.json').write_text(json.dumps(sets,indent=2));pd.DataFrame(rows).to_csv(O/'00_input_audit/coverage_overlap.tsv',sep='\t',index=False)
for name in ['pilot_prespecified.tsv','spatial_tissue_patient_audit.tsv','frozen_signatures.json','section_mapping.tsv']:shutil.copy2(P/'00_inputs_and_QC'/name,O/'00_input_audit'/name)
checks=[]
for p in [P/'COMPREHENSIVE_REPORT.md',P/'scripts/score.R',Path('/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/LLPS_77genes.rdata')]+list((P/'01_spatial_scores').glob('*.parquet')):checks.append(dict(path=str(p.resolve()),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
pd.DataFrame(checks).to_csv(O/'00_input_audit/input_hashes.tsv',sep='\t',index=False)
design_path=O/'00_input_audit/DESIGN.md'
if not design_path.exists():design_path.write_text('''# Frozen exploratory design\nSeed 123; serial BLAS/R, maximum four CPU affinity. Original scores read-only. Primary Five UCell vs LLPS77 UCell; AUCell and A26/core sensitivity. Present features only; maxRank1500; AUCell ceil(5% full section feature universe). Common LLPS set is coverage sensitivity, not calibrated LLPS activity. High coverage >=90%. High masks strictly above 95/90/80th percentile (ties reduce selected fraction).\nRank-residual adjustment: Fibroblast, UMI, nFeature, Epithelial. Matching without replacement with per-covariate percentile caliper0.1, before/after pooled-SD SMD; balance requires >=50% high matched, >=10 pairs, all absolute SMD<=0.1. Two- and three-covariate matching reported. Rings0–3 graph steps; epithelial fraction threshold0.3 (0.2/0.4 sensitivity). HH controls Five-high/LLPS-not-high, matching composition.\nSpatial-null sensitivity: random hex-lattice translations of LLPS map, compare only observed overlap >=50% original spots, exclude zero shift; 199 shifts. Preserves translated local values/geometry on overlap, not full irregular tissue boundary; report overlap and use as exploratory sensitivity, not definitive spatial P. No iid spot-based inference. Patient aggregation within dataset/patient, median and1000 patient bootstrap; unknown patients excluded, section descriptions retained.\nSparse-hotspot diagnosis: positive-score-only correlations and top thresholds, exact Five leave-one-gene-out UCell contributions; do not interpret higher Moran as specificity. Formal C5 NO-GO preserved.\n''')
for name,gs in sets.items():(O/'00_input_audit'/f'{name}.txt').write_text('\n'.join(gs)+'\n')
print('LLPS common',len(common),'sets', {k:len(v) for k,v in sets.items()})
