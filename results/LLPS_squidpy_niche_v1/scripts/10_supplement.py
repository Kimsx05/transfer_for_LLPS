from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,scipy.sparse as sp
from scipy import stats
import anndata as ad,squidpy as sq
O=Path(__file__).resolve().parents[1];d=pd.read_parquet(O/'04_niche_characterization/spot_metadata.parquet');native=pd.read_parquet(O/'03_clustering/native_labels.parquet');d=d.merge(native[['spot_id']+[c for c in native if c not in d]],on='spot_id',validate='one_to_one');d.to_parquet(O/'04_niche_characterization/spot_metadata.parquet',index=False,compression='zstd')
# Validate directed observed adjacency counts against actual Squidpy nhood_enrichment on first section.
sid=d.section_id.iloc[0];g=d[d.section_id==sid];A=sp.load_npz(O/'01_spatial_graph'/f'{sid}_R1.npz');z=ad.AnnData(sp.csr_matrix((len(g),1)));z.obs['niche']=pd.Categorical(g.niche.to_numpy());z.obsp['spatial_connectivities']=A;sq.gr.nhood_enrichment(z,cluster_key='niche',n_perms=1000,seed=123,n_jobs=1,show_progress_bar=False);obs=z.uns['niche_nhood_enrichment']['count'];labels=z.obs.niche.cat.codes.to_numpy();u,v=A.nonzero();K=len(z.obs.niche.cat.categories);ours=np.bincount(labels[u]*K+labels[v],minlength=K*K).reshape(K,K);assert np.array_equal(obs,ours)
(O/'07_neighborhood_enrichment/squidpy_validation.json').write_text(json.dumps(dict(section_id=sid,observed_counts_exact_match=True,n_perms=1000,seed=123,n_jobs=1,zero_variance_z='NA downstream',permutation_rng='Our moments use numpy default_rng; Squidpy uses version-specific RNG; observed counts must match exactly, null moments are Monte Carlo estimates.'),indent=2))
c=pd.read_csv(O/'05_C5_LLPS_enrichment/section_level_contrasts.tsv',sep='\t');keys=['dataset','niche','feature','model','threshold','min_spots'];c.groupby(keys,dropna=False).agg(mean_effect=('effect','mean'),median_effect=('effect','median'),n_sections=('effect','count'),positive_fraction=('effect',lambda x:np.mean(x.dropna()>0) if x.notna().any() else np.nan)).reset_index().to_csv(O/'05_C5_LLPS_enrichment/dataset_enrichment_direction.tsv',sep='\t',index=False)
# Input report presence/provenance audit across all requested historical result directories.
names=['C5_spatial_identifiability_v1','C5_fivegene_LLPS77_spatial_v1','C5_spatial_colocalization_esda_v1','C5_Endothelial_NicheNet_spatial_v1','C5_Endothelial_GSEA_nominalNicheNet_v1'];rows=[]
for name in names:
 path=O.parent/name;reports=[p for p in path.glob('*.md')];rows.append(dict(result=name,path=str(path),exists=path.exists(),reports=';'.join(p.name for p in reports),usage='provenance/context only; no communication or score enters clustering'))
pd.DataFrame(rows).to_csv(O/'00_input_audit/historical_result_audit.tsv',sep='\t',index=False)
# Dominant composition mixture-loss and model associations.
d.groupby('section_id').agg(median_dominant_fraction=('dominant_fraction','median'),fraction_dominant_below_half=('dominant_fraction',lambda x:np.mean(x<.5)),fraction_dominant_below_075=('dominant_fraction',lambda x:np.mean(x<.75))).to_csv(O/'03_clustering/mixing_information_loss.tsv',sep='\t')
# Feature characterization at all niches; both central raw and adjusted scores included.
features=['five_UCell','LLPS77_UCell','A26_UCell','Fibroblast_abundance','Endothelial_abundance','five_UCell__adjusted','LLPS77_UCell__adjusted'];d.groupby(['section_id','niche'])[features].mean().groupby('niche').mean().to_csv(O/'06_C5_endothelial/niche_feature_summary.tsv',sep='\t')
clinical_cols=['section_id','dataset','patient_uid','tissue','source_sample_type','niche']
cs=d.groupby(clinical_cols,dropna=False)[features].mean().reset_index();cs.to_csv(O/'08_clinical_association/clinical_section_niche_scores.tsv',sep='\t',index=False)
cs.groupby(['dataset','tissue','source_sample_type','niche'],dropna=False)[features].agg(['mean','count']).to_csv(O/'08_clinical_association/clinical_score_descriptions.tsv',sep='\t')
from sklearn.metrics import silhouette_score,calinski_harabasz_score,davies_bouldin_score
pilot=np.load(O/'03_clustering/pilot_indices.npy');native_quality=[]
for radius in [1,2,3]:
 X=np.load(O/'02_neighborhood_composition'/f'dominant_R{radius}.npy')[pilot];y=d[f'native_R{radius}'].to_numpy()[pilot];k=len(set(y));native_quality.append(dict(radius=radius,n_niches_full=d[f'native_R{radius}'].nunique(),n_niches_pilot=k,silhouette=silhouette_score(X,y,sample_size=2000,random_state=123),calinski_harabasz=calinski_harabasz_score(X,y),davies_bouldin=davies_bouldin_score(X,y),scope='full pooled native labels evaluated on balanced pilot'))
pd.DataFrame(native_quality).to_csv(O/'03_clustering/native_quality_metrics.tsv',sep='\t',index=False)
import subprocess
subprocess.run([str(O/'scripts/run_python.sh'),str(O/'scripts/12_convergence_diagnostic.py')],check=True)
