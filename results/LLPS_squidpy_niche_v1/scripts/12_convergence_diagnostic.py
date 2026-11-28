from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
from sklearn.cluster import KMeans
O=Path(__file__).resolve().parents[1];D=O/'03_clustering';rec=pd.read_csv(D/'recommended_models.tsv',sep='\t');pilot=np.load(D/'pilot_indices.npy');original=np.load(D/'pilot_labels.npz',allow_pickle=True);out=[]
for role in ['major_best','subtype_best']:
 p=rec[rec.role==role].iloc[0];x=np.load(O/'02_neighborhood_composition'/f'{p.kind}_R{p.radius}.npy')[pilot]
 if p['transform']=='hellinger':x=np.sqrt(x)
 z=ad.AnnData(x);sc.pp.neighbors(z,n_neighbors=int(p.knn),use_rep='X',random_state=123,transformer='pynndescent');base=None
 for seed in [123,17,321,2026]:
  sc.tl.leiden(z,resolution=p.resolution,random_state=seed,flavor='igraph',n_iterations=-1,key_added='niche');y=z.obs.niche.astype(str).to_numpy()
  if base is None:base=y.copy()
  out.append(dict(model=p.model,diagnostic='Leiden_until_convergence',seed=seed,n_niches=len(set(y)),ARI_vs_two_iterations_seed123=adjusted_rand_score(original[p.model],y),ARI_vs_converged_seed123=adjusted_rand_score(base,y),NMI_vs_converged_seed123=normalized_mutual_info_score(base,y)))
 # K-means seed sensitivity at requested control's chosenK for main feature space only.
 if role=='major_best':
  k=int(pd.read_csv(D/'kmeans_quality.tsv',sep='\t').sort_values('silhouette',ascending=False).iloc[0].k);base=None
  for seed in [123,17,321,2026]:
   y=KMeans(n_clusters=k,n_init=10,random_state=seed).fit_predict(x)
   if base is None:base=y.copy()
   out.append(dict(model=p.model,diagnostic='Kmeans_seed_sensitivity',seed=seed,n_niches=k,ARI_vs_converged_seed123=adjusted_rand_score(base,y),NMI_vs_converged_seed123=normalized_mutual_info_score(base,y)))
pd.DataFrame(out).to_csv(D/'convergence_diagnostic.tsv',sep='\t',index=False)
(D/'CONVERGENCE_DIAGNOSTIC.md').write_text('Added after the initial composition-only seed/balance audit revealed instability, to distinguish incomplete optimization from intrinsic partition ambiguity. No C5/LLPS/clinical data used. Runs Leiden until convergence on unchanged balanced pilot/KNN and compares seeds; K-means seed variation included. This diagnostic does not replace frozen main labels or select niches using biological associations.\n')
