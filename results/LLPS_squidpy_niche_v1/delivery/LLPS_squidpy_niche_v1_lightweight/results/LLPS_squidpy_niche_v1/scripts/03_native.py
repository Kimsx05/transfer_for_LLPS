from pathlib import Path
import numpy as np,pandas as pd,scipy.sparse as sp,scanpy as sc,anndata as ad,squidpy as sq
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
O=Path(__file__).resolve().parents[1];d=pd.read_parquet(O/'00_input_audit/audited_spots.parquet');cols=[c for c in d if c.startswith('major_mean__')];x=d[cols].to_numpy();p=x/x.sum(1)[:,None];dominant=np.array(cols)[p.argmax(1)];out=d[['spot_id','section_id','dataset']].copy();out['dominant_type']=dominant;out['dominant_fraction']=p.max(1)
# Squidpy 1.6.5 does not forward random_state into neighborhood flavor.
# Explicit adapters enforce seed123 and consistent Leiden settings at actual Scanpy calls.
# Installed1.6.5 pads variable neighborhoods with final-node labels (NaN/-1 mask mismatch).
# Replace profile construction with validated sparse count normalization; retain native niche clustering.
import importlib
niche_module=importlib.import_module('squidpy.gr._niche')
original_profile=niche_module._calculate_neighborhood_profile
def exact_profile(adata,groups,matrix,abs_nhood):
 categories=np.unique(adata.obs[groups].astype(str));codes=pd.Categorical(adata.obs[groups].astype(str),categories=categories).codes
 onehot=sp.csr_matrix((np.ones(len(codes)),(np.arange(len(codes)),codes)),shape=(len(codes),len(categories)))
 mat=matrix.tocsr();counts=(mat@onehot).toarray();den=np.asarray(mat.sum(1)).ravel()
 values=counts if abs_nhood else np.divide(counts,den[:,None],out=np.zeros_like(counts),where=den[:,None]>0)
 return pd.DataFrame(values,index=adata.obs_names)
toy=ad.AnnData(sp.csr_matrix((4,1)));toy.obs['type']=pd.Categorical(['A','B','C','D']);toygraph=sp.eye(4,format='lil');toygraph[0,1]=toygraph[1,0]=1;toygraph=toygraph.tocoo()
bad=original_profile(toy,'type',toygraph,False).to_numpy();good=exact_profile(toy,'type',toygraph,False).to_numpy();assert np.allclose(good[2],[0,0,1,0]);assert not np.allclose(bad,good)
np.savetxt(O/'03_clustering/native_padding_bug_original.tsv',bad,delimiter='\t');np.savetxt(O/'03_clustering/native_padding_bug_corrected.tsv',good,delimiter='\t')
niche_module._calculate_neighborhood_profile=exact_profile
orig_n=sc.pp.neighbors;orig_l=sc.tl.leiden
sc.pp.neighbors=lambda *a,**kw:orig_n(*a,**(kw|dict(random_state=123,transformer="pynndescent")))
sc.tl.leiden=lambda *a,**kw:orig_l(*a,**(kw|dict(random_state=123,flavor='igraph',n_iterations=2)))
As=[sp.load_npz(O/'01_spatial_graph'/f'{sid}_R1.npz') for sid in d.section_id.drop_duplicates()];A=sp.block_diag(As,format='csr');C=sp.eye(len(d),format='csr',dtype='float32')
for r in [1,2,3]:
 C=(C+C@A).astype(bool).astype('float32');z=ad.AnnData(sp.csr_matrix((len(d),1),dtype='float32'));z.obs['dominant']=pd.Categorical(dominant);z.obsp['spatial_connectivities']=C
 # Supply cumulative binary radius graph, with diagonal. distance=1 avoids native path-count weighting at higher hops.
 profile=exact_profile(z,'dominant',C.tocoo(),False).to_numpy();assert np.allclose(profile,np.load(O/'02_neighborhood_composition'/f'dominant_R{r}.npy'),atol=1e-6)
 sq.gr.calculate_niche(z,flavor='neighborhood',groups='dominant',n_neighbors=30,resolutions=.4,scale=False,abs_nhood=False,distance=1,random_state=123)
 out[f'native_R{r}']=z.obs['nhood_niche_res=0.4'].astype(str).to_numpy();print('native',r,flush=True)
out.to_parquet(O/'03_clustering/native_labels.parquet',index=False)
out.groupby('section_id').dominant_fraction.agg(['mean','median','min','max']).to_csv(O/'03_clustering/dominant_information_loss.tsv',sep='\t')
(O/'03_clustering/native_API_notes.md').write_text('Squidpy 1.6.5 calculate_niche(neighborhood), scale=False, resolution0.4, feature KNN30, all spots pooled. Exact cumulative binary R1/R2/R3 graphs include centers, distance=1. Avoid native higher-hop path weighting: each observed spot contributes once. Block-diagonal spatial graph prevents cross-section edges. Native internal Scanpy random defaults explicitly overridden to123; igraph Leiden2 iterations matches continuous method. Installed1.6.5 profile padding bug corrected through an explicitly replaced sparse profile helper, validated on a four-node toy and against full actual dominant composition; native clustering remains unchanged. Dominant profiles discard mixture magnitude; dominant_fraction quantifies retained share. Native comparison is pooled whereas main method is section-balanced mapped; both representation and training differences must be considered.\n')
