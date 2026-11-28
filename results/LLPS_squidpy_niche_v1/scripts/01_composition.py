from pathlib import Path
import inspect,json
import numpy as np,pandas as pd,scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
import anndata as ad,squidpy as sq
O=Path(__file__).resolve().parents[1]
d=pd.read_parquet(O/'00_input_audit/audited_spots.parquet');maj=[c for c in d if c.startswith('major_mean__')];sub=[c for c in d if c.startswith('sub_mean__')]
qc=[];sizes=[];cover=[];arrays={};covs=[]
for r in [1,2,3]:
 for kind in ['major','subtype','major_equal_spot','subtype_equal_spot','dominant']:arrays[f'{kind}_R{r}']=[]
for sid,g in d.groupby('section_id',sort=False):
 n=len(g);rawxy=g[['array_row','array_col']].to_numpy();assert np.equal(rawxy,np.floor(rawxy)).all();xy=rawxy.astype(int);look={tuple(v):i for i,v in enumerate(xy)}
 edges=np.array([(i,look[(r+dr,c+dc)]) for i,(r,c) in enumerate(xy) for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)] if (r+dr,c+dc) in look]);A=sp.csr_matrix((np.ones(len(edges),dtype='float32'),(edges[:,0],edges[:,1])),shape=(n,n))
 z=ad.AnnData(sp.csr_matrix((n,1),dtype='float32'));z.obsm['spatial']=np.column_stack([xy[:,1]/2,xy[:,0]*np.sqrt(3)/2])
 sq.gr.spatial_neighbors(z,coord_type='grid',n_neighs=6,n_rings=1,set_diag=False)
 native=z.obsp['spatial_connectivities'].astype(bool).astype('float32');extra=native-native.multiply(A);missing=A-A.multiply(native)
 assert (A!=A.T).nnz==0 and A.diagonal().sum()==0
 nc,lab=connected_components(A);deg=np.asarray(A.sum(1)).ravel();qc.append(dict(section_id=sid,n_spots=n,edges=A.nnz//2,isolates=int((deg==0).sum()),boundary_spots=int((deg<6).sum()),components=nc,squidpy_extra_directed_edges=extra.nnz,squidpy_missing_directed_edges=missing.nnz,cross_section_edges=0,graph='Squidpy grid validated/corrected by exact array adjacency'))
 sp.save_npz(O/'01_spatial_graph'/f'{sid}_R1.npz',A)
 M=g[maj].to_numpy('float32');S=g[sub].to_numpy('float32');tot=M.sum(1);assert (tot>0).all();mp=M/tot[:,None];ss=S/S.sum(1)[:,None];dom=np.eye(len(maj),dtype='float32')[mp.argmax(1)]
 C=sp.eye(n,format='csr',dtype='float32');newcov=g[['spot_id']].copy()
 for r,theory in [(1,7),(2,19),(3,37)]:
  C=(C+C@A).astype(bool).astype('float32');cnt=np.asarray(C.sum(1)).ravel();assert cnt.max()<=theory
  newcov[f'R{r}_n_spots']=cnt.astype('int16');newcov[f'R{r}_coverage']=(cnt/theory).astype('float32')
  for kind,V in [('major',M),('subtype',S)]:
   V=C@V;arrays[f'{kind}_R{r}'].append((V/V.sum(1)[:,None]).astype('float32'))
  for kind,V in [('major_equal_spot',mp),('subtype_equal_spot',ss),('dominant',dom)]:arrays[f'{kind}_R{r}'].append(((C@V)/cnt[:,None]).astype('float32'))
  sizes.append(dict(section_id=sid,radius=r,theoretical=theory,min=int(cnt.min()),median=float(np.median(cnt)),mean=float(cnt.mean()),max=int(cnt.max())))
  cover.append(dict(section_id=sid,radius=r,mean_coverage=float(np.mean(cnt/theory)),complete_fraction=float(np.mean(cnt==theory)),under_half_fraction=float(np.mean(cnt<theory/2))))
 covs.append(newcov);print(sid,flush=True)
for name,vals in arrays.items():np.save(O/'02_neighborhood_composition'/f'{name}.npy',np.concatenate(vals))
pd.concat(covs).to_parquet(O/'02_neighborhood_composition/coverage.parquet',index=False)
for name,rows in [('spatial_graph_QC',qc),('neighborhood_size_summary',sizes),('neighborhood_coverage_summary',cover)]:pd.DataFrame(rows).to_csv(O/'01_spatial_graph'/f'{name}.tsv',sep='\t',index=False)
(O/'02_neighborhood_composition/features.json').write_text(json.dumps(dict(major=maj,subtype=sub),indent=2))
(O/'01_spatial_graph/API_audit.txt').write_text(str(inspect.signature(sq.gr.spatial_neighbors))+'\n'+str(inspect.signature(sq.gr.calculate_niche_neighborhood)) if hasattr(sq.gr,'calculate_niche_neighborhood') else 'Native niche API unavailable; composition equivalent required')
