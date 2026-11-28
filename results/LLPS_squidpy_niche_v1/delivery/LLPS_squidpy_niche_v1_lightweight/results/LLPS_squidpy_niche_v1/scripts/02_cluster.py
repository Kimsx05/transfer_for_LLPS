from pathlib import Path
import json,itertools,time
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from sklearn.metrics import silhouette_score,calinski_harabasz_score,davies_bouldin_score,adjusted_rand_score,normalized_mutual_info_score
from sklearn.cluster import KMeans
from sklearn.neighbors import KNeighborsClassifier
O=Path(__file__).resolve().parents[1];D=O/'03_clustering';rng=np.random.default_rng(123)
# Deliberate leakage barrier: load identities only; scores/clinical data never enter selection.
m=pd.read_parquet(O/'00_input_audit/audited_spots.parquet',columns=['spot_id','section_id','dataset']);pilot=np.sort(np.concatenate([rng.choice(g.index,min(200,len(g)),replace=False) for _,g in m.groupby('section_id',sort=False)]));np.save(D/'pilot_indices.npy',pilot)
rows=[];labels={};params={};cache={}
def feature(kind,r,t):
 x=np.load(O/'02_neighborhood_composition'/f'{kind}_R{r}.npy',mmap_mode='r');return np.sqrt(x) if t=='hellinger' else x

def quality(x,y):
 k=len(np.unique(y));counts=np.unique(y,return_counts=True)[1]
 if k<2 or k>=len(y):return dict(n_niches=k,silhouette=np.nan,ch=np.nan,db=np.nan,min_fraction=counts.min()/len(y))
 return dict(n_niches=k,silhouette=float(silhouette_score(x,y,sample_size=min(2000,len(y)),random_state=123)),ch=float(calinski_harabasz_score(x,y)),db=float(davies_bouldin_score(x,y)),min_fraction=counts.min()/len(y))
def fitgrid(kind,r,t,ks):
 x=feature(kind,r,t)[pilot]
 for k in ks:
  z=ad.AnnData(np.array(x));sc.pp.neighbors(z,n_neighbors=k,use_rep='X',random_state=123,transformer='pynndescent')
  for res in [.2,.4,.6,.8,1.]:
   name=f'{kind}_R{r}_{t}_k{k}_res{res}';sc.tl.leiden(z,resolution=res,random_state=123,flavor='igraph',n_iterations=2,key_added='niche');y=z.obs.niche.astype(str).to_numpy();p=dict(model=name,kind=kind,radius=r,transform=t,knn=k,resolution=res,method='leiden');params[name]=p;labels[name]=y;rows.append(p|quality(x,y));print(name,rows[-1]['n_niches'],rows[-1]['silhouette'],flush=True)
for kind,r,t in itertools.product(['major','subtype'],[1,2,3],['raw','hellinger']):fitgrid(kind,r,t,[30])
q=pd.DataFrame(rows)
# Refine k only for each annotation level's best R2 transform/resolution family.
for kind in ['major','subtype']:
 a=q[(q.kind==kind)&(q.radius==2)].copy();a=a[a.n_niches.between(3,20)];best=a.sort_values('silhouette',ascending=False).iloc[0];fitgrid(kind,2,best['transform'],[15,50])
q=pd.DataFrame(rows);q['selection_score']=q.silhouette-.10*(q.min_fraction<.01)-.03*np.maximum(q.n_niches-12,0);q.to_csv(D/'clustering_quality_metrics.tsv',sep='\t',index=False);q.drop(columns=['silhouette','ch','db','min_fraction','selection_score']).to_csv(D/'clustering_parameter_grid.tsv',sep='\t',index=False)
def best_for(sub):return sub[sub.n_niches.between(3,20)].sort_values(['selection_score','model'],ascending=[False,True]).iloc[0]
chosen={};chosen['major_best']=best_for(q[q.kind=='major']);chosen['subtype_best']=best_for(q[q.kind=='subtype']);best=best_for(q);chosen['recommended']=best
for r in [1,2,3]:chosen[f'R{r}_representative']=best_for(q[q.radius==r])
# Recommendation fixed BEFORE reading any score, stage or clinical variable.
rec=pd.DataFrame([dict(role=role,**row.to_dict()) for role,row in chosen.items()]);rec.to_csv(D/'recommended_models.tsv',sep='\t',index=False)
full=m.copy();conf=[];stability=[];agreements=[];fitted={}
for name in rec.model.unique():
 p=params[name];x=feature(p['kind'],p['radius'],p['transform']);y=labels[name];model=KNeighborsClassifier(n_neighbors=15,weights='distance',n_jobs=4).fit(x[pilot],y);pred=model.predict(x);prob=model.predict_proba(x).max(1);full[name]=pred;full[name+'__confidence']=prob.astype('float32');fitted[name]=pred
 for sid,g in m.groupby('section_id',sort=False):conf.append(dict(model=name,section_id=sid,mean_confidence=float(prob[g.index].mean()),below_0_6=float(np.mean(prob[g.index]<.6))))
 z=ad.AnnData(np.array(x[pilot]));sc.pp.neighbors(z,n_neighbors=int(p['knn']),use_rep='X',random_state=123,transformer='pynndescent')
 for seed in [17,321,2026]:
  sc.tl.leiden(z,resolution=p['resolution'],random_state=seed,flavor='igraph',n_iterations=2,key_added='niche');v=z.obs.niche.astype(str).to_numpy();stability.append(dict(model=name,seed=seed,ARI=adjusted_rand_score(y,v),NMI=normalized_mutual_info_score(y,v),scope='same_balanced_pilot'))
 # Density-normalized sensitivity, same settings and training spots.
 equal=feature(p['kind']+'_equal_spot',p['radius'],p['transform']);ze=ad.AnnData(np.array(equal[pilot]));sc.pp.neighbors(ze,n_neighbors=int(p['knn']),use_rep='X',random_state=123,transformer='pynndescent');sc.tl.leiden(ze,resolution=p['resolution'],random_state=123,flavor='igraph',n_iterations=2,key_added='niche');ye=ze.obs.niche.astype(str).to_numpy();full[name+'__equal_spot']=KNeighborsClassifier(n_neighbors=15,weights='distance',n_jobs=4).fit(equal[pilot],ye).predict(equal);agreements.append(dict(model_a=name,model_b=name+'__equal_spot',scope='all_mapped_spots',ARI=adjusted_rand_score(pred,full[name+'__equal_spot']),NMI=normalized_mutual_info_score(pred,full[name+'__equal_spot'])))
 print('MAPPED',name,flush=True)
for role,row in chosen.items():full[role+'_niche']=full[row.model]
# Independent K-means, all requested K, selected using silhouette only.
p=params[best.model];x=feature(p['kind'],p['radius'],p['transform']);kmrows=[];kms={}
for k in range(4,13):
 km=KMeans(n_clusters=k,n_init=10,random_state=123).fit(x[pilot]);v=km.predict(x[pilot]);kmrows.append(dict(k=k,**quality(x[pilot],v)));kms[k]=km
kbest=int(pd.DataFrame(kmrows).sort_values('silhouette',ascending=False).iloc[0].k);full['kmeans_niche']=kms[kbest].predict(x).astype(str);pd.DataFrame(kmrows).to_csv(D/'kmeans_quality.tsv',sep='\t',index=False)
# Full pooled Leiden for major/subtype winners: validate mapping against full-data clustering.
for role in ['major_best','subtype_best']:
 p=params[chosen[role].model];x=feature(p['kind'],p['radius'],p['transform']);z=ad.AnnData(np.array(x));sc.pp.neighbors(z,n_neighbors=int(p['knn']),use_rep='X',random_state=123,transformer='pynndescent');sc.tl.leiden(z,resolution=p['resolution'],random_state=123,flavor='igraph',n_iterations=2,key_added='niche');full[role+'__pooled']=z.obs.niche.astype(str).to_numpy();sc.tl.umap(z,random_state=123);np.save(D/(role+'_umap.npy'),z.obsm['X_umap']);agreements.append(dict(model_a=chosen[role].model,model_b=role+'__pooled',scope='all_spots_balanced_vs_pooled',ARI=adjusted_rand_score(full[chosen[role].model],full[role+'__pooled']),NMI=normalized_mutual_info_score(full[chosen[role].model],full[role+'__pooled'])));print('POOLED',role,flush=True)
cols=[c for c in full if c.endswith('_niche')]
for a,b in itertools.combinations(cols,2):agreements.append(dict(model_a=a,model_b=b,scope='all_spots',ARI=adjusted_rand_score(full[a],full[b]),NMI=normalized_mutual_info_score(full[a],full[b])))
full.to_parquet(D/'spot_labels.parquet',index=False,compression='zstd');pd.DataFrame(conf).to_csv(D/'label_transfer_confidence.tsv',sep='\t',index=False);pd.DataFrame(stability).to_csv(D/'clustering_stability.tsv',sep='\t',index=False);pd.DataFrame(agreements).to_csv(D/'clustering_ARI_NMI.tsv',sep='\t',index=False)
np.savez_compressed(D/'pilot_labels.npz',**labels)
(O/'03_clustering/model_selection_report.md').write_text('''# Composition-only model selection\n\nBalanced pilot: up to 200 spots per section, seed 123. All 12 radius/annotation/transformation combinations evaluated at k=30 and resolutions 0.2–1.0. k=15/50 refined for each annotation level at R2 using its highest-silhouette transform. Selection score = silhouette -0.10 if any niche <1% of pilot -0.03 per niche above12. 3–20 niches allowed; target5–10 not enforced. This prespecified numerical rule does not inspect signatures, clinical variables, coordinates or communications. Spatial diagnostics and reproducibility must still inform reliability ratings.\n\nUnified training labels mapped via 15-neighbor distance-weighted classifier in the same composition feature space; maximum class probability saved per spot. Full pooled Leiden validates both annotation winners. UMAP is only visualization and is fitted on composition KNN, never used as clustering input. Igraph Leiden uses 2 iterations consistently for tractability; seed sensitivity is saved.\n\nNo Z-score or batch correction. All candidate pilot labels saved; every chosen model labels all spots. K-means K4–12 is an independent comparison in recommended composition feature space.\n''')

# Freeze complete model artifacts before downstream score association in future replays.
import hashlib,datetime
paths=['03_clustering/recommended_models.tsv','03_clustering/spot_labels.parquet','03_clustering/clustering_parameter_grid.tsv','03_clustering/pilot_labels.npz']
(O/'03_clustering/clustering_lock.json').write_text(json.dumps(dict(locked_at_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),before_signature_enrichment=True,selection_inputs=['continuous_cell_composition','section_balancing'],forbidden_variables_used=False,files={p:hashlib.sha256((O/p).read_bytes()).hexdigest() for p in paths}),indent=2))
