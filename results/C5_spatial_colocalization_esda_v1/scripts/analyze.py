"""Frozen per-section PySAL analysis. No expression/scoring/model re-fitting."""
import os,time,json,hashlib,inspect,sys,resource,shutil
from pathlib import Path
import numpy as np,pandas as pd,pyarrow as pa
from scipy import sparse,stats
from scipy.sparse.csgraph import connected_components
import esda,libpysal
from libpysal.weights import W,lag_spatial
pa.set_cpu_count(1);pa.set_io_thread_count(1)
O=Path('results/C5_spatial_colocalization_esda_v1');S=Path('results/C5_spatial_exploratory_all63_v1');PARTNERS=['Endothelial','Epithelial','Myeloid','T_NK'];SEED=123;PERM=999
META=pd.read_csv(S/'00_inputs_and_QC/spatial_tissue_patient_audit.tsv',sep='\t').set_index('section_id');PILOT=pd.read_csv(S/'00_inputs_and_QC/pilot_prespecified.tsv',sep='\t').section_id.tolist()
def save(d,p):pd.DataFrame(d).to_csv(p,sep='\t',index=False)
def bh(p):
 p=np.asarray(p,float);out=np.full(len(p),np.nan);ok=np.isfinite(p);a=p[ok];ix=np.argsort(a);q=np.minimum.accumulate((a[ix]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];v=np.empty(len(a));v[ix]=np.minimum(q,1);out[ok]=v;return out
def graph(d):
 rc=d[['array_row','array_col']].to_numpy();assert np.isfinite(rc).all() and np.all(rc==rc.astype(int));lookup={tuple(v):i for i,v in enumerate(rc)};assert len(lookup)==len(d)
 nbs={i:[lookup[(r+dr,c+dc)] for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)] if (r+dr,c+dc) in lookup] for i,(r,c) in enumerate(rc)}
 return nbs
def weights(nbs):
 w=W(nbs,id_order=list(range(len(nbs))),silence_warnings=True);w.transform='r';return w
def rings(nbs,k):
 out={}
 for i in nbs:
  seen={i};front={i}
  for _ in range(k):front={j for v in front for j in nbs[v]}-seen;seen|=front
  out[i]=sorted(front)
 return out
def resid(v,controls):
 a=np.column_stack([np.ones(len(v))]+[stats.rankdata(z)/len(v) for z in controls]);r=stats.rankdata(v)/len(v);return r-a@np.linalg.lstsq(a,r,rcond=None)[0]
def moran(x,y,w):
 if np.std(x)<1e-12 or np.std(y)<1e-12:return dict(I=np.nan,p_sim=np.nan,p_two_sided=np.nan,status='CONSTANT')
 np.random.seed(SEED);m=esda.Moran_BV(x,y,w,transformation='r',permutations=PERM)
 return dict(I=float(m.I),p_sim=float(m.p_sim),p_two_sided=min(1,2*float(m.p_sim)),status='OK')
def audit():
 for f in ['spatial_tissue_patient_audit.tsv','pilot_prespecified.tsv','frozen_signatures.json']:shutil.copyfile(S/'00_inputs_and_QC'/f,O/'00_input_audit'/f)
 assert json.load(open(S/'00_inputs_and_QC/frozen_signatures.json'))['five']==['MMP11','COL11A1','INHBA','COL5A1','CTHRC1']
 imports={m:__import__(m).__version__ for m in ['esda','libpysal','numpy','scipy','pandas','numba','pyarrow']};imports['python']=sys.version;imports['global_API']=str(inspect.signature(esda.Moran_BV));imports['local_API']=str(inspect.signature(esda.Moran_Local_BV));json.dump(imports,open(O/'00_input_audit/python_versions.json','w'),indent=2)
 (O/'00_input_audit/esda_installed_source.txt').write_text(inspect.getsource(esda.Moran_BV)+'\n'+inspect.getsource(esda.Moran_Local_BV))
 h=[]
 for f in sorted((S/'01_spatial_scores').glob('*.parquet')):h.append(dict(path=str(f.resolve()),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 hp=O/'00_input_audit/source_hashes.tsv'
 if hp.exists():
  previous=pd.read_csv(hp,sep='\t');current=pd.DataFrame(h);assert previous.set_index('path').sha256.to_dict()==current.set_index('path').sha256.to_dict(),'Input hashes changed: create a new output version, never silently rebaseline.'
 else:save(h,hp)
 assert len(h)==63
 json.dump(dict(seed=SEED,permutations=PERM,n_jobs=1,keep_simulations=False,geoda_quads=False,quadrants={'1':'HH','2':'LH','3':'LL','4':'HL'},primary='Five UCell -> raw major posterior mean; ring 1',p_value='native folded/minimum-tail p_sim retained; conservative min(1,2*p_sim) used for primary BH',global_BH='252 forward primary tests; each sensitivity and reverse separate 252-test family',local_BH='per section: four partners times all valid graph-connected nodes; doubled p_sim',isolates='all analyses exclude original isolates; no fake edges; missing required values excluded then newly isolated nodes excluded; local output NA for excluded nodes',patient='only nonmissing audited patient_uid; mean per patient then bootstrap median 2000 replicates',rings='exact shortest-path rings, never cumulative; per-ring isolates excluded',primary_frozen_before_pilot=True),open(O/'00_input_audit/frozen_analysis_plan.json','w'),indent=2)
def section(sid):
 done=O/'00_input_audit'/f'{sid}.json'
 if done.exists():print('REUSE',sid,flush=True);return
 start=time.time();d=pd.read_parquet(S/'01_spatial_scores'/f'{sid}.parquet');assert d.spot_id.is_unique and d.section_id.eq(sid).all();n=len(d);original=graph(d)
 majors=[x for x in d if x.startswith('major_mean__')];assert len(majors)==9
 required=['five_UCell','five_AUCell','A26_UCell','core15_UCell','nCount','nFeature']+majors+['major_q05__'+p for p in PARTNERS]
 valid=np.isfinite(d[required].to_numpy()).all(1);missing=int((~valid).sum());oldvalid=valid.copy()
 while True:
  valid=np.array([valid[i] and any(valid[j] for j in original[i]) for i in range(n)])
  if np.array_equal(valid,oldvalid):break
  oldvalid=valid.copy()
 ix=np.flatnonzero(valid);rev={v:i for i,v in enumerate(ix)};nbs={rev[i]:[rev[j] for j in original[i] if j in rev] for i in ix};v=d.iloc[ix];w=weights(nbs);assert w.id_order==list(range(len(v)));assert np.allclose(np.asarray(w.sparse.sum(1)).ravel(),1);assert all(i in nbs[j] for i in nbs for j in nbs[i]);assert max(map(len,nbs.values()))<=6
 assert len(v)>10
 node=d[['spot_id','array_row','array_col']].copy();node['analysis_index']=-1;node.loc[ix,'analysis_index']=np.arange(len(ix));save(node,O/'01_spatial_weights'/f'{sid}_nodes.tsv.gz')
 edges=[(i,j) for i in nbs for j in nbs[i] if i<j];save(pd.DataFrame(edges,columns=['source','target']),O/'01_spatial_weights'/f'{sid}_edges.tsv.gz')
 meta=dict(section_id=sid,**META.loc[sid].to_dict());nc,lab=connected_components(w.sparse,directed=False);qc=dict(section_id=sid,n_total=n,n_usable=len(ix),missing=missing,excluded_isolates=n-missing-len(ix),original_isolates=sum(not a for a in original.values()),components=nc,edges=len(edges),degree_min=min(map(len,nbs.values())),degree_max=max(map(len,nbs.values())),largest_component=int(np.bincount(lab).max()),node_order_sha256=hashlib.sha256('\n'.join(v.spot_id).encode()).hexdigest())
 x=v.five_UCell.to_numpy();xau=v.five_AUCell.to_numpy();fib=v.major_mean__Fibroblast.to_numpy();mural=v.major_mean__Mural.to_numpy();epi=v.major_mean__Epithelial.to_numpy();base=[fib,mural,v.nCount.to_numpy(),v.nFeature.to_numpy()]
 out=d[['spot_id','x','y','array_row','array_col','five_UCell','five_AUCell','nCount','nFeature']+majors+['major_q05__'+p for p in PARTNERS]].copy();out['analysis_included']=valid
 results=[];loc=[];dist=[];distributions=[];neighborhood=[]
 for partner in PARTNERS:
  y=v['major_mean__'+partner].to_numpy();ctrl=base+([] if partner=='Epithelial' else [epi]);xr=resid(x,ctrl);yr=resid(y,ctrl);xb=resid(x,base);yb=resid(y,base)
  methods={'primary':(x,y),'reverse':(y,x),'AUCell':(xau,y),'log1p_abundance':(x,np.log1p(y)),'major_fraction':(x,y/v[majors].sum(axis=1).to_numpy()),'q05':(x,v['major_q05__'+partner].to_numpy()),'adjusted_base':(xb,yb),'adjusted_full':(xr,yr),'A26':(v.A26_UCell.to_numpy(),y),'Core15':(v.core15_UCell.to_numpy(),y)}
  for method,(xx,yy) in methods.items():results.append(dict(**meta,partner=partner,method=method,n=len(v),edges=len(edges),same_spot_spearman=stats.spearmanr(xx,yy).statistic,controls='Fibroblast,Mural,nCount,nFeature'+(',Epithelial' if partner!='Epithelial' else '') if method=='adjusted_full' else '',**moran(xx,yy,w)))
  ml=esda.Moran_Local_BV(x,y,w,transformation='r',permutations=PERM,geoda_quads=False,n_jobs=1,keep_simulations=False,seed=SEED)
  expected=np.where(ml.zx>0,np.where(lag_spatial(w,ml.zy)>0,1,4),np.where(lag_spatial(w,ml.zy)>0,2,3));assert np.array_equal(ml.q,expected),'Quadrant mapping mismatch'
  assert ml.sim is None or np.size(ml.sim)<=1,'Unexpected retained simulations'
  for name,arr in [('I',ml.Is),('quadrant',ml.q),('p_sim',ml.p_sim),('p_two_sided',np.minimum(1,2*ml.p_sim)),('lag_abundance',lag_spatial(w,y)),('lag_standardized_abundance',lag_spatial(w,ml.zy)),('adjusted_five',xr),('adjusted_abundance',yr)]:out[f'{partner}__{name}']=np.nan;out.loc[ix,f'{partner}__{name}']=arr
  for label,z in [('Five_UCell',x),(partner,y),(partner+'_q05',v['major_q05__'+partner].to_numpy())]:distributions.append(dict(section_id=sid,variable=label,n=len(z),min=z.min(),p25=np.quantile(z,.25),median=np.median(z),p75=np.quantile(z,.75),max=z.max(),zero_fraction=np.mean(z==0)))
  high=x>np.quantile(x,.9);zy=(y-y.mean())/y.std();neighbor=lag_spatial(w,zy);neighborhood.append(dict(section_id=sid,partner=partner,n_high=high.sum(),n_five_high_same_low_neighbor_high=int(np.sum(high&(zy<0)&(neighbor>0))),fraction_high_same_low_neighbor_high=np.mean((zy<0)[high]&(neighbor>0)[high]),same_low_neighbor_high_all=int(np.sum((zy<0)&(neighbor>0)))))
  loc.append(dict(section_id=sid,partner=partner,raw_p_min=float(np.min(ml.p_sim)),n_tests=len(ix),HH_descriptive=int(np.sum(ml.q==1))))
 # Local BH combines ALL four partners, not separate partner families.
 p=np.concatenate([out.loc[ix,p+'__p_two_sided'].to_numpy() for p in PARTNERS]);qs=bh(p).reshape(4,len(ix))
 for j,partner in enumerate(PARTNERS):
  out[partner+'__FDR']=np.nan;out.loc[ix,partner+'__FDR']=qs[j];out[partner+'__class']='Excluded';q=out.loc[ix,partner+'__quadrant'].map({1:'HH',2:'LH',3:'LL',4:'HL'}).to_numpy();out.loc[ix,partner+'__class']=np.where(qs[j]<.05,q,'Not significant');loc[j].update(HH_FDR=int(np.sum((q=='HH')&(qs[j]<.05))),HL_FDR=int(np.sum((q=='HL')&(qs[j]<.05))),LH_FDR=int(np.sum((q=='LH')&(qs[j]<.05))),LL_FDR=int(np.sum((q=='LL')&(qs[j]<.05))),FDR_min=float(np.min(qs[j])))
 for k in [1,2,3]:
  rn=nbs if k==1 else rings(nbs,k);keep=np.array([bool(rn[i]) for i in range(len(v))]);ki=np.flatnonzero(keep);rm={j:i for i,j in enumerate(ki)};nn={rm[i]:[rm[j] for j in rn[i] if j in rm] for i in ki};rw=weights(nn)
  for partner in PARTNERS:
   yy=v['major_mean__'+partner].to_numpy();ctrl=base+([] if partner=='Epithelial' else [epi]);xxr=resid(x,ctrl);yyr=resid(yy,ctrl)
   for method,xx,yv in [('raw',x,yy),('adjusted',xxr,yyr)]:dist.append(dict(**meta,partner=partner,method=method,ring=k,n=len(ki),edges=sum(map(len,nn.values()))//2,excluded_ring_isolates=int((~keep).sum()),**moran(xx[ki],yv[ki],rw)))
 # R and plotting share small analysis-only tables, never expression matrices.
 out.to_parquet(O/'03_local_bivariate_moran'/f'{sid}.parquet',index=False);save(out,O/'figures'/f'{sid}_plot_data.tsv.gz')
 for rows,folder,suffix in [(results,'02_global_bivariate_moran','results'),(loc,'03_local_bivariate_moran','summary'),(dist,'05_distance_sensitivity','rings'),(distributions,'00_input_audit','distributions'),(neighborhood,'05_distance_sensitivity','peripheral')]:save(rows,O/folder/f'{sid}_{suffix}.tsv')
 qc.update(seconds=time.time()-start,max_RSS_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,quadrant_verified=True,simulations_retained=False,seed=SEED,permutations=PERM);json.dump(qc,open(done,'w'),indent=2);print('DONE',sid,round(qc['seconds'],2),'sec',qc['max_RSS_MB'],'MB',flush=True)
def aggregate():
 for folder,suffix,dest in [('02_global_bivariate_moran','results','global_results'),('03_local_bivariate_moran','summary','local_summary'),('05_distance_sensitivity','rings','ring_results'),('00_input_audit','distributions','input_distributions'),('05_distance_sensitivity','peripheral','peripheral_patterns')]:
  z=pd.concat([pd.read_csv(f,sep='\t') for f in sorted((O/folder).glob('*_'+suffix+'.tsv'))],ignore_index=True)
  if 'p_two_sided' in z:
   keys=['method']+(['ring'] if 'ring' in z else []);z['FDR']=z.groupby(keys).p_two_sided.transform(lambda x:bh(x));z['FDR_native_p_sim']=z.groupby(keys).p_sim.transform(lambda x:bh(x))
  save(z,O/'tables'/f'{dest}.tsv')
 save([json.load(open(O/'00_input_audit'/f'{sid}.json')) for sid in META.index if (O/'00_input_audit'/f'{sid}.json').exists()],O/'01_spatial_weights/graph_QC.tsv')
if __name__=='__main__':
 mode=sys.argv[1]
 if mode=='audit':audit()
 elif mode=='aggregate':aggregate()
 else:
  sections=PILOT if mode=='pilot' else sorted(META.index) if mode=='all' else [mode]
  for s in sections:section(s)
  aggregate()
