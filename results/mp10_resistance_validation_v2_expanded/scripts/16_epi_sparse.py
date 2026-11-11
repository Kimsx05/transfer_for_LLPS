from env import *
from scipy.sparse import csc_matrix
from scipy.stats import rankdata
import time,importlib
C=P/'neighborhood_mp10_calibration_v5/01_score_fidelity/counts_cache';genes=C.joinpath('genes.txt').read_text().splitlines();cells=C.joinpath('cells.txt').read_text().splitlines();K=4501;assert len(genes)==90011
f=pd.read_parquet(P/'neighborhood_mp10_calibration_v5/01_score_fidelity/cell_level_fidelity.parquet').set_index('cell_id').loc[cells].reset_index();f['sample']=f.dataset+'|'+f['orig.ident'];pat=read('00_plan_audit/epi_patient_identity.tsv');f=f.merge(pat[['sample','patient','verified']],on='sample',how='left',validate='many_to_one');assert f.cell_id.tolist()==cells
base=Path('/home/data/t070721/Metabolistic_scRNA/Multi2');cm=pd.read_csv(base/'3.3 recognize malignant/rough_cluster_infercnv_score_cell_metadata.csv',usecols=['cell','rough_cluster']);st=pd.read_csv(base/'3.4 malignant cell process/NMF1/rough_cluster_malignant_like_status_summary.csv');cm=cm.merge(st[['rough_cluster','malignant_like_status']],on='rough_cluster',validate='many_to_one');f=f.merge(cm.rename(columns={'cell':'cell_id'}),on='cell_id',how='left',validate='one_to_one')
ind=np.memmap(C/'indices.bin',dtype='<i4',mode='r');dat=np.memmap(C/'data.bin',dtype='<f8',mode='r');ptr=np.memmap(C/'indptr.bin',dtype='<i4',mode='r');assert ptr[-1]==len(dat)==len(ind)
# Map signatures to exact symbols in protected feature background. Unique aliases only if exact unavailable.
mod=importlib.import_module('06_prepare');gm=[mod.mapgene(g,g)[0] for g in genes];pos={}
for i,g in enumerate(gm):pos.setdefault(g,[]).append(i)
m=read('01_signatures/members.tsv');mf=read('01_signatures/signature_manifest.tsv');sets={};cov=[];mapping=[]
for sn in mf.loc[mf.method=='ssGSEA','signature']:
 r=m[m.signature==sn];ii=[];missing=[]
 for _,q in r.iterrows():
  candidates=[genes.index(q.gene)] if q.gene in genes else pos.get(q.gene,[]);hit=candidates[0] if len(candidates)==1 else None
  mapping.append(dict(signature=sn,original_member=q.original_member,symbol=q.gene,input_gene=genes[hit] if hit is not None else None,reason='exact/unique' if hit is not None else 'absent/ambiguous'))
  if hit is None:missing.append(str(q.original_member))
  else:ii.append(hit)
 ok=len(r)>0 and len(ii)/len(r)>=.8
 cov.append(dict(dataset='Epi',signature=sn,n_original=len(r),n_present=len(ii),coverage=len(ii)/len(r) if len(r) else 0,missing=';'.join(missing),status='READY' if ok else 'BLOCKED_COVERAGE'))
 if ok:sets[sn]=np.array(ii)
save(pd.DataFrame(mapping),'00_plan_audit/epi_signature_feature_mapping.tsv');save(pd.DataFrame(cov),'00_plan_audit/epi_coverage.tsv')
# Pilot inputs for official methods: bounded 256 cells, original90,011 background.
tmp=B/'02_epi/pilot';tmp.mkdir(exist_ok=True,parents=True);n=256;a=int(ptr[n]);np.asarray(dat[:a]).tofile(tmp/'data.bin');np.asarray(ind[:a]).tofile(tmp/'indices.bin');np.asarray(ptr[:n+1]).tofile(tmp/'indptr.bin');(tmp/'genes.txt').write_text('\n'.join(genes));(tmp/'cells.txt').write_text('\n'.join(cells[:n]));save(f.iloc[:n][['cell_id','sample']],'02_epi/pilot/meta.tsv')
if '--pilot-only' in sys.argv:print('PILOT_READY');exit()
started=time.time();levels=sorted(f['sample'].unique());codes=pd.Categorical(f['sample'],categories=levels).codes;pb=np.zeros((len(genes),len(levels)),dtype=float);scores=np.zeros((len(f),len(sets)));dets=np.zeros_like(scores);union=np.unique(np.concatenate(list(sets.values())));rev=np.full(len(genes),-1,dtype=int);rev[union]=np.arange(len(union));sidx={k:np.searchsorted(union,v) for k,v in sets.items()};detect=np.zeros((len(levels),len(union)),dtype=np.int32)
for i in range(len(f)):
 a,b=ptr[i:i+2];v=np.asarray(dat[a:b]);gi=np.asarray(ind[a:b]);pb[gi,codes[i]]+=v;vals=np.zeros(len(union));j=rev[gi];ok=j>=0;vals[j[ok]]=v[ok];sv=np.sort(v);left=np.searchsorted(sv,vals,'left');right=np.searchsorted(sv,vals,'right');eq=right-left;eq[vals==0]=len(genes)-np.count_nonzero(v);r=1+len(sv)-right+(eq-1)/2;detect[codes[i]]+=(vals>0)
 for k,(sn,si) in enumerate(sidx.items()):scores[i,k]=np.maximum(0,K-r[si]).sum()/sum(K-z for z in range(1,len(si)+1));dets[i,k]=np.count_nonzero(vals[si])
 if i%20000==0:print(i,'seconds',time.time()-started,flush=True)
err=np.max(np.abs(scores[:,list(sets).index('MP10')]-f.Top20_score));assert err<1e-12
for k,sn in enumerate(sets):f[sn]=scores[:,k];f[sn+'_detected']=dets[:,k]
save(f,'02_epi/cell_average_tie_scores.tsv.gz');save(pd.DataFrame(detect,index=levels,columns=np.array(genes)[union]).reset_index(names='sample'),'02_epi/gene_detection_counts.tsv.gz')
# Preserve distinct unmapped backgrounds, true same symbols summed only for pseudobulk.
pbf=pd.DataFrame(pb,index=gm,columns=levels).groupby(level=0,sort=False).sum();pbf.index.name='gene';save(pbf.reset_index(),'02_epi/pseudobulk_counts.tsv.gz');sm=f.groupby('sample',sort=False).agg(dataset=('dataset','first'),patient=('patient','first'),verified=('verified','first'),Type=('Type','first'),n_cells=('cell_id','size')).reset_index();save(sm,'02_epi/pseudobulk_samples.tsv');save(pd.DataFrame({'input_feature':genes,'mapped_gene':gm}),'00_plan_audit/epi_background_mapping.tsv')
(B/'02_epi/average_tie_execution.json').write_text(json.dumps(dict(n_cells=len(f),n_samples=len(sm),n_background=len(genes),n_sets=len(sets),seconds=time.time()-started,v1_MP10_max_error=err),indent=2));print('COMPLETE',time.time()-started,flush=True)
