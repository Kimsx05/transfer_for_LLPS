import os,sys,json
from pathlib import Path
import numpy as np,pandas as pd,h5py,pyarrow as pa
from scipy import sparse,stats
from scipy.spatial import cKDTree
pa.set_cpu_count(1);pa.set_io_thread_count(1)
O=Path('results/C5_Endothelial_NicheNet_spatial_v1');S=Path('results/C5_spatial_exploratory_all63_v1');R=Path('/home/jinshengxi/LLPS/spatial_data/LLPS_spatial_transfer_manual');sys.path.insert(0,str(R/'scripts'));from data_io import index,blocks
M=pd.read_csv(S/'00_inputs_and_QC/spatial_tissue_patient_audit.tsv',sep='\t').set_index('section_id');mp=pd.read_csv(S/'00_inputs_and_QC/section_mapping.tsv',sep='\t').set_index('section_id');axes=pd.read_csv(O/'07_spatial_expression_validation/frozen_axes_subunits.tsv',sep='\t');program=pd.read_csv(O/'02_receiver_geneset/frozen_program_genes.tsv',sep='\t');five=['MMP11','COL11A1','INHBA','COL5A1','CTHRC1'];ligands=axes.ligand.unique().tolist();program=program[program.program.isin(pd.read_csv(O/'04_NicheNet_ligand_activity/all_program_activities.tsv',sep='\t').program.unique()) & ~program.gene.isin(five+ligands)]
requested=pd.read_csv(O/'06_CellChat_integration/requested_axes_receptor_definitions.tsv',sep='\t');extra=set(requested.ligand)|set(';'.join(requested.required_subunits.dropna()).split(';'))
genes=sorted(extra|set(axes.ligand)|set(';'.join(axes.required_subunits).split(';'))|set(program.gene)|{'RAMP2','RAMP3'})
def save(x,p):pd.DataFrame(x).to_csv(p,sep='\t',index=False)
def rho(x,y):return float(stats.spearmanr(x,y).statistic) if len(x)>5 and np.std(x)>0 and np.std(y)>0 else np.nan
def residual(x,C):z=stats.rankdata(x)/len(x);return z-C@np.linalg.lstsq(C,z,rcond=None)[0]
def graph(s):
 rc=s[['array_row','array_col']].to_numpy().astype(int);lookup={tuple(v):i for i,v in enumerate(rc)};rows=[];cols=[]
 for i,(r,c) in enumerate(rc):
  for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)]:
   j=lookup.get((r+dr,c+dc));
   if j is not None:rows.append(i);cols.append(j)
 A=sparse.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(s),len(s)));B=(A@A)>0;B=B.astype(float)-B.multiply((A+sparse.eye(len(s)))>0);B.eliminate_zeros();return rc,lookup,[sparse.eye(len(s),format='csr'),A,B]
def smd(x,a,b):
 
 if len(a)<2 or len(b)<2:return np.nan
 den=np.sqrt((np.var(x[a],ddof=1)+np.var(x[b],ddof=1))/2);delta=np.mean(x[a])-np.mean(x[b]);return delta/den if den>0 else (0.0 if abs(delta)<1e-12 else np.inf)
pilot=pd.read_csv(S/'00_inputs_and_QC/pilot_prespecified.tsv',sep='\t').section_id.tolist();order=pilot+[x for x in M.index if x not in pilot]
for sid in order:
 dest=O/'07_spatial_expression_validation'/sid;dest.mkdir(exist_ok=True)
 if (dest/'DONE').exists():continue
 print('SPATIAL',sid,flush=True);s=pd.read_parquet(S/'01_spatial_scores'/f'{sid}.parquet');n=len(s);info=dict(section_id=sid,**M.loc[sid,['dataset','patient_uid','patient_id_status','tissue']].to_dict())
 with h5py.File(R/mp.loc[sid,'expression_relative_path']) as f:
  gn=list(index(f['var']));obs=list(index(f['obs']));assert obs==s.barcode.tolist();lookupgene={g:i for i,g in enumerate(gn)};present=[g for g in genes if g in lookupgene];ix=[lookupgene[g] for g in present];chunks=[];tot=[]
  for _,x in blocks(f):chunks.append(x[:,ix].toarray());tot.extend(np.asarray(x.sum(1)).ravel())
  raw=np.vstack(chunks);assert np.array_equal(tot,s.nCount.to_numpy())
 X=pd.DataFrame(np.log1p(raw/np.maximum(s.nCount.to_numpy()[:,None],1)*1e4),columns=present);D=pd.DataFrame(raw>0,columns=present)
 save([dict(**info,gene=g,present=g in present,detection=D[g].mean() if g in present else np.nan) for g in genes],dest/'coverage.tsv')
 rc,loc,Ws=graph(s)
 oldroot=Path('results/C5_spatial_colocalization_esda_v1/01_spatial_weights');nodes=pd.read_csv(oldroot/f'{sid}_nodes.tsv.gz',sep='\t');edges=pd.read_csv(oldroot/f'{sid}_edges.tsv.gz',sep='\t');nodes=nodes[nodes.analysis_index>=0].set_index('analysis_index');position=dict(zip(s.spot_id,range(n)));ii=nodes.spot_id.map(position);expected={(int(ii.loc[v.source]),int(ii.loc[v.target])) for v in edges.itertuples()};expected={tuple(sorted(v)) for v in expected};actual=set(zip(*sparse.triu(Ws[1],k=1).nonzero()));assert expected==actual,'Frozen graph mismatch'
 cov=['major_mean__'+c for c in ['Fibroblast','Endothelial','Epithelial','Mural']]+['nCount','nFeature'];C0=np.column_stack([stats.rankdata(s[c])/n for c in cov]+[np.asarray(Ws[1].sum(1)).ravel()/6]);C=np.column_stack([np.ones(n),C0]);Z=(C0-C0.mean(0))/np.maximum(C0.std(0),1e-8)
 endo=s.major_mean__Endothelial.to_numpy();er=endo>=np.quantile(endo,.75);fib=s.major_mean__Fibroblast.to_numpy();fr=fib>=np.quantile(fib,.75);score=s.five_UCell.to_numpy();ranks=stats.rankdata(score)/n
 plotted=s[['spot_id','x','y','five_UCell','major_mean__Endothelial','major_mean__Fibroblast','nCount']].copy()
 for g in set(axes.ligand)|set(';'.join(axes.required_subunits).split(';'))|{'RAMP2','RAMP3'}:
  plotted[g]=X[g] if g in X else np.nan
 for name,p in program.groupby('program'):
  gs=[g for g in p.gene if g in X];plotted[name]=X[gs].mean(1) if len(gs)>=max(5,len(p)*.5) else np.nan
 rows=[];bal=[];matchrows=[];rng=np.random.default_rng(123)
 # Rigid non-wrapping array translations preserve geometry within retained overlap.
 shifts=[]
 for dr in range(-24,25,4):
  for dc in range(-48,49,8):
   if dr==dc==0:continue
   src=[];dst=[]
   for i,(r,c) in enumerate(rc):
    j=loc.get((r+dr,c+dc))
    if j is not None:src.append(i);dst.append(j)
   if len(src)>=n*.5:shifts.append((np.array(src),np.array(dst),dr,dc))
 if len(shifts)>99:shifts=[shifts[i] for i in rng.choice(len(shifts),99,False)]
 for threshold in [.05,.1,.2]:
  high=(score>np.quantile(score,1-threshold))&(score>0)&fr;low=(ranks<=.5)&fr;aa=np.flatnonzero(high);pool=np.flatnonzero(low);a=[];b=[]
  if len(aa) and len(pool):
   tree=cKDTree(Z[pool]);dist,nn=tree.query(Z[aa],k=min(10,len(pool)));dist=np.atleast_2d(dist) if len(aa)==1 else dist;nn=np.atleast_2d(nn) if len(aa)==1 else nn
   if dist.ndim==1:dist=dist[:,None];nn=nn[:,None]
   used=set()
   for i,dd,jj in zip(aa,dist,nn):
    for dd0,j in zip(np.atleast_1d(dd),np.atleast_1d(jj)):
     j=int(pool[j])
     if dd0<=1 and j not in used:a.append(i);b.append(j);used.add(j);break
  a=np.asarray(a,int);b=np.asarray(b,int);smds=[]
  for k,c in enumerate(cov+['degree']):
   v=C0[:,k];val=smd(v,a,b);smds.append(val);bal.append(dict(**info,top_fraction=threshold,covariate=c,SMD_before=smd(v,aa,pool),SMD_after=val,n_high=len(aa),n_pairs=len(a),match_rate=len(a)/max(len(aa),1)))
  balanced=len(a)>=10 and np.all(np.isfinite(smds)) and max(abs(np.asarray(smds)))<.2
  for i,j in zip(a,b):matchrows.append(dict(top_fraction=threshold,high_spot=s.spot_id.iloc[i],control_spot=s.spot_id.iloc[j]))
  for ax in axes.itertuples():
   sub=ax.required_subunits.split(';');available=ax.ligand in X and all(g in X for g in sub)
   if not available:
    rows.append(dict(**info,ligand=ax.ligand,receptor=ax.receptor,role=ax.role,top_fraction=threshold,status='MISSING_FEATURE',ring=1));continue
   l=X[ax.ligand].to_numpy();rec=X[sub].min(1).to_numpy();joint=D[sub].all(1).to_numpy();lres=residual(l,C);rres=residual(rec,C)
   for ring,W in enumerate(Ws):
    WE=W.multiply(er[None,:]).tocsr();den=np.asarray(WE.sum(1)).ravel();valid=den>0;lag=np.divide(WE@rec,den,out=np.zeros(n),where=valid);lagres=np.divide(WE@rres,den,out=np.zeros(n),where=valid);support=np.divide(WE@joint.astype(float),den,out=np.zeros(n),where=valid);mask=high&valid;ix0=np.flatnonzero(valid)
    pairvalid=valid[a]&valid[b] if len(a) else np.array([],bool);ap=a[pairvalid];bp=b[pairvalid];product=l*lag;effect=np.mean(product[ap]-product[bp]) if len(ap)>=5 else np.nan
    trans=[]
    if threshold==.1 and ring==1:
     observed=rho(lres[ix0],lagres[ix0])
     for src,dst,dr,dc in shifts:
      keep=valid[src];trans.append(rho(lres[dst[keep]],lagres[src[keep]]))
    trans=np.array(trans);trans=trans[np.isfinite(trans)];observed=rho(lres[ix0],lagres[ix0]);pv=(1+np.sum(np.abs(trans)>=abs(observed)))/(len(trans)+1) if len(trans) else np.nan
    rows.append(dict(**info,ligand=ax.ligand,receptor=ax.receptor,role=ax.role,top_fraction=threshold,ring=ring,status='OK',n_senders=mask.sum(),n_receiver_spots=er.sum(),n_pairs=len(ap),balanced=balanced,ligand_detection_sender=np.mean(l[mask]>0) if mask.any() else np.nan,receptor_joint_detection_receiver=joint[er].mean(),sender_q05_positive=np.mean(s.major_q05__Fibroblast.to_numpy()[high]>0) if high.any() else np.nan,receiver_q05_positive=np.mean(s.major_q05__Endothelial.to_numpy()[er]>0),adjacent_joint_support=np.mean((l[mask]>0)*support[mask]) if mask.any() else np.nan,mean_ligand_sender=np.mean(l[mask]) if mask.any() else np.nan,mean_receptor_neighbor=np.mean(lag[mask]) if mask.any() else np.nan,raw_rho=rho(l[ix0],lag[ix0]),adjusted_rho=observed,matched_product_effect=effect,detection_matched_n=int(np.sum(((l[ap]>0)==(l[bp]>0))&(np.abs(l[ap]-l[bp])<.5))),detection_matched_receptor_effect=float(np.mean((lag[ap]-lag[bp])[((l[ap]>0)==(l[bp]>0))&(np.abs(l[ap]-l[bp])<.5)])) if np.sum(((l[ap]>0)==(l[bp]>0))&(np.abs(l[ap]-l[bp])<.5))>=5 else np.nan,translation_p=pv,n_translations=len(trans)))
    if threshold==.1:
     for name,p in program.groupby('program'):
      v=plotted[name].to_numpy();pr=np.divide(WE@v,den,out=np.zeros(n),where=valid);rows.append(dict(**info,ligand=ax.ligand,receptor=name,role='curated_program',top_fraction=threshold,ring=ring,status='OK',n_senders=mask.sum(),n_pairs=len(ap),balanced=balanced,raw_rho=rho(l[ix0],pr[ix0]),adjusted_rho=rho(lres[ix0],residual(pr,C)[ix0]),matched_product_effect=np.mean(pr[ap]-pr[bp]) if len(ap)>=5 else np.nan))
 save(rows,dest/'effects.tsv');save(bal,dest/'balance.tsv');save(matchrows,dest/'matched_spots.tsv');plotted.to_parquet(dest/'expression.parquet',index=False);(dest/'DONE').write_text('QC: exact barcode and count sum; six-neighbor graph; no missing gene imputation\n')
print('SPATIAL_COMPLETE',flush=True)
