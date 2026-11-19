from pathlib import Path
import importlib,sys
import numpy as np,pandas as pd,h5py
from scipy.stats import rankdata
smod=importlib.import_module('04_spatial');O=smod.O;S=smod.S;R=smod.R
# Standard GO membership, not genes selected from spatial associations.
gmt=Path('results/llps_C0_functional_state_v1/resources/c5.go.bp.v2025.1.Hs.symbols.gmt');terms={}
for line in gmt.read_text().splitlines():
 z=line.split('\t')
 if z[0] in ['GOBP_EXTRACELLULAR_MATRIX_ORGANIZATION','GOBP_COLLAGEN_FIBRIL_ORGANIZATION']:terms[z[0]]=z[2:]
ecm=set(sum(terms.values(),[]));sets={name:[g for g in z.gene if g not in ecm] for name,z in smod.program.groupby('program')};defs=[]
for name,z in smod.program.groupby('program'):
 for g in z.gene:defs.append(dict(program=name,gene=g,excluded_standard_ECM=g in ecm))
pd.DataFrame(defs).to_csv(O/'02_receiver_geneset/ECM_depleted_program_definition.tsv',sep='\t',index=False);rows=[];genes=sorted(set(sum(sets.values(),[])))
for sid in smod.M.index:
 print('ECM_SENSITIVITY',sid,flush=True);s=pd.read_parquet(S/'01_spatial_scores'/f'{sid}.parquet');d=pd.read_parquet(O/'07_spatial_expression_validation'/sid/'expression.parquet');assert s.spot_id.tolist()==d.spot_id.tolist();n=len(s)
 with h5py.File(R/smod.mp.loc[sid,'expression_relative_path']) as f:
  gn=list(smod.index(f['var']));look={g:i for i,g in enumerate(gn)};present=[g for g in genes if g in look];ix=[look[g] for g in present];chunks=[x[:,ix].toarray() for _,x in smod.blocks(f)]
 X=pd.DataFrame(np.log1p(np.vstack(chunks)/np.maximum(s.nCount.to_numpy()[:,None],1)*1e4),columns=present);rc,loc,W=smod.graph(s);cov=['major_mean__'+c for c in ['Fibroblast','Endothelial','Epithelial','Mural']]+['nCount','nFeature'];C=np.column_stack([np.ones(n)]+[rankdata(s[c])/n for c in cov]+[np.asarray(W[1].sum(1)).ravel()/6]);er=s.major_mean__Endothelial.to_numpy()>=s.major_mean__Endothelial.quantile(.75);WE=W[1].multiply(er.to_numpy()[None,:] if hasattr(er,'to_numpy') else er[None,:]).tocsr();den=np.asarray(WE.sum(1)).ravel();valid=den>0
 for name,gs in sets.items():
  have=[g for g in gs if g in X];v=X[have].mean(1).to_numpy() if len(have)>=max(5,len(gs)*.5) else np.full(n,np.nan);lag=np.divide(WE@v,den,out=np.zeros(n),where=valid);r=smod.residual(lag,C)
  for ax in smod.axes.itertuples():
   l=d[ax.ligand].to_numpy();rows.append(dict(section_id=sid,dataset=smod.M.loc[sid,'dataset'],patient_uid=smod.M.loc[sid,'patient_uid'],ligand=ax.ligand,program=name,n_program_genes=len(gs),n_covered=len(have),adjusted_rho=smod.rho(smod.residual(l,C)[valid],r[valid])))
pd.DataFrame(rows).to_csv(O/'08_patient_dataset_sensitivity/ECM_depleted_program_effects.tsv',sep='\t',index=False);print('ECM_SENSITIVITY_COMPLETE')
