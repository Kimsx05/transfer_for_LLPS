import os,sys
from pathlib import Path
import numpy as np,pandas as pd,h5py,pyarrow as pa
from scipy import sparse
pa.set_cpu_count(1);pa.set_io_thread_count(1)
O=Path('results/C5_Endothelial_GSEA_nominalNicheNet_v1');S=Path('results/C5_spatial_exploratory_all63_v1');R=Path('/home/jinshengxi/LLPS/spatial_data/LLPS_spatial_transfer_manual');sys.path.insert(0,str(R/'scripts'));from data_io import index,blocks
axes=pd.read_csv(O/'06_CellChat_spatial_integration/frozen_new_spatial_axes.tsv',sep='\t');assert len(axes)<=3
M=pd.read_csv(S/'00_inputs_and_QC/spatial_tissue_patient_audit.tsv',sep='\t').set_index('section_id');mp=pd.read_csv(S/'00_inputs_and_QC/section_mapping.tsv',sep='\t').set_index('section_id');genes=sorted(set(axes.test_ligand)|set(';'.join(axes.required_subunits).split(';')))
dest=O/'06_CellChat_spatial_integration/light_spatial';dest.mkdir(exist_ok=True);allq=[];alla=[]
for sid in M.index:
 fn=dest/f'{sid}_axes.tsv';cf=dest/f'{sid}_coverage.tsv'
 if fn.exists() and cf.exists():alla.append(pd.read_csv(fn,sep='\t'));allq.append(pd.read_csv(cf,sep='\t'));continue
 s=pd.read_parquet(S/'01_spatial_scores'/f'{sid}.parquet');info=dict(section_id=sid,**M.loc[sid,['dataset','patient_uid','patient_id_status','tissue']].to_dict());n=len(s)
 with h5py.File(R/mp.loc[sid,'expression_relative_path']) as f:
  gn=list(index(f['var']));assert list(index(f['obs']))==s.barcode.tolist();gi={g:i for i,g in enumerate(gn)};present=[g for g in genes if g in gi];chunks=[];total=[]
  for _,x in blocks(f):chunks.append(x[:,[gi[g] for g in present]].toarray());total.extend(np.asarray(x.sum(1)).ravel())
  raw=np.vstack(chunks);assert np.array_equal(total,s.nCount)
 X=pd.DataFrame(np.log1p(raw/np.maximum(s.nCount.to_numpy()[:,None],1)*1e4),columns=present);D=X>0
 cov=pd.DataFrame([dict(**info,gene=g,present=g in present,detection=D[g].mean() if g in present else np.nan,n_spots=n) for g in genes]);cov.to_csv(cf,sep='\t',index=False);allq.append(cov)
 rc=s[['array_row','array_col']].to_numpy().astype(int);lookup={tuple(v):i for i,v in enumerate(rc)};ii=[];jj=[]
 for i,(r,c) in enumerate(rc):
  for dr,dc in [(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1)]:
   j=lookup.get((r+dr,c+dc))
   if j is not None:ii.append(i);jj.append(j)
 ii=np.array(ii);jj=np.array(jj);score=s.five_UCell.to_numpy();fr=s.major_mean__Fibroblast>=s.major_mean__Fibroblast.quantile(.75);er=s.major_mean__Endothelial>=s.major_mean__Endothelial.quantile(.75);high=(score>np.quantile(score,.9))&(score>0)&fr
 rows=[]
 for ax in axes.itertuples():
  sub=ax.required_subunits.split(';');base=dict(**info,ligand=ax.test_ligand,receptor=ax.receptor,required_subunits=ax.required_subunits,ligand_in_Fivegene=ax.test_ligand in ['MMP11','COL11A1','INHBA','COL5A1','CTHRC1'],n_senders=int(high.sum()),n_receivers=int(er.sum()),scope='LIGHTWEIGHT_EXPRESSION_AUDIT_NOT_SPATIAL_INFERENCE')
  if ax.test_ligand not in X or not all(g in X for g in sub):rows.append(dict(**base,status='MISSING_GENE'));continue
  joint=D[sub].all(axis=1).to_numpy();rec=X[sub].min(axis=1).to_numpy();lig=X[ax.test_ligand].to_numpy();edge=high.to_numpy()[ii]&er.to_numpy()[jj];src=ii[edge];tar=jj[edge]
  rows.append(dict(**base,status='EVALUABLE',ligand_detection_sender=float(np.mean(lig[high]>0)) if high.sum() else np.nan,ligand_mean_logCP10k_sender=float(np.mean(lig[high])) if high.sum() else np.nan,receptor_joint_detection_receiver=float(np.mean(joint[er])),receptor_min_subunit_mean_logCP10k_receiver=float(np.mean(rec[er])),n_directed_sender_receiver_edges=len(src),n_unique_sender_edges=len(np.unique(src)),n_unique_receiver_edges=len(np.unique(tar)),edge_joint_detection=float(np.mean((lig[src]>0)&joint[tar])) if len(src) else np.nan,edge_mean_ligand_minreceptor_product=float(np.mean(lig[src]*rec[tar])) if len(src) else np.nan,sender_fibro_q05_positive=float(np.mean(s.major_q05__Fibroblast[high]>0)) if high.sum() else np.nan,receiver_endo_q05_positive=float(np.mean(s.major_q05__Endothelial[er]>0))))
 d=pd.DataFrame(rows);d.to_csv(fn,sep='\t',index=False);alla.append(d);print('AUDIT',sid,flush=True)
q=pd.concat(allq);a=pd.concat(alla);q.to_csv(O/'06_CellChat_spatial_integration/new_axes_coverage.tsv',sep='\t',index=False);a.to_csv(O/'06_CellChat_spatial_integration/new_axes_expression_audit.tsv',sep='\t',index=False)
a.groupby(['ligand','dataset']).agg(n_sections=('section_id','size'),median_ligand_detection=('ligand_detection_sender','median'),median_joint_receptor_detection=('receptor_joint_detection_receiver','median'),median_edge_detection=('edge_joint_detection','median')).reset_index().to_csv(O/'06_CellChat_spatial_integration/new_axes_dataset_description.tsv',sep='\t',index=False)
print('LIGHT_SPATIAL_COMPLETE')
