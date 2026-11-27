from pathlib import Path
import pandas as pd,numpy as np,hashlib,subprocess,struct
O=Path(__file__).resolve().parents[1];V=O.parent/'llps_C5_cellchat_v2'; checks=[]
def ck(name,ok,detail=''):
 checks.append(dict(check=name,passed=bool(ok),detail=str(detail)));print(name,bool(ok),detail,flush=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
m=pd.read_csv(O/'00_input_audit/patient_sample_mapping.csv');ck('mapping_59_specimens_42_patients',len(m)==59 and m.patient_id.nunique()==42);ck('dataset_prefixed_patient_ids',all(str(r.patient_id).startswith(str(r.dataset)) for r in m.itertuples()))
for t,folder in [(10,'01_threshold10'),(20,'02_threshold20'),(30,'03_threshold30')]:
 s=pd.read_csv(V/folder/'network_all_specimens.csv');ck(f't{t}_strict_cell_boundary',(s.eligible==((s.sender_cells>t)&(s.receiver_cells>t))).all())
 ck(f't{t}_ineligible_NA',s.loc[~s.eligible,'weight'].isna().all());s=s[s.clinical_group!='Normal']
 ck(f't{t}_one_tumor_condition_per_patient',s.groupby('patient_id').clinical_group.nunique().max()==1)
 keys=['patient_id','dataset','clinical_group','sender','receiver'];p=s.groupby(keys,dropna=False).weight.median().reset_index();old=pd.read_csv(O/f'tables/patient_network_t{t}.csv');j=p.merge(old,on=keys,suffixes=('_new','_old'));ck(f't{t}_independent_patient_median_reconstruction',len(j)==len(old) and np.allclose(j.weight_new,j.weight_old,equal_nan=True))
 a=p.groupby(['sender','receiver']).weight.median();n=pd.read_csv(O/f'tables/network_summary_t{t}.csv').set_index(['sender','receiver']);ck(f't{t}_network_medians',np.allclose(a.reindex(n.index),n.patient_condition_median,equal_nan=True))
 # Reconstruct every C5-Endothelial detected/unobserved LR using the native specimen universe.
 lr=pd.read_csv(V/folder/'LR_detected_all_specimens.csv.gz');lr=lr[(lr.clinical_group!='Normal')&(((lr.sender=='Fibro_C5')&(lr.receiver=='Endothelial'))|((lr.receiver=='Fibro_C5')&(lr.sender=='Endothelial')))];fs=pd.read_csv(O/f'tables/C5_LR_t{t}.csv.gz');fs=fs[(fs.stratification=='tumor_only')&(((fs.sender=='Fibro_C5')&(fs.receiver=='Endothelial'))|((fs.receiver=='Fibro_C5')&(fs.sender=='Endothelial')))];bad=[]
 for (se,re),z in fs.groupby(['sender','receiver']):
  base=s[(s.sender==se)&(s.receiver==re)&s.eligible][['specimen_id','patient_id']];obs=lr[(lr.sender==se)&(lr.receiver==re)]
  wide=obs.pivot_table(index='specimen_id',columns='interaction_id',values='weight',aggfunc='sum').reindex(index=base.specimen_id,columns=z.feature).fillna(0);wide.index=base.patient_id.values
  pm=wide.groupby(level=0).median();pos=(wide>0).groupby(level=0).any().sum();zz=z.set_index('feature');ok=np.allclose(pm.median().reindex(zz.index),zz.patient_condition_median,equal_nan=True) and np.array_equal(pos.reindex(zz.index),zz.positive_patients) and (zz.eligible_patients==base.patient_id.nunique()).all()
  if not ok:bad.append((se,re))
 ck(f't{t}_all_Endothelial_LR_median_positive_denominator',not bad,bad)
c=pd.read_csv(O/'00_input_audit/category_sum_check.csv');ck('category_sum_matches_native_network',np.allclose(c.sum_category,c.network,equal_nan=True,atol=1e-8))
idx=pd.read_csv(O/'FIGURE_INDEX.csv');ck('unique_figure_names',not idx.figure.duplicated().any(),len(idx));missing=[str(O/r[c]) for _,r in idx.iterrows() for c in ['path','pdf_path','data_path','script'] if not (O/r[c]).is_file()];ck('all_figure_data_scripts_exist',not missing,missing)
baddpi=[];raster=[]
for _,r in idx.iterrows():
 p=O/r.path;b=p.read_bytes();off=8;dpi=None
 while off+12<=len(b):
  n=struct.unpack('>I',b[off:off+4])[0];typ=b[off+4:off+8]
  if typ==b'pHYs':
   x,y,u=struct.unpack('>IIB',b[off+8:off+17]);dpi=(x*.0254,y*.0254) if u==1 else None
  off+=12+n
 if dpi is None or min(dpi)<299:baddpi.append(r.figure)
 if r.figure_class!='Exploratory legacy figures':
  out=subprocess.run(['pdfimages','-list',str(O/r.pdf_path)],capture_output=True,text=True)
  if out.returncode or len(out.stdout.strip().splitlines())>2:raster.append(r.figure)
ck('PNG_300dpi_including_legacy',not baddpi,baddpi);ck('new_PDF_vector_no_embedded_rasters',not raster,raster)
l=pd.read_csv(O/'legacy_figure_mapping.csv');ck('all_28_legacy_figures_mapped',len(l)==28 and all(sha(O/r.archive_path)==r.source_sha256 for r in l.itertuples()))
inv=pd.read_csv(O/'00_input_audit/source_inventory.csv');changed=[]
for r in inv.itertuples():
 p=Path(r.path)
 if not p.exists() or p.stat().st_size!=r.bytes or sha(p)!=r.sha256:changed.append(str(p))
ck('historical_source_inventory_unchanged',not changed,f'{len(inv)} files checked; changed={changed}')
ck('new_main_data_no_highlow_labels',not any('C5-high' in (O/r.data_path).read_text() or 'C5-low' in (O/r.data_path).read_text() for r in idx.itertuples() if r.figure_class=='Main figure candidates'))
res=pd.DataFrame(checks);res.to_csv(O/'00_input_audit/final_QC_checks.csv',index=False)
text='# 最终 QC 报告\n\n'+f'{len(res)} 项检查，{int(res.passed.sum())} 项通过。\n\n'+'\n'.join(f'- {"PASS" if r.passed else "FAIL"} {r.check}: {r.detail}' for r in res.itertuples())+'\n\n独立复核覆盖三个阈值的完整网络和 C5–Endothelial 全部 LR 的患者中位数、阳性患者数与分母。绘图数量不表示独立实验数量。源文件哈希复核从本轮清单建立时起有效；原始对象未加载或写入。资源采用单线程串行包装器；没有声称测量完整流程峰值内存。\n'
(O/'FINAL_QC_REPORT.md').write_text(text)
if not res.passed.all():raise SystemExit('QC failed: review before delivery')
