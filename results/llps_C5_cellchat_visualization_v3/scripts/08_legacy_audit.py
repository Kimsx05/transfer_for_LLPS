import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import pandas as pd,hashlib,shutil,json
B=Path('/home/jinshengxi/bladder_metabolism/results');O=B/'llps_C5_cellchat_visualization_v3';H=B/'llps_C5_highlow_cellchat_ligand_v1';dest=O/'10_legacy_highlow_exploratory';rows=[];ix=[]
for p in sorted((H/'figures').iterdir()):
 if p.is_file() and p.suffix in ['.png','.pdf','.csv']:shutil.copy2(p,dest/p.name)
mapping={'CellChat_circle_high':'01_circle_plot/C5_TME_outgoing.png','CellChat_circle_low':'01_circle_plot/C5_TME_incoming.png','communication_dual_triangle':'02_dual_triangle_heatmap/fibro_TME_numbers.png','communication_difference_forest':'03_communication_strength/C5_strength_IQR_forest.png','pathway_high_low_heatmap':'04_pathway_heatmap/C5_outgoing_pathway_top15.png','LR_high_low_bubble':'05_LR_bubble/C5_outgoing_data_driven.png','integration_CellChat':'08_multi_evidence_integration/evidence_CellChat_strength.png','integration_expression':'08_multi_evidence_integration/evidence_Expression_difference.png','integration_correlation':'08_multi_evidence_integration/evidence_Exploratory_correlation.png','integration_NicheNet_rank':'08_multi_evidence_integration/evidence_NicheNet_rank.png','integration_spatial':'08_multi_evidence_integration/evidence_Spatial.png','multi_evidence_integration':'08_multi_evidence_integration/multi_evidence_complete_C5.png'}
for p in sorted((H/'figures').glob('*.png')):
 name=p.stem;cont=name=='LLPS_ligand_correlation_heatmap' or name.startswith('ligand_');tech=name=='technical_covariate_paired_QC';new=mapping.get(name,str((dest/p.name).relative_to(O)))
 action='Retain continuous LLPS association as exploratory; unchanged statistics' if cont else 'Archive high/low-specific technical QC; no full-C5 reinterpretation' if tech else 'Archive original and redesign using formal full-C5 v2 data'
 rows.append(dict(legacy_figure=name,legacy_source=str(p),archive_path=str((dest/p.name).relative_to(O)),source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),treatment=action,new_corresponding_figure=new,semantics='Original statistical meaning retained; no transplanted high/low network values'))
 ix.append(dict(figure='legacy_'+name,path=str((dest/p.name).relative_to(O)),pdf_path=str((dest/(name+'.pdf')).relative_to(O)),data_path=str((dest/(name+'_data.csv')).relative_to(O)),script='scripts/08_legacy_audit.py',data_source=str(H/'figures'/p.name),threshold='Historical, see original report',statistical_unit='Historical cell/patient or paired high/low',purpose=action,limitations='Exploratory legacy; not part of full C5 main network',figure_class='Exploratory legacy figures'))
pd.DataFrame(rows).to_csv(O/'legacy_figure_mapping.csv',index=False);pd.DataFrame(ix).to_csv(O/'00_input_audit/index_08_legacy.csv',index=False)
# Reproducible source inventory, excluding large objects. All entries are read-only inputs.
roots=['llps_C5_cellchat_v2','llps_C5_highlow_cellchat_ligand_v1','C5_Endothelial_NicheNet_spatial_v1','C5_Endothelial_GSEA_nominalNicheNet_v1','C5_fivegene_LLPS77_spatial_v1','C5_spatial_colocalization_esda_v1','llps_C5_parallel_v1','llps_C5_distribution_audit_20261110_v1/v2'];sources=[]
for root in roots:
 d=B/root
 if not d.exists():raise FileNotFoundError(d)
 for p in sorted(d.rglob('*')):
  if not p.is_file() or p.suffix.lower() in ['.rds','.qs','.qs2','.png','.pdf','.zip','.pyc'] or p.stat().st_size>30*2**20:continue
  if p.name.endswith('.log') or 'provisional' in str(p) or '/delivery/' in str(p):continue
  sources.append(dict(path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
pd.DataFrame(sources).to_csv(O/'00_input_audit/source_inventory.csv',index=False)
(O/'10_legacy_highlow_exploratory/README.md').write_text('本目录完整保留旧版28张图及原始PDF、绘图数据。连续相关仅作探索；high/low通信及技术QC仅归档，不作为本轮完整C5主分析。逐图去向见 ../legacy_figure_mapping.csv。未重新计算高低分组或模型。\n')
print('LEGACY MAPPED',len(rows),'SOURCE FILES',len(sources))
