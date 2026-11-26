import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']: os.environ[k]='1'
from pathlib import Path
import pandas as pd,numpy as np
O=Path('/home/jinshengxi/bladder_metabolism/results/llps_C5_highlow_cellchat_ligand_v1')
s=pd.read_csv(O/'04_pathway_LR_analysis/within_median_t10_statistics.csv.gz');m=s[s.stratum=='Tumor_ALL'];n=m[m.kind=='network']
def fmt(x):return 'NA' if pd.isna(x) else f'{x:.4g}'
def describe(r):
 return f"{r.direction} {r.partner}: n={int(r.n_patients)}, high={fmt(r.high_median)}, low={fmt(r.low_median)}, 配对中位差={fmt(r.difference_median)}, 95% CI [{fmt(r.CI_low)}, {fmt(r.CI_high)}], FDR={fmt(r.FDR)}, high>low患者比例={r.high_greater_fraction:.1%}。"
def table(d,cols):
 return '| '+' | '.join(cols)+' |\n| '+' | '.join(['---']*len(cols))+' |\n'+'\n'.join('| '+' | '.join(fmt(v) if isinstance(v,(float,np.floating)) else str(v) for v in row)+' |' for row in d[cols].itertuples(index=False,name=None))
t=['## 9. 十个研究问题的直接回答','### 1. 可评估患者数','主展示有12名肿瘤患者 (7 MIBC、5 NMIBC)，Normal另有1名。各边因partner数量要求有不同分母；Endothelial 12，Epithelial 10，Myeloid 12，T/NK 11，B 11，Plasma 10，Mast 8，Other fibroblasts 12，Pericyte 12。阈值5/20/30的肿瘤覆盖为14/4/3；后两者不作患者级检验。','### 2. 总体通讯是否增强']
t += [describe(r) for _,r in n[n.partner=='TOTAL'].iterrows()]
t += ['主方案总outgoing与incoming均得到患者配对FDR支持 (约0.021与0.020)，但不能称每位患者或每条轴都增强。上述检验是观察性患者配对证据；即使显著也不能排除RNA复杂度、分组依赖性或共同细胞状态造成的关联。','### 3. 哪些伙伴变化最明显']
t += [describe(r) for _,r in n[n.partner!='TOTAL'].sort_values('difference_median',key=abs,ascending=False).head(6).iterrows()]
t += ['按绝对配对强度变化展示，不代表跨不同伙伴分母的正式效应比较。','### 4. Endothelial是否稳定']
t += [describe(r) for _,r in n[n.partner=='Endothelial'].iterrows()]
z=s[(s.kind=='network')&(s.partner=='Endothelial')&s.stratum.str.startswith(('dataset:','exclude:'))]
t += [table(z,['stratum','direction','n_patients','difference_median','CI_low','CI_high','FDR']),'两个队列效应中位数同为正，但Endothelial双向总体差异仅InhouseData分别通过FDR约0.040和0.045，GSE222315分别约0.281和0.468；不能声称已获得两个队列的显著重复。需同时参考分层方向及患者数；小队列CI/P为NA不代表没有效应。患者间未拆分C5连续分析与本轮细胞内high/low并非同一问题，不能互相当作独立验证。','### 5. 通路、LR和ECM类别']
z=m[(m.kind=='category')&m.partner.isin(['TOTAL','Endothelial'])]
t += [table(z,['direction','partner','feature','n_patients','difference_median','CI_low','CI_high','FDR'])]
for kind in ['pathway','LR']:
 z=m[(m.kind==kind)&(m.partner=='Endothelial')].sort_values('difference_median',key=abs,ascending=False).head(12)
 t += [table(z,['direction','feature','n_patients','positive_high_patients','positive_low_patients','difference_median','FDR'])]
t += ['主方案没有单个LR或pathway通过全量BH-FDR<0.05，以下仅为效应和患者支持线索。ECM类别占较大的绝对强度变化，但只有incoming ECM类别在该类别检验集合中获FDR支持；outgoing ECM类别FDR约0.074，不能概括成所有ECM通讯显著增强。具体LR及通路的全量阴性/反向结果保留于统计表。新候选指本轮新的观察性线索，不声称新发现配体结合关系。','### 6. 哪些配体有患者级关联','患者C5至少10细胞的主模型包含15名患者。dataset调整后，HSPG2、FN1、COL1A1、POSTN、COL1A2、COL3A1、LAMA1呈正相关的重点候选FDR支持，WNT5A呈负相关。再加入患者平均RNA深度和基因复杂度后，仅HSPG2保持22个预设候选的FDR支持；在全部606个配体相关基因的探索性校正中，没有基因通过该深度调整模型的FDR<0.05。患者内部depth-adjusted关联支持POSTN、COL1A1、COL1A2、COL3A1、SEMA3C；各队列单独检验支持较弱，不能称为已跨队列验证。','### 7. 既有候选表现','POSTN的患者内调整相关中位数约0.164 (FDR约0.033)，患者间dataset调整关联FDR约0.008，但额外深度调整后不再支持。FN1患者间关联强，患者内调整相关未获FDR支持，paired pseudobulk high-low表达有重点候选FDR支持。INHBA、TGFB1、ANGPTL4没有建立重点候选FDR支持的连续LLPS关联；不能因旧NicheNet排名而称它们随LLPS普遍增强。POSTN–ITGAV/ITGB5在high和low均为7/12阳性，仅6/12患者high更高；HSPG2–DAG1也是两组均8/12阳性、仅6/12患者high更高。因此配体表达关联不能直接推出这两条内皮轴的功能增强。下表保留相关内皮LR的实际结果。']
focus=['POSTN','INHBA','FN1','TGFB1','ANGPTL4'];z=m[(m.kind=='LR')&(m.partner=='Endothelial')&(m.direction=='outgoing')&m.ligand.isin(focus)]
t += [table(z,['ligand','feature','n_patients','positive_high_patients','positive_low_patients','difference_median','FDR']),'### 8. 队列、细胞数与少数患者影响','主方案的12名肿瘤患者完全分布为InhouseData的7名MIBC和GSE222315的5名NMIBC，临床组与队列完全混杂，因此不能据此分离临床分期效应和队列效应。主方案两组细胞数最多相差1个，未额外进行不必要的等量下采样；全局中位数组成不同且可能不平衡，作为独立敏感性方案。high组平均nCount和nFeature偏高，不能将未经调整的表达或网络差异等同于LLPS独立效应。各阈值/分组及LODO见tables，通信逐患者留出范围见07_sensitivity_analysis/CellChat_network_leave_patient_out.csv。','### 9. 与既有NicheNet/空间证据一致和不一致的部分','ECM相关表达线索与旧基质程序证据在主题上相容，但不能证明来自同一C5细胞。POSTN在名义NicheNet Set B及多数LOPO中排名靠前，但该排名并不等同于FDR显著的实测内皮响应；也缺乏已验证C5-specific空间定位。INHBA/ANGPTL4的旧候选身份没有获得本轮连续LLPS关联支持。HSPG2的患者级关联值得关注，但旧空间调整后关联并非稳定显著。这些不一致均保留，不通过加权评分掩盖。','### 10. 优先验证的机制假说','优先研究“C5的LLPS相关转录状态与基质/配体表达程序是否共同变化，以及这些变化是否伴随特定内皮LR响应”。HSPG2适合验证深度调整后仍存在的患者间关联；POSTN及COL1A1/COL1A2/COL3A1/FN1适合以基质程序为背景检验患者内表达及具体受体轴，但必须控制RNA复杂度，COL1A1还需持续使用去基因评分。SEMA3C也适合作为补充的患者内状态候选：其调整后患者内表达关联有支持，SEMA3C–PLXND1在9/12患者呈high更高，但具体LR未通过全量FDR，只能探索验证。INHBA、TGFB1、ANGPTL4保留为旧候选对照，不因先验而提升优先级。最终需身份明确的原位测量、配体/受体阻断及LLPS相关扰动；现有数据不支持直接建立LLPS促进血管生成或EndMT的因果机制。']
p=O/'COMPREHENSIVE_REPORT.md';base=p.read_text().split('## 9. 十个研究问题的直接回答')[0];p.write_text(base+'\n\n'+'\n\n'.join(t)+'\n')
print('TEN QUESTIONS APPENDED')

with (O/'COMPREHENSIVE_REPORT.md').open('a') as f:
 f.write("\n\n## 10. 稳健性与技术限制补充\n\n排除中位数ties后主阈值结果方向与支持基本一致。阈值5覆盖14名患者，总outgoing和Endothelial outgoing仍获FDR支持，但两项incoming的FDR均约0.077；Top/Bottom30在阈值10仅7名患者，Endothelial outgoing/incoming FDR约0.075/0.187；全局中位数阈值10仅6名患者，对应约0.040/0.118。因此不能概括为所有分组和阈值都显著。20/30患者不足5名，仅描述。\n\n主方案Endothelial outgoing配对差与组间log(nCount)及log(nFeature)差的相关系数分别约0.617、0.742；incoming约0.246、0.400。RNA复杂度混杂尤其影响outgoing解释。深度调整模型的零深度差参考截距属于小样本参考点估计，不是因果调整后的CellChat效应，不据此宣称效应消失或得到独立验证。逐患者留出后Endothelial outgoing中位差范围0.223–0.238，incoming 0.251–0.284，方向均为正；这降低单一患者驱动中位数的担忧，但不消除队列和技术混杂。\n\n55个原生网络核验全部通过；最终1670项QC、132项独立pseudobulk核验及110项独立患者内关联核验通过。执行中断与恢复见00_input_audit/execution_recovery.md。原始对象、基因集与历史结果保持不变，未重跑NicheNet或空间分析。\n")
p=O/'COMPREHENSIVE_REPORT.md';s=p.read_text();title,body=s.split('\n',1);summary='''
本轮本地分析完成：55个联合网络，主阈值10有12名配对肿瘤患者。LLPS-high总体outgoing/incoming较高 (FDR约0.021/0.020)，但没有单个LR或pathway通过全量校正。Endothelial差异方向在两个队列一致，显著支持主要来自InhouseData；RNA复杂度及队列与临床组完全混杂限制机制解释。深度调整后HSPG2保留预设22候选的患者级关联支持，但全606基因筛选没有FDR支持。具体10个问题见第9节。

交付状态：本地报告、统计表、图片、代码和QC已完成。GitHub写入接口返回403，尚未上传。轻量结果包及上传范围见delivery/UPLOAD_STATUS.md；不含表达矩阵和重复对象。
''';p.write_text(title+'\n'+summary+'\n'+body)
