from pathlib import Path
import pandas as pd,numpy as np
O=Path('results/C5_fivegene_LLPS77_spatial_v1');p=O/'COMPREHENSIVE_REPORT.md';s=p.read_text();needle='## A. Five-gene 表现'
text='''## 主要判断 (按证据强弱)

- **可靠描述**：Five-gene 在62/63张出现正向空间自相关；GSE171351_D1未获得聚集支持 (I=0.0138，FDR=0.245)。Five的零分率中位35.3%，高于A26的5.8%和Core15的1.6%；其Moran I及最大连通分量比例反而更低。因此不支持“Five整体比A26/Core15更连贯、更集中”的泛化判断；更稀疏、较局部的显示确实存在，但不能等同更好的生物学特异性。控制检测数量后仍有残余空间结构，也不能据此排除全部技术影响。
- **具有跨样本支持的探索性关联**：Five–LLPS原始相关患者中位0.318，调整后0.107 (95% CI0.090–0.205)。七来源的调整相关切片中位数均为正 (0.029–0.220)，但并非每张切片都为正。仅三个来源具有已确认患者；逐一删除这些来源后的患者效应仍为正，其中删除GSE319536后中位降至0.062 (CI0.031–0.205)。这是小幅条件关联，不能归因于真实C5或细胞内LLPS活性。
- **覆盖和基因贡献敏感性**：共同36基因与原始LLPS评分的切片内相关中位约0.960 (UCell)/0.968 (AUCell)，但36基因只保留46.8%的原列表，不是等效校正。去掉LLPS中的COL1A1后Five的患者调整关联仍为0.116 (CI0.074–0.182)；Five本身不含COL1A1，与LLPS77无重叠。逐一删除Five中的任一基因后，患者调整关联中位仍为0.068–0.108，CI均未跨零；删除COL5A1的衰减最大。全部九大类组成调整约0.127、q05组成调整约0.103、Fib-high内部调整约0.082，均保留正向患者区间。残余未测量混杂仍无法排除。
- **受基因/组织异质性限制**：GSE285715_T1/T2/T3的Five-high秩权重分别约89.4%/95.2%/93.4%来自COL5A1。这三张不能作为稳定多基因共同表达的代表，也不宜据此做C5正常/肿瘤结论。LLPS评分本身同样是若干可检测基因的组合，而不是77个基因等量贡献。
- **空间共定位有限支持**：top10双高Jaccard的来源中位数约0.065–0.396。调整后的平移敏感性在36/63张达到近似BH<0.05，但平移对不规则边界仅近似保留结构，不能当作严格空间因果检验。双高区域跨来源可被描述，位置及强度并不统一，也没有C5真值标签。
- **微环境不能概括为促血管/免疫抑制**：控制上皮丰度后，上皮距离的患者相关为−0.044 (CI−0.104–0.007；FDR0.158)，没有稳健的独立远近关系。Five-high同spot Endothelial在主阈值下较低 (患者log1p差−0.035，CI−0.056至−0.013；FDR0.00095)，六个来源中位同向，但top20%效应区间跨零，不能称为阈值无关的普遍特征。Five-high Myeloid主结果区间跨零。双高相对Five-high/LLPS非高的Myeloid降低在top10%仅边界支持 (FDR约0.050)、top20%较强、top5%不明确，属于阈值依赖候选关联。
- **与A26的一致和差异**：同spot Endothelial/Myeloid效应方向分别在约84.1%/76.2%的切片与历史A26一致；不能据此认为相同细胞身份。历史A26全匹配患者CI跨零，而Five的Endothelial主比较较明确，但来源和阈值依赖仍在。ring1–3及平衡子集见完整表，不选择单一有利ring作为结论。

**下一步最值得独立验证的是“Five相关基质转录程序与LLPS-associated程序的小幅条件共现”，而不是预设的促血管生成轴。** 可优先选择有多基因贡献、匹配平衡、跨来源同向的切片进行独立多重空间RNA/蛋白验证；不能把本轮筛出的切片重新称作独立验证集。是否开展C5-specific解卷积仍需改善参考可辨识性，正式NO-GO不变。

'''
s=s.replace(needle,text+needle);s=s.replace('所有不足平衡切片见 balance.tsv。','所有不足平衡切片见 `03_Fivegene_controls/insufficient_balance_sections.tsv`，包括未能匹配的情况。')
s += '\n## 附加组成/覆盖与逐基因敏感性汇总\n\n```text\n'+pd.read_csv(O/'07_dataset_heterogeneity/additional_sensitivity_summary.tsv',sep='\t').round(4).to_csv(sep='\t',index=False)+'```\n\n全部主要/敏感性患者检验集合：\n\n```text\n'+pd.read_csv(O/'00_input_audit/patient_BH_families.tsv',sep='\t').to_csv(sep='\t',index=False)+'```\n\n说明：高覆盖子集仍含29位可确认患者，并不构成低覆盖未知患者来源的独立验证。切片内百分位对零分ties也赋予正的中位排名，不能将非零百分位当作基因检出。\n'
p.write_text(s)
# Concrete stage-specific entry reports.
stages={
'01_LLPS77_spatial_scores': '63张全部计算完成；原始77、共同36、去COL1A1的76基因方案。UCell/AUCell中位相关0.942；每spot检测数的切片中位范围4–31。共同36不是等效校正；无方向基因列表不代表LLPS活性。详见qc.tsv、coverage_overlap.tsv及综合报告B部分。',
'02_Fivegene_spatial_distribution':'62/63张正向聚集；GSE171351_D1未支持。Five零分率35.3%、Moran中位0.584。GSE285715_T1/T2/T3由COL5A1主导。见clustering.tsv、gene_contribution.tsv和全部genes空间图。',
'03_Fivegene_controls':'Fib/UMI/Epi主匹配成功率中位93.9%；35/63张达到全部平衡标准。完整SMD前后、失败切片、主/阈值敏感性配对ID均保留。Five仍与Fib、深度及竞争程序相关，残余不等于细胞特异表达。',
'04_Fivegene_LLPS77_association':'患者主调整关联0.107 (95% CI0.090–0.205，29患者)，较原始0.318明显减弱。无Five/LLPS重叠；共同36、去COL1A1、Five逐基因删除和组成敏感性方向保留。见association.tsv及07目录患者/LODO表。',
'05_spatial_colocalization':'连续邻域、top5/10/20重叠及199次目标平移敏感性全部保存；部分切片合格平移少于199。36/63张调整关联的近似空间敏感性BH<0.05，但不规则边界使其不是精确零模型。统称Five-gene-high / LLPS77-high spatial regions。',
'06_epithelial_and_TME_neighborhood':'调整Epithelial后上皮距离关系不明确。Five-high Endothelial主same-spot差为负，top20不明确；Myeloid主结果不稳健。双高Myeloid降低依赖阈值。完整ring0–3、全部匹配和平衡子集均保留；不解释为通讯/因果。',
'07_dataset_heterogeneity':'七来源描述完整，只有三个来源可做29位已确认患者汇总。七来源调整相关中位均正；已知来源LODO仍正。未确认患者不补全；Normal3/Tumor17/UNRESOLVED43，Normal患者未知。不能将未知来源LODO中不变的患者集合视为新的独立验证。'}
for name,body in stages.items():(O/name/'REPORT.md').write_text('# '+name+'\n\n'+body+'\n\n完整方法、数值与限制见[综合报告](../COMPREHENSIVE_REPORT.md)。正式C5-specific空间定位仍为NO-GO。\n')
print('Evidence interpretation and stage reports updated')
