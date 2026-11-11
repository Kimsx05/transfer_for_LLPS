import csv,gzip,json,subprocess,struct
from pathlib import Path
O=Path('/home/jinshengxi/bladder_metabolism/results/llps_C5_cellchat_v2')
def read(path):
 p=O/path
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f))
def table(rows,cols):
 return '| '+' | '.join(cols)+' |\n| '+' | '.join(['---']*len(cols))+' |\n'+'\n'.join('| '+' | '.join(str(r.get(c,'')) for c in cols)+' |' for r in rows)+'\n'
def fmt(rows,cols):
 for r in rows:
  for c in cols:
   if r.get(c) not in (None,'NA',''):
    r[c]=f'{float(r[c]):.4g}'
 return rows
candidate=read('05_threshold_comparison/C5_LR_candidate_priorities_t20.csv')
network=read('05_threshold_comparison/C5_partner_network_all_thresholds.csv')
status=[]
for t,folder in [(10,'01_threshold10'),(20,'02_threshold20'),(30,'03_threshold30')]:
 r=read(folder+'/sample_status.csv');net=read(folder+'/network_all_specimens.csv');s={x['specimen_id'] for x in r if x['status']=='SUCCESS'};pat={x['patient_id'] for x in net if x['specimen_id'] in s};c5=[x for x in net if x['sender']=='Fibro_C5' and int(x['sender_cells'])>t];status.append(dict(min_cells=t,successful_specimens=len(s),successful_patients=len(pat),C5_eligible_specimens=len({x['specimen_id'] for x in c5}),C5_eligible_patients=len({x['patient_id'] for x in c5})))
mainnet=[x for x in network if x['min_cells']=='20' and x['stratification']=='tumor_only']
outgoing=sorted([x for x in mainnet if x['sender']=='Fibro_C5'],key=lambda x:-float(x['patient_condition_median']) if x['patient_condition_median']!='NA' else 0)
incoming=sorted([x for x in mainnet if x['receiver']=='Fibro_C5'],key=lambda x:-float(x['patient_condition_median']) if x['patient_condition_median']!='NA' else 0)
fmt(outgoing+incoming,['patient_condition_median','patient_prevalence'])
focus={}
for other in ['Endothelial','Epithelial','Myeloid','T_NK']:
 focus[other]=[x for x in candidate if x['sender']=='Fibro_C5' and x['receiver']==other][:5]
pathways=read('05_threshold_comparison/pathway_rankings_all_thresholds.csv.gz')
pathfocus=[x for x in pathways if x['min_cells']=='20' and x['sender']=='Fibro_C5' and x['receiver'] in ['Endothelial','Epithelial'] and int(x['positive_patients'])>0]
pwshow=[]
for g in ['Endothelial','Epithelial']:pwshow += [x for x in pathfocus if x['receiver']==g][:5]
core=read('05_threshold_comparison/core_LR_support_all_thresholds.csv')
support={}
for r in core:
 if r['stratification']=='tumor_only':support.setdefault((r['sender'],r['receiver'],r['feature']),{})[int(r['min_cells'])]=r
stable=[]
for key,vals in support.items():
 if len(vals)==3 and all(int(r['positive_patients'])>=3 for r in vals.values()):
  stable.append(dict(sender=key[0],receiver=key[1],LR=key[2],t10=f"{vals[10]['positive_patients']}/{vals[10]['eligible_patients']}",t20=f"{vals[20]['positive_patients']}/{vals[20]['eligible_patients']}",t30=f"{vals[30]['positive_patients']}/{vals[30]['eligible_patients']}",positive_datasets_t20=vals[20]['positive_datasets']))
validation=read('00_audit/final_validation.csv')
# Verify vector export, image resolution and source-sized artifacts.
figchecks=[]
for p in sorted((O/'06_figures').glob('*.pdf')):
 lines=subprocess.check_output(['pdfimages','-list',str(p)],text=True).splitlines();assert len(lines)==2,p;figchecks.append(dict(path=str(p),check='PDF_vector_no_raster',pass_=True))
for p in sorted((O/'06_figures').glob('*.png')):
 b=p.read_bytes();i=b.index(b'pHYs');x,y,u=struct.unpack('>IIB',b[i+4:i+13]);assert u==1 and abs(x*.0254-300)<1,p;figchecks.append(dict(path=str(p),check='PNG_300dpi',pass_=True))
with (O/'00_audit/figure_validation.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['path','check','pass_']);w.writeheader();w.writerows(figchecks)
readme=f'''# LLPS C5 CellChat v2

日期：2026-11-11。seed=123。本轮完成 CellChat，不运行 NicheNet。原始对象、前期审计及旧 CellChat 结果未修改。

## 完成情况

{table(status,['min_cells','successful_specimens','successful_patients','C5_eligible_specimens','C5_eligible_patients'])}

成功推断样本数与 C5 可评估样本数不同：全部样本仍输出完整 14 群网络，C5 不符合阈值的通信必须为 NA。本轮最终验证 {len(validation)} 项通过；图形另检查 PDF 无栅格图像及 PNG 300 dpi。失败尝试和导出恢复日志保留于 00_audit，最终样本状态见各阈值 sample_status.csv。

## 输入、数据库与资源约束

最终 atlas：`/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/results/objects/major_annotated_with_LLPS77_AUC.qs`。

正式 fibro 对象：`/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order/global_res012_object.qs`。

患者映射来自 `/home/jinshengxi/bladder_metabolism/results/llps_C5_distribution_audit_20261110_v1/v2/patient_sample_mapping_final.csv`。59 样本、42 患者；716 个 C5 分布于 37 样本、31 患者。InhouseData 13 样本、12 患者，a/b 保持独立取材输入。

本轮重新加载原对象验证 QC、barcode、orig.ident、dataset、大类与正式 C0–C6 身份。RNA Assay5 layers 为 counts/data/scale.data；使用 data，LogNormalize scale.factor=10000。每样本首细胞全基因 normalized/counts 校验最大误差 0；旧 per-sample RNA 缓存的全部表达值和 barcode/14 群标签与本轮 atlas 精确一致后只读复用，未复制大型表达对象。最终有效 atlas=353623，通信输入=344377；7324 个 major Cycling 和 1922 个历史排除 stromal 不进入通信，正式 C6 保留，QC 标准不改变。

实际 CellChat 2.1.2。完整 CellChatDB.human=3233 条：Secreted Signaling 1280、ECM-Receptor 424、Cell-Cell Contact 535、Non-protein Signaling 994。旧流程曾排除非蛋白类且保存的是 min.cells=30 过滤对象，因此本轮没有复用旧通信概率，全部重新推断，不与旧版总强度作直接数值比较。数据库提取基因 1446，表达矩阵包含 1408 (97.37%)，其中 1407 实际表达；表达矩阵共 71683 基因，70813 非零表达。非蛋白类条目的 ligand/receptor 字段可能表示代谢合成或相关基因，不应全部解读为字面意义的分泌蛋白结合。

脚本通过 taskset 将推断工作固定于 CPU 0、1、2，每进程一个样本、future sequential、data.table/BLAS/OpenMP 等单线程。汇总/导出检查使用 CPU 3，整个分析限制在四个 CPU 内。环境变量详见 scripts/run_limited.sh 和 launch_workers.sh；无嵌套并行。resource_observations/resource_thread_details 还记录辅助进程的 OS 睡眠维护线程；它们继承同一单 CPU 亲和性，不是额外数据并行 worker。上游只保存 compact net/LR/meta 缓存，不保存重复 Seurat 或完整表达型 CellChat 对象。

## 参数和检出定义

每 specimen 独立：identifyOverExpressedGenes(do.fast=FALSE)，其他选择阈值为软件默认；computeCommunProb(type="truncatedMean",trim=0.1,population.size=FALSE,nboot=100,seed.use=123,raw.use=TRUE)。没有合并患者或样本来推断，没有 projection/downsampling。每样本共享一份未过滤 net，然后分别执行 filterCommunication、computeCommunProbPathway、aggregateNet。

实际 filterCommunication 使用 n<=min.cells 剔除，10/20/30 分别要求至少 11/21/31 个 sender 和 receiver 细胞。LR 数以数据库 interaction ID 计，另区分唯一 LR/唯一 pathway 与跨伙伴通信事件数，见 C5_LR_pathway_counts_all_thresholds.csv 及逐样本总强度表。LR 和 aggregateNet 使用原生 p<0.05；原生 pathway 计算采用 p<=0.05，二者边界差异明确保留。100 次置换使 p-value 分辨率有限，p=0 不是数学意义的零概率。

原始 CellChat p-value 是样本内模型置换值，不是独立患者显著性。C5 candidate CSV 另提供每样本每方向的候选 LR BH 校正值，仅作诊断，不替换本轮原生筛选，不声称解决候选选择或跨患者推断问题。没有 Normal vs Tumor 正式显著性检验，也没有将重复阈值当作独立验证。

## NA、零与完整 LR/pathway 输出

每阈值每样本有完整 14×14 network CSV、compact_networks.rds、全群 LR_detected CSV、全群 pathway_detected CSV 和 C5_LR_candidates.csv.gz。紧凑 RDS 保留全部候选 LR 数组、p-values、完整 pathway 数组、counts、eligibility 与最新样本信息。未过滤缓存位于 00_audit/unfiltered；其表达矩阵未保存。

不合格边：CSV weight/probability/cellchat_p_value=NA。原生过滤数组内部置零是 CellChat 实现，使用 compact RDS 时必须同时应用 eligible_groups 掩码。C5 candidate 表 raw_probability/p_value 是共享上游诊断值，不能越过 eligible 标志解释。

eligible 且无显著检出：权重 0。未入过表达候选集的数据库 LR 也属于本流程未检出，但没有有效置换 p-value，不伪造 p=1。每阈值 complete_universe.csv.gz 提供完整 C5-其他群数据库条目和检出状态；汇总表包括从未阳性的条目并给出正确分母/零值。patient-condition 明细只物化至少一份样本阳性的 feature；其余由完整 universe 和网络 eligibility 明确定义为 eligible 零或 ineligible NA。

## 患者汇总和队列异质性

先按 patient_id × clinical_group × dataset 汇总 specimen 强度中位数，保留 Normal/NMIBC/MIBC 分离。患者阳性定义为至少一份 eligible specimen 阳性；患者 prevalence=阳性独立患者/eligible 独立患者。overall 表的强度是 patient-condition 单元分布，不假称这些单元全部独立；独立患者计数仍按 patient_id 去重。分临床组和 tumor-only 表用于解释，核心图使用 min.cells=20 的 tumor-only 作为展示参考。没有据此将 20 固定为以后唯一阈值。

表中同时给 eligible/positive specimens、eligible/positive patients、prevalence、specimen median、patient-condition median/Q25/Q75/min/max、positive/eligible datasets。分层包括 overall、clinical_group、dataset、dataset_clinical、tumor_only、tumor_dataset。同一患者若跨 Normal/Tumor，不得将分组患者数相加作为总体独立 N。

候选优先级按患者支持、跨队列重复、prevalence、强度排序，无 POSTN 或其他通路先验。Tier A 要求至少 3 名阳性患者且至少两个 dataset 各有至少 2 名阳性患者；Tier B 至少 3 名患者但队列支持受限；Tier C 为稀疏描述。阈值敏感性讨论患者组成，不把共同样本概率相同解释为验证。`LR_dataset_heterogeneity_all_thresholds.csv.gz` 对至少 2 名 eligible 患者的队列报告 prevalence 范围、中位强度范围/比值和 dominant dataset 阳性患者占比；中位数含 0 时不人为添加伪计数计算比值。

## C5 发送和接收的主要伙伴

以下为 t20 肿瘤患者中位数强度排序；分母较小的条目须结合 eligible N 阅读。

Outgoing：

{table(outgoing[:8],['receiver','eligible_patients','positive_patients','patient_prevalence','patient_condition_median','positive_datasets'])}

Incoming：

{table(incoming[:8],['sender','eligible_patients','positive_patients','patient_prevalence','patient_condition_median','positive_datasets'])}

## 重点 pathway 和 LR

C5→Endothelial/Epithelial pathway (优先患者支持)：

{table(fmt(pwshow,['patient_prevalence','patient_condition_median']),['receiver','feature','eligible_patients','positive_patients','patient_prevalence','positive_datasets','patient_condition_median'])}
'''
for other,rows in focus.items():readme+=f'\nC5→{other}：\n\n'+table(fmt(rows,['patient_prevalence','patient_condition_median']),['ligand','receptor','pathway','annotation','eligible_patients','positive_patients','patient_prevalence','positive_datasets','evidence_tier'])
readme+='''
## 跨阈值患者支持

以下候选在三个阈值均有至少 3 名肿瘤患者检出；数字为阳性/eligible 患者。这仅表示患者范围改变后仍有支持，并不是独立重复实验。

'''+table(stable[:20],['sender','receiver','LR','t10','t20','t30','positive_datasets_t20'])
readme+='''
`core_LR_added_specimen_support.csv` 提供低阈值新增 specimen 的 C5 和伙伴细胞数、LR 是否阳性、强度以及是否新增患者。`shared_upstream_probability_check.csv` 逐样本验证共同 eligible 群的概率完全相同。重点对照患者支持变化而非重复概率相关性。同一 C5-eligible specimen 的总 outgoing/incoming 强度也可能因伙伴群新增 eligibility 而改变；总强度表同时保留 eligible_partners，区分样本/患者组成变化与伙伴组成变化。

## C5 vs C0 和其他 fibro 亚群

每阈值的 `C5_vs_fibro_common_specimen_strength.csv` 对每个亚群、伙伴、方向仅保留 C5、比较亚群和伙伴三者均 eligible 的共同 specimen。`C5_vs_fibro_common_partner_total_strength.csv` 对同一 specimen 的共同伙伴求和并报告伙伴数。患者 LR/pathway 表给共同患者的中位数、配对差、composition 和阳性数；不从不相同的总体样本覆盖推断 C5 强于 C0。此比较为探索性，没有跨患者假设检验。比较伙伴排除 C5 和对应比较亚群本身，保证比较的是同一第三方伙伴；C5 核心 outgoing/incoming 总强度排除自分泌，完整网络及候选表仍保留自分泌。

主要比较结论见下面的数据驱动摘要；所有比较亚群均保留，核心图只展示 C5 vs C0 的四个重点伙伴，图中点为患者-临床条件单元。

'''
cmp=read('02_threshold20/C5_vs_fibro_common_specimen_strength.csv')
# Report patient-level paired differences for focal tumor comparisons without patient-level p-values.
from statistics import median
agg={}
for r in cmp:
 if r['other']=='Fibro_C0' and r['partner'] in focus and r['clinical_group']!='Normal':agg.setdefault((r['partner'],r['direction'],r['patient_id'],r['clinical_group']),[]).append(float(r['difference']))
grp={}
for (partner,direction,patient,clinical),v in agg.items():grp.setdefault((partner,direction),[]).append(median(v))
cr=[dict(partner=k[0],direction=k[1],patient_condition_units=len(v),paired_median_difference=f'{median(v):.4g}',C5_greater_units=sum(x>0 for x in v)) for k,v in grp.items()]
readme+=table(cr,['partner','direction','patient_condition_units','paired_median_difference','C5_greater_units'])
readme+='''
## 研究限制和后续候选

C5 的 92.88% 细胞来自 InhouseData/GSE222315。Normal eligible specimen 在 10/20/30 仅 2/1/0；20 下 MIBC 主要来自 InhouseData，NMIBC 来自 GSE222315。队列与临床分组混杂，单队列高 prevalence 不能视为跨疾病阶段规律。低细胞数、文库深度、候选基因选择及复杂体表达都可能影响推断；仅做了细胞数阈值组成敏感性，没有新增细胞重采样稳定性分析。

后续 NicheNet 可优先审查 Tier A 中患者重复且在三个阈值仍受支持的、具有明确 ligand 和 receiver 转录反应假设的条目，再审查 Tier B 的队列限定轴。ECM/contact 和非蛋白条目不能不加区分地等同于可直接进入 NicheNet ligand activity 的候选；需再核对资源覆盖和 receiver geneset。以下清单是本轮真实结果筛选，不是已验证机制，也没有运行 NicheNet。

'''+table(candidate[:12],['sender','receiver','ligand','receptor','pathway','annotation','evidence_tier','positive_patients','support_by_dataset'])
readme+='\n按同一患者支持规则筛选的重点方向分泌信号假设 (用于后续资源核对，不代表已运行 NicheNet)：\n\n'+table(read('05_threshold_comparison/secreted_LR_followup_hypotheses.csv')[:16],['sender','receiver','ligand','receptor','pathway','positive_patients','eligible_patients','evidence_tier','support_by_dataset'])
readme+='''
## 图形、复现与断点

核心图共 7 张：t20 tumor-only outgoing/incoming circle、LR bubble、pathway strength、共同患者 C5 vs C0、跨阈值患者支持、跨 dataset 重复性。圈图显示至少 3 名 eligible 患者且中位强度>0 的边，所有节点保留；两个圈图共享线宽标度，0.4 mm 为可读性下限。未画出不意味着生物学不存在。图形遵循已安装 jins-routine-plotting 公共主题和导出器，白底 Arial 12 pt、线宽至少 0.4 mm、矢量 PDF、300 dpi PNG，有图例时输出 full/panel_only/legend_only PDF。图的确切数值见 06_figures/*display_data.csv。

顺序：01_preflight.R → launch_workers.sh → 02_run_specimens.R --export-only (如需恢复导出) → 03_summarize.R → 04_compare_plot.R → 05_validation.R → 06_report.py。R/图形任务用 run_limited.sh；推断期间其他计算固定 CPU 3，避免额外并行。缓存存在且配置一致时不重复 computeCommunProb；每阈值 SUCCESS.csv 存在时不重复阈值导出。中断后沿用同一目录恢复；输入/参数变更必须新目录，不套用旧缓存。

首批导出遇到 data.table 单逻辑列索引兼容性错误；已改为显式 ==TRUE，只恢复导出，不重算成功的概率推断。原始失败状态归档为 *.attempt_*，worker 与 export_repair 日志保留。启动器收尾的脚本读取冲突还中断了最后一份尚未保存概率的 T4469479 (S057)；仅这份未完成推断按原参数重启，其他 58 份复用缓存，记录见 worker_2_resume.log。最终状态和 validation 才代表交付完成状态。

全部路径、文件大小见 output_manifest.csv。大表采用 CSV.gz；compact RDS 用于保留原生完整数组。未保存新的大型 Seurat 或重复表达矩阵。
'''
readme+='''
## 研究结论与候选优先级

以下患者支持均指 tumor-only，避免 Normal/Tumor 重复患者混入解释。C5 与其他 fibro 亚群之间的总强度较高，但 C6/C4 等分母很小。在非 fibro 伙伴中，Endothelial 的双向强度和覆盖较充分；C5 outgoing 还涉及 T/NK、Myeloid 和 Epithelial。不同边的 eligible 分母不同，不能将强度排名直接视为偏好或细胞特异性。

C5→Endothelial 的核心为 COLLAGEN/FN1，包括 COL1A1/2–ITGA1_ITGB1 和 FN1–ITGA5_ITGB1，t20 均为 12/12 患者，来自 GSE222315 的 5 名及 InhouseData 的 7 名。C5→Epithelial 的胶原/FN1–SDC1 在 t20 为 9/9 患者，两队列分别 2/7 名。接收端不同且上皮覆盖较少，不能据此作两方向强弱的正式统计推断。

可优先讨论后续 ligand activity 验证的轴：

| 方向及候选 LR | t10 阳性/eligible 患者 | t20 | t30 |
| --- | --- | --- | --- |
| Myeloid→C5 SPP1–ITGAV_ITGB1 / ITGAV_ITGB5 | 15/15 | 13/13 | 10/10 |
| C5→Endothelial MDK–ITGA6_ITGB1 | 15/15 | 12/12 | 10/10 |
| C5→Endothelial ADM–CALCRL | 15/15 | 12/12 | 10/10 |
| C5→Epithelial WNT5A–FZD6 | 11/13 | 8/9 | 6/6 |
| C5→Epithelial TNFSF12–TNFRSF12A | 10/13 | 8/9 | 6/6 |

前 3 类在扩展患者范围后仍有普遍支持，可作为第一优先级的资源/receiver 反应核查对象。上皮方向 WNT5A/TNFSF12 为第二优先级：高阈值的 100% prevalence 在扩大患者范围后下降，提示患者异质性，不能仅展示 t30。这里的优先级是后续研究建议，没有运行 NicheNet，也没有验证 ligand-target 因果关系。APP–CD74 (C5→Myeloid，t20 13/13) 和 ECM 轴适合平行的接触/基质机制研究，不应仅因强度高就默认最适合 NicheNet。

单队列例子：C5→T/NK COL4A5–ITGA1_ITGB1 在 InhouseData 为 8/8，而 GSE222315 为 0/4；T/NK→C5 MIF–ACKR3 分别 7/8 与 0/4。这些是队列限定候选，不能从 pooled 支持外推普遍性。完整异质性表同时提供分母和各队列强度。

C5 vs C0 共同患者对照中，四个重点伙伴 outgoing 的配对差中位数均大于 0；Endothelial 为 +0.4419 (5 个患者-条件单元)，Epithelial 为 +0.779 (4 个)，Myeloid 为 +0.4638 (6 个)，T/NK 为 +0.3593 (6 个)。incoming 差较小且个别患者相反。MIBC 的 4 名共同患者中，C5→Endothelial 的 FN1 composition 中位数约 10.38%，C0 约 5.27%；C5→Epithelial 分别 9.97% 与 5.62%。COLLAGEN 的相对组成则相近，说明总强度提升不等同于全新的通路组成。共同患者主要来自 InhouseData，此结果仍不能认定为跨队列 C5 特异机制。

阈值从 30 降至 20 增加 4 个 C5-eligible specimen、3 名患者；20 降至 10 再增加 3 个 specimen、2 名患者。增益还包括伙伴群达到门槛的变化。建议以 20 作清晰展示的参考，同时保留 10 的覆盖扩展和 30 的严格子集结果；本次不以任何一个阈值取代其余结果。共同样本概率完全相同是共享上游推断的预期，不构成独立稳定性验证。

独立 Python 复核直接从逐样本 LR 输出重新计算了三个阈值下所有阳性 tumor-only LR 的 eligible/positive 患者数与中位数，均与 R 汇总一致，见 independent_patient_check_final.log。源对象修改时间的最后检查采用数值时间比较以避免字符串微秒精度差异，大小及修改时间均未改变。图形已人工检查，患者气泡面积以 0 为起点并使用整数图例，避免夸大 12 与 13 名患者的差别。
'''
(O/'README.md').write_text(readme)
files=sorted(p for p in O.rglob('*') if p.is_file())
with (O/'output_manifest.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['path','bytes']);w.writeheader();w.writerows(dict(path=str(p),bytes=p.stat().st_size) for p in files)
print(table(status,['min_cells','successful_specimens','successful_patients','C5_eligible_specimens','C5_eligible_patients']))
print('TOP CANDIDATES');print(table(candidate[:10],['sender','receiver','ligand','receptor','pathway','evidence_tier','positive_patients','support_by_dataset']))
print('COMMON C5 vs C0');print(table(cr,['partner','direction','patient_condition_units','paired_median_difference']))
print('REPORT COMPLETE',O/'README.md')
