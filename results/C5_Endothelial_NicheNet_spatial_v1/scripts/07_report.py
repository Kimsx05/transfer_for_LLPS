from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
O=Path('results/C5_Endothelial_NicheNet_spatial_v1')
def read(p):return pd.read_csv(O/p,sep='\t')
def texttable(d):return d.to_string(index=False,float_format=lambda x:f'{x:.4g}')
a=read('04_NicheNet_ligand_activity/all_program_activities.tsv');named=read('04_NicheNet_ligand_activity/requested_ligand_ranks.tsv');sp=read('tables/spatial_effects.tsv');ci=read('08_patient_dataset_sensitivity/patient_bootstrap_LODO.tsv');gate=read('02_receiver_geneset/gate.tsv');ax=read('07_spatial_expression_validation/frozen_axes_subunits.tsv');cov=read('tables/spatial_gene_coverage.tsv');bal=read('tables/matching_balance.tsv');spe=read('03_candidate_ligands/specificity_summary.tsv');program=read('02_receiver_geneset/GO_program_status.tsv')
primary=ci[(ci.omitted_dataset=='NONE')&(ci.top_fraction==.1)&(ci.ring==1)&(ci.role!='curated_program')];sliceprimary=sp[(sp.top_fraction==.1)&(sp.ring==1)&(sp.role!='curated_program')];pools=read('03_candidate_ligands/expression_summary.tsv');selected=spe[spe.gene.isin(ax.ligand)]
valid=[]
for p in (O/'07_spatial_expression_validation').glob('*/DONE'):valid.append(dict(check=p.parent.name,passed=True))
for x in ['01_receiver_pseudobulk/DE_primary.tsv','04_NicheNet_ligand_activity/all_program_activities.tsv','06_CellChat_integration/full_evidence.tsv','08_patient_dataset_sensitivity/patient_bootstrap_LODO.tsv']:valid.append(dict(check=x,passed=(O/x).exists()))
pd.DataFrame(valid).to_csv(O/'00_input_audit/final_validation.tsv',sep='\t',index=False)
coverage=cov[cov.gene.isin(set(ax.ligand)|set(';'.join(ax.required_subunits).split(';')))].groupby('gene').agg(sections_present=('present','sum'),median_detection=('detection','median')).reset_index()
summary=sliceprimary.groupby('ligand').agg(n_sections=('section_id','nunique'),median_raw_rho=('raw_rho','median'),median_adjusted_rho=('adjusted_rho','median'),median_matched_effect=('matched_product_effect','median'),balanced_sections=('balanced','sum'),technically_evaluable_sections=('technical_evaluable','sum'),median_sender_detection=('ligand_detection_sender','median'),median_receiver_joint_detection=('receptor_joint_detection_receiver','median')).reset_index()
summary.to_csv(O/'tables/final_axis_summary.tsv',sep='\t',index=False)
report=f'''# C5 → Endothelial NicheNet 与 Visium 空间验证 v1

## 执行与证据边界
实际执行方案 B (`CURATED_EXPLORATORY`)。完整完成固定候选筛选、官方 nichenetr 2.2.1.1 human v2 配体排序、CellChat 整合及 {sp.section_id.nunique()} 张切片的针对性表达分析。Five-gene 固定且直接复用，仍是 C5-associated program，不是真实 C5 丰度。没有重跑 CellChat、cell2location、全局 Moran 或原 Five-gene 评分。

## 1. 患者级内皮响应能否建立？为什么选择方案 B？
冻结门槛 Fibro_C0–C6 >=50、Endothelial >=20，合并同患者的肿瘤取材，不混入 Normal。21 名患者合格，其中 InhouseData 9、GSE222315 8，另两来源共4。原始 counts 缓存逐样本/细胞类群数量与审计一致，患者身份重新映射。连续 C5 fraction，TMM/voom/limma，模型包含 dataset，效应每个队列分析内 fraction 的1 SD。
主模型11445个基因，最小BH-FDR=0.420835，没有达到FDR<0.10的正向基因。GSE222315单独存在6个正向FDR<0.10基因，但InhouseData与合并分析不支持构建预设响应。调整采样数、严格采样与leave-one-dataset-out完整表已保存。此结果是当前设计下缺乏可重复响应证据，不能证明没有真实影响。
Gate：
```
{texttable(gate)}
```
没有通过筛选的响应集，因此不进行无意义的响应集LOPO优化。临床阶段敏感性模型经审查为满秩(21患者、6参数、15残差自由度)，单独保存；主模型仍可能残留临床阶段与取材混杂。dataset解释约38.3%的C5 fraction总变异，但两主要来源内部均有变异。

## 2. 实际功能程序与 NicheNet 排名
标准 MSigDB human v2025.1，Angiogenesis、TGF-beta signaling 和 blood-vessel remodeling 分别分析，不合并为虚构DEG。Activation原始仅11基因，按冻结>=15背景匹配基因规则NOT_EVALUABLE。全部来源、SHA256及准确基因列表见02_receiver_geneset。
```
{texttable(program[['program','original_genes','matched_genes','status']])}
```
Endothelial背景7070基因；C5表达且受体有患者重复支持的246个broad ligands。 pooled detection>=10%，C5>=5细胞、Endo>=20细胞的至少2患者支持。三套程序的前10排名如下；分数是对先验网络/指定程序的优先排序，不是实际活性或显著性：
```
{texttable(a[a['rank']<=10][['program','test_ligand','rank','aupr_corrected','CellChat_supported']])}
```

## 3. 指定配体的实际排名
```
{texttable(named[['program','test_ligand','rank','aupr_corrected','CellChat_supported']])}
```
未出现在表中的配体应查 broad_LR_all.tsv 的表达/数据库门槛，不能解释为零活性。所有NicheNet-only、CellChat-only、共同支持和受体条目保留于完整整合表。CellChat使用tumor-only t20，t10/30患者分母完整保留；配体排名未按空间显著性调整。

## 4. 表达特异性与预测 targets
共同患者中 C5 vs C0/C1/C4/其余Fibro比较使用相同标准化RNA表达的患者均值差，完整连续效应与方向支持保存；不以逐细胞P值代替患者证据。
```
{texttable(selected)}
```
均值差不等于C5-specific；大量CAF共享的配体只能称shared CAF ligands。前10配体target表使用官方get_weighted_ligand_target_links(n=250)，展示先验top250 targets与指定程序交集；空白不是零调控潜力。target属于标准内皮相关功能程序，但不是在本队列观察到的C5诱导target。

## 5. 空间候选、覆盖与表达
候选在空间分析前按各程序标准化排名均值和CellChat支持确定，再取每配体CellChat中位强度最高的受体轴；FN1/COL1A1/COL1A2为固定广泛ECM对照。该聚合只用于候选优先级，单独程序排名仍完整保留。
```
{texttable(ax[['ligand','receptor','required_subunits','role']])}
```
受体复合体从CellChatDB展开必要亚基；复合体连续表达取最小亚基log1p(CP10k)，检测要求所有必要亚基检出。CALCRL单基因定义不能证明完整ADM受体，RAMP2/3另外审查。不同资源定义不能互相替代。
```
{texttable(coverage)}
```
所有63张按切片读取counts，barcode和UMI总数与旧评分精确对齐。缺feature标记NA/不可评估，不填成零。

## 6. 空间邻接是否支持、组成是否影响？
发送区域Five切片内top10%且Fib上四分位；接收为Endothelial上四分位，q05支持比例保留。原始阵列六邻接，Ring1主要、Ring0/2对照，无微米换算。连续配体与邻近Endo-rich受体表达相关，以及高分发送区域与对照的表达乘积差，均是混合spot的探索性关联。
```
{texttable(summary)}
```
rank-residual控制Fib、Endo、Epi、Mural、nCount、nFeature和图度数。匹配同切片Fib-rich低Five区，七维标准化秩距离<=1、无放回，平衡要求所有|SMD|<0.2且>=10对。完整匹配和不平衡结果都保留。另有配体检出一致且log表达差<0.5的匹配敏感性。
空间功能程序使用排除Five和候选ligand后的实际覆盖基因平均log1p(CP10k)，至少5基因且覆盖>=50%；不等同于内皮细胞内表达。
非环绕阵列平移最多99次，至少50%重叠，保持平移区域的局部结构但受边界和不规则组织限制。其P值仅为敏感性。BH在全部有平移检验的候选/对照×切片集合内校正。患者推断不使用spot数作为n。

## 7. 独立患者、数据集与论文可用性
患者级主汇总先应用检测充分性规则：发送区域至少10个ligand检出spot、内皮富集区至少10个必要受体亚基共同检出spot；匹配效应至少5对。规则在技术执行期间、跨切片聚合前固定；所有未通过切片仍完整保存描述结果，不替换候选。已确认患者先平均其合格多张切片效应；未知患者仅切片描述。患者内聚合后分dataset bootstrap 2000次，seed123；所有数据集/leave-one-dataset-out见08_patient_dataset_sensitivity。置信区间是探索性的未同时校正区间。各metric与LODO范围分别在全部轴、阈值、ring和程序检验内BH调整（较保守）。
```
{texttable(primary[['ligand','metric','n_patients','n_datasets','median_effect','ci_low','ci_high','BH_FDR']])}
```
可靠可写：实际表达、可测覆盖、固定程序的计算排名、按患者汇总的空间效应及不确定性。具有患者CI支持且跨来源/平移/匹配方向一致的结果才可作为有限探索性邻接证据；不能仅凭未调整相关或CellChat与NicheNet交集宣布通讯。
仅先验支持：所有方案B配体–功能程序关联本质上是候选生成，无观察到的C5相关内皮DEG支撑。空间检测低、组成调整后效应衰减/反向或LODO不稳的轴只能留作探索候选。没有任何轴在本分析中证明C5身份、细胞内target激活或因果通讯。

## 8. 最终判读与小文章取舍

**可靠表达/计算结果**：资源SHA256验证通过，246个配体、三套独立程序、63切片/265002spots完整分析。原始Angiogenesis 36基因中20匹配背景，TGF-beta 54中47，vascular-remodeling 54中21；后者仅39%进入内皮表达背景，排名含义已受到覆盖限制。ADM/MDK/SEMA3C/WNT2及FN1/COLLAGEN均保留真实排名，没有强制成为阳性。上述主要基因均在63张存在，但WNT2中位spot检测约9.8%，CALCRL/RAMP2/RAMP3约16.0%/13.6%/13.2%，feature存在不等于复合体稳定检出。

**有限患者支持的探索关联**：FN1–ITGA5/ITGB1调整后患者中位rho=0.074 (95%CI 0.035–0.180，29患者)，七来源切片中位方向均正，剔除各个有患者身份的数据集后CI仍在零以上。它是预先指定的共享ECM对照，最适合写为基质ECM与内皮富集区的描述性邻接证据，不能改写成C5特异机制。所有504个空间平移检验均未通过整体BH-FDR，因此严格保留空间结构零模型的不确定性。

TGFB1–ACVR1/TGFBR1/TGFBR2的调整后患者中位rho=0.046 (CI 0.018–0.078，26患者)，但移除GSE319536后CI跨零；C5相对C0/C4表达未稳定富集。可作为后续实验的**非C5特异、队列依赖**候选，不宜作为文章核心因果轴。CellChat复合体中ACVR1与NicheNet直接受体映射不完全一致，不能只写成所有资源完全一致的单受体轴。

**先验/功能程序支持为主**：HSPG2–DAG1的调整后CI跨零。COL3A1相对C0/C4表达偏高，但COL3A1–ADGRG1无NicheNet直接LR边，原始空间rho偏负、调整后接近零；只能区分配体target先验与具体受体轴证据。候选target包括JAG1、EPAS1、MEF2C、SERPINE1及SMAD3/SMAD7等标准程序成员，它们是先验预测而非本队列已观察到的C5诱导响应。

**空间支持不足/不支持普遍正向邻接**：LAMA1与WNT4在发送区域的中位检测率约10%；技术上可评估分别35/36张，患者推断分别18/16人。LAMA1调整后方向偏负，WNT4调整后CI跨零。不能将原始弱正相关作为稳定正向通讯证据。

**匹配的关键限制**：top10%主方案仅2/63切片满足所有|SMD|<0.2和>=10对，其中只有1名已确认患者。已核查两组相同常量图度数的SMD应为0，并保留修正审计；其余主要失衡仍来自Fibroblast组成。全匹配样本的乘积差与其P值仅作残余失衡下的描述，不可宣称已充分控制组成。平衡子集无法提供可靠患者级CI，不能据此证明任何轴。

论文可使用真实排名、表达覆盖、阴性/不稳定结果、FN1共享ECM的描述性空间证据。若要突出C5→Endothelial，应先改善独立患者参考和内皮响应设计，并用细胞身份明确的空间原位检测/功能干预验证；本轮没有建立可归因于C5的配体–受体–target因果链。保留现有正式C5定位NO-GO，不能用方案B推翻它。

## 9. 复现与文件导航
- 00_input_audit/PRESPECIFIED.md：查看DEG/排名前冻结的选择规则。
- 01_receiver_pseudobulk：原始counts小型汇总、全部模型与SE/P/FDR。
- 03_candidate_ligands：逐患者/样本表达、受体筛选、C5特异性。
- 04_NicheNet_ligand_activity：全部246配体×可评估程序、20次大小/检测分箱匹配随机程序（诊断而非显著性检验）。
- 06_CellChat_integration：三个阈值与所有交并集；complete_expression_activity_LR_target_evidence.tsv整合C5表达、内皮必要亚基、患者分母和targets。
- 07_spatial_expression_validation：每切片coverage、spot_id表达Parquet、效果/匹配表。
- 08_patient_dataset_sensitivity：表达门槛5/10/20%、C0/C4 sender、LODO候选保留、空间患者CI与LODO。
- figures：A–E核心图及原始绘图数据；D使用预先确定pilot和预先排名前两候选，不按空间显著性选图。
- scripts：顺序01_response.R →02_expression.R →03_nichenet.R →03d_GO.R →03b_subunits.R →03c_support_sensitivity.py →04_spatial.py →05_summary.py →图数据导出→06_plots.R→07_report.py。

资源最多4核，实际R/Python/BLAS串行，seed123；没有下载模型或修改原数据。正式C5空间定位NO-GO保持不变。
'''
(O/'COMPREHENSIVE_REPORT.md').write_text(report)
(O/'README.md').write_text('# C5–Endothelial NicheNet spatial v1\n\nResponse mode: CURATED_EXPLORATORY. See COMPREHENSIVE_REPORT.md and frozen protocol in 00_input_audit/PRESPECIFIED.md. All historical inputs read-only; no CellChat/deconvolution/Five scoring rerun. Reproduce with scripts/run_R.sh or scripts/run_python.sh in project root.\n')
print('REPORT_COMPLETE')
