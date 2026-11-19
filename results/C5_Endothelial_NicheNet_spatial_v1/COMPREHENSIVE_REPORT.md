# C5 → Endothelial NicheNet 与 Visium 空间验证 v1

## 执行与证据边界
实际执行方案 B (`CURATED_EXPLORATORY`)。完整完成固定候选筛选、官方 nichenetr 2.2.1.1 human v2 配体排序、CellChat 整合及 63 张切片的针对性表达分析。Five-gene 固定且直接复用，仍是 C5-associated program，不是真实 C5 丰度。没有重跑 CellChat、cell2location、全局 Moran 或原 Five-gene 评分。

## 1. 患者级内皮响应能否建立？为什么选择方案 B？
冻结门槛 Fibro_C0–C6 >=50、Endothelial >=20，合并同患者的肿瘤取材，不混入 Normal。21 名患者合格，其中 InhouseData 9、GSE222315 8，另两来源共4。原始 counts 缓存逐样本/细胞类群数量与审计一致，患者身份重新映射。连续 C5 fraction，TMM/voom/limma，模型包含 dataset，效应每个队列分析内 fraction 的1 SD。
主模型11445个基因，最小BH-FDR=0.420835，没有达到FDR<0.10的正向基因。GSE222315单独存在6个正向FDR<0.10基因，但InhouseData与合并分析不支持构建预设响应。调整采样数、严格采样与leave-one-dataset-out完整表已保存。此结果是当前设计下缺乏可重复响应证据，不能证明没有真实影响。
Gate：
```
 n_eligible  datasets_ge4  positive_response_genes  gate
         21             2                        0 False
```
没有通过筛选的响应集，因此不进行无意义的响应集LOPO优化。临床阶段敏感性模型经审查为满秩(21患者、6参数、15残差自由度)，单独保存；主模型仍可能残留临床阶段与取材混杂。dataset解释约38.3%的C5 fraction总变异，但两主要来源内部均有变异。

## 2. 实际功能程序与 NicheNet 排名
标准 MSigDB human v2025.1，Angiogenesis、TGF-beta signaling 和 blood-vessel remodeling 分别分析，不合并为虚构DEG。Activation原始仅11基因，按冻结>=15背景匹配基因规则NOT_EVALUABLE。全部来源、SHA256及准确基因列表见02_receiver_geneset。
```
                         program  original_genes  matched_genes                  status
    GOBP_BLOOD_VESSEL_REMODELING              54             21               EVALUABLE
GOBP_ENDOTHELIAL_CELL_ACTIVATION              11              5 NOT_EVALUABLE_TOO_SMALL
```
Endothelial背景7070基因；C5表达且受体有患者重复支持的246个broad ligands。 pooled detection>=10%，C5>=5细胞、Endo>=20细胞的至少2患者支持。三套程序的前10排名如下；分数是对先验网络/指定程序的优先排序，不是实际活性或显著性：
```
                     program test_ligand  rank  aupr_corrected  CellChat_supported
       HALLMARK_ANGIOGENESIS         B2M     1         0.06993               False
       HALLMARK_ANGIOGENESIS        SDC2     2         0.06555                True
       HALLMARK_ANGIOGENESIS       SPON1     3         0.06334               False
       HALLMARK_ANGIOGENESIS       VEGFA     4         0.05979                True
       HALLMARK_ANGIOGENESIS        FGF1     5         0.05718                True
       HALLMARK_ANGIOGENESIS       TGFB2     6         0.03863                True
       HALLMARK_ANGIOGENESIS      COL5A1     7         0.03853               False
       HALLMARK_ANGIOGENESIS        FGF2     8         0.03852                True
       HALLMARK_ANGIOGENESIS      COL4A1     9          0.0345                True
       HALLMARK_ANGIOGENESIS      COL5A2    10         0.03424               False
 HALLMARK_TGF_BETA_SIGNALING       TGFB1     1          0.1532                True
 HALLMARK_TGF_BETA_SIGNALING       BMP8A     2          0.1343                True
 HALLMARK_TGF_BETA_SIGNALING         FST     3          0.1341               False
 HALLMARK_TGF_BETA_SIGNALING        BMP1     4          0.1135               False
 HALLMARK_TGF_BETA_SIGNALING        BMP2     5           0.108               False
 HALLMARK_TGF_BETA_SIGNALING        BMP4     6          0.1068               False
 HALLMARK_TGF_BETA_SIGNALING        BMP5     7          0.1013                True
 HALLMARK_TGF_BETA_SIGNALING       INHBA     8         0.09029                True
 HALLMARK_TGF_BETA_SIGNALING        PSAP     9          0.0883               False
 HALLMARK_TGF_BETA_SIGNALING         ENG    10         0.08166               False
GOBP_BLOOD_VESSEL_REMODELING        CTSD     1         0.06585               False
GOBP_BLOOD_VESSEL_REMODELING       TIMP1     2         0.06275               False
GOBP_BLOOD_VESSEL_REMODELING       MFGE8     3         0.06056               False
GOBP_BLOOD_VESSEL_REMODELING    SERPINE1     4         0.05813               False
GOBP_BLOOD_VESSEL_REMODELING        WNT4     5         0.05658                True
GOBP_BLOOD_VESSEL_REMODELING       ALCAM     6         0.05622               False
GOBP_BLOOD_VESSEL_REMODELING       SPARC     7         0.04376               False
GOBP_BLOOD_VESSEL_REMODELING        JAG1     8         0.04281                True
GOBP_BLOOD_VESSEL_REMODELING        FGF2     9         0.03631                True
GOBP_BLOOD_VESSEL_REMODELING       PTHLH    10         0.03457               False
```

## 3. 指定配体的实际排名
```
                     program test_ligand  rank  aupr_corrected  CellChat_supported
       HALLMARK_ANGIOGENESIS         FN1    12         0.03231                True
       HALLMARK_ANGIOGENESIS        WNT2    15         0.02888                True
       HALLMARK_ANGIOGENESIS      COL1A1    16         0.02776                True
       HALLMARK_ANGIOGENESIS         MDK    41          0.0214                True
       HALLMARK_ANGIOGENESIS      SEMA3C   115         0.01761                True
       HALLMARK_ANGIOGENESIS      COL1A2   168         0.01517                True
       HALLMARK_ANGIOGENESIS         ADM   222         0.01078                True
 HALLMARK_TGF_BETA_SIGNALING      COL1A2    28         0.04852                True
 HALLMARK_TGF_BETA_SIGNALING      SEMA3C    58          0.0434                True
 HALLMARK_TGF_BETA_SIGNALING         FN1   108         0.03857                True
 HALLMARK_TGF_BETA_SIGNALING         MDK   125          0.0371                True
 HALLMARK_TGF_BETA_SIGNALING      COL1A1   170         0.03279                True
 HALLMARK_TGF_BETA_SIGNALING        WNT2   218         0.01912                True
 HALLMARK_TGF_BETA_SIGNALING         ADM   230         0.01559                True
GOBP_BLOOD_VESSEL_REMODELING      COL1A2    64         0.01581                True
GOBP_BLOOD_VESSEL_REMODELING         FN1    72         0.01539                True
GOBP_BLOOD_VESSEL_REMODELING         MDK    84         0.01472                True
GOBP_BLOOD_VESSEL_REMODELING      COL1A1   157         0.01238                True
GOBP_BLOOD_VESSEL_REMODELING         ADM   203        0.008932                True
GOBP_BLOOD_VESSEL_REMODELING      SEMA3C   208        0.008148                True
GOBP_BLOOD_VESSEL_REMODELING        WNT2   226        0.006212                True
```
未出现在表中的配体应查 broad_LR_all.tsv 的表达/数据库门槛，不能解释为零活性。所有NicheNet-only、CellChat-only、共同支持和受体条目保留于完整整合表。CellChat使用tumor-only t20，t10/30患者分母完整保留；配体排名未按空间显著性调整。

## 4. 表达特异性与预测 targets
共同患者中 C5 vs C0/C1/C4/其余Fibro比较使用相同标准化RNA表达的患者均值差，完整连续效应与方向支持保存；不以逐细胞P值代替患者证据。
```
  gene comparison  n_patients  median_difference  positive_fraction
COL1A1   Fibro_C0          13              1.199                  1
COL1A2   Fibro_C0          13             0.8665                  1
COL3A1   Fibro_C0          13              1.103                  1
   FN1   Fibro_C0          13              1.034             0.9231
 HSPG2   Fibro_C0          13            0.09495             0.6923
 LAMA1   Fibro_C0          13                  0             0.4615
 TGFB1   Fibro_C0          13           -0.03583             0.4615
  WNT4   Fibro_C0          13           0.007573             0.5385
COL1A1   Fibro_C1          20              1.522               0.85
COL1A2   Fibro_C1          20             0.9893               0.95
COL3A1   Fibro_C1          20              1.052               0.95
   FN1   Fibro_C1          20             0.5063               0.65
 HSPG2   Fibro_C1          20             0.5758               0.95
 LAMA1   Fibro_C1          20            0.07651                0.6
 TGFB1   Fibro_C1          20            0.02303               0.55
  WNT4   Fibro_C1          20            0.01979               0.55
COL1A1   Fibro_C4          14              1.024             0.9286
COL1A2   Fibro_C4          14              0.878                  1
COL3A1   Fibro_C4          14              1.011             0.7857
   FN1   Fibro_C4          14              1.235                  1
 HSPG2   Fibro_C4          14             0.3228             0.9286
 LAMA1   Fibro_C4          14           -0.02622             0.3571
 TGFB1   Fibro_C4          14          -0.006303                0.5
  WNT4   Fibro_C4          14           -0.08644             0.2143
COL1A1       REST          20              1.258                0.9
COL1A2       REST          20             0.9168               0.95
COL3A1       REST          20             0.9815               0.95
   FN1       REST          20             0.6227                0.7
 HSPG2       REST          20             0.3996               0.95
 LAMA1       REST          20             0.0388                0.6
 TGFB1       REST          20            0.04176               0.55
  WNT4       REST          20            0.00779                0.5
```
均值差不等于C5-specific；大量CAF共享的配体只能称shared CAF ligands。前10配体target表使用官方get_weighted_ligand_target_links(n=250)，展示先验top250 targets与指定程序交集；空白不是零调控潜力。target属于标准内皮相关功能程序，但不是在本队列观察到的C5诱导target。

## 5. 空间候选、覆盖与表达
候选在空间分析前按各程序标准化排名均值和CellChat支持确定，再取每配体CellChat中位强度最高的受体轴；FN1/COL1A1/COL1A2为固定广泛ECM对照。该聚合只用于候选优先级，单独程序排名仍完整保留。
```
ligand    receptor   required_subunits        role
COL1A1 ITGA1_ITGB1         ITGA1;ITGB1 ECM_control
COL1A2 ITGA1_ITGB1         ITGA1;ITGB1 ECM_control
COL3A1      ADGRG1              ADGRG1   candidate
   FN1 ITGA5_ITGB1         ITGA5;ITGB1 ECM_control
 HSPG2        DAG1                DAG1   candidate
 LAMA1 ITGA6_ITGB1         ITGA6;ITGB1   candidate
 TGFB1 ACVR1_TGFbR ACVR1;TGFBR1;TGFBR2   candidate
  WNT4   FZD4_LRP5           FZD4;LRP5   candidate
```
受体复合体从CellChatDB展开必要亚基；复合体连续表达取最小亚基log1p(CP10k)，检测要求所有必要亚基检出。CALCRL单基因定义不能证明完整ADM受体，RAMP2/3另外审查。不同资源定义不能互相替代。
```
  gene  sections_present  median_detection
 ACVR1                63            0.2798
ADGRG1                63            0.4654
COL1A1                63              0.96
COL1A2                63            0.9295
COL3A1                63            0.9609
  DAG1                63            0.5029
   FN1                63             0.845
  FZD4                63            0.1879
 HSPG2                63            0.4448
 ITGA1                63            0.2819
 ITGA5                63            0.4457
 ITGA6                63            0.3871
 ITGB1                63            0.6155
 LAMA1                63           0.06635
  LRP5                63            0.4679
 TGFB1                63            0.2715
TGFBR1                63             0.423
TGFBR2                63            0.5395
  WNT4                63            0.1159
```
所有63张按切片读取counts，barcode和UMI总数与旧评分精确对齐。缺feature标记NA/不可评估，不填成零。

## 6. 空间邻接是否支持、组成是否影响？
发送区域Five切片内top10%且Fib上四分位；接收为Endothelial上四分位，q05支持比例保留。原始阵列六邻接，Ring1主要、Ring0/2对照，无微米换算。连续配体与邻近Endo-rich受体表达相关，以及高分发送区域与对照的表达乘积差，均是混合spot的探索性关联。
```
ligand  n_sections  median_raw_rho  median_adjusted_rho  median_matched_effect  balanced_sections  technically_evaluable_sections  median_sender_detection  median_receiver_joint_detection
COL1A1          63          0.1276              0.01913                 0.3799                  2                              59                        1                           0.4136
COL1A2          63           0.138              0.02302                 0.3383                  2                              59                        1                           0.4136
COL3A1          63         -0.1538            -0.002795                 0.5151                  2                              59                        1                           0.8225
   FN1          63          0.1873              0.06745                 0.3929                  2                              59                        1                           0.5557
 HSPG2          63        0.008281              0.01503                 0.0745                  2                              55                   0.7107                           0.7831
 LAMA1          63         0.02479              -0.0148               0.001353                  2                              35                  0.09924                           0.5641
 TGFB1          63         0.06921              0.02938               0.004064                  2                              51                   0.4648                           0.2219
  WNT4          63         0.04644              0.01587                      0                  2                              36                   0.0981                           0.2161
```
rank-residual控制Fib、Endo、Epi、Mural、nCount、nFeature和图度数。匹配同切片Fib-rich低Five区，七维标准化秩距离<=1、无放回，平衡要求所有|SMD|<0.2且>=10对。完整匹配和不平衡结果都保留。另有配体检出一致且log表达差<0.5的匹配敏感性。
空间功能程序使用排除Five和候选ligand后的实际覆盖基因平均log1p(CP10k)，至少5基因且覆盖>=50%；不等同于内皮细胞内表达。
非环绕阵列平移最多99次，至少50%重叠，保持平移区域的局部结构但受边界和不规则组织限制。其P值仅为敏感性。BH在全部有平移检验的候选/对照×切片集合内校正。患者推断不使用spot数作为n。

## 7. 独立患者、数据集与论文可用性
患者级主汇总先应用检测充分性规则：发送区域至少10个ligand检出spot、内皮富集区至少10个必要受体亚基共同检出spot；匹配效应至少5对。规则在技术执行期间、跨切片聚合前固定；所有未通过切片仍完整保存描述结果，不替换候选。已确认患者先平均其合格多张切片效应；未知患者仅切片描述。患者内聚合后分dataset bootstrap 2000次，seed123；所有数据集/leave-one-dataset-out见08_patient_dataset_sensitivity。置信区间是探索性的未同时校正区间。各metric与LODO范围分别在全部轴、阈值、ring和程序检验内BH调整（较保守）。
```
ligand                 metric  n_patients  n_datasets  median_effect    ci_low    ci_high    BH_FDR
COL1A1           adjusted_rho          29           3       0.001336  -0.01839    0.06999    0.2821
COL1A1 matched_product_effect          29           3         0.3842    0.3307     0.8839 1.788e-07
COL1A1                raw_rho          29           3        0.06792 -0.009297     0.1686   0.01761
COL1A2           adjusted_rho          29           3        0.01402 -0.001241    0.06343   0.06623
COL1A2 matched_product_effect          29           3         0.3604    0.2747     0.7579 1.788e-07
COL1A2                raw_rho          29           3         0.1041   0.02733     0.1723  0.002168
COL3A1           adjusted_rho          29           3       -0.01894  -0.07899    0.01328    0.3264
COL3A1 matched_product_effect          29           3         0.5838    0.2562     0.9204 5.364e-07
COL3A1                raw_rho          29           3        -0.2015   -0.2817     -0.091 2.411e-05
   FN1           adjusted_rho          29           3        0.07418   0.03532     0.1802 5.692e-05
   FN1 matched_product_effect          29           3         0.3579    0.2147     0.8135 0.0006802
   FN1                raw_rho          29           3          0.182    0.1076     0.3423 2.891e-05
 HSPG2           adjusted_rho          28           3        0.01533 -0.004001    0.03441      0.25
 HSPG2 matched_product_effect          28           3         0.0694  0.002475     0.1314   0.05042
 HSPG2                raw_rho          28           3        0.02644  0.002418    0.06936    0.1918
 LAMA1           adjusted_rho          18           2       -0.06548  -0.09217   -0.01853    0.1179
 LAMA1 matched_product_effect          18           2        0.01189 -0.006598    0.04129    0.2354
 LAMA1                raw_rho          18           2        0.05229   0.02699    0.08411 0.0008392
 TGFB1           adjusted_rho          26           3        0.04556   0.01834    0.07781   0.00551
 TGFB1 matched_product_effect          25           3       0.006952 -0.001691    0.02262   0.05042
 TGFB1                raw_rho          26           3        0.08236   0.06757     0.1305 1.431e-06
  WNT4           adjusted_rho          16           2        0.05278  -0.03689    0.09211     0.175
  WNT4 matched_product_effect          16           2      -0.001587 -0.006374 -3.655e-05   0.05616
  WNT4                raw_rho          16           2        0.06062    0.0226     0.1024 0.0004962
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


## 补充：功能程序、来源稳定性与解释限制

在排除Five及入选ligand后，TGFB1与邻近内皮富集spots的标准TGF-beta程序仍有弱条件关联：患者中位rho约0.067 (95%CI 0.022–0.094)。这不是内皮细胞内SMAD响应的测量，也不是方案A失败后的替代实测DEG。
追加标准GO ECM-organization/collagen-fibril-organization去除敏感性后，TGFB1–TGF-beta程序rho约0.067 (CI 0.025–0.094)，FN1–Angiogenesis约0.228 (CI 0.161–0.304)。但本去除规则实际仅从已去重叠程序中再去除2个Angiogenesis、1个TGF-beta、6个vascular-remodeling基因，不能称已彻底排除ECM状态。完整定义、每切片结果与患者LODO见ECM_depleted*文件。这是预先要求的混杂诊断扩展，不是独立验证，也未反向改变候选排序。

WNT4在移除InhouseData后不再通过原表达候选门槛，进一步提示来源依赖；其余所列主候选在两主要来源的逐一留出表达门槛下保留。表达门槛保留不等于NicheNet响应的跨队列验证，因为本轮没有患者关联response。

单细胞患者与Visium患者没有进行未经证实的一一配对。只有三个空间来源有确认患者身份，其余来源的LODO不会提供新的患者层面独立证据。全部空间组织标签原样保留，UNRESOLVED没有改为Tumor。

### 证据分级
- 可靠表达与计算：固定资源/程序的全部配体排名、表达检测、63切片覆盖、受体亚基审查。
- 有限患者支持的探索关联：FN1共享ECM对照最稳定；TGFB1弱且依赖GSE319536。原位细胞身份和空间零模型支持尚不足。
- 先验网络和功能程序候选：HSPG2、TGFB1等配体–target关联；不能称C5诱导。
- 空间支持不足：COL3A1–ADGRG1有网络边缺口且邻接非稳定正向；LAMA1、WNT4检出/方向/来源不稳。ADM/MDK/SEMA3C/WNT2保留真实排序与覆盖，并未因不是主候选而记作阴性。

可用于小文章的是审慎的表达与描述性空间证据组合。若后续只推进一个机制，建议先补充分辨C5与其他CAF的原位证据和独立内皮响应数据，而非直接将任一当前轴写为已确立的C5特异通讯。
