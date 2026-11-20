# C5–Endothelial GSEA 与 nominal-DEG NicheNet v1

完成日期：2026-11-20。全部分析使用现有冻结数据；随机种子 123，BLAS/OpenMP/R 单线程、CPU affinity 限于 4 核。没有修改历史结果、重新注释或重新运行 CellChat/cell2location。

## 核心结论

**没有获得主模型 FDR 支持的 Endothelial EMT、Angiogenesis 或 TGF-β 响应。** 主模型 EMT NES=1.183，BH-FDR=0.320；采样规模调整后 EMT 得到支持，但不能用敏感性模型替换主模型。标准 EMT-NicheNet 为 **REFERENCE_ONLY**，不等于观察到 EndMT。

名义 Set A/B/C 经双队列正向及固定背景筛选，分别为 **28/123/302 genes**。Set B Top10 为：POSTN, INHBA, ANGPTL4, IGF2, TGFB1, PSAP, TGFB2, BMP2, HMGB1, SFRP1。这些不是 FDR 显著 DEG，也不是已验证的 C5-induced targets。完整 246 配体均输出，不只保留与 CellChat 的交集。

最值得保留的是“真实内皮转录变化下的候选排序及其不稳定性”，而非确定的 C5→内皮机制。TGFB1 在多个阈值/标准程序中排名较高，但其来源并非 C5 特异，既有空间支持仍受数据集及组成限制。

## 1. 输入、排名、资源和可重复性

- 21 名符合原标准的独立肿瘤患者：InhouseData 9、GSE222315 8、CNP0000460 3、GSE135337 1；Fibro C0–C6≥50、Endothelial≥20。同患者多个肿瘤 specimen 已聚合。
- 主模型 `~ dataset + standardized C5 fraction`，edgeR TMM → limma-voom → robust eBayes。主要 DEG 复用，不从大型单细胞对象重建。LOPO 用原小型 pseudobulk 重拟合主模型与两个主要队列模型。
- 重拟合原模型与原 t 值逐基因一致，见 `00_input_audit/model_reproduction.tsv`。
- GSEA 使用完整有限 moderated t 排名，不按 DEG 显著性截断。各模型检测背景不同，未补入缺失基因。混合 ENSG 标识原样保留在完整排名中，只对符号匹配的 Hallmark 基因计命中；未猜测别名或将其移至末尾。重复符号为 0、非有限 t 为 0；详见审查表。
- MSigDB **2026.1.Hs** Hallmark 50，与 C5 既有 GSEA 同版；fgseaMultilevel，minSize=10、maxSize=500、eps=0、nPermSimple=10000，SerialParam(RNGseed=123)。每模型全部 50 通路 BH 校正。少量 t ties 按 gene symbol 确定输入顺序，日志记录原 fgsea 警告；不存在结果导向重排。
- C5 旧 GSEA 使用 avg_log2FC、不同 detection 过滤，内皮使用患者级 moderated t，故不能将两者 NES 数值当作直接可比效应量。
- 配体/受体表达与背景沿用上一轮的全部合格肿瘤表达审查，和本轮 21 位 pseudobulk 患者的统计样本不是完全相同的集合；因此表达支持人数不能冒充本轮响应的患者分母。随机匹配也沿用原背景 pooled detection，而不是重新选择背景。
- NicheNet 固定 human v2、nichenetr 2.2.1.1、7,070 背景与 246 broad ligands。官方 `predict_ligand_activities()`；以 corrected AUPR、Pearson、gene symbol 排序。资源 SHA256 已复核，未下载或重新筛选 sender/background。
- Set A/B/C：主模型 logFC>0，P<0.01/0.05/0.10，InhouseData 与 GSE222315 logFC 均>0，且属于固定背景及 LT target universe。未在某队列通过表达过滤的基因标为未检验，不能算阴性或强行满足方向条件。
- 每套集合的精确列表见 `02_nominal_response_genesets/frozen_genesets.tsv`；原始阈值、排除清单和流程见 `00_input_audit/prespecified_parameters.json`。

## 2. Endothelial Hallmark GSEA

表中为 NES (BH-FDR)，所有 250 项检验见 `01_Endothelial_Hallmark_GSEA/all_models.tsv`。

| pathway | GSE222315 | InhouseData | primary | sampling_adjusted | stage_adjusted |
| --- | --- | --- | --- | --- | --- |
| HALLMARK_ANGIOGENESIS | +1.501 (0.0941) | +1.021 (0.594) | +1.168 (0.453) | +1.606 (0.0963) | +1.162 (0.447) |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | +0.947 (0.684) | +1.497 (0.0417) | +1.183 (0.32) | +1.423 (0.0341) | +1.177 (0.345) |
| HALLMARK_GLYCOLYSIS | +0.754 (0.942) | +1.165 (0.362) | -1.203 (0.297) | -1.103 (0.533) | -1.200 (0.312) |
| HALLMARK_HYPOXIA | -1.237 (0.126) | +1.443 (0.0574) | -0.859 (0.913) | +1.056 (0.548) | +0.836 (0.941) |
| HALLMARK_TGF_BETA_SIGNALING | +1.648 (0.0301) | +0.994 (0.609) | +1.220 (0.342) | +1.190 (0.346) | +1.178 (0.407) |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | -1.103 (0.312) | +1.243 (0.238) | -0.781 (0.999) | +1.027 (0.637) | +0.818 (0.941) |

主模型 FDR<0.10 有 12 个通路，其中 11 个与两个主要队列同向。以下是“主模型 FDR 支持 + 队列方向一致”，不是所有队列各自显著，也不是独立验证：

| pathway | NES | padj |
| --- | --- | --- |
| HALLMARK_ADIPOGENESIS | -1.911 | 1.593e-05 |
| HALLMARK_APICAL_SURFACE | 1.842 | 0.01314 |
| HALLMARK_DNA_REPAIR | -2.46 | 2.634e-10 |
| HALLMARK_FATTY_ACID_METABOLISM | -1.437 | 0.06689 |
| HALLMARK_MITOTIC_SPINDLE | 1.641 | 0.008344 |
| HALLMARK_MYC_TARGETS_V1 | -3.019 | 2.107e-26 |
| HALLMARK_MYC_TARGETS_V2 | -2.237 | 4.486e-05 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | -3.109 | 1.693e-28 |
| HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | -1.819 | 0.008738 |
| HALLMARK_UV_RESPONSE_UP | -1.559 | 0.03637 |
| HALLMARK_XENOBIOTIC_METABOLISM | -1.55 | 0.02586 |

E2F targets 主模型负向，而 GSE222315 正向，明确属于来源不一致。EMT 在 InhouseData BH=0.0417，但 GSE222315 BH=0.684；采样调整模型 BH=0.0341。Angiogenesis 在 GSE222315 BH=0.0941、采样调整 BH=0.0963，主模型不支持；TGF-β 仅 GSE222315 BH=0.0301。不能选择这些更有利的模型替代预定主结果。

C5 自身 EMT 在旧 5%/10%/25% 检测阈值下均为正且 FDR<0.05，是用户重点程序中最稳定的上调特征。Glycolysis 在 25% 阈值 FDR≈0.0545，稳定性较弱。C5 的 EMT leading-edge 未用于任何内皮 geneset_oi。

## 3. 名义基因集与来源稳定性

| set | p_threshold | positive_primary | positive_both_cohorts | matched | ECM |
| --- | --- | --- | --- | --- | --- |
| A | 0.01 | 102 | 40 | 28 | 1 |
| B | 0.05 | 351 | 167 | 123 | 1 |
| C | 0.1 | 748 | 400 | 302 | 1 |

| set | n | min_FDR | median_primary_P | median_logFC | cohort_logFC_spearman | stage_positive_fraction | sampling_positive_fraction | ECM_n | ECM_genes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A | 28 | 0.4208 | 0.006116 | 0.4788 | 0.0208 | 1 | 1 | 1 | FAP |
| B | 123 | 0.4208 | 0.02313 | 0.3761 | -0.1204 | 1 | 1 | 1 | FAP |
| C | 302 | 0.4208 | 0.06163 | 0.3269 | -0.1563 | 1 | 1 | 1 | FAP |

三套固定间质清单均只命中 **FAP**；保留原集，同时另存 noECM 版本 (27/122/301 genes)。这只是预定义标志排除敏感性，不能说明没有广泛 ECM 信号或污染。双队列方向条件利用了参与合并模型的患者，不构成独立复制。

21 次 LOPO 均可拟合，名义 Set B 数量 96–180，原集 Jaccard 中位数 0.775 (范围 0.335–0.913)。所有 21 次均达到预定至少 15 基因门槛；没有补齐基因。原 Set B 基因在留出中的正向比例、入选比例、最大 P 值和效应范围已逐基因保存；28/123 个在全部留出中仍入选。LOPO 中队列内方向也重新估计，未使用被留出患者。POSTN 最差留出为 InhouseData::T4295396，名次由 1 降至 114；这是单患者影响的具体例子，不能只报告中位名次仍为 1。

## 4. Set B NicheNet 与阈值/患者敏感性

| test_ligand | rank | aupr_corrected | auroc | pearson | CellChat_supported | LOPO_rank_min | LOPO_rank_max | LOPO_top10_fraction | empirical_exceedance | empirical_BH_across246 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| POSTN | 1 | 0.008527 | 0.5189 | 0.004145 | True | 1 | 114 | 0.9524 | 0.01961 | 0.5063 |
| INHBA | 2 | 0.008312 | 0.5198 | -0.001263 | True | 1 | 6 | 1 | 0.01961 | 0.5063 |
| ANGPTL4 | 3 | 0.005714 | 0.5092 | 0.05914 | True | 3 | 29 | 0.9048 | 0.01961 | 0.5063 |
| IGF2 | 4 | 0.00469 | 0.5401 | 0.03995 | True | 3 | 13 | 0.8571 | 0.03922 | 0.5063 |
| TGFB1 | 5 | 0.004381 | 0.5533 | 0.02353 | True | 1 | 12 | 0.9048 | 0.01961 | 0.5063 |
| PSAP | 6 | 0.003752 | 0.5402 | 0.04111 | False | 1 | 39 | 0.8571 | 0.01961 | 0.5063 |
| TGFB2 | 7 | 0.002444 | 0.523 | 0.007357 | True | 2 | 34 | 0.6667 | 0.05882 | 0.5063 |
| BMP2 | 8 | 0.002363 | 0.5207 | 0.01003 | False | 3 | 67 | 0.4762 | 0.07843 | 0.5063 |
| HMGB1 | 9 | 0.002317 | 0.5483 | 0.005815 | False | 6 | 76 | 0.619 | 0.07843 | 0.5063 |
| SFRP1 | 10 | 0.001987 | 0.5489 | 0.01028 | False | 4 | 58 | 0.381 | 0.09804 | 0.5063 |

ANGPTL4 在 Set B 的官方 top250-target 交集为空，未通过改阈值补造链接；activity 可计算与可展示 top250 link 是不同结果。

Top 配体 corrected AUPR 较小，且部分 AUROC 接近 0.5、Pearson 接近 0 或负值；不能只凭“第 1 名”写成配体实际活性增强。gene set 大小不同，跨 A/B/C/EMT 比较仅讨论排序，不将 AUPR 差值当作方案优劣。

| program1 | program2 | rank_spearman | top10_intersection | top20_intersection |
| --- | --- | --- | --- | --- |
| Set_A | Set_B | 0.5821 | 4 | 8 |
| Set_A | Set_C | 0.4914 | 2 | 7 |
| Set_A | Set_B_noECM | 0.5641 | 4 | 8 |
| Set_A | EMT | 0.09128 | 2 | 3 |
| Set_B | Set_C | 0.6962 | 7 | 12 |
| Set_B | Set_B_noECM | 0.9953 | 10 | 19 |
| Set_B | EMT | -0.03023 | 2 | 2 |
| Set_C | Set_B_noECM | 0.6858 | 7 | 13 |
| Set_C | EMT | 0.009413 | 1 | 2 |

患者留出排名相关中位数 0.865 (范围 0.449–0.967)，Top10 与完整 Set B 交集为 3–10。这是影响诊断，不是独立测试性能或患者级置信区间。逐配体中位/最小/最大名次见 `LOPO_ligand_stability.tsv`。

重点配体真实名次：

| test_ligand | EMT | Set_A | Set_B | Set_C |
| --- | --- | --- | --- | --- |
| ADM | 221 | 219 | 197 | 8 |
| COL1A1 | 111 | 65 | 149 | 91 |
| COL1A2 | 90 | 62 | 95 | 122 |
| FN1 | 21 | 108 | 116 | 180 |
| MDK | 170 | 85 | 222 | 188 |
| SEMA3C | 166 | 218 | 89 | 130 |
| TGFB1 | 1 | 9 | 5 | 4 |
| VEGFA | 217 | 19 | 42 | 66 |
| WNT2 | 230 | 66 | 12 | 18 |

尤其 ADM 从 A/B 的低排名升至 Set C 第 8，不能据此挑选宽松阈值作为唯一主分析。TGFB1 的 A/B/C 名次为 9/5/4；WNT2 为 66/12/18。MDK、SEMA3C、FN1/COLLAGEN 未一致进入名义 Top10。

50 次随机对照严格匹配 Set B 数量和内皮 pooled detection 十分位分布，每次 123 genes，无放回抽样，允许与原 Set B 重叠。`empirical_exceedance=(1+null>=observed)/(51)`；另给全部 246 配体 BH，最低经验 BH 约 0.506。分辨率最低 1/51，这些是网络/程序匹配经验校准，不是患者级显著性，更没有模拟完整名义基因筛选过程。0/246 配体通过此经验 BH<0.05。

预定义 FAP 删除后前三名不变。另按冻结 Set B Top10 各配体的最大先验权重 target，逐一删除唯一 target 重新评分；不是选择使结果更好的删除方案：

| removed_gene | rank_spearman | top10_overlap | max_top10_rank_change |
| --- | --- | --- | --- |
| AOPEP | 0.9907 | 9 | 80 |
| ATXN1 | 0.9968 | 9 | 1 |
| FAP | 0.9953 | 10 | 1 |
| NPEPPS | 0.9942 | 9 | 13 |
| PARP12 | 0.9964 | 9 | 3 |
| SLC16A14 | 0.9952 | 9 | 2 |
| SMAD7 | 0.9791 | 7 | 47 |

影响最大的名次变化 (完整表保留全部配体)：

| test_ligand | removed_gene | rank_original | rank_drop | rank_change |
| --- | --- | --- | --- | --- |
| POSTN | AOPEP | 1 | 81 | 80 |
| INHBA | SMAD7 | 2 | 49 | 47 |
| BMP2 | SMAD7 | 8 | 26 | 18 |
| TGFB2 | SMAD7 | 7 | 25 | 18 |
| PSAP | NPEPPS | 6 | 19 | 13 |
| HMGB1 | SMAD7 | 9 | 6 | -3 |
| HMGB1 | PARP12 | 9 | 12 | 3 |
| SFRP1 | SMAD7 | 10 | 7 | -3 |
| PSAP | SMAD7 | 6 | 4 | -2 |
| TGFB2 | SLC16A14 | 7 | 9 | 2 |
| ANGPTL4 | AOPEP | 3 | 2 | -1 |
| INHBA | AOPEP | 2 | 1 | -1 |

这说明排名是否依赖个别 target，而不能把 LT 权重份额等同评分的精确贡献或真实靶基因激活。

## 5. 标准 EMT 程序及污染/EndMT 风险

EMT 原始 200 genes，固定背景匹配 98；状态 **REFERENCE_ONLY**。官方 NicheNet 可计算，但没有通过患者响应门槛。Top10：TGFB1, TGFB2, TNC, B2M, ENG, VCAN, ICAM1, LTBP1, BMP1, CTSD。

EMT 与 nominal Set B 的 Top10 交集为 2 个；共同候选是 TGFB1, TGFB2。既有 Angiogenesis/TGF-β 使用 2025.1，本轮 EMT 用 2026.1，比较图显式标注版本差异，仅比较优先级，不当作同版数值验证。

冻结内皮标注中，合格患者共有 7167 个内皮细胞；至少检出 1/2 个 PECAM1/VWF/CDH5/KDR 的比例分别为 96.1%/88.0%。至少检出 2 个 FAP/COL1A1/COL1A2/POSTN/DCN/LUM 的比例为 8.5%，同时≥2 内皮与≥2间质标志为 7.5%；≥2间质但无上述内皮标志为 0.32%。这些细胞层比例仅为身份诊断，不作独立患者推断。共检出使用冻结 LogNormalize10000 稀疏矩阵中的非零值，与原始计数非零对应，不是新增打分。

原冻结流程已经进行 QC/doublet 过滤；本次未新聚类、未改变标签。共同检出可来自状态变化、残余 doublet 或环境 RNA，无法据此确认或排除污染；检测不足也会漏掉真正内皮标志。患者级 CPM、检测率与 specimen/患者级共同检出表均保存。适当名称是 Endothelial EMT-associated transcriptional signature，**不能宣称 EndMT 或 C5 诱导 EndMT**。

## 6. CellChat、C5 特异性与已有空间证据

Set B Top10 中有 t20 阳性 CellChat 的配体及每配体强度最高的既有轴如下；完整多轴、t10/t20/t30 和 eligible 分母均保留：

| test_ligand | rank | receptor | positive_patients | eligible_patients | patient_condition_median | C5_pooled_pct | Fibro_C0_median_difference | Fibro_C4_median_difference | required_subunits | Endo_min_subunit_pct |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| POSTN | 1 | ITGAV_ITGB5 | 6 | 12 | 0.001004 | 0.7638 | 0.9604 | 0.7624 | ITGAV;ITGB5 | 0.1396 |
| INHBA | 2 | ACVR1B_ACVR2A | 3 | 12 | 0 | 0.6786 | 0.5612 | 0.5257 | ACVR1B;ACVR2A | 0.06444 |
| ANGPTL4 | 3 | CDH5 | 5 | 12 | 0 | 0.1555 | 0 | -0.02185 | CDH5 | 0.6618 |
| IGF2 | 4 | IGF1R | 1 | 12 | 0 | 0.1988 | 0 | 0.003989 | IGF1R | 0.1814 |
| TGFB1 | 5 | ACVR1_TGFbR | 2 | 12 | 0 | 0.4664 | -0.03583 | -0.006303 | ACVR1;TGFBR1;TGFBR2 | 0.1121 |
| TGFB2 | 7 | ACVR1B_TGFbR2 | 1 | 12 | 0 | 0.1435 | -0.09022 | 5.442e-05 | ACVR1B;TGFBR2 | 0.06444 |

POSTN 相对 C0/C4 的差值中位数约 +0.960/+0.762，患者正向比例 12/13、13/14；INHBA 约 +0.561/+0.526，正向比例 13/13、11/14。这两者相对富集较清楚，但分别存在患者/target 依赖及 Five-gene 共享基因循环性。HMGB1 也有相对表达富集，但没有本次 t20 C5→Endothelial CellChat 支持；ANGPTL4、TGFB1 则不能称为 C5 富集配体。

表达差异是共同可评估患者中 mean-log-expression 的差值中位数，不是患者级 DE FDR。只有相对 C0/C4 更高及患者方向支持的配体可称相对富集，不能推广为 C5-specific。TGFB1 并非 C5 特异；NicheNet 的 ligand-target 先验不验证某一 CellChat receptor complex。

原空间分析复用：FN1–ITGA5/ITGB1 为共享 ECM 对照，调整 Ring1 rho≈0.074 (29 位已确认患者)，TGFB1–ACVR1/TGFBR1/TGFBR2 rho≈0.046 (26 位患者)，但 TGFB1 leave-GSE319536-out CI 跨 0。原匹配平衡子集仅 1 位已确认患者，504 项平移检验无 BH<0.05；因此不能用原正相关宣称 C5 特异通讯。其他原轴结果不因本轮排序而重算或隐藏。

按冻结规则新增且仅新增 POSTN–ITGAV/ITGB5、INHBA–ACVR1B/ACVR2A、ANGPTL4–CDH5 的轻量表达审查 (63 切片)。以下为逐切片中位数，无新的空间显著性推断：

| ligand | n_sections | median_sender_detection | median_joint_receptor_detection | median_edge_detection |
| --- | --- | --- | --- | --- |
| ANGPTL4 | 63 | 0.1503 | 0.3229 | 0.0715 |
| INHBA | 63 | 0.8288 | 0.2643 | 0.1058 |
| POSTN | 63 | 0.7931 | 0.5685 | 0.4788 |

发送区仍为原 Five-high (>切片 90% 分位、排除零分) 且 Fibro-rich (≥75% 分位)，接收区为 Endothelial-rich (≥75% 分位)，仅实际六邻接边。受体复杂体以全部必要亚基共同检出审查，不用单个亚基代替；表达用 raw counts→log1p(CP10k)。有边的 joint detection 只是表达相容性，没有匹配或空间零模型，**不构成新轴的邻接富集证据**。INHBA 本身是 Five-gene 成员，其发送区高表达部分是定义造成的；不能视为独立空间验证。INHBA 的必要亚基单细胞最小检测比例低于 10%，完整复合体支持需谨慎。

本轮不重新运行 63 切片的完整空间统计。已确认患者按照原身份审查的两个合法状态聚合 (30 切片/29 患者)，其余切片仅描述；未将带有 UNVERIFIED 的标签误作 verified。新增轻量审查的患者描述有 29 位，不能与单细胞患者擅自配对。

## 7. 小文章中的证据分层与十二问回答

1. **C5 自身最稳健的重点上调程序是 EMT**，旧多检测阈值均支持；内皮主模型没有复制该显著性。
2. 内皮 EMT/Angiogenesis/TGF-β 主模型均未获 BH<0.10 支持。队列/采样调整模型出现支持，应完整呈现而不替换主分析。
3. 两主要队列对上述三程序均正向，但显著性不一致；Hypoxia、TNFA 等方向不一致。全 50 通路比较已保存。
4. A/B/C 为 28/123/302，方向筛选前后的全部数量见上表；都是名义候选。
5. Set B Top10 见第 4 节；不强制任何旧候选进入。
6. 阈值与 LOPO 稳定性为连续诊断，不能概括为全部稳定；具体相关、交集及逐配体范围已给出。
7. 重点候选完整名次见上表。ADM 的宽松阈值提升及 MDK/SEMA3C 的低/变化排名必须保留。
8. Top10 中 CellChat 支持见第 6 节；预测强度与 positive/eligible 分母不能混为普遍发生率。
9. C5 相对 C0/C1/C4/其余 Fibro 的表达差值和方向比例均在整合表，shared CAF ligand 不能称 C5-specific。
10. EMT/Set B 只有部分共同候选，主要体现 TGF 家族先验一致性；EMT 是 reference-only。
11. 间质共同检出和模型敏感性构成 EndMT 解释风险；没有充分证据定义 EndMT，也不能仅凭共同检出判定污染。
12. 正文可以展示透明的患者级 GSEA、固定 Set B 排名及关键敏感性，结论限于探索性候选；完整 A/C、LOPO、随机对照、参考 EMT 和新增空间表达审查应作为补充，避免用这些弱证据形成确定机制主线。

### FDR 支持且跨队列方向一致的观察

以第 2 节列出的主模型通路为依据，尤其 OXPHOS/MYC/DNA repair 负向及 apical surface/mitotic spindle 正向。它们是观察性患者关联，未证明由 C5 因果驱动。不同通路的表达相关结构仍会影响 GSEA，通路 FDR 不等于 DEG FDR。

### 名义 DEG 支持的探索性结果

固定 Set B 的候选排序、A/C 敏感性、21 次 LOPO、target 删除及随机程序校准；无 FDR 显著主模型 DEG。TGFB1/POSTN/INHBA/ANGPTL4 等仅按各自完整稳定性和来源支持描述，不统一宣称阳性。

### 标准功能程序 NicheNet 先验预测

完整 Hallmark EMT 及旧 Angiogenesis/TGF-β 的网络排序。EMT 不通过主门槛，参考程序不能替代观察到的内皮响应。

### CellChat 与空间表达补充证据

保留 t20 主分析、t10/t30、受体亚基、患者分母、共享 CAF 表达；旧空间混杂限制不变。新三轴只证明技术表达可评估程度，不证明空间富集、真实通讯或细胞内靶反应。

### 不支持或无法确认

不支持主模型显著 C5-associated EndMT/血管生成/TGF-β 响应，不支持真实 C5 空间身份、C5-specific ligand 或配体因果效应。正式 C5-specific 空间定位 NO-GO 不变。最有价值后续是独立患者、内皮纯化/多重原位及针对候选的扰动验证，而不是继续放宽 P 阈值。

## 文件与图说明

- `01_Endothelial_Hallmark_GSEA/`：5 模型、250 通路、完整排名、leading edges、方向比较和客观选曲线。
- `02_nominal_response_genesets/`：完整基因审查、固定列表、去间质版本、身份/污染诊断。
- `03_nominal_NicheNet/`：所有 activity、全部 Top10 links、checkpoint；`04_rank_sensitivity/`：21 LOPO、50 random、target influence。
- `05_EMT_program_NicheNet/`：REFERENCE_ONLY、原始/匹配数量、标准程序排名对比。
- `06_CellChat_spatial_integration/`：完整证据、多阈值 CellChat、冻结三轴与表达审查。
- Figures A: all Hallmark/focus，星号 * BH<0.10、** BH<0.05；B: 主模型 FDR<0.10 中最大的两个正 NES 与两个负 NES；C: 名义筛选/队列效应；D: Set B、阈值、LOPO、随机；E: Top10 prior links 的最高 max-weight 35 targets (完整 top10 links 表不作展示截断；链接定义仍为每配体官方全 target top250 与 geneset 交集，空白不等于零调控潜力)；F: 标准程序排名 (版本明确)。
- 最终自动审查见 `00_input_audit/final_validation.tsv`，核对原输入哈希、250项GSEA的BH、固定集合/配体数量、21次LOPO、50次检测率匹配随机集合、三轴63切片、患者身份、PDF矢量和PNG分辨率。
- 所有图带原始绘图数据；PDF 为矢量、PNG 300 dpi、Arial 12pt、白底；有图例时导出主体与独立图例 PDF。LOPO range 不是 CI，随机 95% 区间不是患者 CI。
