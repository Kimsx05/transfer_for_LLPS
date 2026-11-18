# C5-associated 全量空间探索报告

已完成 **63/63 张切片、265,002 spots**。全量结果用于描述性探索；既有正式 C5-specific 定位 **NO-GO 保持不变**。本轮没有重训 cell2location、没有改变 A26、没有用真实空间图反向选基因。

## 交付与数据身份

- [63张缩略图总览](figures/all63_percentile_atlas.png)、[可浏览图册](GALLERY.html)、[逐切片汇总](tables/section_summary.tsv)。每张切片10个连续值panel及辅助HE视图，PDF矢量点图、PNG 300 dpi。HE本身是原始栅格图像。
- 所有spot级原始评分、百分位、检测基因数、nCount/nFeature、30亚型和9大类mean/q05均保存在`01_spatial_scores/*.parquet`，保留原spot_id及barcode，未重复保存表达矩阵。
- 原7张pilot的五基因/三基因两种评分与历史结果逐spot完全一致，最大绝对误差0。技术通过后扩展全部63张，没有模型AUROC门槛。
- 已知患者身份切片30张、29位患者，来自3个dataset；其余切片不虚构患者身份。可靠Normal 3、Tumor 17、UNRESOLVED 43。CanDis/Inhouse/NPJ来源标签沿用UNVERIFIED限制。
- 所有63张HE只有来源/尺寸/边界几何检查，未独立组织学核验。主图不叠加HE；辅助图明确标为provisional，并使用原始上方原点像素方向，不旋转、不镜像、不转换微米。

## 主要观察

A26呈现可复现的空间程序结构，并常见于Fibroblast较高、上皮较低的区域；但微环境关联在进一步控制上皮组成后明显减弱。当前最有支持的是**C5-associated stromal/ECM程序的描述性空间结构**，而不是C5特异细胞定位或促血管/免疫机制。不能把高Moran I或漂亮图像当成新的辨识性验证。

## 九个问题

1. **A26呈现什么分布？** 全部切片连续分数均已输出。切片A26中位数范围0.000–0.280，这是未校准的描述值，不构成跨dataset绝对丰度差异。原始score统一[0,1]显示；百分位panel与总览只比较切片内部位置；并列零分也会得到中间排名，不能将有颜色理解为已检出强表达。CanDis_2026_Sample10约59.5%的spot A26 UCell为0，完整零分比例表已保存。MMP11、丰度和UMI采用全63共享99%上限，剪裁数逐图保存。
2. **哪些有聚集？** 按预先记录的操作性定义，63/63张呈正空间自相关且BH<0.05。A26 Moran I中位数0.713，范围0.109–0.932。未满足定义的切片：无。其中较弱的4张(描述性I<0.3，非预设生物学门槛)为：GSE171351_C1, GSE171351_D1, GSE285715_T3, GSE319536_BC10。统计可检测不等于视觉上形成强聚集。199次置换导致P值分辨率为0.005，不能区分所有达到该下限的切片谁更显著，应看I与连通性。

这表明混合spot程序有空间结构，不能证明离散C5 niche。严格top10%集合最大连通分量占选中spot比例的切片中位数为0.321，最小0.033，因此不少高分集合仍是多个分散区域。逐切片5%/10%/20%高分连通性见[组件表](03_spatial_clustering/components.tsv)。
3. **signature是否一致？** A26与core15、五基因、AUCell的逐切片相关分布见下表，包含差异而非只选最好切片。共享基因会产生部分一致性，不能当作独立验证。
4. **是否主要由Fibroblast驱动？** A26-high落在Fibroblast最高四分位的比例，切片中位数0.690；落在最低四分位为0.004。基线四分位各约25%，但丰度并列会改变比例。Fibroblast相关及调整UMI/上皮后的偏相关均保存；仍有5张呈负相关：CanDis_2026_Sample12, GSE238145_B11, GSE285715_N2, Inhousedata_T4471438, NPJ_sample5，说明不能把所有切片解释为同一CAF模式。匹配成功率切片中位数0.957，Fibroblast匹配后绝对SMD中位数0.067、UMI为0.022；但Fibroblast绝对SMD最大仍为0.863，不能将所有切片视为充分匹配。匹配不是细胞类型特异表达去卷积；不能凭残余结构认定C5真实存在。
5. **深度伪影？** UMI/nFeature关系见下表和逐dataset表。调整Fibroblast、UMI、nFeature和上皮后的rank残差，63/63张仍呈正自相关，残差Moran I中位数0.300。这排除了“全部空间结构只是这些变量线性rank关系”的简单解释，不能排除非线性、遗漏组织因素或共同程序影响。
6. **C2/C4竞争重叠？** 单细胞预定义联合C2/C4程序与A26高分集合重叠，切片中位数0.194，范围0.000–0.631。这是组合程序诊断，不能分别识别真实C2/C4，也不能用低重叠推翻既有C4困难负对照结果；对比程序仅覆盖所选20个基因。
7. **上皮/内皮/髓系空间关系？** 上皮mean占9大类总量≥0.3作为主定义，0.2/0.4敏感性；距离为观察spot图步数，非微米。主阈值可估计的切片61张，A26与距离相关中位数0.281。距离相关为正表示远离上皮富集区的score较高，为负则相反。内皮/髓系提供同spot匹配对比及ring0–3邻域曲线；下表效应正值表示高分区域周围高于匹配对照，负值表示低于对照。Fibroblast/UMI匹配的匹配率和SMD公开；相邻ring与两条件可能重叠，不能当独立重复，更不能解释为通讯/因果。
8. **跨来源/患者重复？** 下列逐dataset诊断和患者聚合CI展示方向是否一致。所有患者内多切片先聚合；未知患者只作切片描述。患者CI是探索性bootstrap，没有使用spot独立性；只重采患者，不包含dataset级抽样不确定性，已知患者仅3个来源，不能替代可靠的dataset层级推断。Normal切片患者ID未确认，无法完成可靠的患者配对Normal/Tumor检验。Normal仅来自GSE285715，不能把合并组差异解释为C5丰度变化；同来源GSE285715的3张Normal切片score中位数的中位数为0.099、3张Tumor为0.022，属于未配对描述，不能据此得出真实C5丰度在Normal更高，也不支持预设的Tumor必然更高。各组表见`07_dataset_heterogeneity/dataset_tissue_summary.tsv`。不比较固定top10%比例。留一dataset结果见[患者效应区间](07_dataset_heterogeneity/patient_effect_intervals.tsv)。
9. **值得独立验证的pattern？** 更值得验证的是多来源重复、调整后仍保留的A26-associated基质程序结构及其C5/C4细胞来源；与上皮成分分隔可作为描述性取样线索，而不是预设的血管/髓系富集机制。所有7个来源组的原始上皮距离相关中位数为正，但条件分析明显衰减。 可优先从同时具有A26连续空间结构、调整后结构仍保留、不同signature一致且多个来源重复的区域开展外部验证，但不能据此挑选/重训signature。本轮提供全部切片供人工评估，不重新给“正式C5定位PASS”。后续价值在于独立C5/C4/C2参考与多标记原位验证，检验这些区域是C5、其他CAF或共享ECM状态。

## 相关性与dataset异质性

以下均为逐切片相关的汇总，未将265,002 spots作为独立患者。

| variable | median | min | max |
| --- | --- | --- | --- |
| A26_AUCell | 0.959 | 0.844 | 0.993 |
| C2C4_program_UCell | 0.223 | -0.236 | 0.757 |
| core15_UCell | 0.811 | 0.239 | 0.967 |
| five_UCell | 0.764 | 0.303 | 0.974 |
| major_mean__Epithelial | -0.180 | -0.871 | 0.706 |
| major_mean__Fibroblast | 0.522 | -0.418 | 0.909 |
| major_mean__Mural | 0.298 | -0.360 | 0.792 |
| nCount | 0.023 | -0.671 | 0.759 |
| nFeature | 0.011 | -0.627 | 0.749 |

| dataset | n_sections | median_moran | A26_high_in_fib_high | A26_high_in_fib_low | C2C4_overlap | rho_Fib | rho_UMI | rho_features | rho_core | rho_five | rho_C2C4 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CanDis__UNVERIFIED | 13 | 0.701 | 0.617 | 0.003 | 0.221 | 0.434 | -0.213 | -0.081 | 0.825 | 0.810 | 0.279 |
| GSE171351 | 4 | 0.280 | 0.481 | 0.117 | 0.138 | 0.242 | 0.087 | 0.092 | 0.586 | 0.573 | 0.114 |
| GSE238145 | 10 | 0.758 | 0.581 | 0.001 | 0.061 | 0.607 | 0.445 | 0.441 | 0.861 | 0.849 | 0.067 |
| GSE285715 | 6 | 0.414 | 0.727 | 0.009 | 0.238 | 0.434 | -0.114 | -0.047 | 0.435 | 0.546 | 0.135 |
| GSE319536 | 16 | 0.780 | 0.882 | 0.000 | 0.203 | 0.600 | 0.198 | 0.164 | 0.743 | 0.756 | 0.192 |
| Inhousedata__UNVERIFIED | 7 | 0.713 | 0.577 | 0.006 | 0.194 | 0.528 | 0.025 | -0.053 | 0.883 | 0.861 | 0.477 |
| NPJ__UNVERIFIED | 7 | 0.646 | 0.720 | 0.005 | 0.383 | 0.525 | -0.400 | -0.403 | 0.766 | 0.740 | 0.565 |

## 匹配邻域的逐切片描述汇总

| variable | ring | n_sections | median_effect | n_positive |
| --- | --- | --- | --- | --- |
| Endothelial | 0 | 63 | -0.094 | 12 |
| Endothelial | 1 | 63 | -0.012 | 21 |
| Endothelial | 2 | 63 | 0.012 | 37 |
| Endothelial | 3 | 63 | 0.043 | 45 |
| Myeloid | 0 | 63 | -0.130 | 11 |
| Myeloid | 1 | 63 | -0.039 | 19 |
| Myeloid | 2 | 63 | 0.000 | 32 |
| Myeloid | 3 | 63 | 0.046 | 40 |

## 已知患者聚合的探索性效应

同spot丰度使用log1p差，距离使用Spearman系数；不能混用单位。每位患者先汇总重复切片，bootstrap 1000次。CI没有校正跨终点多重比较，因此不作为确认性显著发现；本轮仅对空间置换P值按score内63切片BH。小于3位患者不计算CI。

| analysis | variable | n_patients | median_effect | ci_low | ci_high | status |
| --- | --- | --- | --- | --- | --- | --- |
| epithelial_distance | 0.2 | 28 | 0.159 | 0.034 | 0.223 | EXPLORATORY_PATIENT_BOOTSTRAP |
| epithelial_distance | 0.3 | 27 | 0.136 | 0.065 | 0.291 | EXPLORATORY_PATIENT_BOOTSTRAP |
| epithelial_distance | 0.4 | 26 | 0.182 | 0.079 | 0.297 | EXPLORATORY_PATIENT_BOOTSTRAP |
| matched_same_spot | A26_UCell | 29 | 0.138 | 0.125 | 0.151 | EXPLORATORY_PATIENT_BOOTSTRAP |
| matched_same_spot | C2C4_program_UCell | 29 | -0.005 | -0.007 | -0.001 | EXPLORATORY_PATIENT_BOOTSTRAP |
| matched_same_spot | major_mean__Endothelial | 29 | -0.114 | -0.168 | -0.032 | EXPLORATORY_PATIENT_BOOTSTRAP |
| matched_same_spot | major_mean__Epithelial | 29 | -0.329 | -0.465 | -0.249 | EXPLORATORY_PATIENT_BOOTSTRAP |
| matched_same_spot | major_mean__Myeloid | 29 | -0.136 | -0.234 | -0.061 | EXPLORATORY_PATIENT_BOOTSTRAP |

## 上皮组成敏感性 (诊断后补充)

由于主匹配仍有上皮组成不平衡，追加了Fibroblast+UMI+Epithelial三变量匹配，未改signature或挑选空间图。下面既给全部可匹配切片，也给match rate≥0.5且三项绝对SMD≤0.1的平衡子集；后者是探索性敏感性，不能选取其中更有利的CI替代主分析。

| variable | subset | dataset | n_sections | median_section_effect | n_patients | median_patient_effect | ci_low | ci_high |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Endothelial | all_matched | ALL | 63 | -0.029 | 29 | -0.026 | -0.057 | 0.002 |
| Endothelial | balanced_50pct | ALL | 27 | -0.074 | 14 | -0.059 | -0.089 | -0.010 |
| Myeloid | all_matched | ALL | 63 | -0.013 | 29 | -0.032 | -0.053 | 0.006 |
| Myeloid | balanced_50pct | ALL | 27 | -0.015 | 14 | -0.028 | -0.041 | 0.008 |

上皮距离的逐切片相关中位数，控制Fibroblast/UMI/nFeature后为0.126，再控制同spot上皮丰度后为-0.021。这提示原距离关联很大程度可由局部组成解释；不能声称发现独立的上皮边界C5 niche。内皮/髓系负关联在额外上皮匹配后也明显减弱，全部匹配已知患者CI跨0。远端ring的正差异还受两组邻域重叠及残余组成影响，不支持C5特异促血管/免疫因果机制。

## 基因覆盖与少数基因贡献

签名-基因-切片记录中缺失特征5条，完整表记录每个基因的实际检测率。评分仅使用present genes；如缺失会改变分母，必须结合覆盖表解释。A26在63张均为26/26特征覆盖，但每切片每spot检出基因数的中位数范围为2.000–19.000；特征覆盖不等于实际表达检出。具体例外：GSE238145_B5竞争程序缺ADH1B/MT1M (18/20)，GSE285715_N2及T3的marker3缺COL10A1 (2/3)，GSE285715_T3 core15缺SFRP4 (14/15)。单基因MMP11只是诊断，不和多基因证据等价。

完整的UCell排名权重重建还用于逐切片核对官方评分实现，误差表见`00_inputs_and_QC/UCell_rank_reconstruction.tsv`。

在A26-high区域，按UCell正向rank权重描述，最主要的中位贡献基因为：MMP2 (11.0%), COL5A1 (10.7%), LUM (9.5%), MMP14 (8.5%), CTHRC1 (7.5%), THBS2 (7.5%)。每切片前三个基因贡献合计中位数34.2%。这是排名权重分解，不是基因因果效应，也不等于RNA贡献；未因此修改A26。详见[A26贡献表](04_fibroblast_and_depth_controls/A26_gene_rank_contributions.tsv)。

## 方法边界与结论层级

**可靠的描述性事实**：输入完整、分数可复现、所显示的连续程序和观察图上的聚集/相关/重叠可直接由保存的数据复算。颜色共用尺度，但不同dataset绝对score仍未校准。

**探索性关联**：Fibroblast/测序深度条件下残余结构、上皮距离及内皮/髓系邻域差异、组织标签分组与患者聚合方向。受空间自相关、匹配不完全、组织区域未知、非线性组成混杂、来源标签及患者数量限制。

**无法确认**：真实C5细胞丰度、C5独有空间niche、恶性侵袭前沿、细胞通讯、C5导致血管生成/免疫抑制。已有低比例混合spot及跨dataset泛化限制不因本轮63张图而消失。正式定位NO-GO保持；描述性程序分析已完成。

复现参数、原始版本来源、图像坐标、检验集合见[README](README.md)、冻结基因列表及脚本。输入与历史文件哈希复核通过。
