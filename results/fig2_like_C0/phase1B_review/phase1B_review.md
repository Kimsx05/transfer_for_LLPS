# C0 Phase 1B 补充核查

完成日期：2026-10-25；seed=123。唯一分母为冻结的 C0–C6 全部 Fibro+Pericyte。未改动注释、UMAP、LLPS77_AUC 或原始 Phase 1 图表；SHA256 核验见 sources/phase1_preservation_check.txt。绘图使用原 SCP Paired 配色，PDF 与 350 dpi PNG 独立输出。

## 1. 比例变化：方向较稳健，但不是绝对数量或保护作用

HRA003620 的七对身份由 Chen et al. 原文 Table 1 明确验证，并有 NGDC 个体标识交叉支持，不是依据名称猜测。原始来源：https://pmc.ncbi.nlm.nih.gov/articles/PMC9893503/#T1 。其中两个肿瘤标本没有保留的 stromal 细胞，分母为零，比例保留 NA，因此新增五对可用于比例分析。旧 Phase 1 的“独立样本”检验没有利用新证据，其 P 值仅作为历史结果，不再作为本轮主要推断；历史图不覆盖。

排除同一患者的额外标本 T4463209a（沿用事先规定保留 b），共 53 个有效标本、40 个患者单元。dataset-adjusted 线性模型估计 Tumor−Normal 为 −0.522（约 −52.2 百分点），患者聚类稳健 95% CI [−0.689, −0.355]，P=3.06e−7。协方差采用有限样本修正的 patient-cluster sandwich，t 自由度为患者簇数减设计矩阵列数。线性模型用于估计比例差，不用于产生越界预测。逐个移除 dataset 后方向均为负，P=2.39e−7 至 3.23e−4；要求分母至少 50 细胞的敏感性分析也保持方向。

13 对有效配对中 12 对下降、1 对不变，paired Wilcoxon P=0.00253，配对差中位数 −0.632。分层数据中 CNP0000460、GSE222315、HRA003620 均为负方向，但各自只有 5、3、5 对；各 dataset 的配对检验均未达到 0.05。Forest 图的 CI 是配对均值差的患者 bootstrap 区间，与 Wilcoxon 检验并非同一估计方法，小样本时不能将两者视为等价证据。GSE135337 只有一个 Normal，仅描述；InhouseData 无 Normal，不能提供独立的组织差异对照。

这一结果支持 C0 在所捕获 Fibro+Pericyte 中的相对构成下降，其他亚群扩张也可导致该现象。单细胞捕获不是单位组织绝对计数，不能推出 C0 绝对减少或保护作用。

## 2. C0 内部 LLPS：未发现稳定组织差异

按标本汇总原始 AUCell 中位数。主分析要求每标本至少 10 个 C0：15 Normal + 15 Tumor，26 个患者单元。组中位数分别 0.1420 与 0.1370；dataset-adjusted 差 −0.0103，95% CI [−0.0259, 0.00524]，P=0.182。四对满足该细胞数条件的标本 paired Wilcoxon P=0.584。至少 1/20 个 C0、逐 dataset 排除的敏感性分析均未支持稳定差异。

标本级 LLPS 与 median nFeature/nCount 强相关（Spearman rho 0.863 / 0.830），加入这两个技术协变量分别得到 P=0.662 / 0.469。技术调整不是证明生物学信号全由技术产生，但提示不能忽略测序复杂度。因有效 C0 较少和可用配对仅四对，阴性结果不等于两种组织状态完全相同。C0 整体 LLPS 中位数偏高的描述可以保留；不能宣称它显著高于 C5，也不能将比例下降解释为 LLPS 降低。

## 3. 功能特异性：炎症/TNF/缺氧较稳健；EMT/TGF-beta 需降调

使用既有注释确定 C4、C5 为主要 fibro-like comparator；C1/C3 为 mural、C6 为 cycling mural-enriched，排除于该参照。C2 既有审查仍为身份不明确，因此主敏感性参照不包含 C2，并另做 C2+C4+C5。没有重新分群或重命名。

沿用原样本筛选与 pseudobulk edgeR robust quasi-likelihood、样本阻断设计。每一侧至少 20 细胞：C4+C5 参照得到 11 个标本（8 Tumor、3 Normal），6208 个受检基因、240 个上调基因（FDR<0.05、log2FC>0.25）；含 C2 参照为 13 标本、8441 基因、565 上调。另用同样 11 标本重新估计 all-others 对比，以区分样本集合变化与参照变化。患者在每个 DEG 模型中只出现一次。GO BP/KEGG 使用各模型实际受检且可映射的基因作为背景；Hallmark 对完整受检排序表做 GSEA，BH 校正。

C0 vs C4+C5：TNF/NF-kB NES=2.425、FDR=8.28e−14；Hypoxia NES=1.766、FDR=0.000841；Inflammatory response NES=1.540、FDR=0.0328。三者方向和显著性保留。EMT FDR=0.123，而同 11 标本 all-others 为 0.000408，提示细胞身份参照影响较大。TGF-beta FDR=0.140，同 11 标本 all-others 也不显著（0.117），因此其变化还与样本集合及统计功效有关。含 C2 参照时 EMT/TGF-beta 又达到 FDR<0.05，应报告参照依赖性，不能据此重新定义 C2。

GO 的发育、形态发生、迁移等主题部分保留；在 fibro-like 参照下 KEGG 仅 complement/coagulation 达 FDR<0.05（0.00209），Wnt 和 TGF-beta 未保留显著性。不能将原始 all-others 的全部富集均解释为 C0 特异机制。

原始 leading edge 中 EMT 有 36/81（44.4%）基因为 GO ECM 注释，其中 COL*/FN1 三个；炎症 3/44、TNF 5/84、缺氧 9/62、TGF-beta 3/23 为 ECM。故炎症/TNF 并非主要由胶原驱动，但这些通路也并非独立证据：SERPINE1 等跨多个通路重复，完整交集/Jaccard 表已保存。ECM 定义为 org.Hs.eg.db 的 GO:0031012/GO:0062023 及子项，通过 GOALL 获取，不等同于所有可能的基质功能基因。

## 4. 完成文件与解释边界

本目录含八张独立正式图，每张 PDF + PNG：

- FigS_C0_fraction_dataset：dataset 分层构成；FigS_C0_fraction_forest：分层效应量；FigS_C0_fraction_paired_updated：更新配对。
- FigS1_C0_LLPS_TumorNormal、FigS2_C0_LLPS_paired、FigS3_C0_LLPS_dataset：标本级 LLPS。
- FigS_leading_edge_gene_pathway：基因通路关系；FigS_functional_comparator_GSEA：参照敏感性。

完整绝对路径清单见 file_index.tsv。输入、处理和检验可在 scripts/01_inputs_and_pairing.R、02_sample_statistics.R、03_functional_specificity.R 与 common.R 复现。必要稀疏 sample×cluster counts 仅保留一份；没有复制单细胞大对象。原始身份审查路径保存在 fibro_like_comparator_definition.tsv 的来源说明及代码 FIBROOT 中。

C0 更符合 PI16/SFRP2-associated、保留基质/发育程序且可响应炎症与缺氧的 fibroblast 状态；当前不支持把它命名为“肿瘤特异 LLPS 增强”或“保护性”亚群。比例复核没有破坏 C0 定义或 pseudobulk 基础，因此可以进入 Phase 2；外部丰度结论仍必须等待患者独立参考验证。
