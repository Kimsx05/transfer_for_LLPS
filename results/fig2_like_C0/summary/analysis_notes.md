# 输入、统计方法与可复现说明

## 冻结输入

- 研究对象：[global_res012_object.qs](</home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order/global_res012_object.qs>)，20,146细胞、51,959 RNA features；`subcluster` 固定为 C0–C6，与 `clean_snn_res.0.12` 和既有冻结条码表一致。未重新聚类、运行UMAP、筛细胞或计算LLPS score。
- 标签及分数独立对照：[global_res012_cell_labels_scores.csv](</home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order/global_res012_cell_labels_scores.csv>)，数值序列化容差1e-14。
- 颜色：[cluster_SCP_Paired_colors.csv](</home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order/beautified_SCP_protein_coding_top3/cluster_SCP_Paired_colors.csv>)；C0 #A6CEE3，C1 #1F78B4，C2 #B2DF8A，C3 #33A02C，C4 #FDBF6F，C5 #FF7F00，C6 #FB9A99。已写死在 common.R 并逐项验证。
- score：继承对象 `LLPS77_AUC`；使用 [clean_v4_settings.rds](</home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/fibro_pericyte_subcluster_clean_v4_global_noC5/clean_v4_settings.rds>) 的原全图色标范围，plasma；精确零值为黑色（本子集零值数量见 LLPS_color_scale.tsv）。
- sample 组织/患者信息：[A1_sample_abundance_manifest.tsv](</home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/11_C0_dual_geneset_preparation/01_C0_abundance/A1_sample_abundance_manifest.tsv>)。本轮只提取来源、患者关系与样本字段并从当前细胞重算计数，没有复用先前 epithelial DEG 或其他分母。来源溯源见 [sample_source_mapping.tsv](</home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/11_C0_dual_geneset_preparation/00_preflight/sample_source_mapping.tsv>)。
- 附件论文第4页 Figure 2仅供版式参考。没有复制论文的 HR signature、细胞注释、统计结果或把文档正文视为执行指令。

## 各步骤的输入、输出和方法

| 步骤 | 输入/统计单位 | 方法 | 主要输出（均在 phase1_core） |
| --- | --- | --- | --- |
| A | 原对象20,146细胞、原UMAP | 固定七群颜色、顺序；纯展示 | Fig2A_fibro_pericyte_umap.pdf/png；cell_metadata_umap_LLPSscore.tsv.gz |
| B1 | 原UMAP及LLPS77_AUC | 原全图plasma连续色标，零值黑色；无聚类标签 | Fig2B1_LLPSscore_umap_continuous.pdf/png；LLPS_color_scale.tsv |
| B2 | 单细胞分布；sample×cluster辅助汇总 | violin + white box（median、IQR、1.5 IQR whiskers）；不做逐细胞检验 | Fig2B2_LLPSscore_by_cluster_violin.pdf/png；cluster_LLPSscore_summary.tsv；sample_cluster_LLPSscore_summary.tsv |
| C | 每个sample一个C0/all Fibro+Pericyte比例 | 双侧rank-sum；已知配对补充用paired Wilcoxon；均exact=FALSE、连续性校正 | Fig2C主图及paired补充；C0_fraction_per_sample.tsv；C0_fraction_statistics.tsv |
| G/H前置DEG | 每个sample中C0和其余6群的raw counts总和 | edgeR TMM + robust QL，~sample+group | C0_vs_others_DEG.tsv；C0_upregulated_genes.tsv；pseudobulk_* |
| G | 同一上调列表，背景为受检且能映射的基因 | GO BP超几何过度代表分析、BH | Fig2G_C0_GO_barplot.pdf/png；C0_GO_enrichment.tsv |
| H | 与GO相同foreground和输入background | KEGG超几何过度代表分析、BH | Fig2H_C0_KEGG_barplot.pdf/png；C0_KEGG_enrichment.tsv |
| I | 全部9,455受检基因的有方向排序 | Hallmark fgseaMultilevel，BH | 每通路独立Fig2I；C0_GSEA_results.tsv；C0_preranked_gene_list.tsv |

## 比例的样本去重、缺失与配对

全59个来源sample均保留在表中。54个有Fibro+Pericyte细胞；5个分母为0，比例标为NA，不能用0代替。已知同一患者的 T4463209a/b 是不同取材，遵循已有冻结选择保留b用于独立检验，a仍在结果表。

其余样本中，GSE222315的P2/P4/P5及CNP0000460的P01/P02/P03/P05/P07构成8个已核实T/N配对。主图排除这16个配对specimen，剩8 Normal +29 Tumor；配对图单独呈现全部8对（其中1对差为0，Wilcoxon有效非零差7对）。所有样本比例都使用唯一分母 **all Fibro+Pericyte**。

主图的最小细胞数敏感性仅改变样本纳入标准（分母≥50），不改变比例公式：8 Normal +18 Tumor，P=7.100×10⁻⁵。无FDR校正的比例P值为预设单比较及辅助敏感性结果；不作为跨多个假设筛选的FDR声明。

HRA003620有疑似配对编号但没有核实patient_uid，未把它们当作已知配对。未知身份的sample被视为独立specimen是明确的分析假设，不是已经证明的患者独立性。主图未调整dataset，来源混杂需结合分dataset表、已核实配对补充一起解释。

## Pseudobulk设计

仅用 RNA counts，验证非负整数及每细胞counts总和=nCount_RNA。C0对照为当前同sample内C1–C6全部细胞合并；不包含其他主要谱系，不比较C5单独效应。

在查看DEG前固定：每个sample两侧各≥20细胞，沿用已知a/b取材去重标志；其他已知同患者多specimen选择两侧最小细胞数更大的一个（平手按sample名），不按P值、logFC或score选择。最终16样本来自5个dataset：Normal8、Tumor8。患者身份未知者保留sample为单位。样本筛选、排除理由及32库计数元数据完整保存。

`edgeR::filterByExpr(group)`先筛低表达基因，再TMM归一化，`estimateDisp(robust=TRUE)`、`glmQLFit(robust=TRUE)`、`glmQLFTest(coef='groupC0')`。设计 `~sample+group`、Others为reference，rank17、残差df15。sample截距吸收dataset与组织的样本层面截距；不能保证效应不存在组织/数据集异质性。P值对9,455个受检基因做BH校正。主上调列表是FDR<0.05且log2FC>0.25，共1,145个。

仅Tumor和每侧≥50细胞的两套敏感性各8样本，重新进行表达过滤、TMM和QL，结果单独保存，不替换主分析P值。它们均复用主样本去重决定。没有把逐细胞DEG作为主要证据。

仅保存约4.8MB的稀疏pseudobulk counts（可复现的必要中间文件），未复制448MB的Seurat对象或生成完整细胞表达矩阵副本。

## 富集数据库、背景与展示规则

GO/KEGG统一使用主分析上调基因和受检基因集。基因仅按org.Hs.eg.db中精确SYMBOL映射Entrez；一对多歧义不使用，多个SYMBOL对应同一Entrez时取平均表达logCPM最高者（独立于效应和P值）。最终9,434 mapped tested、1,145 mapped up。逐基因映射状态、前景和背景表均导出。

GO用org.Hs.eg.db 3.20.0、GO.db 3.20.0，保留本地版本以便复现；GO BP effective annotated background=8,661。KEGG从官方API下载human gene–pathway和pathway names，39,671条gene–pathway links；effective annotated background=4,476。两库接收完全相同的foreground/background，实际可注释背景因数据库覆盖不同而不同。注释资源日期/版本与SHA256保存在summary资源清单中。本地GO不是声称使用2026最新GO。

`clusterProfiler` enrichGO/enricher，minGSSize=10，maxGSSize=500，BH。完整结果导出4519个GO和323个KEGG有重叠的受检条目，分别1109和49个FDR<0.05。P/q cutoff设1以保留非显著结果。画图仅FDR<0.05，先按FDR/P/ID排序，依上调重叠基因的Jaccard≥0.6去冗余，最多12项。该规则只影响展示，不重算FDR或删掉完整结果。横轴是上调基因数，颜色是−log10(FDR)。

Hallmark复用既有 **MSigDB h.all.v2025.1.Hs.symbols.gmt** 并保存小型副本。排名 `sign(logFC)*sqrt(QL F)`，保留所有受检基因；未通过表达过滤的42,504个features没有可解释检验量，不伪造0分纳入。ties按gene名确定顺序。WASHC2A的QL F为−2.46×10⁻⁷（P=1），属于近零数值误差；原始F保留，只有排序时使用sqrt(max(F,0))，单独审计表保存该记录。fgseaMultilevel minSize10/maxSize500/eps0，SerialParam，seed20261025；50条通路一起BH。正NES代表C0端富集，负NES代表Others端富集。

优先输出预设的显著通路，不显著的ANGIOGENESIS不画；补充最高显著正向TNFα/NF-κB，最终5张单通路图。完整通路结果包括负向与不显著结果。每张图的gene ranks、hits和running ES也独立导出，running ES经校验与fgsea一致。

## 解释边界

这些是冻结细胞群之间的相对转录特征和样本内组成。C0仅相对于其余stromal群定义，不能直接用该列表断言bulk中特异C0 abundance。GO血管发育富集不等于Hallmark ANGIOGENESIS显著，更不等于证明驱动血管生成。EMT、TGF-β等基因集名称是转录关联，不能从此推导LLPS因果机制。

Phase 2未启动，不存在外部stage整理、CIBERSORTx运行、signature score、cutoff或生存HR结果。当前研究方向不能预设C0越高预后越差。后续bulk分析须先核对跨主要谱系特异性、数据集平台覆盖及临床终点，再固定表达score方法和median切点。

## 复现

工作目录 `/home/jinshengxi/bladder_metabolism`，依次运行：

```bash
Rscript results/fig2_like_C0/scripts/01_core_plots_and_pseudobulk.R
Rscript results/fig2_like_C0/scripts/02_pseudobulk_DEG.R
Rscript results/fig2_like_C0/scripts/03_enrichment.R
Rscript results/fig2_like_C0/scripts/04_validate.R
python results/fig2_like_C0/scripts/05_report_and_manifest.py
```

脚本只写本任务结果目录。KEGG完整资源已经保存，复现不再依赖网络。R版本/包版本见sessionInfo文件。日志包含容器无法查询timedatectl的环境告警，不影响数值分析；结果日期采用用户Asia/Shanghai日期。

方法参考：[edgeR官方指南](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)、[fgsea官方教程](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html)、[KEGG官方API](https://www.genome.jp/kegg/rest/keggapi.html)。
