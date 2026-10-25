# C0 Fig.2-like：Phase 1 结果汇报

完成日期：2026-10-25（Asia/Shanghai）。**Phase 1 已完成；按照“先汇报，再决定 Phase 2”的要求，外部 bulk 分期/生存尚未启动。**

输出目录：[phase1_core](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core>)。共 **12 张独立图，每张 PDF + 350 dpi PNG**，没有拼图。脚本位于 [scripts](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/scripts>)；详细输入、统计方法和限制见 [analysis_notes.md](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/summary/analysis_notes.md>)。

## 关键结果

- **C0 为 LLPS-high**：10,717 个 C0 细胞，LLPS77 AUCell 中位数 **0.14134**，七群最高；按 sample 汇总后的中位数也居首。B2 为描述性分布图，没有逐细胞大规模标星。
- **C0 proportion 在 Tumor 中下降，未观察到升高。** 主图 Normal n=8、Tumor n=29，中位数分别 **97.39% 与 9.68%**；双侧 Wilcoxon rank-sum P=**2.652×10⁻⁵**。已确认的 T/N 配对单独分析：8 对，中位数 73.23% 与 10.11%，配对 Wilcoxon P=**0.02249**。未知患者关系及来源混杂限制主图 P 值的解释。
- **pseudobulk DEG**：16 个样本（8 Normal、8 Tumor；5 个数据集），32 个 C0/Others 聚合库，9,455 个受检基因。`~sample + group` 的 edgeR QL 模型识别 **1,145 个上调基因**（BH FDR<0.05，log2FC>0.25）。
- **GO**：血管发育、上皮细胞增殖、细胞外结构组织、生长因子响应等；血管发育 FDR=1.78×10⁻¹⁴。**KEGG**：Wnt、补体/凝血、细胞因子相互作用及 IL-17/TNF 等信号。
- **GSEA**：TNFα/NF-κB、炎症反应、缺氧、EMT 和 TGF-β 为已输出的 5 条显著正向通路。50 条 Hallmark 全部保留于结果表，21 条 FDR<0.05。ANGIOGENESIS FDR=0.455、胆固醇稳态 FDR=0.439、脂肪酸代谢 FDR=0.777，未强行出图。氧化磷酸化为显著负向通路（NES=-2.972，FDR=1.08×10⁻²⁵），完整结果亦保留。
- **稳健性辅助结果**：仅 Tumor 的8样本分析中，主分析上调基因里1,110个可检验，方向全部仍为正，607个仍满足上调门槛；每侧≥50细胞的8样本分析中，1,137个可检验、方向均正，583个仍满足门槛。不能据此声称所有功能通路已在独立队列验证。

## 单独图文件

| 图 | PDF | PNG |
| --- | --- | --- |
| Fig2A_fibro_pericyte_umap | [PDF](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2A_fibro_pericyte_umap.pdf>) | [PNG](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2A_fibro_pericyte_umap.png>) |
| Fig2B1_LLPSscore_umap_continuous | [PDF](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2B1_LLPSscore_umap_continuous.pdf>) | [PNG](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2B1_LLPSscore_umap_continuous.png>) |
| Fig2B2_LLPSscore_by_cluster_violin | [PDF](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2B2_LLPSscore_by_cluster_violin.pdf>) | [PNG](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2B2_LLPSscore_by_cluster_violin.png>) |
| Fig2C_C0_fraction_tumor_vs_normal | [PDF](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2C_C0_fraction_tumor_vs_normal.pdf>) | [PNG](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2C_C0_fraction_tumor_vs_normal.png>) |
| Fig2C_C0_fraction_tumor_vs_normal_paired | [PDF](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2C_C0_fraction_tumor_vs_normal_paired.pdf>) | [PNG](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2C_C0_fraction_tumor_vs_normal_paired.png>) |
| Fig2G_C0_GO_barplot | [PDF](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2G_C0_GO_barplot.pdf>) | [PNG](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2G_C0_GO_barplot.png>) |
| Fig2H_C0_KEGG_barplot | [PDF](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2H_C0_KEGG_barplot.pdf>) | [PNG](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2H_C0_KEGG_barplot.png>) |
| Fig2I_GSEA_EMT | [PDF](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2I_GSEA_EMT.pdf>) | [PNG](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2I_GSEA_EMT.png>) |
| Fig2I_GSEA_HYPOXIA | [PDF](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2I_GSEA_HYPOXIA.pdf>) | [PNG](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2I_GSEA_HYPOXIA.png>) |
| Fig2I_GSEA_INFLAMMATORY_RESPONSE | [PDF](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2I_GSEA_INFLAMMATORY_RESPONSE.pdf>) | [PNG](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2I_GSEA_INFLAMMATORY_RESPONSE.png>) |
| Fig2I_GSEA_TGFB | [PDF](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2I_GSEA_TGFB.pdf>) | [PNG](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2I_GSEA_TGFB.png>) |
| Fig2I_GSEA_TNFA_SIGNALING_VIA_NFKB | [PDF](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2I_GSEA_TNFA_SIGNALING_VIA_NFKB.pdf>) | [PNG](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/Fig2I_GSEA_TNFA_SIGNALING_VIA_NFKB.png>) |

## GSEA 关注通路

| pathway | NES | padj | exported |
| --- | --- | --- | --- |
| HALLMARK_ANGIOGENESIS | -1.07689499730248 | 0.455082331689015 | FALSE |
| HALLMARK_CHOLESTEROL_HOMEOSTASIS | 1.10141751557747 | 0.438834600072385 | FALSE |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 1.72784902533808 | 0.000420634741612692 | TRUE |
| HALLMARK_FATTY_ACID_METABOLISM | 0.878378067324375 | 0.777301046197249 | FALSE |
| HALLMARK_HYPOXIA | 1.78967229263135 | 0.000137783992246076 | TRUE |
| HALLMARK_INFLAMMATORY_RESPONSE | 2.10946094711882 | 9.90137878965115e-07 | TRUE |
| HALLMARK_TGF_BETA_SIGNALING | 1.57545109917505 | 0.0375165307480429 | TRUE |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | 2.54431937996606 | 2.1015027258704e-17 | TRUE |

## 核心数据表

- [cluster_LLPSscore_summary.tsv](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/cluster_LLPSscore_summary.tsv>)
- [C0_fraction_per_sample.tsv](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/C0_fraction_per_sample.tsv>)
- [C0_fraction_statistics.tsv](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/C0_fraction_statistics.tsv>)
- [C0_vs_others_DEG.tsv](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/C0_vs_others_DEG.tsv>)
- [C0_upregulated_genes.tsv](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/C0_upregulated_genes.tsv>)
- [C0_GO_enrichment.tsv](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/C0_GO_enrichment.tsv>)
- [C0_KEGG_enrichment.tsv](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/C0_KEGG_enrichment.tsv>)
- [C0_GSEA_results.tsv](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/C0_GSEA_results.tsv>)
- [C0_preranked_gene_list.tsv](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/C0_preranked_gene_list.tsv>)
- [pseudobulk_model_summary.tsv](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/pseudobulk_model_summary.tsv>)
- [C0_DEG_sensitivity_comparison.tsv](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/C0_DEG_sensitivity_comparison.tsv>)
- [validation_checks.tsv](</home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/validation_checks.tsv>)

## 需要人工核对的事项

1. **HRA003620 患者配对身份**：既有映射仅含候选配对编号，patient_uid 未核实。本轮没有根据相似名字构造配对；主分析暂以 sample/specimen 为单位。因此 rank-sum 的独立性仍是工作假设，需要来源方核对；8 对补充分析仅使用已核实 study-patient IDs。
2. **临床含义**：Normal 包括邻近正常组织；组织来源和数据集组成不平衡。比例下降是该冻结子集内的结果，不能直接解释为肿瘤整体绝对细胞数下降。
3. **颜色已解决，无需重新分配**：复用后续 publication-style 图明确保存的 SCP Paired mapping；更早的探索图曾用 Dark 3，本轮没有采用旧的探索版颜色。
4. **富集不等于机制验证**：GO 血管相关条目显著，而 Hallmark ANGIOGENESIS 不显著，两者必须分别汇报。EMT signature 在既有间质细胞中的富集不能证明它们发生上皮到间质转化，亦不能证明 LLPS 导致这些功能。
5. **Phase 2 尚未执行**：尚无外部 stage、HR、KM 或 abundance 结论；C0-vs-stromal markers 还需检查相对其他主要谱系的特异性后再用于 bulk。若采用 expression-based score，必须明确区别于严格去卷积丰度。
