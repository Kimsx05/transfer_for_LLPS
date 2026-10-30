# C5 parallel v1

完成日期：2026-10-30。六模块均已执行。结果解释见 [result_summary.md](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/result_summary.md)，证据对比见 [C0_vs_C5_evidence_comparison.md](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/C0_vs_C5_evidence_comparison.md)。

## 复用代码与修改范围

项目根：`/home/jinshengxi/bladder_metabolism`。所有C0源代码和结果只读使用；新代码及输出仅写入本目录。没有调用新聚类、注释、LLPS/HR打分、UMAP/Harmony或CIBERSORTx。没有下载Bulk或基因集，没有创建环境冻结或文件哈希。

| 模块 | C0源代码（相对项目根） | 本轮脚本 | 修改 |
| --- | --- | --- | --- |
| 1 Normal/Tumor | results/fig2_like_C0/phase1B_review/scripts/common.R、02_sample_statistics.R；scripts/01_core_plots_and_pseudobulk.R | scripts/01_composition.R | C0→C5；采用reviewed配对；新输出路径；仅截取丰度区段；补100%堆叠图和全atlas分母；300dpi透明背景 |
| 2 分期 | results/fig2_like_C0/music_optimization/scripts/06_clinical.R；既有工作簿解析结果及样本映射 | scripts/02_stage.R | 复用stage色表、box+jitter、KW/Spearman/Wilcoxon逻辑；将已映射工作簿T成分与C5比例连接；加均值柱图；稀疏组不做成对推断 |
| 3 Bulk OS | results/cibersortx_phase2_clinical_v1/scripts/01_integrate.py、02_analyze.R | scripts/03_bulk_inputs.py、03_bulk_survival.R、03b_reference_comparison.R | 直接复用临床整合表，不重跑整合下载；动态提取fit_cox、choose_covars、KM和forest代码；目标C0/C5；仅OS；增加PNG和横向比较表 |
| 4 DEG | results/llps_C0_functional_signature_v2/scripts/01_markers.R、common.R；Fig2 scripts/02_pseudobulk_DEG.R；llps_C0_functional_state_v1/scripts/04_DEG_LLPS.R火山图段 | scripts/04_markers.R、04b_pseudobulk.R、04c_volcano.R | C5 vs setdiff(C0–C6,C5)；焦点基因换POSTN/COL11A1/COL10A1；对称down标准；原pseudobulk函数及过滤；不足4样本记录NOT_FIT；火山图沿用几何但按v2阈值与Bonferroni轴名 |
| 5 ORA | results/llps_C0_functional_signature_v2/scripts/02_ORA.R；v3/scripts/01_curated_ORA.R绘图表达式 | scripts/05_ORA.R、05b_ORA_figures.R | 新C5前景及背景；up/down输出；保留实际enrichGO/enrichKEGG/简化；用新显著ID替换C0固定10条；10条GO、实际4条KEGG |
| 6 GSEA | results/llps_C0_functional_signature_v3/scripts/02_GSEA.R、common.R | scripts/06_GSEA.R | 指向C5排名；保留C0 v2同版Hallmark资源；去掉不适用的C0 v2 NES对照；曲线由C5真实结果选择；阈值/重复symbol/fgsea/布局不变 |
| 资源与图检 | results/llps_C0_functional_signature_v3/scripts/run_limited.sh、08_figure_QC.py | scripts/run_limited.sh、08_figure_QC.py | 运行限制原样复用；图检仅去掉C0专用的固定17张断言，动态列出本轮38张 |

原脚本均先阅读后复用。动态提取代码片段的脚本读取原始R文件并在当前新输出路径执行，未编辑原文件。模块2没有找到同口径的单细胞C0分期图，因而复用原Fig2临床图配色和统计结构，并使用已有工作簿映射生成C5单细胞丰度分析，未将Bulk signature分期结论直接移植。

实际参数差异：C0 Fig2早期common的seed为20261025，本轮统一123；早期Fig2图为白底350dpi，按本轮明确要求改为透明300dpi。C0 v2 marker的`p_val_adj`是**Bonferroni**，并非BH FDR，本轮保留原判定；另存`FDR_BH`但不用于改变marker筛选。C0 v3 GSEA参数与任务描述完全一致。

## 输入及分析对象

- 冻结Fibro/Pericyte对象：[global_res012_object.qs](/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order/global_res012_object.qs)
- 冻结细胞标签/已有分数：[cell_metadata_umap_LLPSscore.tsv.gz](/home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1_core/cell_metadata_umap_LLPSscore.tsv.gz)
- 复核样本/患者/配对：[sample_metadata_reviewed.tsv](/home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase1B_review/sample_metadata_reviewed.tsv)
- 精确全atlas分母及工作簿行映射：[A1_sample_abundance_manifest.tsv](/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/11_C0_dual_geneset_preparation/01_C0_abundance/A1_sample_abundance_manifest.tsv)
- 已有工作簿临床行（F列为TNM）：[user_workbook_relevant_scRNA_rows.tsv](/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7/11_C0_dual_geneset_preparation/00_preflight/user_workbook_relevant_scRNA_rows.tsv)
- Full14/Focused参考细胞映射：[reference_cells_LOCAL_ONLY.tsv](/home/jinshengxi/bladder_metabolism/fig2_C0_multimethod_deconvolution/03_methods/CIBERSORTx_web/manifest/reference_cells_LOCAL_ONLY.tsv)
- Full14原始输出与队列清单：[input_inventory.tsv](/home/jinshengxi/bladder_metabolism/results/cibersortx_phase2_clinical_v1/input_inventory.tsv)
- 已整理Bulk表达和临床位置：[01_bulk_cohorts](/home/jinshengxi/bladder_metabolism/fig2_C0_multimethod_deconvolution/01_bulk_cohorts)
- 原OS定义及来源：[endpoint_definitions.tsv](/home/jinshengxi/bladder_metabolism/results/cibersortx_phase2_clinical_v1/endpoint_definitions.tsv)
- Hallmark同版资源：[gene_sets.rds](/home/jinshengxi/bladder_metabolism/results/llps_C0_functional_signature_v2/resources/gene_sets.rds)
- Hallmark版本记录：[gene_set_sources.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C0_functional_signature_v2/resources/gene_set_sources.tsv)
- KEGG人类快照：[KEGG_hsa_gene_pathway.tsv](/home/jinshengxi/bladder_metabolism/results/fig2_like_C0/resources/KEGG_hsa_gene_pathway.tsv)
- KEGG通路名称快照：[KEGG_hsa_pathway_names.tsv](/home/jinshengxi/bladder_metabolism/results/fig2_like_C0/resources/KEGG_hsa_pathway_names.tsv)
- gene symbol映射：[gencode_v36_gene_map.tsv](/home/jinshengxi/bladder_metabolism/fig2_C0_multimethod_deconvolution/03_methods/CIBERSORTx_web/resources/gencode_v36_gene_map.tsv)

Bulk四队列OS及协变量直接读取原 `results/cibersortx_phase2_clinical_v1/by_dataset/<cohort>/clinical_merged.tsv`，其具体原始去卷积文件、样本数、临床来源与映射核对记录在 [input_inventory.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/input_inventory.tsv)。Full14的7个队列及现有focused输出列清单见 [existing_deconvolution_inventory.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/existing_deconvolution_inventory.tsv)。没有找到单独保存的既有Full14/focused最终一致性统计表；因此仅对已有TCGA和UROMOL输出新增描述性相关性核对，未重新运行去卷积。

- 冻结对象：51,959 RNA features ×20,146 cells；C5=716，Rest=C0+C1+C2+C3+C4+C6=19,430。RNA已有LogNormalize data（scale.factor=10000），未重归一化。
- 组成：保留原59个来源样本清单；分母为0时比例NA，C5=0且分母有效时保留。复核主要分析53样本、40 patient_key；16 Normal/37 Tumor；13对有有效分母的已验证配对。
- 分期：37个肿瘤样本/37 patient_key；NMIBC/MIBC=17/20；Ta/T1/T2/T3/T4=3/14/11/5/4。只解析工作簿明确T字段，保留clinical/pathologic/unspecified类型。采用原样本独立性标记，不从样本名新推断患者关系。
- OS：TCGA 403/176事件、GSE31684 93/65、GSE13507 165/69、GSE48075 73/45；调整样本/事件分别399/175、93/65、165/69、72/44。UROMOL、GSE32894和GSE87304没有对应有效OS，未纳入OS推断。
- Pseudobulk：原始RNA counts；每样本C5和Rest至少各20细胞；已知患者保留一个样本；13样本/26文库、两dataset、均Tumor；7,377基因可检验。patient_key中未知患者以既有sample单位表示，不宣称获得了新增患者身份验证。

## 方法

1. 组成使用原phase1B配对Wilcoxon、分dataset比较/5,000次bootstrap，以及`cluster_lm`的dataset固定效应与患者聚类协方差；沿用leave-dataset-out和min50敏感性。独立样本图排除已验证配对。堆叠图同时给出pooled-cell和sample-equal版本，两者仅为描述性图；颜色C0–C6为原SCP Paired映射。
2. 分期使用样本比例，NMIBC/MIBC Wilcoxon，T组Kruskal–Wallis及探索性Spearman；不作稀疏组pairwise显著性筛选。散点仅横向jitter，避免把零比例抖到坐标范围外。
3. Cox直接复用`fit_cox`（Efron ties、cox.zph transform=km）和`choose_covars`：分期优先，其次age/sex/grade，基于缺失/变异及原10事件/参数条件，绝不按P选变量。TCGA使用原病理总体AJCC分期I/II、III、IV调整，与单细胞T分期分开。每1 SD由**所有原队列去卷积样本**确定；KM亦使用该总体中位数，High>median。C0与C5完全同口径；跨队列BH另存，不混合端点、不作合并HR。
4. DEG复用`FindMarkers(test.use='wilcox', min.pct=0, logfc.threshold=0, only.pos=FALSE, random.seed=123)`生成完整表；marker候选原参数min.pct=.25/logfc.threshold=.75/only.pos=TRUE，最终Bonferroni P<.05且logFC>.75。down对称。51,959行均保留；ranking再排除无表达/无效统计，symbol重复按检出细胞数、整体均值、原名决定，非按P或方向择优。样本方向检查复用原top30+指定基因逻辑。Pseudobulk原edgeR TMM、filterByExpr(group)、robust estimateDisp/glmQLFit、`~sample+group`，原BH FDR<.05且logFC>.25用于该验证表。
5. ORA原enrichGO(org.Hs.eg.db/BP)/enrichKEGG(hsa)，BH，p/q cutoff=.05，min/max gene set=10/500；GO Wang simplify cutoff=.7。背景为任一组检出≥25%的可映射基因；GO/KEGG按各数据库注释交集形成有效背景。GO原库版本及KEGG实际快照记录见ORA_annotation_sources.tsv。前景分up/down；选择的ID/FDR/Count/GeneRatio与完整结果逐一核对。绘图横轴Count，浅黄`#FFF3B0`至橙红`#D84120`表示−log10(FDR)。
6. GSEA使用C0 v2已缓存的MSigDB **2026.1.Hs Hallmark50**（不是该资源中KEGG_LEGACY来做ORA）。排名为新C5 avg_log2FC，所有符合检测阈值且有效的基因；主10%，敏感5%/25%；fgseaMultilevel eps=0、nPermSimple=10000、minSize10/maxSize500、SerialParam(RNGseed123)。阈值内分别BH；运行曲线逐一核对ES。四条实际代表曲线按原4正/2负上限客观选择。

## 直接复用与新计算

直接复用：冻结对象及已有表达层、身份和颜色、sample/patient映射、工作簿临床记录、全atlas分母、Full14/Focused去卷积数值、OS/临床整合表、GO库/KEGG快照/Hallmark资源、C0已完成结果。既有LLPS与HR signature只在比较总结中引用，未重新计分。

新计算：C5样本组成/分期统计与图、C5完整DEG与样本验证、C5 up/down ORA、C5三阈值GSEA、同口径C0/C5 OS模型和KM、既有Full14/focused的C0相关性核对。没有新增signature、风险模型或用户范围之外的分析。

## 完成状态与边界

六模块主要分析均成功。预先规定的min50 pseudobulk敏感性仅3样本，记录NOT_FIT；这是遵守原至少4样本规则，未调低门槛。COL10A1在pseudobulk表达过滤后未检验；非所有逐细胞marker都获样本层支持。KEGG-up仅4条显著，不补足10条。25% GSEA有1个基因集尺寸不足，保留状态。三个没有OS的Bulk队列不作OS；focused无独立C5。未选做C5/全部“Fibro”去卷积比值：Full14名称Fibro_C0–C6中包含mural/pericyte类别，不将其总和误称全部可识别成纤维细胞。

C5来源分布与参考来源偏倚、稀疏样本分母、GSE13507零中位数和调整模型全局PH偏离均在中文总结中报告。C0/C5各自Rest不同，不能直接把NES差异当作C0 vs C5对比检验。

## 运行与验证

从项目根运行：

```bash
bash results/llps_C5_parallel_v1/scripts/run_all.sh
```

R任务依次运行，复用run_limited.sh限制CPU亲和性最多8核，BLAS/OpenMP各1线程，future sequential，data.table单线程，seed123。没有并行高内存任务。DEG已有完整表时脚本可读其缓存；GSEA也沿用原结果缓存行为。若主动更改输入/算法，应另建版本目录，避免将缓存误作新结果。

[validation.txt](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/validation.txt)：核对716细胞、分母/零值、完整DEG及up/down规则、选定ORA统计精确对应、GSEA阈值与BH；本轮重算的8个C0 Cox HR/N与既有结果一致。[figure_QC.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/figure_QC.tsv)：**38对独立矢量PDF/透明300dpi PNG**；所有PDF单页且无栅格嵌入，PNG alpha正常。已目视核对代表性组成、分期、KM、火山、ORA及GSEA图，必要的零点显示修正已完成。

## 主要数据表

- [01_Normal_Tumor/C5_fraction_per_sample_reviewed.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/C5_fraction_per_sample_reviewed.tsv)
- [01_Normal_Tumor/all_cluster_sample_composition.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/all_cluster_sample_composition.tsv)
- [01_Normal_Tumor/C5_fraction_adjusted_and_leave_dataset_out.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/C5_fraction_adjusted_and_leave_dataset_out.tsv)
- [01_Normal_Tumor/C5_fraction_paired_stats_updated.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/C5_fraction_paired_stats_updated.tsv)
- [01_Normal_Tumor/C5_independent_TN_statistics.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/C5_independent_TN_statistics.tsv)
- [01_Normal_Tumor/C5_fraction_dataset_stats.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/C5_fraction_dataset_stats.tsv)
- [01_Normal_Tumor/C5_all_cells_dataset_statistics.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/C5_all_cells_dataset_statistics.tsv)
- [02_Clinical_Stage/sample_stage_mapping.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/02_Clinical_Stage/sample_stage_mapping.tsv)
- [02_Clinical_Stage/stage_statistics.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/02_Clinical_Stage/stage_statistics.tsv)
- [02_Clinical_Stage/NMIBC_MIBC_group_summary.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/02_Clinical_Stage/NMIBC_MIBC_group_summary.tsv)
- [02_Clinical_Stage/T_stage_group_summary.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/02_Clinical_Stage/T_stage_group_summary.tsv)
- [02_Clinical_Stage/T_stage_dataset_counts.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/02_Clinical_Stage/T_stage_dataset_counts.tsv)
- [03_Bulk_Survival/input_inventory.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/input_inventory.tsv)
- [03_Bulk_Survival/C5_reference_mapping.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/C5_reference_mapping.tsv)
- [03_Bulk_Survival/Full14_focused_consistency.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/Full14_focused_consistency.tsv)
- [03_Bulk_Survival/C0_C5_KM_comparison.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/C0_C5_KM_comparison.tsv)
- [03_Bulk_Survival/C0_C5_Cox_comparison.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/C0_C5_Cox_comparison.tsv)
- [04_C5_vs_Rest_DEG/DEG_all.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/04_C5_vs_Rest_DEG/DEG_all.tsv)
- [04_C5_vs_Rest_DEG/markers_C5_up.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/04_C5_vs_Rest_DEG/markers_C5_up.tsv)
- [04_C5_vs_Rest_DEG/markers_C5_down.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/04_C5_vs_Rest_DEG/markers_C5_down.tsv)
- [04_C5_vs_Rest_DEG/marker_unit_consistency.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/04_C5_vs_Rest_DEG/marker_unit_consistency.tsv)
- [04_C5_vs_Rest_DEG/marker_unit_contrasts.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/04_C5_vs_Rest_DEG/marker_unit_contrasts.tsv)
- [04_C5_vs_Rest_DEG/pseudobulk/model_summary.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/04_C5_vs_Rest_DEG/pseudobulk/model_summary.tsv)
- [04_C5_vs_Rest_DEG/pseudobulk/C5_vs_Rest_primary.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/04_C5_vs_Rest_DEG/pseudobulk/C5_vs_Rest_primary.tsv)
- [04_C5_vs_Rest_DEG/prespecified_gene_pseudobulk_validation.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/04_C5_vs_Rest_DEG/prespecified_gene_pseudobulk_validation.tsv)
- [04_C5_vs_Rest_DEG/cell_vs_pseudobulk_comparison.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/04_C5_vs_Rest_DEG/cell_vs_pseudobulk_comparison.tsv)
- [05_GO_KEGG/GO_BP_all.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/05_GO_KEGG/GO_BP_all.tsv)
- [05_GO_KEGG/KEGG_all.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/05_GO_KEGG/KEGG_all.tsv)
- [05_GO_KEGG/selected_GO_BP.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/05_GO_KEGG/selected_GO_BP.tsv)
- [05_GO_KEGG/selected_KEGG.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/05_GO_KEGG/selected_KEGG.tsv)
- [05_GO_KEGG/ORA_summary.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/05_GO_KEGG/ORA_summary.tsv)
- [05_GO_KEGG/down/GO_BP_all.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/05_GO_KEGG/down/GO_BP_all.tsv)
- [05_GO_KEGG/down/KEGG_all.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/05_GO_KEGG/down/KEGG_all.tsv)
- [05_GO_KEGG/requested_theme_review_all.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/05_GO_KEGG/requested_theme_review_all.tsv)
- [06_Hallmark_GSEA/GSEA_primary_10pct.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/GSEA_primary_10pct.tsv)
- [06_Hallmark_GSEA/GSEA_all_thresholds.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/GSEA_all_thresholds.tsv)
- [06_Hallmark_GSEA/threshold_direction_comparison.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/threshold_direction_comparison.tsv)
- [06_Hallmark_GSEA/leading_edge.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/leading_edge.tsv)
- [06_Hallmark_GSEA/preranked_5pct.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/preranked_5pct.tsv)
- [06_Hallmark_GSEA/preranked_10pct.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/preranked_10pct.tsv)
- [06_Hallmark_GSEA/preranked_25pct.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/preranked_25pct.tsv)
- [C0_vs_C5_evidence_comparison.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/C0_vs_C5_evidence_comparison.tsv)

各队列各亚群的`cox_PH_tests.tsv`、`cox_covariate_decisions.tsv`、`cox_sample_inclusion.tsv`、`KM_risk_tables.tsv`位于`03_Bulk_Survival/by_dataset/<cohort>/<C0或C5>/`。所有表、图及脚本的绝对路径另见 [output_index.tsv](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/output_index.tsv)。

## 全部图片

| 图 | 矢量PDF | 300dpi PNG |
| --- | --- | --- |
| C5_Normal_Tumor | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/C5_Normal_Tumor.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/C5_Normal_Tumor.png) |
| C5_all_cells_dataset | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/C5_all_cells_dataset.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/C5_all_cells_dataset.png) |
| FigS_C5_fraction_dataset | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/FigS_C5_fraction_dataset.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/FigS_C5_fraction_dataset.png) |
| FigS_C5_fraction_forest | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/FigS_C5_fraction_forest.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/FigS_C5_fraction_forest.png) |
| FigS_C5_fraction_paired_updated | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/FigS_C5_fraction_paired_updated.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/FigS_C5_fraction_paired_updated.png) |
| Normal_Tumor_stacked_pooled_cells | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/Normal_Tumor_stacked_pooled_cells.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/Normal_Tumor_stacked_pooled_cells.png) |
| Normal_Tumor_stacked_sample_equal | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/Normal_Tumor_stacked_sample_equal.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/01_Normal_Tumor/Normal_Tumor_stacked_sample_equal.png) |
| NMIBC_MIBC_barplot | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/02_Clinical_Stage/NMIBC_MIBC_barplot.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/02_Clinical_Stage/NMIBC_MIBC_barplot.png) |
| NMIBC_MIBC_boxplot | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/02_Clinical_Stage/NMIBC_MIBC_boxplot.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/02_Clinical_Stage/NMIBC_MIBC_boxplot.png) |
| T_stage_barplot | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/02_Clinical_Stage/T_stage_barplot.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/02_Clinical_Stage/T_stage_barplot.png) |
| T_stage_boxplot | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/02_Clinical_Stage/T_stage_boxplot.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/02_Clinical_Stage/T_stage_boxplot.png) |
| C5_Cox_forest | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/C5_Cox_forest.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/C5_Cox_forest.png) |
| C0_Cox_forest | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE13507/C0/figures/C0_Cox_forest.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE13507/C0/figures/C0_Cox_forest.png) |
| C0_KM_OS | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE13507/C0/figures/C0_KM_OS.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE13507/C0/figures/C0_KM_OS.png) |
| C5_Cox_forest | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE13507/C5/figures/C5_Cox_forest.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE13507/C5/figures/C5_Cox_forest.png) |
| C5_KM_OS | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE13507/C5/figures/C5_KM_OS.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE13507/C5/figures/C5_KM_OS.png) |
| C0_Cox_forest | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE31684/C0/figures/C0_Cox_forest.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE31684/C0/figures/C0_Cox_forest.png) |
| C0_KM_OS | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE31684/C0/figures/C0_KM_OS.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE31684/C0/figures/C0_KM_OS.png) |
| C5_Cox_forest | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE31684/C5/figures/C5_Cox_forest.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE31684/C5/figures/C5_Cox_forest.png) |
| C5_KM_OS | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE31684/C5/figures/C5_KM_OS.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE31684/C5/figures/C5_KM_OS.png) |
| C0_Cox_forest | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE48075/C0/figures/C0_Cox_forest.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE48075/C0/figures/C0_Cox_forest.png) |
| C0_KM_OS | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE48075/C0/figures/C0_KM_OS.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE48075/C0/figures/C0_KM_OS.png) |
| C5_Cox_forest | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE48075/C5/figures/C5_Cox_forest.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE48075/C5/figures/C5_Cox_forest.png) |
| C5_KM_OS | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE48075/C5/figures/C5_KM_OS.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/GSE48075/C5/figures/C5_KM_OS.png) |
| C0_Cox_forest | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/TCGA_BLCA/C0/figures/C0_Cox_forest.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/TCGA_BLCA/C0/figures/C0_Cox_forest.png) |
| C0_KM_OS | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/TCGA_BLCA/C0/figures/C0_KM_OS.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/TCGA_BLCA/C0/figures/C0_KM_OS.png) |
| C5_Cox_forest | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/TCGA_BLCA/C5/figures/C5_Cox_forest.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/TCGA_BLCA/C5/figures/C5_Cox_forest.png) |
| C5_KM_OS | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/TCGA_BLCA/C5/figures/C5_KM_OS.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/by_dataset/TCGA_BLCA/C5/figures/C5_KM_OS.png) |
| C5_top20_marker_dotplot | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/04_C5_vs_Rest_DEG/C5_top20_marker_dotplot.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/04_C5_vs_Rest_DEG/C5_top20_marker_dotplot.png) |
| C5_vs_Rest_volcano | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/04_C5_vs_Rest_DEG/C5_vs_Rest_volcano.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/04_C5_vs_Rest_DEG/C5_vs_Rest_volcano.png) |
| C5_GO_BP_barplot | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/05_GO_KEGG/C5_GO_BP_barplot.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/05_GO_KEGG/C5_GO_BP_barplot.png) |
| C5_KEGG_barplot | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/05_GO_KEGG/C5_KEGG_barplot.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/05_GO_KEGG/C5_KEGG_barplot.png) |
| GSEA_EPITHELIAL_MESENCHYMAL_TRANSITION | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/figures/GSEA_EPITHELIAL_MESENCHYMAL_TRANSITION.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/figures/GSEA_EPITHELIAL_MESENCHYMAL_TRANSITION.png) |
| GSEA_GLYCOLYSIS | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/figures/GSEA_GLYCOLYSIS.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/figures/GSEA_GLYCOLYSIS.png) |
| GSEA_MYOGENESIS | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/figures/GSEA_MYOGENESIS.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/figures/GSEA_MYOGENESIS.png) |
| GSEA_TNFA_SIGNALING_VIA_NFKB | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/figures/GSEA_TNFA_SIGNALING_VIA_NFKB.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/figures/GSEA_TNFA_SIGNALING_VIA_NFKB.png) |
| Hallmark_GSEA_summary | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/figures/Hallmark_GSEA_summary.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/figures/Hallmark_GSEA_summary.png) |
| NES_threshold_comparison | [PDF](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/figures/NES_threshold_comparison.pdf) | [PNG](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/06_Hallmark_GSEA/figures/NES_threshold_comparison.png) |
