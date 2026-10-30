from pathlib import Path
import csv
B=Path('/home/jinshengxi/bladder_metabolism');O=B/'results/llps_C5_parallel_v1'
def read(p):
    return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rows):
    with p.open('w') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def link(p,label=None):
    p=Path(p);return f'[{label or p.name}]({p})'
def fmt(v,d=3):
    try:return f'{float(v):.{d}g}'
    except:return str(v)
def table(rows,cols):
    return '| '+' | '.join(cols)+' |\n| '+' | '.join('---' for _ in cols)+' |\n'+'\n'.join('| '+' | '.join(str(r.get(c,'')) for c in cols)+' |' for r in rows)
figs=read(O/'figure_QC.tsv')
cx=read(O/'03_Bulk_Survival/C0_C5_Cox_comparison.tsv');km=read(O/'03_Bulk_Survival/C0_C5_KM_comparison.tsv')
comparison=[
 {'证据维度':'细胞与来源覆盖','C0':'10,717 cells；43 个有检出样本','C5':'716 cells；37 个有检出样本；InhouseData 382 (53.4%)，GSE222315 283 (39.5%)','解释与限制':'C5 92.9% 来自两个 dataset；有检出样本数不等于丰度分析总样本数。'},
 {'证据维度':'Normal/Tumor 丰度','C0':'dataset 调整后 Tumor−Normal −52.20 个百分点；P=3.06e−7；13 对配对 P=0.00253','C5':'调整后 +11.61 个百分点；P=0.0378；13 对配对 P=0.0528；min50 P=0.235','解释与限制':'C5 更偏向肿瘤富集，但对小分母样本及来源敏感；C0 的下降不能解释为生物学无意义。'},
 {'证据维度':'临床分期','C0':'未发现可直接复用的同口径单细胞 C0 分期检验结果；既有 Bulk signature 分期结果不替代本项','C5':'NMIBC/MIBC：17/20 样本，P=0.0643；T 分组总体 P=0.0628；Ta/T4 仅3/4样本','解释与限制':'尚不足以断言随临床分期单调升高；T 分类不等同于总体病理分期。'},
 {'证据维度':'Full14 连续丰度与 OS','C0':'同四队列中仅 GSE13507 名义支持：调整 HR=1.35，P=0.0245；跨队列 BH q=0.0981','C5':'TCGA/GSE48075 未调整 HR=1.19/1.27；TCGA 跨队列 BH q=0.0144；所有调整模型 P>0.05','解释与限制':'C5 有未调整关联线索，尚无稳定独立预后证据；C0/C5 使用相同队列、OS、SD 和中位数规则。'},
 {'证据维度':'marker 与样本层验证','C0':'沿用 C0 v2 Wilcoxon marker 流程作为方法基准','C5':'POSTN、COL11A1、COL10A1、MMP11 等；POSTN/COL11A1 获13个样本 pseudobulk 支持；COL10A1 被表达过滤','解释与限制':'560个可检验 C5-up 中501个方向一致，293个 pseudobulk BH FDR<0.05且正向；样本来自两个 dataset。'},
 {'证据维度':'GO/KEGG','C0':'既有 v3 显著代表通路覆盖 ECM、迁移/黏附、血管、炎症与脂质调节','C5':'ECM组织、胶原代谢/降解/合成、黏附与血管发育；KEGG ECM-receptor/integrin 等','解释与限制':'C5 的基质重塑证据集中且有表达支持；不存在显著通路的主题不强行补齐。'},
 {'证据维度':'Hallmark 主分析及稳定性','C0':'EMT NES=1.65；TNFα/NF-κB NES=2.75；炎症正向，OXPHOS负向','C5':'EMT NES=1.93；TNFα/NF-κB NES=−2.38；两者三阈值均显著；糖酵解正向但25% FDR=0.0545','解释与限制':'提示不同功能状态；两者各自 Rest 不同，NES差异不是直接 C0 vs C5 检验，也不证明 EMT 转分化或 LLPS 机制。'},
 {'证据维度':'既有 LLPS 与 HR signature（仅引用）','C0':'LLPS细胞中位数0.1413；HR细胞中位数0.2140；患者等权HR中位数0.1959','C5':'LLPS细胞中位数0.1260；HR细胞中位数0.2429；患者等权HR中位数0.1974','解释与限制':'逐细胞 HR 优势在患者等权后明显缩小；未重算分数，未用TCGA筛出的694基因建立新模型。'},
 {'证据维度':'去卷积可识别性','C0':'Full14 和 focused 均有独立C0；TCGA/UROMOL 两方案 C0 Spearman=0.787/0.759','C5':'仅Full14有独立C5；150参考细胞来自GSE222315的7个patient_key','解释与限制':'两方案一致性是C0信息，不能验证C5准确性；参考来源单一限制C5的跨队列解释。'}]
write(O/'C0_vs_C5_evidence_comparison.tsv',comparison)
(O/'C0_vs_C5_evidence_comparison.md').write_text('# C0 vs C5 证据对比\n\n'+table(comparison,list(comparison[0]))+'\n')
bulk=[]
for co in ['TCGA_BLCA','GSE31684','GSE13507','GSE48075']:
    u=next(r for r in cx if r['dataset']==co and r['cluster']=='C5' and r['model']=='Univariable')
    a=next(r for r in cx if r['dataset']==co and r['cluster']=='C5' and r['model']=='Multivariable')
    k=next(r for r in km if r['dataset']==co and r['cluster']=='C5')
    bulk.append({'队列':'TCGA' if co=='TCGA_BLCA' else co,'OS N/事件':f"{u['N']}/{u['Event_N']}",'未调整 HR (95%CI)':f"{float(u['HR']):.2f} ({float(u['lower95']):.2f}–{float(u['upper95']):.2f})",'P / 跨队列BH q':f"{fmt(u['p_value'])} / {fmt(u['BH_across_cohorts'])}",'调整 HR (95%CI), P':f"{float(a['HR']):.2f} ({float(a['lower95']):.2f}–{float(a['upper95']):.2f}), {fmt(a['p_value'])}",'KM log-rank P':fmt(k['logrank_p'])})
gs=read(O/'06_Hallmark_GSEA/GSEA_primary_10pct.tsv');sens={r['pathway']:r for r in read(O/'06_Hallmark_GSEA/threshold_direction_comparison.tsv')}
grows=[]
for r in gs:
    if float(r['padj'])<.05:
        q=sens[r['pathway']];grows.append({'Hallmark':r['pathway'].removeprefix('HALLMARK_'),'主分析NES':fmt(r['NES']),'5% FDR':fmt(q['padj_0.05']),'10% FDR':fmt(q['padj_0.1']),'25% FDR':fmt(q['padj_0.25'])})
summary=f'''# C5 六模块平行分析结果

完成日期：2026-10-30。六模块均已执行；未重新聚类、注释、打分或去卷积。全部新结果位于 {link(O,'llps_C5_parallel_v1')}。

**主要结论：C5 最有说服力的证据是 ECM/胶原重塑与 EMT 基因集富集。肿瘤富集有提示，但对样本分母与 dataset 敏感；目前不足以认定为稳定的独立预后亚群。**

## 1. Normal vs Tumor

使用复核后的患者映射，主要模型包含 **53 个有效样本（Normal 16、Tumor 37），40 个 patient_key**；其中 Normal 10 个、Tumor 7 个样本 C5=0，均保留。样本丰度中位数分别为 0 和 5.88%。

dataset 调整、患者聚类标准误模型的 Tumor−Normal 差为 **+11.61 个百分点**（95%CI 0.69–22.53，P=0.0378）。13 对已验证配对中8对升高、2对降低、3对相同，中位配对差+2.84个百分点，配对Wilcoxon P=0.0528。剔除已知配对后，独立比较仅有3个Normal及24个Tumor，P=0.0329，不能忽略Normal样本过少。

稳定性有限：要求每样本至少50个Fibro/Pericyte后，调整效应降为+2.61个百分点（P=0.235）；去掉HRA003620后为+2.63个百分点（P=0.103）。HRA003620仅贡献7个C5细胞，但一些样本分母很小。因此应表述为“**提示Tumor偏高，尚不稳健**”。另已用现有精确样本映射补充 C5/全部atlas细胞结果，二者分母分开报告。

## 2. 临床分期

37个有效肿瘤样本均有本轮所用分组记录。NMIBC 17例、MIBC 20例，C5比例中位数分别4.11%、10.43%，Wilcoxon P=0.0643。Ta/T1/T2/T3/T4分别 **3/14/11/5/4例**；总体Kruskal–Wallis P=0.0628。未经多重校正的探索性Spearman rho=0.332、P=0.0449，不能据此断言单调进展关系：Ta组仅3例，且其中位数高于T1/T2，T4也仅4例。没有进行稀疏分期的成对显著性推断。

T来自既有工作簿“TNM stage”中明确的T成分；记录中同时有pT和未注明临床/病理属性的T，本轮保留其来源类型。**总体病理AJCC分期未被替换成T分期。** 来源构成和分期构成仍可能混杂。

## 3. Full14 CIBERSORTx estimated abundance 与 OS

已核对Full14参考细胞：150个`Fibro_C5`参考条形码均对应冻结C5；均来自GSE222315的7个patient_key。7个现有Bulk队列均有有效独立C5列，其中4个有OS。UROMOL仅PFS、GSE32894仅DSS、GSE87304无可用生存端点，未冒充OS。

以下HR均按**每1个队列总体SD**增加。C0/C5采用相同原始Full14样本、OS定义和标准化方式。

{table(bulk,list(bulk[0]))}

TCGA和GSE48075支持未调整的风险升高，其中仅TCGA在四队列BH校正后仍显著。**四队列的调整模型均未支持独立C5关联**。TCGA调整分期、年龄、性别及原流程允许的grade，N=399/175事件；其他队列按原有可用字段政策调整，见完整模型表。TCGA的C5项PH检验P=0.498/0.328（未调整/调整），调整模型全局PH P=0.268；GSE13507调整模型全局PH P=0.00465，存在PH偏离，C5项P=0.160，需谨慎解释模型。

KM固定采用**全部去卷积样本的队列中位数**，High > median、Low ≤ median，不优化cutoff。GSE13507有61.8%零值，中位数=0，分组为Low 102/High 63，等同zero/nonzero，不能称作平衡的上下半组。GSE48075 cutoff来自142个原始样本，只有73个有效OS样本，因此OS分组为18/55；这是沿用C0的固定规则。C0-focused没有独立C5。使用现存结果补充的C0两方案相关性不能替代C5准确性验证。

## 4. DEG与样本层验证

冻结RNA/data：20,146 cells，C5 716、Rest 19,430；完整表 **51,959基因**。完全沿用C0 v2的Seurat Wilcoxon；marker标准为max检测率≥25%、|avg_log2FC|>0.75、**Bonferroni-adjusted P<0.05**，得到 **564个上调、405个下调**。完整表另外提供`FDR_BH`，仅供透明报告，未用它替换原筛选标准。极小P在原软件中数值下溢为0，保留原值；火山图仅对显示的−log10作数值下限保护。

主要上调包括 **COL10A1 (avg_log2FC=6.51)、COL11A1 (6.24)、POSTN (3.59)、MMP11 (4.56)、PLPP4 (5.39)、GJB2 (4.50)**。POSTN/COL11A1/COL10A1在5个dataset的均值方向均为正，但患者/样本层面并非普遍检出：37个同时含两组的样本中，POSTN有29个正向，COL11A1和COL10A1各17个正向；小样本捕获与零表达不能忽略。这是C5 vs Rest的相对表达，不证明全atlas范围专一性。

复用原Fig2 edgeR pseudobulk代码：每臂≥20细胞、已知患者去重，**13个患者/样本单位、26个文库，全部Tumor、仅两个dataset**；7,377基因可检验。560个通过pseudobulk过滤的C5-up中501个方向一致（89.5%），293个正向且BH FDR<0.05。POSTN：pseudobulk logFC=2.61、FDR=2.34e−5；COL11A1：logFC=4.37、FDR=5.08e−5。**COL10A1未通过原有filterByExpr，状态为未检验，不解释为验证失败或无效。** ≥50细胞/臂仅剩3个样本，按原规则不拟合；Tumor-only与主模型使用相同13个样本，不能当成独立复现。

## 5. GO BP / KEGG ORA

重新运行C0 v2实际富集代码，没有把C0 v3固定通路重绘当作C5富集。GO有效背景4,210、前景518个Entrez；KEGG背景2,402、前景272个Entrez；同一检测率定义和项目KEGG快照。GO有70条显著通路，展示10条显著代表通路；**KEGG仅4条显著，全部展示，不补齐10条**。下调ORA完整表亦已保存。

GO最突出的是ECM组织（Count=61，FDR=5.63e−17）、胶原代谢（27，3.99e−9）、细胞黏附（97，1.04e−5）、血管发育（2.25e−5）及血管生成（0.00613），并包括胶原降解/合成、结缔组织发育、细胞骨架及Wnt调节。未据此宣称炎症或脂质/胆固醇过程在C5-up中显著富集；完整主题检索表保留阴性结果。

KEGG为Protein digestion and absorption、Cytoskeleton in muscle cells、ECM-receptor interaction、Integrin signaling。前一条的前景以胶原基因为主，应结合基因构成解读为基质相关信号，不能直接称为消化功能增强。富集不等于通路活性或因果机制证明。

## 6. Hallmark GSEA

C0 v3的实际实现与指定参数一致：MSigDB 2026.1.Hs，fgseaMultilevel，minSize=10/maxSize=500，eps=0，nPermSimple=10000，seed=123，各阈值内BH。使用**新C5全DEG**，按avg_log2FC降序，5%/10%/25%分别11,262/8,855/4,481个基因；分别测试50/50/49个Hallmark。25%时尺寸不足的集合明确标为NOT_TESTED，未伪造NES。

{table(grows,list(grows[0]))}

**EMT正向与TNFα/NF-κB负向最稳定，三阈值均显著。** 糖酵解、myogenesis在25%时FDR略高于0.05；MYC_TARGETS_V2仅主阈值显著。TGF-β、angiogenesis和hypoxia等Hallmark没有达到主分析FDR<0.05；不把GO血管相关富集外推成Hallmark也显著。按C0 v3客观展示规则（最多4条正向+2条负向），本轮绘制实际满足条件的2正+2负曲线，所有50条均在NES总览及完整表中。

## C0 vs C5证据对比

{table(comparison,list(comparison[0]))}

C5的潜在优势是肿瘤基质重塑相关表达与稳定EMT富集，而非“P值更小”。C0保有更高的既有LLPS分数、不同的炎症状态以及不同的丰度/生存关联。C5尚不能被认定为优于C0的独立预后标志物，更不能由这些结果证明LLPS因果机制。

主要文件：{link(O/'README.md','复用方法、输入和全部图表索引')}；{link(O/'03_Bulk_Survival/C0_C5_Cox_comparison.tsv','同队列C0/C5 Cox对照')}；{link(O/'C0_vs_C5_evidence_comparison.tsv','证据对比表TSV')}。
'''
(O/'result_summary.md').write_text(summary)
reuse=[
 {'模块':'1 Normal/Tumor','原始实际代码':'results/fig2_like_C0/phase1B_review/scripts/common.R、02_sample_statistics.R；scripts/01_core_plots_and_pseudobulk.R','本轮代码':'scripts/01_composition.R','必要改动':'C0→C5；采用reviewed配对；新输出路径；仅截取丰度区段；补100%堆叠图和全atlas分母；300dpi透明背景'},
 {'模块':'2 分期','原始实际代码':'results/fig2_like_C0/music_optimization/scripts/06_clinical.R；既有工作簿解析结果及样本映射','本轮代码':'scripts/02_stage.R','必要改动':'复用stage色表、box+jitter、KW/Spearman/Wilcoxon逻辑；将已映射工作簿T成分与C5比例连接；加均值柱图；稀疏组不做成对推断'},
 {'模块':'3 Bulk OS','原始实际代码':'results/cibersortx_phase2_clinical_v1/scripts/01_integrate.py、02_analyze.R','本轮代码':'scripts/03_bulk_inputs.py、03_bulk_survival.R、03b_reference_comparison.R','必要改动':'直接复用临床整合表，不重跑整合下载；动态提取fit_cox、choose_covars、KM和forest代码；目标C0/C5；仅OS；增加PNG和横向比较表'},
 {'模块':'4 DEG','原始实际代码':'results/llps_C0_functional_signature_v2/scripts/01_markers.R、common.R；Fig2 scripts/02_pseudobulk_DEG.R；llps_C0_functional_state_v1/scripts/04_DEG_LLPS.R火山图段','本轮代码':'scripts/04_markers.R、04b_pseudobulk.R、04c_volcano.R','必要改动':'C5 vs setdiff(C0–C6,C5)；焦点基因换POSTN/COL11A1/COL10A1；对称down标准；原pseudobulk函数及过滤；不足4样本记录NOT_FIT；火山图沿用几何但按v2阈值与Bonferroni轴名'},
 {'模块':'5 ORA','原始实际代码':'results/llps_C0_functional_signature_v2/scripts/02_ORA.R；v3/scripts/01_curated_ORA.R绘图表达式','本轮代码':'scripts/05_ORA.R、05b_ORA_figures.R','必要改动':'新C5前景及背景；up/down输出；保留实际enrichGO/enrichKEGG/简化；用新显著ID替换C0固定10条；10条GO、实际4条KEGG'},
 {'模块':'6 GSEA','原始实际代码':'results/llps_C0_functional_signature_v3/scripts/02_GSEA.R、common.R','本轮代码':'scripts/06_GSEA.R','必要改动':'指向C5排名；保留C0 v2同版Hallmark资源；去掉不适用的C0 v2 NES对照；曲线由C5真实结果选择；阈值/重复symbol/fgsea/布局不变'},
 {'模块':'资源与图检','原始实际代码':'results/llps_C0_functional_signature_v3/scripts/run_limited.sh、08_figure_QC.py','本轮代码':'scripts/run_limited.sh、08_figure_QC.py','必要改动':'运行限制原样复用；图检仅去掉C0专用的固定17张断言，动态列出本轮38张'}]
write(O/'script_reuse.tsv',reuse)
inputs=[
 ('冻结Fibro/Pericyte对象','/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes/fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order/global_res012_object.qs'),
 ('冻结细胞标签/已有分数',str(B/'results/fig2_like_C0/phase1_core/cell_metadata_umap_LLPSscore.tsv.gz')),
 ('复核样本/患者/配对',str(B/'results/fig2_like_C0/phase1B_review/sample_metadata_reviewed.tsv')),
 ('精确全atlas分母及工作簿行映射',str(B/'cellchat_nichechat_fibro7/11_C0_dual_geneset_preparation/01_C0_abundance/A1_sample_abundance_manifest.tsv')),
 ('已有工作簿临床行（F列为TNM）',str(B/'cellchat_nichechat_fibro7/11_C0_dual_geneset_preparation/00_preflight/user_workbook_relevant_scRNA_rows.tsv')),
 ('Full14/Focused参考细胞映射',str(B/'fig2_C0_multimethod_deconvolution/03_methods/CIBERSORTx_web/manifest/reference_cells_LOCAL_ONLY.tsv')),
 ('Full14原始输出与队列清单',str(B/'results/cibersortx_phase2_clinical_v1/input_inventory.tsv')),
 ('已整理Bulk表达和临床位置',str(B/'fig2_C0_multimethod_deconvolution/01_bulk_cohorts')),
 ('原OS定义及来源',str(B/'results/cibersortx_phase2_clinical_v1/endpoint_definitions.tsv')),
 ('Hallmark同版资源',str(B/'results/llps_C0_functional_signature_v2/resources/gene_sets.rds')),
 ('Hallmark版本记录',str(B/'results/llps_C0_functional_signature_v2/resources/gene_set_sources.tsv')),
 ('KEGG人类快照',str(B/'results/fig2_like_C0/resources/KEGG_hsa_gene_pathway.tsv')),
 ('KEGG通路名称快照',str(B/'results/fig2_like_C0/resources/KEGG_hsa_pathway_names.tsv')),
 ('gene symbol映射',str(B/'fig2_C0_multimethod_deconvolution/03_methods/CIBERSORTx_web/resources/gencode_v36_gene_map.tsv'))]
write(O/'input_sources.tsv',[{'用途':x,'路径':p} for x,p in inputs])
main_tables=[
 '01_Normal_Tumor/C5_fraction_per_sample_reviewed.tsv','01_Normal_Tumor/all_cluster_sample_composition.tsv','01_Normal_Tumor/C5_fraction_adjusted_and_leave_dataset_out.tsv','01_Normal_Tumor/C5_fraction_paired_stats_updated.tsv','01_Normal_Tumor/C5_independent_TN_statistics.tsv','01_Normal_Tumor/C5_fraction_dataset_stats.tsv','01_Normal_Tumor/C5_all_cells_dataset_statistics.tsv',
 '02_Clinical_Stage/sample_stage_mapping.tsv','02_Clinical_Stage/stage_statistics.tsv','02_Clinical_Stage/NMIBC_MIBC_group_summary.tsv','02_Clinical_Stage/T_stage_group_summary.tsv','02_Clinical_Stage/T_stage_dataset_counts.tsv',
 '03_Bulk_Survival/input_inventory.tsv','03_Bulk_Survival/C5_reference_mapping.tsv','03_Bulk_Survival/Full14_focused_consistency.tsv','03_Bulk_Survival/C0_C5_KM_comparison.tsv','03_Bulk_Survival/C0_C5_Cox_comparison.tsv',
 '04_C5_vs_Rest_DEG/DEG_all.tsv','04_C5_vs_Rest_DEG/markers_C5_up.tsv','04_C5_vs_Rest_DEG/markers_C5_down.tsv','04_C5_vs_Rest_DEG/marker_unit_consistency.tsv','04_C5_vs_Rest_DEG/marker_unit_contrasts.tsv','04_C5_vs_Rest_DEG/pseudobulk/model_summary.tsv','04_C5_vs_Rest_DEG/pseudobulk/C5_vs_Rest_primary.tsv','04_C5_vs_Rest_DEG/prespecified_gene_pseudobulk_validation.tsv','04_C5_vs_Rest_DEG/cell_vs_pseudobulk_comparison.tsv',
 '05_GO_KEGG/GO_BP_all.tsv','05_GO_KEGG/KEGG_all.tsv','05_GO_KEGG/selected_GO_BP.tsv','05_GO_KEGG/selected_KEGG.tsv','05_GO_KEGG/ORA_summary.tsv','05_GO_KEGG/down/GO_BP_all.tsv','05_GO_KEGG/down/KEGG_all.tsv','05_GO_KEGG/requested_theme_review_all.tsv',
 '06_Hallmark_GSEA/GSEA_primary_10pct.tsv','06_Hallmark_GSEA/GSEA_all_thresholds.tsv','06_Hallmark_GSEA/threshold_direction_comparison.tsv','06_Hallmark_GSEA/leading_edge.tsv','06_Hallmark_GSEA/preranked_5pct.tsv','06_Hallmark_GSEA/preranked_10pct.tsv','06_Hallmark_GSEA/preranked_25pct.tsv','C0_vs_C5_evidence_comparison.tsv']
for p in main_tables: assert (O/p).exists(),p
reuse_md=table([{'模块':r['模块'],'C0源代码（相对项目根）':r['原始实际代码'],'本轮脚本':r['本轮代码'],'修改':r['必要改动']} for r in reuse],['模块','C0源代码（相对项目根）','本轮脚本','修改'])
readme=f'''# C5 parallel v1

完成日期：2026-10-30。六模块均已执行。结果解释见 {link(O/'result_summary.md')}，证据对比见 {link(O/'C0_vs_C5_evidence_comparison.md')}。

## 复用代码与修改范围

项目根：`{B}`。所有C0源代码和结果只读使用；新代码及输出仅写入本目录。没有调用新聚类、注释、LLPS/HR打分、UMAP/Harmony或CIBERSORTx。没有下载Bulk或基因集，没有创建环境冻结或文件哈希。

{reuse_md}

原脚本均先阅读后复用。动态提取代码片段的脚本读取原始R文件并在当前新输出路径执行，未编辑原文件。模块2没有找到同口径的单细胞C0分期图，因而复用原Fig2临床图配色和统计结构，并使用已有工作簿映射生成C5单细胞丰度分析，未将Bulk signature分期结论直接移植。

实际参数差异：C0 Fig2早期common的seed为20261025，本轮统一123；早期Fig2图为白底350dpi，按本轮明确要求改为透明300dpi。C0 v2 marker的`p_val_adj`是**Bonferroni**，并非BH FDR，本轮保留原判定；另存`FDR_BH`但不用于改变marker筛选。C0 v3 GSEA参数与任务描述完全一致。

## 输入及分析对象

'''+ '\n'.join(f'- {label}：{link(path)}' for label,path in inputs)+f'''

Bulk四队列OS及协变量直接读取原 `results/cibersortx_phase2_clinical_v1/by_dataset/<cohort>/clinical_merged.tsv`，其具体原始去卷积文件、样本数、临床来源与映射核对记录在 {link(O/'03_Bulk_Survival/input_inventory.tsv')}。Full14的7个队列及现有focused输出列清单见 {link(O/'03_Bulk_Survival/existing_deconvolution_inventory.tsv')}。没有找到单独保存的既有Full14/focused最终一致性统计表；因此仅对已有TCGA和UROMOL输出新增描述性相关性核对，未重新运行去卷积。

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

{link(O/'validation.txt')}：核对716细胞、分母/零值、完整DEG及up/down规则、选定ORA统计精确对应、GSEA阈值与BH；本轮重算的8个C0 Cox HR/N与既有结果一致。{link(O/'figure_QC.tsv')}：**{len(figs)}对独立矢量PDF/透明300dpi PNG**；所有PDF单页且无栅格嵌入，PNG alpha正常。已目视核对代表性组成、分期、KM、火山、ORA及GSEA图，必要的零点显示修正已完成。

## 主要数据表

'''+ '\n'.join('- '+link(O/p,p) for p in main_tables)+f'''

各队列各亚群的`cox_PH_tests.tsv`、`cox_covariate_decisions.tsv`、`cox_sample_inclusion.tsv`、`KM_risk_tables.tsv`位于`03_Bulk_Survival/by_dataset/<cohort>/<C0或C5>/`。所有表、图及脚本的绝对路径另见 {link(O/'output_index.tsv')}。

## 全部图片

| 图 | 矢量PDF | 300dpi PNG |
| --- | --- | --- |
'''+ '\n'.join(f"| {Path(r['path']).stem} | {link(O/r['path'],'PDF')} | {link((O/r['path']).with_suffix('.png'),'PNG')} |" for r in figs)+'\n'
(O/'README.md').write_text(readme)
status=[{'module':str(i),'status':'COMPLETE','limitation':s} for i,s in enumerate(['Sparse denominators and dataset sensitivity','Sparse Ta/T4; exploratory unadjusted trend','Four OS cohorts; focused C5 unavailable; adjusted associations not significant','min50 pseudobulk NOT_FIT (3 samples); COL10A1 not tested by pseudobulk filter','4 significant KEGG up pathways only','25% threshold: one set not tested by size'],1)]
write(O/'module_status.tsv',status)
rows=[]
for p in sorted(O.rglob('*')):
    if p.is_file() and p.name!='output_index.tsv': rows.append({'path':str(p),'relative_path':str(p.relative_to(O)),'type':p.suffix,'bytes':p.stat().st_size})
write(O/'output_index.tsv',rows)
print(f'Written Chinese summary, README, evidence comparison and index ({len(figs)} figures).')
