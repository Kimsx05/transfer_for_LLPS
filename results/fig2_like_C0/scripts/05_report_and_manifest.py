"""Build a compact, data-linked Phase 1 report and reproducibility manifests."""
from pathlib import Path
import csv
import hashlib
import json

BASE = Path('/home/jinshengxi/bladder_metabolism')
OUT = BASE / 'results/fig2_like_C0'
CORE = OUT / 'phase1_core'
SUMMARY = OUT / 'summary'
ATLAS = Path('/home/jinshengxi/LLPS/2.1 LLPS atlas/LLPS_77_genes')
FIBROOT = ATLAS / 'fibro_pericyte_subcluster_clean_v4_global_noC5/global_resolution_0.12_numeric_order'
OLD = BASE / 'cellchat_nichechat_fibro7/11_C0_dual_geneset_preparation'

def read(name):
    with (CORE / name).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def link(path, label=None):
    return f'[{label or path.name}](<{path}>)'

def table(rows, columns):
    return '\n'.join(['| '+' | '.join(columns)+' |', '| '+' | '.join(['---']*len(columns))+' |']+
                     ['| '+' | '.join(str(row[k]) for k in columns)+' |' for row in rows])

def checksum(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(4*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

def write_tsv(path, rows, columns):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(rows)

stats=read('C0_fraction_statistics.tsv')
gsea=read('GSEA_pathway_plot_selection.tsv')
figs=sorted(CORE.glob('Fig2*.pdf'))
assert len(figs)==12
assert all(x['pass']=='TRUE' for x in read('validation_checks.tsv'))

report = f'''# C0 Fig.2-like：Phase 1 结果汇报

完成日期：2026-10-25（Asia/Shanghai）。**Phase 1 已完成；按照“先汇报，再决定 Phase 2”的要求，外部 bulk 分期/生存尚未启动。**

输出目录：{link(CORE, 'phase1_core')}。共 **12 张独立图，每张 PDF + 350 dpi PNG**，没有拼图。脚本位于 {link(OUT/'scripts', 'scripts')}；详细输入、统计方法和限制见 {link(SUMMARY/'analysis_notes.md', 'analysis_notes.md')}。

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
'''
for f in figs:
    report += f'| {f.stem} | {link(f,"PDF")} | {link(f.with_suffix(".png"),"PNG")} |\n'
report += '\n## GSEA 关注通路\n\n'
report += table(gsea,['pathway','NES','padj','exported'])
report += '\n\n## 核心数据表\n\n'
for name in ['cluster_LLPSscore_summary.tsv','C0_fraction_per_sample.tsv','C0_fraction_statistics.tsv',
             'C0_vs_others_DEG.tsv','C0_upregulated_genes.tsv','C0_GO_enrichment.tsv',
             'C0_KEGG_enrichment.tsv','C0_GSEA_results.tsv','C0_preranked_gene_list.tsv',
             'pseudobulk_model_summary.tsv','C0_DEG_sensitivity_comparison.tsv','validation_checks.tsv']:
    report += '- '+link(CORE/name)+'\n'
report += '''
## 需要人工核对的事项

1. **HRA003620 患者配对身份**：既有映射仅含候选配对编号，patient_uid 未核实。本轮没有根据相似名字构造配对；主分析暂以 sample/specimen 为单位。因此 rank-sum 的独立性仍是工作假设，需要来源方核对；8 对补充分析仅使用已核实 study-patient IDs。
2. **临床含义**：Normal 包括邻近正常组织；组织来源和数据集组成不平衡。比例下降是该冻结子集内的结果，不能直接解释为肿瘤整体绝对细胞数下降。
3. **颜色已解决，无需重新分配**：复用后续 publication-style 图明确保存的 SCP Paired mapping；更早的探索图曾用 Dark 3，本轮没有采用旧的探索版颜色。
4. **富集不等于机制验证**：GO 血管相关条目显著，而 Hallmark ANGIOGENESIS 不显著，两者必须分别汇报。EMT signature 在既有间质细胞中的富集不能证明它们发生上皮到间质转化，亦不能证明 LLPS 导致这些功能。
5. **Phase 2 尚未执行**：尚无外部 stage、HR、KM 或 abundance 结论；C0-vs-stromal markers 还需检查相对其他主要谱系的特异性后再用于 bulk。若采用 expression-based score，必须明确区别于严格去卷积丰度。
'''
(SUMMARY/'README.md').write_text(report)

notes=f'''# 输入、统计方法与可复现说明

## 冻结输入

- 研究对象：{link(FIBROOT/'global_res012_object.qs')}，20,146细胞、51,959 RNA features；`subcluster` 固定为 C0–C6，与 `clean_snn_res.0.12` 和既有冻结条码表一致。未重新聚类、运行UMAP、筛细胞或计算LLPS score。
- 标签及分数独立对照：{link(FIBROOT/'global_res012_cell_labels_scores.csv')}，数值序列化容差1e-14。
- 颜色：{link(FIBROOT/'beautified_SCP_protein_coding_top3/cluster_SCP_Paired_colors.csv')}；C0 #A6CEE3，C1 #1F78B4，C2 #B2DF8A，C3 #33A02C，C4 #FDBF6F，C5 #FF7F00，C6 #FB9A99。已写死在 common.R 并逐项验证。
- score：继承对象 `LLPS77_AUC`；使用 {link(FIBROOT.parent/'clean_v4_settings.rds')} 的原全图色标范围，plasma；精确零值为黑色（本子集零值数量见 LLPS_color_scale.tsv）。
- sample 组织/患者信息：{link(OLD/'01_C0_abundance/A1_sample_abundance_manifest.tsv')}。本轮只提取来源、患者关系与样本字段并从当前细胞重算计数，没有复用先前 epithelial DEG 或其他分母。来源溯源见 {link(OLD/'00_preflight/sample_source_mapping.tsv')}。
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
'''
(SUMMARY/'analysis_notes.md').write_text(notes)

resources = [
 (FIBROOT/'global_res012_object.qs','frozen Seurat object'),
 (FIBROOT/'global_res012_cell_labels_scores.csv','frozen barcode labels and LLPS scores'),
 (FIBROOT/'beautified_SCP_protein_coding_top3/cluster_SCP_Paired_colors.csv','saved publication palette'),
 (FIBROOT.parent/'clean_v4_settings.rds','inherited atlas-wide color limits'),
 (OLD/'01_C0_abundance/A1_sample_abundance_manifest.tsv','existing source/relationship manifest'),
 (OUT/'resources/KEGG_hsa_gene_pathway.tsv','https://rest.kegg.jp/link/pathway/hsa'),
 (OUT/'resources/KEGG_hsa_pathway_names.tsv','https://rest.kegg.jp/list/pathway/hsa'),
 (OUT/'resources/h.all.v2025.1.Hs.symbols.gmt',str(OLD/'00_preflight/h.all.v2025.1.Hs.symbols.gmt')),
]
rows=[dict(path=str(p),source=source,sha256=checksum(p),bytes=p.stat().st_size,date='2026-10-25') for p,source in resources]
write_tsv(SUMMARY/'input_resource_manifest.tsv',rows,['path','source','sha256','bytes','date'])
historical=json.loads((OLD/'00_preflight/protected_input_snapshot_before.json').read_text())
index={x['path']:x for x in historical}
comparisons=[]
for row in rows:
    if row['path'] in index:
        comparisons.append(dict(path=row['path'],matches_historical_sha256=row['sha256']==index[row['path']]['sha256']))
assert all(r['matches_historical_sha256'] for r in comparisons)
write_tsv(SUMMARY/'historical_input_preservation.tsv',comparisons,['path','matches_historical_sha256'])

artifacts=[dict(path=str(p),bytes=p.stat().st_size,sha256=checksum(p)) for p in sorted(OUT.rglob('*'))
           if p.is_file() and p.name!='artifact_manifest.tsv']
write_tsv(SUMMARY/'artifact_manifest.tsv',artifacts,['path','bytes','sha256'])
print(f'Phase 1 complete: {len(figs)} independent PDF/PNG pairs. Reports and manifests written.')
