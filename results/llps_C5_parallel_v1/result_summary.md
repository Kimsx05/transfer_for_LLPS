# C5 六模块平行分析结果

完成日期：2026-10-30。六模块均已执行；未重新聚类、注释、打分或去卷积。全部新结果位于 [llps_C5_parallel_v1](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1)。

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

| 队列 | OS N/事件 | 未调整 HR (95%CI) | P / 跨队列BH q | 调整 HR (95%CI), P | KM log-rank P |
| --- | --- | --- | --- | --- | --- |
| TCGA | 403/176 | 1.19 (1.06–1.35) | 0.0036 / 0.0144 | 1.10 (0.97–1.25), 0.151 | 0.013 |
| GSE31684 | 93/65 | 1.21 (0.96–1.52) | 0.112 / 0.149 | 0.98 (0.74–1.29), 0.874 | 0.121 |
| GSE13507 | 165/69 | 1.15 (0.92–1.43) | 0.231 / 0.231 | 0.87 (0.67–1.12), 0.273 | 0.134 |
| GSE48075 | 73/45 | 1.27 (1.02–1.59) | 0.0342 / 0.0683 | 1.22 (0.95–1.55), 0.114 | 0.315 |

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

| Hallmark | 主分析NES | 5% FDR | 10% FDR | 25% FDR |
| --- | --- | --- | --- | --- |
| EPITHELIAL_MESENCHYMAL_TRANSITION | 1.93 | 2.05e-06 | 9.47e-06 | 2.72e-06 |
| GLYCOLYSIS | 1.62 | 0.00453 | 0.0142 | 0.0545 |
| MYC_TARGETS_V2 | -1.75 | 0.0637 | 0.0448 | 0.173 |
| MYOGENESIS | -1.79 | 0.0284 | 0.0048 | 0.0545 |
| TNFA_SIGNALING_VIA_NFKB | -2.38 | 2.05e-06 | 2.18e-09 | 3.16e-09 |

**EMT正向与TNFα/NF-κB负向最稳定，三阈值均显著。** 糖酵解、myogenesis在25%时FDR略高于0.05；MYC_TARGETS_V2仅主阈值显著。TGF-β、angiogenesis和hypoxia等Hallmark没有达到主分析FDR<0.05；不把GO血管相关富集外推成Hallmark也显著。按C0 v3客观展示规则（最多4条正向+2条负向），本轮绘制实际满足条件的2正+2负曲线，所有50条均在NES总览及完整表中。

## C0 vs C5证据对比

| 证据维度 | C0 | C5 | 解释与限制 |
| --- | --- | --- | --- |
| 细胞与来源覆盖 | 10,717 cells；43 个有检出样本 | 716 cells；37 个有检出样本；InhouseData 382 (53.4%)，GSE222315 283 (39.5%) | C5 92.9% 来自两个 dataset；有检出样本数不等于丰度分析总样本数。 |
| Normal/Tumor 丰度 | dataset 调整后 Tumor−Normal −52.20 个百分点；P=3.06e−7；13 对配对 P=0.00253 | 调整后 +11.61 个百分点；P=0.0378；13 对配对 P=0.0528；min50 P=0.235 | C5 更偏向肿瘤富集，但对小分母样本及来源敏感；C0 的下降不能解释为生物学无意义。 |
| 临床分期 | 未发现可直接复用的同口径单细胞 C0 分期检验结果；既有 Bulk signature 分期结果不替代本项 | NMIBC/MIBC：17/20 样本，P=0.0643；T 分组总体 P=0.0628；Ta/T4 仅3/4样本 | 尚不足以断言随临床分期单调升高；T 分类不等同于总体病理分期。 |
| Full14 连续丰度与 OS | 同四队列中仅 GSE13507 名义支持：调整 HR=1.35，P=0.0245；跨队列 BH q=0.0981 | TCGA/GSE48075 未调整 HR=1.19/1.27；TCGA 跨队列 BH q=0.0144；所有调整模型 P>0.05 | C5 有未调整关联线索，尚无稳定独立预后证据；C0/C5 使用相同队列、OS、SD 和中位数规则。 |
| marker 与样本层验证 | 沿用 C0 v2 Wilcoxon marker 流程作为方法基准 | POSTN、COL11A1、COL10A1、MMP11 等；POSTN/COL11A1 获13个样本 pseudobulk 支持；COL10A1 被表达过滤 | 560个可检验 C5-up 中501个方向一致，293个 pseudobulk BH FDR<0.05且正向；样本来自两个 dataset。 |
| GO/KEGG | 既有 v3 显著代表通路覆盖 ECM、迁移/黏附、血管、炎症与脂质调节 | ECM组织、胶原代谢/降解/合成、黏附与血管发育；KEGG ECM-receptor/integrin 等 | C5 的基质重塑证据集中且有表达支持；不存在显著通路的主题不强行补齐。 |
| Hallmark 主分析及稳定性 | EMT NES=1.65；TNFα/NF-κB NES=2.75；炎症正向，OXPHOS负向 | EMT NES=1.93；TNFα/NF-κB NES=−2.38；两者三阈值均显著；糖酵解正向但25% FDR=0.0545 | 提示不同功能状态；两者各自 Rest 不同，NES差异不是直接 C0 vs C5 检验，也不证明 EMT 转分化或 LLPS 机制。 |
| 既有 LLPS 与 HR signature（仅引用） | LLPS细胞中位数0.1413；HR细胞中位数0.2140；患者等权HR中位数0.1959 | LLPS细胞中位数0.1260；HR细胞中位数0.2429；患者等权HR中位数0.1974 | 逐细胞 HR 优势在患者等权后明显缩小；未重算分数，未用TCGA筛出的694基因建立新模型。 |
| 去卷积可识别性 | Full14 和 focused 均有独立C0；TCGA/UROMOL 两方案 C0 Spearman=0.787/0.759 | 仅Full14有独立C5；150参考细胞来自GSE222315的7个patient_key | 两方案一致性是C0信息，不能验证C5准确性；参考来源单一限制C5的跨队列解释。 |

C5的潜在优势是肿瘤基质重塑相关表达与稳定EMT富集，而非“P值更小”。C0保有更高的既有LLPS分数、不同的炎症状态以及不同的丰度/生存关联。C5尚不能被认定为优于C0的独立预后标志物，更不能由这些结果证明LLPS因果机制。

主要文件：[复用方法、输入和全部图表索引](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/README.md)；[同队列C0/C5 Cox对照](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/03_Bulk_Survival/C0_C5_Cox_comparison.tsv)；[证据对比表TSV](/home/jinshengxi/bladder_metabolism/results/llps_C5_parallel_v1/C0_vs_C5_evidence_comparison.tsv)。
