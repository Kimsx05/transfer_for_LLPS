# MP10 High 定义敏感性分析结果

本轮已实际完成 21 套 fine CellChat 推断。基线 run_id 为 `c2lfinal_20261108T123140Z_0b5fd64f`，旧结果保持原样。C10 是 38 个标本、37 个已记录患者标识；C30 是 30 个标本、30 个已记录患者标识。患者标识来自原始记录/推导，未验证患者关系，不能当作已确认独立患者。每个 Epi 状态分别为 2,651 与 2,479 个抽样细胞，标本 k_s 不相同。

总体判断：数据支持**部分通讯在高分尾部更强**，尤其 Treg、Tfh/CXCL13-like → Epi 的输入，以及 Epi → Arterial EC/Capillary EC 的输出；不支持把全部通讯统一归结为 top 5% 特有信号。H50 已有同向联系，H10 已出现主要变化，H05 进一步增强的幅度因目标而异。因此更合适的表述是“目标与方向特异的尾部增强”，仍不足以确定活性分界。

## 1. 预设三种定义的主结果

以下均为 C10、ALL、min.cells=10、固定 1,933 LR 域、native 公共背景、严格 `probability>0 & internal_P<0.05` 的 W RR。只列 RR>1 的目标。

| High | 方向 | 背景调整 W top3 (RR) |
| --- | --- | --- |
| H50 | 输出 Epi→TME | Mast (1.459)；Neutrophil (1.142)；pDC (1.118) |
| H50 | 输入 TME→Epi | Tfh/CXCL13-like (1.216)；Treg (1.213)；cDC1 (1.074) |
| H10 | 输出 Epi→TME | Mast (1.450)；vCAF (1.141)；CD4 Helper/Conventional (1.135) |
| H10 | 输入 TME→Epi | Plasma (1.979)；Tfh/CXCL13-like (1.370)；Treg (1.337) |
| H05 | 输出 Epi→TME | Mast (1.432)；Arterial EC (1.162)；vCAF (1.160) |
| H05 | 输入 TME→Epi | Plasma (2.209)；Treg (1.482)；Tfh/CXCL13-like (1.455) |


完整 29 个目标的 W/R、原始 High/Low、差/比、背景分子分母、排名和状态见 `tables/native_effects_rankings.tsv.gz`。输出与输入是两种不同方向，不能互换解释。

## 2. 从 H50 到高分尾部的变化

| 方向 | 模式 | 目标 |
| --- | --- | --- |
| output | H50_only_relative_enhancement | 无 |
| output | tail_relative_enhancement_emerges | Venous EC；iCAF |
| output | mixed_direction | Treg |
| output | at_least_one_ratio_undefined | 无 |
| input | H50_only_relative_enhancement | Capillary EC；Pericyte；vCAF |
| input | tail_relative_enhancement_emerges | C1QC TAM；CD4 Helper/Conventional；CD8 Exhausted；CD8 GZMK/Effector-memory-like；Classical Mono；GC B；Mast；Naive/Memory B；Neutrophil；cDC2 |
| input | mixed_direction | SMC |
| input | at_least_one_ratio_undefined | 无 |


“翻转”仅表示 RR 穿过 1 的描述性变化，不是经过标本间差异检验的显著翻转。接近 1 的数值特别需要结合 R、dataset 和 C30 判断。

## 3. W、raw R、dataset 和 C30 复核

表中 dataset 分母是 RR 可定义的 dataset 数，不把 NE 或零分母算作负证据；dataset 中 RR>1 的数量是方向性描述，不是独立重复的显著性检验。各 dataset 还提供与 ALL 在同一背景上的逐项对照。C30 数值来自 C10/C30 共同背景。

| 方向 | 目标 | W_H50 | W_H10 | W_H05 | R_H50 | R_H10 | R_H05 | C30_H05_W | dataset_H05_W>1 | C30_H05_R | dataset_H05_R>1 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| output | Mast | 1.46 | 1.45 | 1.43 | 1.03 | 1.03 | 1.02 | 1.11 | 2/5 | 1.04 | 2/5 |
| output | Arterial EC | 1.06 | 1.13 | 1.16 | 1.08 | 1.15 | 1.18 | 1.14 | 4/4 | 1.15 | 4/4 |
| output | vCAF | 1.11 | 1.14 | 1.16 | 1.13 | 1.16 | 1.18 | 1.08 | 2/3 | 1.09 | 2/3 |
| output | Capillary EC | 1.1 | 1.11 | 1.14 | 1.06 | 1.11 | 1.14 | 1.17 | 4/5 | 1.13 | 4/5 |
| output | Plasma | 0.957 | 0.958 | 0.957 | 0.948 | 0.927 | 0.915 | 0.934 | 3/5 | 0.918 | 2/5 |
| output | Treg | 1.01 | 1 | 0.991 | 1.02 | 1.02 | 1.01 | 1.01 | 2/5 | 1.02 | 1/5 |
| output | Tfh/CXCL13-like | 1.01 | 1.02 | 1.02 | 1.03 | 1.03 | 1.03 | 1.02 | 2/5 | 1.04 | 2/5 |
| output | cDC1 | 0.917 | 0.89 | 0.889 | 0.938 | 0.912 | 0.903 | 0.896 | 2/4 | 0.908 | 2/4 |
| output | Naive/Memory B | 0.923 | 0.902 | 0.887 | 0.944 | 0.915 | 0.901 | 0.893 | 0/4 | 0.909 | 0/4 |
| output | Neutrophil | 1.14 | 1.13 | 1.11 | 1.09 | 1.07 | 1.05 | 1.1 | 0/2 | 1.03 | 0/2 |
| output | pDC | 1.12 | 1.11 | 1.09 | 0.967 | 0.948 | 0.936 | 1.1 | 1/4 | 0.938 | 1/4 |
| input | Mast | 0.916 | 0.974 | 1.05 | 0.99 | 1.02 | 1.07 | 0.996 | 3/5 | 1.03 | 3/5 |
| input | Arterial EC | 1.03 | 1.04 | 1.05 | 1.03 | 1.04 | 1.05 | 1.04 | 2/4 | 1.04 | 3/4 |
| input | vCAF | 1.01 | 0.999 | 0.988 | 0.999 | 0.989 | 0.98 | 0.995 | 1/3 | 0.987 | 2/3 |
| input | Capillary EC | 1 | 0.998 | 0.991 | 0.999 | 0.998 | 0.993 | 0.994 | 1/5 | 0.995 | 3/5 |
| input | Plasma | 1.04 | 1.98 | 2.21 | 1.25 | 1.17 | 1.15 | 1.83 | 3/5 | 1.13 | 1/5 |
| input | Treg | 1.21 | 1.34 | 1.48 | 1.25 | 1.3 | 1.39 | 1.49 | 4/5 | 1.43 | 4/5 |
| input | Tfh/CXCL13-like | 1.22 | 1.37 | 1.46 | 1.24 | 1.34 | 1.39 | 1.37 | 4/5 | 1.33 | 4/5 |
| input | cDC1 | 1.07 | 1.12 | 1.27 | 1.18 | 1.16 | 1.22 | 1.21 | 4/4 | 1.2 | 3/4 |
| input | Naive/Memory B | 0.931 | 0.998 | 1.12 | 1.04 | 1.06 | 1.13 | 1.07 | 3/4 | 1.1 | 3/4 |
| input | Neutrophil | 0.719 | 0.721 | 1.06 | 0.811 | 0.794 | 1.08 | 1.07 | 0/2 | 1.09 | 0/2 |
| input | pDC | 1.05 | 1.08 | 1.12 | 1.08 | 1.1 | 1.12 | 1.09 | 2/4 | 1.1 | 2/4 |


Treg 和 Tfh/CXCL13-like 输入是较一致的尾部增强：H05 的 W/R 在 4/5 dataset 中均 RR>1，C30 保留增强。Arterial EC 输出在 4 个可评估 dataset 中 W/R 均为正，Capillary EC 为 4/5；但不同 dataset 效应幅度不同，部分数值接近 1。cDC1 输入在 4/4 可评估 dataset 中 W>1，raw R 为 3/4，适合作为次级证据。Naive/Memory B 输入在 H05 中由 W<=1 转为 >1，W/R 均有 3/4 dataset 同向，C30 效应较小；其输出在所有可评估 dataset 的 W/R 均未相对增强。

Plasma 输入虽然是 H10/H05 的 W top1，raw R 只有约 1.15，且仅 1/5 dataset 的 raw R>1；其 H50 raw R 反而更高。C30 保留 W 增强，但不能把 W 的排名跃升直接解释为原始通讯强度集中于尾部。Mast 输出在三方案的 C10 W 都排名第一，而 H05 raw R 约 1.02，只有 2/5 dataset 同向；C30 W RR 降至约 1.11。这两个候选都应保留，但证据弱于只看 W top3 所暗示的程度。

H50 独有的输入 RR>1 主要见于 Capillary EC、Pericyte、vCAF，数值仅约 1.001–1.012；它们不能构成强的“宽泛状态独有”证据。Treg 输出的 W 从略大于 1 变为略小于 1，但 raw R 和 C30 并不保持这一翻转，宜视为边界敏感。ALL 29 个目标三种定义均可评估，主比较没有 NE；dataset 的缺失应按覆盖表解释。

## 4. 节点范围、计算域和支持阈值

ALL 的 29 个 TME 在 C10 和 C30 均超过 30 个抽样细胞。因此 ALL 的 min10/min30 后筛不改变可评估范围；共同背景检查确认 H、L、RR、E 相同。各 dataset 小节点不同，具体为：

| scope | node | N |
| --- | --- | --- |
| CNP0000460 | Neutrophil | 6 |
| CNP0000460 | Matrix_CAF | 29 |
| CNP0000460 | rCAF | 30 |
| CNP0000460 | vCAF | 20 |
| CNP0000460 | Lymphatic EC | 10 |
| GSE326225 | LAMP3 DC | 7 |
| GSE326225 | Neutrophil | 10 |
| GSE326225 | cDC1 | 4 |
| GSE326225 | cDC2 | 26 |
| GSE326225 | pDC | 9 |
| GSE326225 | rCAF | 1 |
| GSE326225 | vCAF | 9 |
| GSE326225 | SMC | 0 |
| GSE326225 | Arterial EC | 1 |
| GSE326225 | Capillary EC | 25 |
| GSE326225 | Lymphatic EC | 1 |
| GSE326225 | Venous EC | 5 |
| GSE326225 | GC B | 2 |
| HRA003620 | Neutrophil | 9 |
| HRA003620 | Matrix_CAF | 8 |
| HRA003620 | rCAF | 0 |
| HRA003620 | vCAF | 0 |
| HRA003620 | Pericyte | 3 |
| HRA003620 | SMC | 4 |
| HRA003620 | Arterial EC | 16 |
| HRA003620 | Capillary EC | 19 |
| HRA003620 | Lymphatic EC | 21 |
| InhouseData | Venous EC | 16 |
| InhouseData | GC B | 24 |
| InhouseData | Naive/Memory B | 8 |
| PRJNA662018 | pDC | 26 |


这些节点中 N<=10 的节点没有建模；10<N<=30 的节点参与推断，但在 min30 后筛为 NE。两种后筛在共同可评估节点及背景上数值一致。它们不是两个独立实验。C30 则真正改变标本构成并重新推断，不能与后筛混为一谈。

数据库保持 CellChat 1.6.1 的 1,939 条原定义。旧表达缓存缺少 IFNLR1、IL2RA，已从同一原对象 RNA:data 补齐。6 条因原始 RNA 缺少必需基因而不可计算，详见 LR coverage 表；未将其补零。新固定域 1,933 条包含旧 ALL fine 实际计算域全部 1,875 条，扩展 58 条。C10 ALL W 的新域与旧域交集 RR 最大绝对差为 0；完整逐目标结果见 `LR_domain_sensitivity.tsv`。

Plasma 输入需要特别区分 W 与 R：W 的高分尾部增强幅度明显大于 raw R。`focus_LR_background_contributions.tsv.gz` 同时保留各 LR 的原始概率、内部 P 和支持状态，可定位阈值作用。不同 High 定义也改变参与置换的细胞，因此即使 Low/TME 成员固定，内部 P 仍可能变化。不能把内部 P 当作标本间 High/Low 差异检验。

## 5. 通路与 LR

重点候选由旧六个候选和本轮 ALL 正向 top3 并集预先规则确定。每个目标/方向的作图通路按 C10 ALL H05 W 的正向 E_path 选最多五条，之后固定展示三种 High。图内给出每种定义总体 RR；总体未增强的目标即使有局部正贡献，也不改写总体结论。

下面列出 H05 的主要正向通路贡献 (总体背景比固定，所有通路可加总为总体 E)：

| direction | target | pathway | E_path | RR |
| --- | --- | --- | --- | --- |
| output | Capillary EC | VEGF | 0.0536 | 1.14 |
| output | Mast | MIF | 0.0406 | 1.43 |
| output | Arterial EC | VEGF | 0.0361 | 1.16 |
| output | Capillary EC | LAMININ | 0.0304 | 1.14 |
| output | Arterial EC | LAMININ | 0.03 | 1.16 |
| output | Arterial EC | VISFATIN | 0.029 | 1.16 |
| output | Capillary EC | VISFATIN | 0.0287 | 1.14 |
| output | pDC | MIF | 0.0241 | 1.09 |
| output | vCAF | LAMININ | 0.0237 | 1.16 |
| input | vCAF | APP | 0.0218 | 0.988 |
| output | Mast | LAMININ | 0.0214 | 1.43 |
| output | Tfh/CXCL13-like | NECTIN | 0.0194 | 1.02 |
| input | Arterial EC | NOTCH | 0.0183 | 1.05 |
| output | Treg | NECTIN | 0.0164 | 0.991 |
| input | Treg | TIGIT | 0.0138 | 1.48 |
| output | vCAF | ANGPTL | 0.0127 | 1.16 |
| output | Treg | LAMININ | 0.0119 | 0.991 |
| output | Tfh/CXCL13-like | LAMININ | 0.0118 | 1.02 |
| output | Naive/Memory B | LAMININ | 0.0114 | 0.887 |
| input | vCAF | NOTCH | 0.0102 | 0.988 |
| input | Tfh/CXCL13-like | TIGIT | 0.00985 | 1.46 |
| input | vCAF | FN1 | 0.00939 | 0.988 |
| input | Arterial EC | EGF | 0.0083 | 1.05 |
| output | Neutrophil | ANNEXIN | 0.00829 | 1.11 |
| input | Arterial EC | EPHA | 0.00751 | 1.05 |
| output | Neutrophil | ADGRE5 | 0.00738 | 1.11 |
| output | pDC | GDF | 0.00734 | 1.09 |
| output | vCAF | EPHA | 0.00667 | 1.16 |
| input | Capillary EC | HSPG | 0.00664 | 0.991 |
| output | Plasma | LAMININ | 0.00654 | 0.957 |
| input | Capillary EC | APP | 0.00637 | 0.991 |
| output | Tfh/CXCL13-like | CLEC | 0.00614 | 1.02 |
| output | Neutrophil | LAMININ | 0.00537 | 1.11 |
| input | Tfh/CXCL13-like | CD6 | 0.00524 | 1.46 |
| input | pDC | APP | 0.00523 | 1.12 |
| input | pDC | EGF | 0.00509 | 1.12 |
| input | Capillary EC | NOTCH | 0.00467 | 0.991 |
| input | Tfh/CXCL13-like | CD96 | 0.00438 | 1.46 |
| input | Mast | EGF | 0.00419 | 1.05 |
| output | pDC | ADGRE5 | 0.0036 | 1.09 |
| output | cDC1 | NECTIN | 0.00329 | 0.889 |
| input | cDC1 | CD226 | 0.0031 | 1.27 |
| output | cDC1 | LAMININ | 0.00295 | 0.889 |
| output | cDC1 | JAM | 0.00288 | 0.889 |
| input | Neutrophil | VISFATIN | 0.00287 | 1.06 |
| output | Naive/Memory B | ADGRE5 | 0.00287 | 0.887 |
| input | Plasma | ADGRE5 | 0.0028 | 2.21 |
| output | Naive/Memory B | GDF | 0.0025 | 0.887 |
| output | Mast | GDF | 0.00208 | 1.43 |
| output | Treg | JAM | 0.00195 | 0.991 |
| input | Treg | CD96 | 0.00154 | 1.48 |
| input | Naive/Memory B | ADGRE5 | 0.00135 | 1.12 |
| input | Mast | ADGRE5 | 0.00129 | 1.05 |
| input | pDC | ADGRE5 | 0.00122 | 1.12 |
| input | cDC1 | EGF | 0.0012 | 1.27 |
| input | Mast | CD46 | 0.00119 | 1.05 |
| output | Plasma | VISFATIN | 0.00107 | 0.957 |
| input | cDC1 | CD46 | 0.00104 | 1.27 |
| input | Treg | CD46 | 0.000988 | 1.48 |
| output | Plasma | ANGPTL | 0.000689 | 0.957 |
| input | Plasma | CD46 | 0.00059 | 2.21 |
| input | Naive/Memory B | VISFATIN | 0.000205 | 1.12 |
| input | Plasma | VISFATIN | 0.000108 | 2.21 |


对应 W/R、High/Low 和主要 LR 见 `focus_top10_positive_LR_contributions.tsv.gz`。关键 H05 项的实际分解如下，H/L probability 是未按内部 P 筛选的值，E_LR 为 supported W 的贡献：

| sender_H | receiver_H | interaction_id | DB_ligand | DB_receptor | H_probability | L_probability | H_internal_P | L_internal_P | E_LR |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Epi High | Mast | MIF_CD74_CD44 | MIF | CD74_CD44 | 0.0406 | 0.0363 | 0.02 | 0.92 | 0.0406 |
| Epi High | Mast | LAMB3_CD44 | LAMB3 | CD44 | 0.0358 | 0.0199 | 0 | 0 | 0.00898 |
| Epi High | Mast | LAMC2_CD44 | LAMC2 | CD44 | 0.00637 | 0 | 0 | 1 | 0.00637 |
| Epi High | Capillary EC | VEGFA_VEGFR1 | VEGFA | FLT1 | 0.0591 | 0.0253 | 0 | 0 | 0.0252 |
| Epi High | Capillary EC | NAMPT_INSR | NAMPT | INSR | 0.0947 | 0.0552 | 0 | 0 | 0.0207 |
| Epi High | Capillary EC | VEGFA_VEGFR1R2 | VEGFA | FLT1_KDR | 0.0416 | 0.0176 | 0 | 0 | 0.018 |
| Epi High | Arterial EC | VEGFA_VEGFR1 | VEGFA | FLT1 | 0.0547 | 0.0233 | 0 | 0 | 0.0234 |
| Epi High | Arterial EC | NAMPT_INSR | NAMPT | INSR | 0.0847 | 0.0491 | 0 | 0 | 0.0189 |
| Epi High | Arterial EC | NAMPT_ITGA5_ITGB1 | NAMPT | ITGA5_ITGB1 | 0.0425 | 0.0242 | 0 | 0 | 0.0101 |
| Treg | Epi High | TIGIT_NECTIN2 | TIGIT | NECTIN2 | 0.0299 | 0.0112 | 0 | 0 | 0.0138 |
| Treg | Epi High | CD96_NECTIN1 | CD96 | NECTIN1 | 0.0054 | 0.00269 | 0 | 0 | 0.00154 |
| Treg | Epi High | CD46_JAG1 | CD46 | JAG1 | 0.000988 | 0 | 0 | 1 | 0.000988 |
| Tfh/CXCL13-like | Epi High | TIGIT_NECTIN2 | TIGIT | NECTIN2 | 0.0212 | 0.00793 | 0 | 0 | 0.00985 |
| Tfh/CXCL13-like | Epi High | CD6_ALCAM | CD6 | ALCAM | 0.0167 | 0.00801 | 0 | 0 | 0.00524 |
| Tfh/CXCL13-like | Epi High | CD96_NECTIN1 | CD96 | NECTIN1 | 0.0156 | 0.0078 | 0 | 0 | 0.00438 |
| cDC1 | Epi High | CD226_NECTIN2 | CD226 | NECTIN2 | 0.00663 | 0.00245 | 0 | 0 | 0.0031 |
| cDC1 | Epi High | CD46_JAG1 | CD46 | JAG1 | 0.00104 | 0 | 0 | 1 | 0.00104 |
| cDC1 | Epi High | AREG_EGFR | AREG | EGFR | 0.00504 | 0.00281 | 0 | 0 | 0.001 |
| Plasma | Epi High | ADGRE5_CD55 | ADGRE5 | CD55 | 0.0028 | 0.00173 | 0 | 1 | 0.0028 |
| Plasma | Epi High | CD46_JAG1 | CD46 | JAG1 | 0.00059 | 0 | 0 | 1 | 0.00059 |
| Plasma | Epi High | NAMPT_INSR | NAMPT | INSR | 0.00143 | 0.000922 | 0 | 0 | 0.000108 |


Plasma 的 ADGRE5_CD55 在 Low 中原始 probability 仍为正，但内部 P=1，故其 supported Low 贡献为零；Mast 输出的 MIF_CD74_CD44 也有正的 Low 原始 probability，而内部 P=0.92、High P=0.02。这直接显示支持掩码如何放大 W 对比，不能把 Low 未支持解释为该联系未计算或无原始表达信号。

本报告输入/输出沿用 CellChat 1.6.1 的 DB 编码角色。例如 TIGIT_NECTIN2 在该 DB 中写为 ligand=TIGIT、receptor=NECTIN2；这里报告的是计算方向，不能据此单独认定 TIGIT 为生物学配体，或 Treg/Tfh 为生物学上游。DB_ligand/DB_receptor、发送细胞、接收细胞已分别给出，机制方向仍需另行验证。

## 6. 可复现性与检查

全部参数、版本、数据库与 feature/LR 域、缺失原因、真实输入顺序/hash、抽样优先级、原始未过滤 probability/P、汇总和直接绘图数据已保存。数值推断对象只有紧凑缓存，不重复保存大型表达或完整 CellChat 对象。`run_all.sh` 可按内容 hash 断点复用；`scripts/04_replot.R` 只读结果表，完全不依赖重跑 CellChat。

本轮输入与结果核验 2476 项，全部通过；另有原生过滤/P 边界、通路加总及共同背景后筛一致性检查。图为 Arial 约 12 pt、白底矢量 PDF、300 dpi PNG；有图例时提供完整、仅主体及仅图例 PDF。直接绘图数据及统一色标见 `plot_data/` 与 `config/figure_manifest.tsv`。热图真实零、NE 和比值未定义分别标识，无伪计数。

建议继续保留三套定义，以 Treg/Tfh 输入和血管内皮输出为后续验证重点，同时将 Plasma、Mast 的支持阈值和标本构成敏感性单独列为限制。这里提供证据优先级，不替用户冻结最终候选或阈值。

本报告保留 H50/H10/H05 全套预设结果，不根据最显著方案挑选最终阈值。H05 应称为高分尾部，不能直接称作真实活跃细胞；其他 95% 也没有被判定为噪声。结果讨论与最终判断需结合上述全部复核，而不能仅凭 pooled ALL 的 top3。
