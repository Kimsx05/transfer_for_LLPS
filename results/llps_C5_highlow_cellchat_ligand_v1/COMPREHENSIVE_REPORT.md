# LLPS C5 High/Low CellChat 与候选配体关联分析 v1

本轮本地分析完成：55个联合网络，主阈值10有12名配对肿瘤患者。LLPS-high总体outgoing/incoming较高 (FDR约0.021/0.020)，但没有单个LR或pathway通过全量校正。Endothelial差异方向在两个队列一致，显著支持主要来自InhouseData；RNA复杂度及队列与临床组完全混杂限制机制解释。深度调整后HSPG2保留预设22候选的患者级关联支持，但全606基因筛选没有FDR支持。具体10个问题见第9节。

交付状态：本地报告、统计表、图片、代码和QC已完成。GitHub写入接口返回403，尚未上传。轻量结果包及上传范围见delivery/UPLOAD_STATUS.md；不含表达矩阵和重复对象。


日期：2026-11-25。观察性分析；seed=123；至多4个CPU；未重跑NicheNet、cell2location或空间共定位。

## 1. 输入与分组可行性

复核59样本、42独立患者；C5为716细胞、37样本、31患者。患者ID沿用最终审计并包含dataset前缀；T4463209a/b及BC4/BCN保持取材层级。Normal独立描述，不与肿瘤合并。

| min_cells | specimens | patients |
| --- | --- | --- |
| 5 | 14 | 14 |
| 10 | 12 | 12 |
| 20 | 4 | 4 |
| 30 | 3 | 3 |

主方案为 specimen 内 C5 score>中位数定义high，其余low，不随机拆同分细胞。无C5的22个样本也按high=low=0保留于审计表并标为ineligible。原分数无零值，单样本恰等于中位数最多1细胞。LLPS77是转录程序评分，不是凝聚体形成或相分离活性的直接测量。主阈值10未依据结果选择；实际CellChat剔除n<=阈值。主阈值配对肿瘤患者12人，20/30仅4/3人，不作患者级正式检验。不同伙伴分母见图和完整表。

四种分组分别推断：within_median、exclude_median_ties、top_bottom30 (严格高于70%/低于30%分位，边界ties排除)、global_tumor_median (仅肿瘤)。全局方案更易反映患者/队列水平差异，不与样本内结果合并。主方案eligible样本两组最多相差1细胞，最大low/high约1.083，因此未增加非必要的等量下采样；全局组不平衡另在计数表保留。

## 2. 联合网络与统计定义

每specimen和分组方法独立联合构建C5-high/C5-low与九个伙伴群：Epithelial、Endothelial、Myeloid、T/NK、B、Plasma、Mast、Other fibroblasts、Pericyte。按用户明确确认，全部正式C5的716细胞保留 (包括原atlas标为Pericyte的22个)；非C5的5895个原Pericyte单列，其余13535个归入Other fibroblasts。两群均不含C5；逐细胞映射见00_input_audit/stromal_identity_final.csv.gz。旧合并网络已隔离，不用于本轮结论。沿用v2的通信分析纳入范围：major Cycling虽通过atlas QC，但作为非明确谱系的状态大类，历史通信输入未纳入；本轮不重新赋予其谱系身份。历史未纳入正式C0–C6的stromal细胞亦不重新加入。Cycling排除不等同于QC-fail，正式C6按照相同优先级分配至Other fibroblasts或Pericyte。

完整CellChatDB.human、RNA LogNormalize data、truncatedMean、trim=0.1、population.size=FALSE、nboot=100、seed.use=123、raw.use=TRUE；实际版本及函数源码见00_input_audit。候选LR在同一联合网络筛选，两组采用同一候选universe。原生meta缺少samples列时CellChat自动赋值sample1；已核验55个对象各自只含一份specimen，样本真实身份保存在job及逐细胞映射中，未因此跨样本混合。四阈值共享同方法同样本未过滤概率；仅当barcode/state标签完全相同时跨方法复用新计算缓存，逐项记录并分别导出。

不满足high、low和partner三方细胞数条件时记NA；满足且未检出记0。数据库未入候选集的LR视为流程未检出，p-value保留NA。原生LR/aggregate使用p<0.05；pathway采用原生p<=0.05。CellChat置换P不是患者high/low差异P。CellChat推断强度不是实际分泌量或分子通讯速率。

总强度使用两组相同eligible partners，排除self及high↔low。先算样本内high-low，再按患者肿瘤specimen取差值中位数；high/low强度亦分别保留患者中位数，差值中位数不一定等于两者相减。至少5患者才计算双侧配对Wilcoxon (差值单样本signed-rank，近似法) 和2000次患者bootstrap中位差percentile 95% CI。BH按分组×阈值×统计类型×方向×分层各自完整检验集合校正，CI为未多重校正探索性区间。dataset与LODO均保留；小于5人的分层只报告效应，P/FDR/CI为NA。

## 3. 总体与伙伴通信变化

| direction | partner | n_patients | high_median | low_median | difference_median | CI_low | CI_high | high_greater_fraction | FDR |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| incoming | B | 11 | 0.3183 | 0.2578 | 0.05339 | -0.007539 | 0.09629 | 0.7273 | 0.1111 |
| incoming | Endothelial | 12 | 2.442 | 2.259 | 0.2674 | 0.1021 | 0.4505 | 0.9167 | 0.01974 |
| incoming | Epithelial | 10 | 0.6811 | 0.6685 | 0.0457 | 0.007351 | 0.2322 | 0.9 | 0.05385 |
| incoming | Mast | 8 | 0.5042 | 0.428 | 0.06993 | 0.009352 | 0.1572 | 0.75 | 0.06045 |
| incoming | Myeloid | 12 | 0.5311 | 0.5321 | 0.07135 | 0.003979 | 0.1096 | 0.75 | 0.1111 |
| incoming | Other_fibroblasts | 12 | 5.062 | 4.042 | 0.7075 | 0.3379 | 1.195 | 0.9167 | 0.02147 |
| incoming | Pericyte | 12 | 4.492 | 3.903 | 0.6252 | 0.312 | 1.018 | 0.9167 | 0.02157 |
| incoming | Plasma | 10 | 0.1883 | 0.1608 | 0.02711 | 0.01056 | 0.03719 | 1 | 0.01974 |
| incoming | T_NK | 11 | 0.3668 | 0.3801 | 0.003252 | -0.01962 | 0.05924 | 0.5455 | 0.4498 |
| outgoing | B | 11 | 2.461 | 2.251 | 0.07151 | 0.004352 | 0.1131 | 0.8182 | 0.02069 |
| outgoing | Endothelial | 12 | 3.828 | 3.756 | 0.2303 | 0.07203 | 0.3361 | 0.75 | 0.02069 |
| outgoing | Epithelial | 10 | 2.974 | 2.828 | 0.1564 | -0.007246 | 0.517 | 0.7 | 0.04545 |
| outgoing | Mast | 8 | 3.202 | 2.974 | 0.2073 | 0.1167 | 0.4146 | 1 | 0.02069 |
| outgoing | Myeloid | 12 | 2.951 | 2.939 | 0.1342 | 0.02827 | 0.2981 | 0.8333 | 0.02069 |
| outgoing | Other_fibroblasts | 12 | 4.121 | 3.672 | 0.2025 | 0.07249 | 0.2865 | 0.8333 | 0.02069 |
| outgoing | Pericyte | 12 | 4.854 | 4.483 | 0.2821 | 0.06658 | 0.4214 | 0.9167 | 0.02069 |
| outgoing | Plasma | 10 | 1.318 | 1.225 | 0.07008 | -0.005831 | 0.1554 | 0.7 | 0.04039 |
| outgoing | T_NK | 11 | 3.29 | 3.262 | 0.09326 | 0.03903 | 0.1675 | 0.8182 | 0.04545 |
| incoming | TOTAL | 12 | 14.64 | 12.77 | 1.762 | 0.8375 | 2.834 | 0.9167 | 0.01974 |
| outgoing | TOTAL | 12 | 28.41 | 25.86 | 1.359 | 0.2818 | 2.133 | 0.75 | 0.02069 |

总体增强必须结合方向、患者比例、CI和FDR，不因某些强边高于low便推断普遍增强。高阈值与低阈值共同样本概率相同是计算复用结果，不是独立验证。

### Endothelial

| direction | partner | n_patients | high_median | low_median | difference_median | CI_low | CI_high | high_greater_fraction | FDR |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| incoming | Endothelial | 12 | 2.442 | 2.259 | 0.2674 | 0.1021 | 0.4505 | 0.9167 | 0.01974 |
| outgoing | Endothelial | 12 | 3.828 | 3.756 | 0.2303 | 0.07203 | 0.3361 | 0.75 | 0.02069 |

### Pathway (按绝对配对效应展示，含反向)

| direction | partner | feature | n_patients | difference_median | CI_low | CI_high | FDR |
| --- | --- | --- | --- | --- | --- | --- | --- |
| incoming | Other_fibroblasts | COLLAGEN | 12 | 0.4665 | 0.3043 | 0.8309 | 1 |
| incoming | Pericyte | COLLAGEN | 12 | 0.4477 | 0.3088 | 0.7571 | 1 |
| incoming | Endothelial | COLLAGEN | 12 | 0.1425 | 0.1013 | 0.2793 | 1 |
| outgoing | Mast | COLLAGEN | 8 | 0.1113 | 0.04669 | 0.1356 | 1 |
| outgoing | Pericyte | COLLAGEN | 12 | 0.09415 | 0.004135 | 0.1645 | 1 |
| outgoing | Other_fibroblasts | COLLAGEN | 12 | 0.08701 | -0.0057 | 0.142 | 1 |
| incoming | Pericyte | LAMININ | 12 | 0.08689 | -0.03299 | 0.1432 | 1 |
| incoming | Other_fibroblasts | LAMININ | 12 | 0.08429 | -0.02644 | 0.1297 | 1 |
| outgoing | Endothelial | COLLAGEN | 12 | 0.06623 | -0.007244 | 0.1054 | 1 |
| outgoing | Pericyte | FN1 | 12 | 0.06066 | -0.0006525 | 0.1185 | 1 |
| incoming | Other_fibroblasts | FN1 | 12 | 0.06025 | 0.0129 | 0.09607 | 1 |
| incoming | Endothelial | LAMININ | 12 | 0.05992 | -0.01628 | 0.07627 | 1 |
| outgoing | T_NK | COLLAGEN | 11 | 0.05881 | 0.007575 | 0.1041 | 1 |
| outgoing | Endothelial | FN1 | 12 | 0.05377 | -0.001379 | 0.06358 | 1 |
| incoming | Pericyte | FN1 | 12 | 0.05343 | 0.01586 | 0.07949 | 1 |
| outgoing | Other_fibroblasts | FN1 | 12 | 0.05074 | -0.001686 | 0.08475 | 1 |
| outgoing | Epithelial | THBS | 10 | 0.05049 | 0 | 0.1267 | 1 |
| outgoing | Mast | FN1 | 8 | 0.05002 | 0.01203 | 0.08355 | 1 |
| outgoing | Epithelial | COLLAGEN | 10 | 0.04623 | -0.008093 | 0.1888 | 1 |
| outgoing | Plasma | COLLAGEN | 10 | 0.0459 | -0.006141 | 0.05557 | 1 |

### LR (按绝对配对效应展示，含反向)

| direction | partner | feature | ligand | receptor | annotation | n_patients | difference_median | FDR |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| incoming | Other_fibroblasts | COL1A2_ITGA11_ITGB1 | COL1A2 | ITGA11_ITGB1 | ECM-Receptor | 12 | 0.03001 | 1 |
| outgoing | Mast | FN1_CD44 | FN1 | CD44 | ECM-Receptor | 8 | 0.02747 | 1 |
| incoming | Other_fibroblasts | COL1A1_ITGA11_ITGB1 | COL1A1 | ITGA11_ITGB1 | ECM-Receptor | 12 | 0.02744 | 1 |
| incoming | Other_fibroblasts | COL6A3_ITGA11_ITGB1 | COL6A3 | ITGA11_ITGB1 | ECM-Receptor | 12 | 0.0254 | 1 |
| incoming | Other_fibroblasts | COL6A2_ITGA11_ITGB1 | COL6A2 | ITGA11_ITGB1 | ECM-Receptor | 12 | 0.02392 | 1 |
| incoming | Pericyte | COL4A2_ITGA11_ITGB1 | COL4A2 | ITGA11_ITGB1 | ECM-Receptor | 12 | 0.02181 | 1 |
| incoming | Pericyte | COL4A1_ITGA11_ITGB1 | COL4A1 | ITGA11_ITGB1 | ECM-Receptor | 12 | 0.02116 | 1 |
| incoming | Pericyte | COL6A2_ITGA11_ITGB1 | COL6A2 | ITGA11_ITGB1 | ECM-Receptor | 12 | 0.02014 | 1 |
| incoming | Other_fibroblasts | COL6A1_ITGA11_ITGB1 | COL6A1 | ITGA11_ITGB1 | ECM-Receptor | 12 | 0.01999 | 1 |
| incoming | Pericyte | COL1A2_ITGA11_ITGB1 | COL1A2 | ITGA11_ITGB1 | ECM-Receptor | 12 | 0.01933 | 1 |
| incoming | Endothelial | COL4A1_ITGA11_ITGB1 | COL4A1 | ITGA11_ITGB1 | ECM-Receptor | 12 | 0.01881 | 1 |
| outgoing | Mast | COL6A3_CD44 | COL6A3 | CD44 | ECM-Receptor | 8 | 0.01877 | 1 |
| outgoing | T_NK | FN1_CD44 | FN1 | CD44 | ECM-Receptor | 11 | 0.01865 | 1 |
| outgoing | Endothelial | FN1_ITGA5_ITGB1 | FN1 | ITGA5_ITGB1 | ECM-Receptor | 12 | 0.01843 | 1 |
| outgoing | B | FN1_CD44 | FN1 | CD44 | ECM-Receptor | 11 | 0.01764 | 1 |
| incoming | Pericyte | COL1A1_ITGA11_ITGB1 | COL1A1 | ITGA11_ITGB1 | ECM-Receptor | 12 | 0.01731 | 1 |
| incoming | Endothelial | COL4A2_ITGA11_ITGB1 | COL4A2 | ITGA11_ITGB1 | ECM-Receptor | 12 | 0.01711 | 1 |
| outgoing | Myeloid | FN1_CD44 | FN1 | CD44 | ECM-Receptor | 12 | 0.01677 | 1 |
| outgoing | Pericyte | FN1_ITGAV_ITGB1 | FN1 | ITGAV_ITGB1 | ECM-Receptor | 12 | 0.01644 | 1 |
| outgoing | Other_fibroblasts | FN1_ITGAV_ITGB1 | FN1 | ITGAV_ITGB1 | ECM-Receptor | 12 | 0.01397 | 1 |

通路命名：POSTN=PERIOSTIN；MDK=MK；TGF-beta=TGFb；WNT包括WNT/ncWNT。完整映射及四类DB注释见00_input_audit/pathway_aliases.csv。新变化条目只是本轮候选，不声称发现新生物学结合关系。

## 4. 配体连续关联与pseudobulk

主评分保持正式LLPS77_AUC；逐细胞Pearson/Spearman仅描述，不将716细胞当716独立重复。各specimen内相关先控制log1p(nCount_RNA)/log1p(nFeature_RNA)，再取患者效应中位数；至少5个C5细胞的specimen用于患者内关联汇总。另有全体肿瘤细胞sample固定效应+深度残差相关，样本固定效应同时吸收patient/dataset，不提供细胞独立P。

肿瘤C5 raw counts按患者合并，TMM库大小归一化后logCPM连续关联 (logCPM保留10位小数以消除零表达值的浮点尾差，防止人为打破Spearman ties)；纳入C5细胞数至少5/10/20，均保留平均/中位LLPS分数。调整模型为score+dataset，审查秩与残差自由度，df<3或患者<5不报告推断P。同患者high/low pseudobulk另用TMM/voom与patient+state设计；保留未通过表达过滤的候选为未评估，不能填零。重点22基因与全部表达DB配体分别BH。低于3名患者检测到的基因不报告患者相关推断P，近零方差按1e-8容差标记无法评估。另提供dataset加患者平均nCount/nFeature调整模型，避免忽略评分与RNA复杂度的混杂。

| gene | n_patients | pearson | spearman | adjusted_beta | adjusted_p | adjusted_p_FDR_focus | residual_df |
| --- | --- | --- | --- | --- | --- | --- | --- |
| HSPG2 | 15 | 0.7927 | 0.8679 | 42.49 | 4.991e-05 | 0.001098 | 11 |
| FN1 | 15 | 0.7501 | 0.75 | 63.18 | 0.0004371 | 0.004161 | 11 |
| COL1A1 | 15 | 0.6892 | 0.7679 | 59.16 | 0.0005674 | 0.004161 | 11 |
| POSTN | 15 | 0.6437 | 0.7107 | 60.48 | 0.001454 | 0.007997 | 11 |
| COL1A2 | 15 | 0.6948 | 0.725 | 45.77 | 0.002086 | 0.00918 | 11 |
| COL3A1 | 15 | 0.6855 | 0.7429 | 47.64 | 0.00339 | 0.01243 | 11 |
| LAMA1 | 15 | 0.738 | 0.6721 | 35.16 | 0.008732 | 0.02744 | 11 |
| WNT5A | 15 | -0.504 | -0.4143 | -55.81 | 0.01432 | 0.03938 | 11 |
| TGFB2 | 15 | -0.1589 | -0.2 | -21.49 | 0.05452 | 0.1333 | 11 |
| VEGFA | 15 | 0.4894 | 0.4607 | 26.28 | 0.06866 | 0.151 | 11 |
| FGF1 | 15 | -0.5364 | -0.5536 | -18.89 | 0.08377 | 0.1675 | 11 |
| ANGPTL4 | 15 | 0.3598 | 0.3601 | 43.49 | 0.09636 | 0.1767 | 11 |
| TGFB1 | 15 | 0.196 | 0.3179 | 16.08 | 0.1152 | 0.195 | 11 |
| MDK | 15 | 0.05909 | 0.03214 | 13.32 | 0.3119 | 0.4828 | 11 |
| TNFSF12 | 15 | 0.1308 | -0.04643 | 14.84 | 0.3292 | 0.4828 | 11 |
| WNT4 | 15 | -0.1622 | -0.06776 | -20.44 | 0.4012 | 0.5224 | 11 |
| WNT2 | 15 | 0.2064 | 0.2107 | 10.86 | 0.4274 | 0.5224 | 11 |
| SEMA3C | 15 | 0.2782 | 0.2607 | -9.059 | 0.4174 | 0.5224 | 11 |
| ADM | 15 | 0.1626 | 0.2143 | 8.237 | 0.4852 | 0.5618 | 11 |
| IGF2 | 15 | -0.3245 | -0.2613 | 11.54 | 0.6121 | 0.6733 | 11 |
| INHBA | 15 | 0.2148 | 0.2357 | 6.719 | 0.7239 | 0.7584 | 11 |
| FGF2 | 15 | 0.3732 | 0.2796 | 1.474 | 0.8935 | 0.8935 | 11 |

### 重点配体的患者内调整关联

| gene | n | median | ci_low | ci_high | FDR_focus |
| --- | --- | --- | --- | --- | --- |
| COL1A1 | 19 | 0.2585 | 0.1932 | 0.4407 | 0.005088 |
| COL3A1 | 19 | 0.2752 | 0.1853 | 0.5091 | 0.005496 |
| COL1A2 | 19 | 0.3287 | 0.0799 | 0.4207 | 0.008781 |
| SEMA3C | 18 | 0.2187 | 0.05004 | 0.3215 | 0.01093 |
| POSTN | 19 | 0.1643 | 0.08558 | 0.272 | 0.03277 |
| LAMA1 | 11 | -0.1116 | -0.1866 | 0.02262 | 0.2506 |
| HSPG2 | 18 | 0.2396 | -0.1283 | 0.359 | 0.3336 |
| ADM | 19 | -0.04599 | -0.1875 | 0.01728 | 0.3336 |
| FN1 | 19 | 0.2194 | -0.02856 | 0.4568 | 0.3369 |
| MDK | 19 | -0.08898 | -0.2001 | 0.0241 | 0.3369 |
| FGF1 | 18 | 0.04329 | -0.04531 | 0.2197 | 0.4132 |
| WNT2 | 19 | 0.03923 | -0.07672 | 0.2589 | 0.4609 |
| TGFB2 | 17 | 0.05394 | -0.04176 | 0.1486 | 0.6671 |
| WNT5A | 18 | -0.05553 | -0.09805 | 0.03102 | 0.7214 |
| FGF2 | 15 | 0.005471 | -0.05714 | 0.1284 | 0.7694 |
| INHBA | 19 | -0.07544 | -0.2427 | 0.1863 | 0.7694 |
| TNFSF12 | 16 | 0.01152 | -0.2081 | 0.09658 | 0.9535 |
| WNT4 | 10 | -0.002139 | -0.1533 | 0.2524 | 0.9972 |
| TGFB1 | 19 | 0.00966 | -0.1256 | 0.1626 | 0.9972 |
| IGF2 | 13 | -0.04838 | -0.08775 | 0.1822 | 0.9972 |
| VEGFA | 19 | -0.04014 | -0.1612 | 0.1849 | 0.9972 |
| ANGPTL4 | 14 | 0.02345 | -0.1065 | 0.1144 | 1 |

### 评分重叠敏感性

重点候选只有COL1A1属于LLPS77；全配体筛选还有FGF10/GAL/TNXB/CDH2/DSC1/L1CAM/NCAM1/NCAM2/LRFN5重叠。已逐基因删除重算AUC，沿用原10%排名门槛。由于原分数使用chunk-wise seed且存在随机ties，新seed123重建完整score与原score相关约0.923，最大差约0.0528。因此敏感性比较的是同一重建排名上的完整/删基因分数，并另列原score，不将两者混作完全同一评分；这一技术限制会降低重叠基因的结论强度。

## 5. 队列、患者与未拆分C5敏感性

测序深度敏感性是重要限制：主模型中的多数基质配体关联在额外加入患者平均log(nCount)/log(nFeature)后不再获FDR支持，HSPG2保留22个预设候选的FDR支持，但全部606基因的深度调整探索性筛选无FDR<0.05结果；这不等于证明其与LLPS有因果联系。具体完整结果见pseudobulk_dataset_depth_adjusted.csv。各队列及leave-one-dataset-out结果位于完整统计表stratum列；重点pseudobulk另有leave-one-patient-out系数。低细胞数相关和小样本Wilcoxon均可能不稳健，需结合效应方向，不能以一次FDR判断机制。

07_sensitivity_analysis/unsplit_C5_continuous_communication.csv复用已有未拆分C5通信，分析患者平均LLPS与强度的连续关系及队列内high/low描述。这回答患者间差异，不能替代患者内部细胞高低比较。不同患者总强度的eligible伙伴组成可能不同，具体伙伴结果比TOTAL更易解释。

## 6. 既有NicheNet和空间证据

综合表逐配体记录C5-vs-C0表达、患者内/患者间LLPS关联、所有C5→Endothelial LR的双向效应、旧NicheNet各program activity/rank、21次LOPO稳定性及空间摘要；不设置人为综合分数。原始每轴数据保留在长表，不用单个最大效应代表所有受体。

既有NicheNet Set B基于名义响应基因，不是FDR显著主模型DEG；旧内皮主模型EMT/Angiogenesis/TGF-beta未通过BH<0.10，不能将标准程序先验或Set B排名解释为已观察到的内皮响应。LOPO范围不是置信区间。

既有63张空间切片的Five-gene分数不是已验证C5丰度；正式C5-specific定位NO-GO不变。部分患者身份未知、细胞组成/深度混杂和匹配失衡仍存在，表达邻接不证明真实细胞通讯。Five与LLPS77的空间关联是program共现，不是本轮细胞内LLPS机制的独立因果验证。

## 7. 证据分类与后续机制问题

可靠部分：正式输入/映射复核、透明阈值与患者分母、保留原始评分、联合网络和全部阴性输出。统计支持的观察性关联仍需独立患者和扰动实验。探索性部分：小患者分层、所有机制解释、名义NicheNet程序及空间邻接。阴性部分按FDR、效应和CI如实展示；无法评估部分以NA及status保存，尤其高阈值和低覆盖受体。

当前最值得验证的假说应由本轮患者内关联、患者间关联及内皮具体LR方向是否一致来决定；优先检验LLPS-associated转录状态与配体表达/内皮响应的关联，不预设促血管生成或EndMT。需用细胞身份明确的原位证据、配体阻断与LLPS相关扰动区分共同细胞状态、基质程序及真正因果链。

## 8. 文件、图形与可复现性

所有图独立保存PDF/PNG300dpi，白底Arial约12pt；有图例另有panel_only/legend_only PDF，每图附CSV。圈图high/low使用同节点和线宽标度、所有九个伙伴均eligible的同一批患者 (完整患者名单见figures/circle_common_patients.csv)，棕色outgoing/蓝色incoming；每条边N见绘图数据。双三角high在左上、low在右下，数值为患者中位强度，灰色NA与0不同。Figure I的不同证据量纲采用独立色标，另输出各证据分图，不设置加权得分。

scripts为完整复现代码；02_CellChat_highlow保留紧凑未过滤/阈值网络checkpoint，不另存大型Seurat或重复样本表达。GitHub只上传轻量报告、图、表和脚本，不上传表达矩阵和RDS网络。所有路径及文件大小见output_manifest.csv，最终QC见00_input_audit/final_QC.csv。



## 9. 十个研究问题的直接回答

### 1. 可评估患者数

主展示有12名肿瘤患者 (7 MIBC、5 NMIBC)，Normal另有1名。各边因partner数量要求有不同分母；Endothelial 12，Epithelial 10，Myeloid 12，T/NK 11，B 11，Plasma 10，Mast 8，Other fibroblasts 12，Pericyte 12。阈值5/20/30的肿瘤覆盖为14/4/3；后两者不作患者级检验。

### 2. 总体通讯是否增强

incoming TOTAL: n=12, high=14.64, low=12.77, 配对中位差=1.762, 95% CI [0.8375, 2.834], FDR=0.01974, high>low患者比例=91.7%。

outgoing TOTAL: n=12, high=28.41, low=25.86, 配对中位差=1.359, 95% CI [0.2818, 2.133], FDR=0.02069, high>low患者比例=75.0%。

主方案总outgoing与incoming均得到患者配对FDR支持 (约0.021与0.020)，但不能称每位患者或每条轴都增强。上述检验是观察性患者配对证据；即使显著也不能排除RNA复杂度、分组依赖性或共同细胞状态造成的关联。

### 3. 哪些伙伴变化最明显

incoming Other_fibroblasts: n=12, high=5.062, low=4.042, 配对中位差=0.7075, 95% CI [0.3379, 1.195], FDR=0.02147, high>low患者比例=91.7%。

incoming Pericyte: n=12, high=4.492, low=3.903, 配对中位差=0.6252, 95% CI [0.312, 1.018], FDR=0.02157, high>low患者比例=91.7%。

outgoing Pericyte: n=12, high=4.854, low=4.483, 配对中位差=0.2821, 95% CI [0.06658, 0.4214], FDR=0.02069, high>low患者比例=91.7%。

incoming Endothelial: n=12, high=2.442, low=2.259, 配对中位差=0.2674, 95% CI [0.1021, 0.4505], FDR=0.01974, high>low患者比例=91.7%。

outgoing Endothelial: n=12, high=3.828, low=3.756, 配对中位差=0.2303, 95% CI [0.07203, 0.3361], FDR=0.02069, high>low患者比例=75.0%。

outgoing Mast: n=8, high=3.202, low=2.974, 配对中位差=0.2073, 95% CI [0.1167, 0.4146], FDR=0.02069, high>low患者比例=100.0%。

按绝对配对强度变化展示，不代表跨不同伙伴分母的正式效应比较。

### 4. Endothelial是否稳定

incoming Endothelial: n=12, high=2.442, low=2.259, 配对中位差=0.2674, 95% CI [0.1021, 0.4505], FDR=0.01974, high>low患者比例=91.7%。

outgoing Endothelial: n=12, high=3.828, low=3.756, 配对中位差=0.2303, 95% CI [0.07203, 0.3361], FDR=0.02069, high>low患者比例=75.0%。

| stratum | direction | n_patients | difference_median | CI_low | CI_high | FDR |
| --- | --- | --- | --- | --- | --- | --- |
| dataset:GSE222315 | incoming | 5 | 0.07613 | -0.1127 | 0.4772 | 0.4679 |
| dataset:GSE222315 | outgoing | 5 | 0.223 | -0.08679 | 0.2899 | 0.2807 |
| dataset:InhouseData | incoming | 7 | 0.2838 | 0.2503 | 0.5298 | 0.04499 |
| dataset:InhouseData | outgoing | 7 | 0.2903 | 0.1527 | 0.478 | 0.04004 |
| exclude:GSE222315 | incoming | 7 | 0.2838 | 0.2503 | 0.5298 | 0.04499 |
| exclude:GSE222315 | outgoing | 7 | 0.2903 | 0.1527 | 0.478 | 0.04004 |
| exclude:InhouseData | incoming | 5 | 0.07613 | -0.1127 | 0.4772 | 0.4679 |
| exclude:InhouseData | outgoing | 5 | 0.223 | -0.08679 | 0.2899 | 0.2807 |

两个队列效应中位数同为正，但Endothelial双向总体差异仅InhouseData分别通过FDR约0.040和0.045，GSE222315分别约0.281和0.468；不能声称已获得两个队列的显著重复。需同时参考分层方向及患者数；小队列CI/P为NA不代表没有效应。患者间未拆分C5连续分析与本轮细胞内high/low并非同一问题，不能互相当作独立验证。

### 5. 通路、LR和ECM类别

| direction | partner | feature | n_patients | difference_median | CI_low | CI_high | FDR |
| --- | --- | --- | --- | --- | --- | --- | --- |
| incoming | Endothelial | Non-protein Signaling | 12 | 0.004579 | 0.0008722 | 0.00954 | 0.05638 |
| outgoing | Endothelial | Non-protein Signaling | 12 | 0 | -0.0002575 | 8.417e-06 | 0.6417 |
| incoming | Endothelial | Secreted Signaling | 12 | 0.02842 | 0.0009277 | 0.05657 | 0.08738 |
| outgoing | Endothelial | Secreted Signaling | 12 | 0.02156 | -0.0191 | 0.07437 | 0.4864 |
| incoming | Endothelial | ECM-Receptor | 12 | 0.1904 | 0.06857 | 0.3909 | 0.0478 |
| outgoing | Endothelial | ECM-Receptor | 12 | 0.1587 | -0.02809 | 0.2794 | 0.1299 |
| incoming | Endothelial | Cell-Cell Contact | 12 | 0.03366 | 0.00674 | 0.04615 | 0.0478 |
| outgoing | Endothelial | Cell-Cell Contact | 12 | 0.03222 | 0.01238 | 0.0438 | 0.07386 |
| incoming | TOTAL | Non-protein Signaling | 12 | 0.03501 | 0.01155 | 0.1105 | 0.04909 |
| outgoing | TOTAL | Non-protein Signaling | 12 | 0.003735 | -0.002043 | 0.01677 | 0.6076 |
| incoming | TOTAL | Secreted Signaling | 12 | 0.1381 | 0.03317 | 0.401 | 0.0478 |
| outgoing | TOTAL | Secreted Signaling | 12 | 0.0245 | -0.08013 | 0.3859 | 0.6076 |
| incoming | TOTAL | ECM-Receptor | 12 | 1.36 | 0.7448 | 2.608 | 0.0478 |
| outgoing | TOTAL | ECM-Receptor | 12 | 1.141 | 0.08839 | 2.131 | 0.07386 |
| incoming | TOTAL | Cell-Cell Contact | 12 | 0.1231 | -0.006793 | 0.1982 | 0.1437 |
| outgoing | TOTAL | Cell-Cell Contact | 12 | 0.1626 | 0.007036 | 0.3721 | 0.1033 |

| direction | feature | n_patients | positive_high_patients | positive_low_patients | difference_median | FDR |
| --- | --- | --- | --- | --- | --- | --- |
| incoming | COLLAGEN | 12 | 12 | 12 | 0.1425 | 1 |
| outgoing | COLLAGEN | 12 | 12 | 12 | 0.06623 | 1 |
| incoming | LAMININ | 12 | 12 | 12 | 0.05992 | 1 |
| outgoing | FN1 | 12 | 12 | 12 | 0.05377 | 1 |
| outgoing | SEMA3 | 12 | 12 | 11 | 0.02175 | 1 |
| outgoing | THBS | 12 | 11 | 11 | 0.01688 | 1 |
| outgoing | LAMININ | 12 | 12 | 12 | 0.01657 | 1 |
| outgoing | VEGF | 12 | 12 | 12 | 0.01646 | 1 |
| outgoing | TENASCIN | 12 | 11 | 9 | 0.01465 | 1 |
| outgoing | ANGPTL | 12 | 10 | 11 | 0.01042 | 1 |
| outgoing | MK | 12 | 11 | 12 | -0.009198 | 1 |
| incoming | JAM | 12 | 11 | 12 | 0.009174 | 1 |

| direction | feature | n_patients | positive_high_patients | positive_low_patients | difference_median | FDR |
| --- | --- | --- | --- | --- | --- | --- |
| incoming | COL4A1_ITGA11_ITGB1 | 12 | 12 | 12 | 0.01881 | 1 |
| outgoing | FN1_ITGA5_ITGB1 | 12 | 12 | 12 | 0.01843 | 1 |
| incoming | COL4A2_ITGA11_ITGB1 | 12 | 12 | 12 | 0.01711 | 1 |
| outgoing | FN1_ITGAV_ITGB1 | 12 | 12 | 12 | 0.01069 | 1 |
| outgoing | FN1_ITGA8_ITGB1 | 12 | 11 | 11 | 0.009673 | 1 |
| incoming | COL4A2_ITGA2_ITGB1 | 12 | 7 | 9 | -0.008579 | 1 |
| incoming | COL4A1_ITGA1_ITGB1 | 12 | 11 | 11 | 0.007214 | 1 |
| incoming | COL4A2_ITGA1_ITGB1 | 12 | 11 | 11 | 0.00699 | 1 |
| outgoing | ADM_CALCRL | 12 | 12 | 12 | -0.006107 | 1 |
| outgoing | COL1A1_ITGA1_ITGB1 | 12 | 12 | 12 | 0.005813 | 1 |
| outgoing | VEGFA_VEGFR1 | 12 | 12 | 10 | 0.005803 | 1 |
| incoming | COL4A1_CD44 | 12 | 10 | 10 | 0.00574 | 1 |

主方案没有单个LR或pathway通过全量BH-FDR<0.05，以下仅为效应和患者支持线索。ECM类别占较大的绝对强度变化，但只有incoming ECM类别在该类别检验集合中获FDR支持；outgoing ECM类别FDR约0.074，不能概括成所有ECM通讯显著增强。具体LR及通路的全量阴性/反向结果保留于统计表。新候选指本轮新的观察性线索，不声称新发现配体结合关系。

### 6. 哪些配体有患者级关联

患者C5至少10细胞的主模型包含15名患者。dataset调整后，HSPG2、FN1、COL1A1、POSTN、COL1A2、COL3A1、LAMA1呈正相关的重点候选FDR支持，WNT5A呈负相关。再加入患者平均RNA深度和基因复杂度后，仅HSPG2保持22个预设候选的FDR支持；在全部606个配体相关基因的探索性校正中，没有基因通过该深度调整模型的FDR<0.05。患者内部depth-adjusted关联支持POSTN、COL1A1、COL1A2、COL3A1、SEMA3C；各队列单独检验支持较弱，不能称为已跨队列验证。

### 7. 既有候选表现

POSTN的患者内调整相关中位数约0.164 (FDR约0.033)，患者间dataset调整关联FDR约0.008，但额外深度调整后不再支持。FN1患者间关联强，患者内调整相关未获FDR支持，paired pseudobulk high-low表达有重点候选FDR支持。INHBA、TGFB1、ANGPTL4没有建立重点候选FDR支持的连续LLPS关联；不能因旧NicheNet排名而称它们随LLPS普遍增强。POSTN–ITGAV/ITGB5在high和low均为7/12阳性，仅6/12患者high更高；HSPG2–DAG1也是两组均8/12阳性、仅6/12患者high更高。因此配体表达关联不能直接推出这两条内皮轴的功能增强。下表保留相关内皮LR的实际结果。

| ligand | feature | n_patients | positive_high_patients | positive_low_patients | difference_median | FDR |
| --- | --- | --- | --- | --- | --- | --- |
| ANGPTL4 | ANGPTL4_CDH11 | 12 | 0 | 0 | 0 | 1 |
| ANGPTL4 | ANGPTL4_CDH5 | 12 | 7 | 5 | 0.0005096 | 1 |
| ANGPTL4 | ANGPTL4_ITGA5_ITGB1 | 12 | 6 | 5 | 0 | 1 |
| ANGPTL4 | ANGPTL4_ITGAV_ITGB3 | 12 | 3 | 3 | 0 | 1 |
| ANGPTL4 | ANGPTL4_SDC1 | 12 | 0 | 0 | 0 | 1 |
| ANGPTL4 | ANGPTL4_SDC2 | 12 | 0 | 0 | 0 | 1 |
| ANGPTL4 | ANGPTL4_SDC3 | 12 | 5 | 3 | 0 | 1 |
| ANGPTL4 | ANGPTL4_SDC4 | 12 | 1 | 1 | 0 | 1 |
| FN1 | FN1_CD44 | 12 | 2 | 2 | 0 | 1 |
| FN1 | FN1_ITGA2B_ITGB3 | 12 | 0 | 0 | 0 | 1 |
| FN1 | FN1_ITGA3_ITGB1 | 12 | 8 | 8 | 0.0001909 | 1 |
| FN1 | FN1_ITGA4_ITGB1 | 12 | 0 | 0 | 0 | 1 |
| FN1 | FN1_ITGA4_ITGB7 | 12 | 0 | 0 | 0 | 1 |
| FN1 | FN1_ITGA5_ITGB1 | 12 | 12 | 12 | 0.01843 | 1 |
| FN1 | FN1_ITGA8_ITGB1 | 12 | 11 | 11 | 0.009673 | 1 |
| FN1 | FN1_ITGAV_ITGB1 | 12 | 12 | 12 | 0.01069 | 1 |
| FN1 | FN1_ITGAV_ITGB3 | 12 | 6 | 6 | 0 | 1 |
| FN1 | FN1_ITGAV_ITGB6 | 12 | 0 | 0 | 0 | 1 |
| FN1 | FN1_ITGAV_ITGB8 | 12 | 0 | 0 | 0 | 1 |
| FN1 | FN1_SDC1 | 12 | 0 | 0 | 0 | 1 |
| FN1 | FN1_SDC4 | 12 | 2 | 2 | 0 | 1 |
| INHBA | INHBA_ACVR1B_ACVR2A | 12 | 4 | 4 | 0 | 1 |
| INHBA | INHBA_ACVR1B_ACVR2B | 12 | 0 | 0 | 0 | 1 |
| POSTN | POSTN_ITGAV_ITGB3 | 12 | 6 | 6 | 0 | 1 |
| POSTN | POSTN_ITGAV_ITGB5 | 12 | 7 | 7 | 0.0001186 | 1 |
| TGFB1 | TGFB1_ACVR1B_TGFBR2 | 12 | 0 | 0 | 0 | 1 |
| TGFB1 | TGFB1_ACVR1C_TGFBR2 | 12 | 0 | 0 | 0 | 1 |
| TGFB1 | TGFB1_ACVR1_TGFBR1 | 12 | 3 | 2 | 0 | 1 |
| TGFB1 | TGFB1_TGFBR1_TGFBR2 | 12 | 0 | 0 | 0 | 1 |

### 8. 队列、细胞数与少数患者影响

主方案的12名肿瘤患者完全分布为InhouseData的7名MIBC和GSE222315的5名NMIBC，临床组与队列完全混杂，因此不能据此分离临床分期效应和队列效应。主方案两组细胞数最多相差1个，未额外进行不必要的等量下采样；全局中位数组成不同且可能不平衡，作为独立敏感性方案。high组平均nCount和nFeature偏高，不能将未经调整的表达或网络差异等同于LLPS独立效应。各阈值/分组及LODO见tables，通信逐患者留出范围见07_sensitivity_analysis/CellChat_network_leave_patient_out.csv。

### 9. 与既有NicheNet/空间证据一致和不一致的部分

ECM相关表达线索与旧基质程序证据在主题上相容，但不能证明来自同一C5细胞。POSTN在名义NicheNet Set B及多数LOPO中排名靠前，但该排名并不等同于FDR显著的实测内皮响应；也缺乏已验证C5-specific空间定位。INHBA/ANGPTL4的旧候选身份没有获得本轮连续LLPS关联支持。HSPG2的患者级关联值得关注，但旧空间调整后关联并非稳定显著。这些不一致均保留，不通过加权评分掩盖。

### 10. 优先验证的机制假说

优先研究“C5的LLPS相关转录状态与基质/配体表达程序是否共同变化，以及这些变化是否伴随特定内皮LR响应”。HSPG2适合验证深度调整后仍存在的患者间关联；POSTN及COL1A1/COL1A2/COL3A1/FN1适合以基质程序为背景检验患者内表达及具体受体轴，但必须控制RNA复杂度，COL1A1还需持续使用去基因评分。SEMA3C也适合作为补充的患者内状态候选：其调整后患者内表达关联有支持，SEMA3C–PLXND1在9/12患者呈high更高，但具体LR未通过全量FDR，只能探索验证。INHBA、TGFB1、ANGPTL4保留为旧候选对照，不因先验而提升优先级。最终需身份明确的原位测量、配体/受体阻断及LLPS相关扰动；现有数据不支持直接建立LLPS促进血管生成或EndMT的因果机制。


## 10. 稳健性与技术限制补充

排除中位数ties后主阈值结果方向与支持基本一致。阈值5覆盖14名患者，总outgoing和Endothelial outgoing仍获FDR支持，但两项incoming的FDR均约0.077；Top/Bottom30在阈值10仅7名患者，Endothelial outgoing/incoming FDR约0.075/0.187；全局中位数阈值10仅6名患者，对应约0.040/0.118。因此不能概括为所有分组和阈值都显著。20/30患者不足5名，仅描述。

主方案Endothelial outgoing配对差与组间log(nCount)及log(nFeature)差的相关系数分别约0.617、0.742；incoming约0.246、0.400。RNA复杂度混杂尤其影响outgoing解释。深度调整模型的零深度差参考截距属于小样本参考点估计，不是因果调整后的CellChat效应，不据此宣称效应消失或得到独立验证。逐患者留出后Endothelial outgoing中位差范围0.223–0.238，incoming 0.251–0.284，方向均为正；这降低单一患者驱动中位数的担忧，但不消除队列和技术混杂。

55个原生网络核验全部通过；最终1670项QC、132项独立pseudobulk核验及110项独立患者内关联核验通过。执行中断与恢复见00_input_audit/execution_recovery.md。原始对象、基因集与历史结果保持不变，未重跑NicheNet或空间分析。
