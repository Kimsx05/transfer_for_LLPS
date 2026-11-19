# C5 Five-gene 空间共定位分析 v1 (PySAL/esda)

**63 张 Visium、265,002 spots 的全量分析完成。Five-gene 为 MMP11、COL11A1、INHBA、COL5A1、CTHRC1；原有评分与解卷积只读复用。正式 C5-specific 定位 NO-GO 不变。**

本轮研究的是 Five-gene-associated program 与四个 TME 大类丰度的空间关联，不是已验证 C5 细胞的定位或通讯。所有主结果、反向、局部和敏感性均完整保存，未根据显著性选择切片。

## 对九个核心问题的直接回答

1. **最强正向原始关联是 Endothelial**：全63张中位I=0.086；若比较各类切片中位效应的绝对值，则Epithelial的负关联更大 (中位−0.171)。不能把“最强正向”与“绝对值最大”混为一谈。
2. **方向依赖调整和来源**：原始Endothelial/T_NK切片中位偏正，Epithelial/Myeloid偏负；full调整后四类切片中位均负。29名已确认患者中，Endothelial原始0.155 [0.079,0.279]，调整−0.041 [−0.082,−0.022]；Epithelial调整−0.103 [−0.189,−0.065]；Myeloid调整−0.029 [−0.063,−0.002]；T_NK调整区间跨零。这些CI是探索性且未作多套CI同时校正。
3–4. **局部HH可以反复观察，但不等于稳健C5 niche**：Endothelial/Epithelial/Myeloid/T_NK分别在63/60/62/63张有至少1个FDR-HH，均覆盖7个来源；连续至少10个HH spots的区域分别出现在49/26/37/40张。局部正向区域可以与全局负向结构并存。未进行63张×全部局部检验的整体atlas级FDR，重复热点数量不是独立患者证据。
5. **不支持普遍的“同spot内皮低、外围内皮高”模式**：同尺度Pearson负而Ring1 Moran正仅1/63张 (使用Spearman为2/63)。Five top10 spots中局部自身低/邻域高的配置比例切片中位约6.6%，存在但不是主导模式；调整后Endothelial在Ring1–3患者中位均偏负。
6. **组成/深度混杂影响实质性**：原始Endothelial正向排名在调整后消失，不能称为独立C5–血管空间富集。Epithelial调整中位在7来源均负，但删除GSE319536后患者CI跨零；Endothelial删除GSE238145或GSE319536后CI也跨零。Myeloid同样受主要来源影响。没有得到跨来源、条件正向邻近的稳健结论。
7. **算法一致性较高但不消除混杂**：UCell/AUCell的逐切片Moran效应相关约0.984–0.993，方向一致约93.7%–98.4%；Moran/Lee相关约0.997–0.999，方向一致98.4%–100%。近似平移零模型的504项结果全部未通过各自252项BH；99次平移分辨率有限，因此这一结果表示稳健支持不足，而非证明空间关联不存在。
8. **CellChat与原始空间结果仅部分一致**：Endothelial在四伙伴中C5 outgoing/incoming中位强度最高 (3.577/2.039，12名eligible患者)，也具有最高原始正Moran，但条件关联反转。Epithelial的C5−C0 outgoing增量最大 (+0.779，4名共同患者)，空间结果却主要为负，不能说其通讯增量已获得正向邻近验证。T_NK强度高不对应稳定条件空间支持；不同伙伴分母不可直接比较偏好。
9. **当前没有伙伴满足“已获得稳健正向空间支持”的LR验证优先标准**。若资源只允许选择一个进行独立验证，Endothelial可作为条件性优先候选：它同时有较充分CellChat支持和原始空间信号，但首先要用独立细胞身份/多重原位检测检验真实邻接及组成混杂，再决定具体LR实验。Epithelial可作为基质–上皮空间分区的次级问题，不能包装成已证实通讯；T_NK不得下推到Treg/CD8。

**证据结论：存在可描述的Five相关空间结构；不足以证明C5特异细胞共定位或真实通讯。后续应优先独立验证细胞身份及组织边界，而不是从本轮相关性直接指定因果配体–受体轴。**

## 核心结果

```text
partner	section_median_Moran	positive_sections	positive_global_FDR_sections	negative_global_FDR_sections	HH_any_sections	HH_ge10_sections	HH_contiguous_ge10_sections	HH_any_datasets	positive_dataset_medians	primary_median	primary_ci_low	primary_ci_high	primary_n_patients	adjusted_full_median	adjusted_full_ci_low	adjusted_full_ci_high	adjusted_full_n_patients	outgoing_CellChat_strength	outgoing_eligible_patients	outgoing_C5_minus_C0	outgoing_paired_patients	incoming_CellChat_strength	incoming_eligible_patients	incoming_C5_minus_C0	incoming_paired_patients
Endothelial	0.08594	38	36	24	63	55	49	7	4	0.1546	0.07874	0.2794	29	-0.04087	-0.08164	-0.02207	29	3.577	12	0.4419	5	2.039	12	0.218	5
Epithelial	-0.1709	15	14	47	60	41	26	7	2	-0.03928	-0.1692	0.105	29	-0.1025	-0.1893	-0.06458	29	2.71	9	0.779	4	0.7518	9	0.2344	4
Myeloid	-0.05138	29	26	33	62	49	37	7	3	0.06659	-0.0561	0.191	29	-0.02857	-0.06264	-0.002191	29	2.8	13	0.4638	6	0.5422	13	0.07861	6
T_NK	0.0159	33	31	25	63	56	40	7	3	0.1303	-0.01148	0.2587	29	-0.0253	-0.06123	0.01753	29	3.231	12	0.3593	6	0.3663	12	0.06174	6
```

## 1. 哪个伙伴最强，正负方向如何？

以下第一组是全部63张切片的描述；第二组只使用29名已确认患者 (30张切片、3个来源)，同患者切片先平均。未知患者的33张不作为独立患者。原始强度排名不能替代调整后稳定性。
```text
partner	median_I	min_I	max_I	positive	negative
Endothelial	0.08594	-0.4277	0.5973	38	25
Epithelial	-0.1709	-0.5272	0.4992	15	48
Myeloid	-0.05138	-0.5825	0.6488	29	34
T_NK	0.0159	-0.5631	0.6247	33	30
```
```text
partner	method	n_patients	median	ci_low	ci_high	positive_patients	negative_patients	n_datasets
Endothelial	adjusted_full	29	-0.04087	-0.08164	-0.02207	9	20	3
Endothelial	primary	29	0.1546	0.07874	0.2794	23	6	3
Epithelial	adjusted_full	29	-0.1025	-0.1893	-0.06458	7	22	3
Epithelial	primary	29	-0.03928	-0.1692	0.105	13	16	3
Myeloid	adjusted_full	29	-0.02857	-0.06264	-0.002191	8	21	3
Myeloid	primary	29	0.06659	-0.0561	0.191	19	10	3
T_NK	adjusted_full	29	-0.0253	-0.06123	0.01753	12	17	3
T_NK	primary	29	0.1303	-0.01148	0.2587	19	10	3
```

七个来源分组如下；UNVERIFIED来源名按原审计保留，不能宣称都是已验证独立研究队列。
```text
dataset	partner	method	median	minimum	maximum	n_sections	positive_sections	negative_sections
CanDis__UNVERIFIED	Endothelial	adjusted_full	-0.04697	-0.1756	0.06096	13	5	8
CanDis__UNVERIFIED	Endothelial	primary	-0.09373	-0.3281	0.5543	13	5	8
CanDis__UNVERIFIED	Epithelial	adjusted_full	-0.1941	-0.4188	0.02468	13	1	12
CanDis__UNVERIFIED	Epithelial	primary	-0.2651	-0.5272	-0.04056	13	0	13
CanDis__UNVERIFIED	Myeloid	adjusted_full	-0.01429	-0.06886	0.1034	13	6	7
CanDis__UNVERIFIED	Myeloid	primary	-0.1949	-0.4223	0.6314	13	5	8
CanDis__UNVERIFIED	T_NK	adjusted_full	-0.003171	-0.1073	0.08541	13	6	7
CanDis__UNVERIFIED	T_NK	primary	-0.0708	-0.4059	0.5972	13	6	7
GSE171351	Endothelial	adjusted_full	0.03449	0.005684	0.05952	4	4	0
GSE171351	Endothelial	primary	-0.0003847	-0.11	0.2307	4	2	2
GSE171351	Epithelial	adjusted_full	-0.01419	-0.08924	0.146	4	2	2
GSE171351	Epithelial	primary	0.03853	-0.4035	0.2151	4	2	2
GSE171351	Myeloid	adjusted_full	-0.002475	-0.02857	0.04728	4	2	2
GSE171351	Myeloid	primary	0.02842	-0.08127	0.1702	4	3	1
GSE171351	T_NK	adjusted_full	-0.02306	-0.06123	0.00797	4	1	3
GSE171351	T_NK	primary	-0.009493	-0.02291	0.1613	4	2	2
GSE238145	Endothelial	adjusted_full	-0.09515	-0.2696	0.06053	10	3	7
GSE238145	Endothelial	primary	0.2432	-0.1794	0.4491	10	7	3
GSE238145	Epithelial	adjusted_full	-0.02726	-0.3336	0.2421	10	5	5
GSE238145	Epithelial	primary	0.1577	-0.2177	0.4992	10	6	4
GSE238145	Myeloid	adjusted_full	-0.03866	-0.2336	0.2331	10	4	6
GSE238145	Myeloid	primary	0.1837	-0.2175	0.6488	10	6	4
GSE238145	T_NK	adjusted_full	0.05946	-0.1614	0.1456	10	8	2
GSE238145	T_NK	primary	0.2967	-0.1165	0.6247	10	7	3
GSE285715	Endothelial	adjusted_full	-0.01072	-0.1205	0.01002	6	2	4
GSE285715	Endothelial	primary	0.1051	-0.1841	0.2962	6	4	2
GSE285715	Epithelial	adjusted_full	-0.1716	-0.3226	-0.03274	6	0	6
GSE285715	Epithelial	primary	-0.1623	-0.485	0.007178	6	1	5
GSE285715	Myeloid	adjusted_full	0.03875	-0.0465	0.1203	6	5	1
GSE285715	Myeloid	primary	-0.04465	-0.4534	0.2726	6	3	3
GSE285715	T_NK	adjusted_full	0.03169	-0.1114	0.08193	6	4	2
GSE285715	T_NK	primary	-0.03232	-0.1731	0.3403	6	3	3
GSE319536	Endothelial	adjusted_full	-0.04137	-0.2615	0.04859	16	3	13
GSE319536	Endothelial	primary	0.1511	-0.1854	0.5973	16	14	2
GSE319536	Epithelial	adjusted_full	-0.2109	-0.4812	0.01768	16	1	15
GSE319536	Epithelial	primary	-0.1356	-0.3181	0.3793	16	5	11
GSE319536	Myeloid	adjusted_full	-0.04525	-0.2447	0.094	16	3	13
GSE319536	Myeloid	primary	0.05746	-0.2246	0.6216	16	10	6
GSE319536	T_NK	adjusted_full	-0.06925	-0.2928	0.07382	16	3	13
GSE319536	T_NK	primary	0.1132	-0.2295	0.6219	16	10	6
Inhousedata__UNVERIFIED	Endothelial	adjusted_full	-0.05377	-0.2364	-0.02727	7	0	7
Inhousedata__UNVERIFIED	Endothelial	primary	0.09746	-0.3026	0.3171	7	5	2
Inhousedata__UNVERIFIED	Epithelial	adjusted_full	-0.2429	-0.5203	-0.09626	7	0	7
Inhousedata__UNVERIFIED	Epithelial	primary	-0.1869	-0.3543	0.2801	7	1	6
Inhousedata__UNVERIFIED	Myeloid	adjusted_full	0.01777	-0.2179	0.1235	7	4	3
Inhousedata__UNVERIFIED	Myeloid	primary	-0.07166	-0.3762	0.3827	7	2	5
Inhousedata__UNVERIFIED	T_NK	adjusted_full	-0.08435	-0.2191	0.2612	7	1	6
Inhousedata__UNVERIFIED	T_NK	primary	0.06397	-0.4056	0.4829	7	4	3
NPJ__UNVERIFIED	Endothelial	adjusted_full	-0.0524	-0.09642	-0.007819	7	0	7
NPJ__UNVERIFIED	Endothelial	primary	-0.1635	-0.4277	0.1973	7	1	6
NPJ__UNVERIFIED	Epithelial	adjusted_full	-0.2771	-0.4674	0.1567	7	1	6
NPJ__UNVERIFIED	Epithelial	primary	-0.3855	-0.5269	-0.08073	7	0	7
NPJ__UNVERIFIED	Myeloid	adjusted_full	-0.03044	-0.1256	0.02987	7	2	5
NPJ__UNVERIFIED	Myeloid	primary	-0.3625	-0.5825	-0.05138	7	0	7
NPJ__UNVERIFIED	T_NK	adjusted_full	-0.04666	-0.1145	0.02756	7	1	6
NPJ__UNVERIFIED	T_NK	primary	-0.3997	-0.5631	0.2804	7	1	6
```

## 2. Local HH 支持与重复性

HH=本切片 Five UCell 高于其均值，邻近目标丰度的标准化空间滞后高于零。不是top10%，也不是统一绝对细胞量阈值。LH/HL含义随X=Five、Y=目标固定；所有局部结果为正向定义。局部BH对每切片四目标×可用spots联合实施，使用2×原生min-tail P (上限1)。
```text
partner	dataset	n_sections	sections_any_HH	sections_ge10_HH	median_HH_fraction	total_HH_spots	sections_contiguous_HH_ge10
Endothelial	CanDis__UNVERIFIED	13	13	10	0.01771	1738	9
Endothelial	GSE171351	4	4	3	0.01824	87	1
Endothelial	GSE238145	10	10	10	0.05515	3334	10
Endothelial	GSE285715	6	6	5	0.01229	207	4
Endothelial	GSE319536	16	16	16	0.04865	3275	16
Endothelial	Inhousedata__UNVERIFIED	7	7	7	0.02601	1166	6
Endothelial	NPJ__UNVERIFIED	7	7	4	0.007229	377	3
Epithelial	CanDis__UNVERIFIED	13	12	9	0.008657	659	4
Epithelial	GSE171351	4	4	2	0.0124	75	1
Epithelial	GSE238145	10	10	10	0.0416	2588	9
Epithelial	GSE285715	6	6	3	0.003678	50	1
Epithelial	GSE319536	16	16	12	0.007193	1115	8
Epithelial	Inhousedata__UNVERIFIED	7	7	3	0.002058	436	2
Epithelial	NPJ__UNVERIFIED	7	5	2	0.0004819	87	1
Myeloid	CanDis__UNVERIFIED	13	12	10	0.01468	1552	8
Myeloid	GSE171351	4	4	2	0.01905	94	2
Myeloid	GSE238145	10	10	10	0.0518	3278	9
Myeloid	GSE285715	6	6	5	0.01113	162	3
Myeloid	GSE319536	16	16	15	0.04231	2322	10
Myeloid	Inhousedata__UNVERIFIED	7	7	6	0.02302	1126	4
Myeloid	NPJ__UNVERIFIED	7	7	1	0.0005928	86	1
T_NK	CanDis__UNVERIFIED	13	13	11	0.02881	1692	8
T_NK	GSE171351	4	4	4	0.01643	90	1
T_NK	GSE238145	10	10	10	0.06322	4075	10
T_NK	GSE285715	6	6	5	0.0103	184	3
T_NK	GSE319536	16	16	16	0.04221	2981	12
T_NK	Inhousedata__UNVERIFIED	7	7	6	0.01207	923	4
T_NK	NPJ__UNVERIFIED	7	7	4	0.00462	294	2
```

“至少1个HH”和“至少10个HH”均为描述性计数，后者不是新增显著性检验。HH存在不保证该切片全局正相关，也不保证控制组织组成后仍然存在。全量FDR地图与未筛选quadrant地图分别保存，未校正quadrant图明确不表示显著区域。

## 3. Endothelial外围邻近是否得到支持？

严格区分同spot与相邻spot。下表统计全切片同spot负相关但Ring1 Moran正向的数量，并给top10 Five spots中“自身Endo低于均值、邻近Endo高于均值”的比例。Pearson与Moran使用同一原始尺度；Spearman方向相反还可能包含变换尺度差异。个别局部配置不等于全体切片稳定模式。
```text
partner	section_Spearman_negative_neighbor_positive	section_Pearson_negative_neighbor_positive	median_high_peripheral_fraction
Endothelial	2	1	0.06589
Epithelial	1	1	0.05124
Myeloid	2	1	0.04834
T_NK	1	2	0.04803
```

精确Ring1–3的患者结果：
```text
partner	method	ring	n_patients	median	ci_low	ci_high	positive_patients	negative_patients	n_datasets
Endothelial	adjusted	1	29	-0.04087	-0.08164	-0.02207	9	20	3
Endothelial	adjusted	2	29	-0.02884	-0.05837	-0.002014	9	20	3
Endothelial	adjusted	3	29	-0.01814	-0.05008	-0.005052	9	20	3
Endothelial	raw	1	29	0.1546	0.07874	0.2794	23	6	3
Endothelial	raw	2	29	0.1332	0.05888	0.2471	22	7	3
Endothelial	raw	3	29	0.1249	0.05833	0.2385	22	7	3
Epithelial	adjusted	1	29	-0.1025	-0.1893	-0.06458	7	22	3
Epithelial	adjusted	2	29	-0.06492	-0.1584	-0.03336	8	21	3
Epithelial	adjusted	3	29	-0.05045	-0.1391	-0.01737	8	21	3
Epithelial	raw	1	29	-0.03928	-0.1692	0.105	13	16	3
Epithelial	raw	2	29	0.005635	-0.1463	0.1167	15	14	3
Epithelial	raw	3	29	0.02271	-0.1191	0.1264	15	14	3
Myeloid	adjusted	1	29	-0.02857	-0.06264	-0.002191	8	21	3
Myeloid	adjusted	2	29	-0.02181	-0.06185	-0.006731	8	21	3
Myeloid	adjusted	3	29	-0.01603	-0.05874	-0.007148	9	20	3
Myeloid	raw	1	29	0.06659	-0.0561	0.191	19	10	3
Myeloid	raw	2	29	0.05113	-0.05108	0.1845	19	10	3
Myeloid	raw	3	29	0.04967	-0.05291	0.1885	19	10	3
T_NK	adjusted	1	29	-0.0253	-0.06123	0.01753	12	17	3
T_NK	adjusted	2	29	-0.00833	-0.06523	0.004069	11	18	3
T_NK	adjusted	3	29	-0.01038	-0.05286	0.01294	13	16	3
T_NK	raw	1	29	0.1303	-0.01148	0.2587	19	10	3
T_NK	raw	2	29	0.1266	-0.01079	0.2356	18	11	3
T_NK	raw	3	29	0.1194	-0.005551	0.2306	19	10	3
```

旧分析的匹配后同spot丰度差与本次原始Moran不是同一估计量，不能用二者符号直接证明“外围内皮”。本轮Ring0 Pearson/Spearman与Ring1–3同尺度结果见tables/ring0_to_ring123_comparison.tsv。

## 4. 组成、深度与方法敏感性

秩残差敏感性：在每切片内，将Five及目标丰度分别对秩转换后的Fibroblast、Mural、nCount、nFeature做含截距OLS；full额外控制Epithelial，但目标是Epithelial时不控制自身。将残差送入同一个esda.Moran_BV。base/full分开保存。reverse仅表示Moran两列变量顺序互换，不能映射CellChat的发送/接收方向。它控制有限线性秩关系，不是空间因果模型；残差负相关也不证明细胞排斥。
```text
partner	n_patients	median	ci_low	ci_high	positive_patients	negative_patients	n_datasets
Endothelial	29	-0.2318	-0.3702	-0.1434	7	22	3
Epithelial	29	-0.06926	-0.1445	-0.02006	6	23	3
Myeloid	29	-0.1749	-0.2307	0.006138	10	19	3
T_NK	29	-0.1704	-0.315	-0.03964	7	22	3
```

方法一致性 (与原始UCell/mean相比)：
```text
partner	method	spearman_vs_primary	sign_agreement	median_difference
Endothelial	reverse	0.9998	1	-7.248e-06
Endothelial	AUCell	0.9859	0.9365	-0.01239
Endothelial	log1p_abundance	0.9878	0.9524	0.01588
Endothelial	major_fraction	0.1454	0.5714	0.03773
Endothelial	q05	0.9954	1	-0.002186
Endothelial	adjusted_base	0.3037	0.4444	-0.2121
Endothelial	adjusted_full	0.04402	0.4444	-0.1639
Endothelial	A26	0.9363	0.9048	0.06231
Endothelial	Core15	0.9281	0.9206	0.0408
Epithelial	reverse	0.9995	1	-6.412e-05
Epithelial	AUCell	0.9844	0.9683	0.001652
Epithelial	log1p_abundance	0.9494	0.9365	-0.01826
Epithelial	major_fraction	0.795	0.9048	-0.09291
Epithelial	q05	0.9999	1	-7.604e-05
Epithelial	adjusted_base	0.3707	0.8254	-0.03608
Epithelial	adjusted_full	0.3707	0.8254	-0.03608
Epithelial	A26	0.9106	0.9524	-0.02223
Epithelial	Core15	0.9226	0.8413	-0.05985
Myeloid	reverse	0.9998	1	0.0002525
Myeloid	AUCell	0.9886	0.9683	-0.00727
Myeloid	log1p_abundance	0.9792	0.9683	0.007553
Myeloid	major_fraction	0.0949	0.5238	0.04813
Myeloid	q05	0.9991	1	0.001261
Myeloid	adjusted_base	0.451	0.6508	-0.08975
Myeloid	adjusted_full	0.07546	0.5397	-0.002723
Myeloid	A26	0.9317	0.9365	0.01951
Myeloid	Core15	0.9349	0.8413	0.01266
T_NK	reverse	0.9998	0.9841	0.0001089
T_NK	AUCell	0.9931	0.9841	-0.009088
T_NK	log1p_abundance	0.9874	0.9048	0.01337
T_NK	major_fraction	-0.1365	0.4286	0.02363
T_NK	q05	0.9986	0.9841	-0.00133
T_NK	adjusted_base	0.5385	0.6349	-0.09178
T_NK	adjusted_full	0.3129	0.6349	-0.03964
T_NK	A26	0.9522	0.9048	0.03929
T_NK	Core15	0.9489	0.873	0.01471
```

q05为原posterior下界辅助值，不能把q05>0等同细胞身份已验证；九大类相对丰度为组成敏感性，闭合约束可引入负关联。Fibroblast与Mural始终分开，T_NK只用大类。

## 5. Moran 与 Lee's L
```text
partner	spearman	sign_agreement
Endothelial	0.9992	0.9841
Epithelial	0.9969	1
Myeloid	0.9988	0.9841
T_NK	0.9993	0.9841
```

Lee.mc使用相同无对角线、row-standardized权重及999次置换。其联合置换(x,y)行的零假设不同于Moran_BV单独置换Y；Local Moran又是条件置换。因此只比较效应方向/模式，不要求相同数值或P值。spdep原生双侧P可在极端秩返回0，保留但不称为零概率。
```text
partner	method	n_patients	median	ci_low	ci_high	positive_patients	negative_patients	n_datasets
Endothelial	adjusted	29	-0.03555	-0.07155	-0.01642	8	21	3
Endothelial	raw	29	0.1477	0.07288	0.2672	22	7	3
Epithelial	adjusted	29	-0.1018	-0.1878	-0.05542	7	22	3
Epithelial	raw	29	-0.02915	-0.1617	0.1036	13	16	3
Myeloid	adjusted	29	-0.02085	-0.06493	-0.00344	9	20	3
Myeloid	raw	29	0.06542	-0.0573	0.1865	19	10	3
T_NK	adjusted	29	-0.02346	-0.06845	0.0144	12	17	3
T_NK	raw	29	0.1282	-0.0129	0.2507	18	11	3
```

## 6. 保留部分空间结构的平移敏感性

每切片使用seed123随机挑选最多99个非零、六角阵列同奇偶的平移，至少保留50%节点重叠；不环绕补边、不跨缺失spot补邻居。每次对重叠诱导图分别调用esda.Moran_BV计算原始与平移I，以同一重叠的差值构成近似双尾概率。边界/可用图改变，因此这不是精确随机化检验；99次的最小双侧P约0.02，BH判断分辨率有限；不能把未通过FDR等同于不存在空间关联。
```text
partner	method	median_translations	min_translations	approx_FDR_supported_sections
Endothelial	adjusted	99	99	0
Endothelial	raw	99	99	0
Epithelial	adjusted	99	99	0
Epithelial	raw	99	99	0
Myeloid	adjusted	99	99	0
Myeloid	raw	99	99	0
T_NK	adjusted	99	99	0
T_NK	raw	99	99	0
```

## 7. 患者、数据集与留一来源

患者CI为按已知dataset分层、患者单位有放回重采样2000次的中位数百分位95% CI，seed123。它不是所有七来源的总体人群CI。CI未对多套敏感性同步校正。删除无已知患者的来源不会改变患者CI，表中removed_confirmed_patients明确为0；另保存全63切片的描述性LODO。
```text
partner	method	omitted_dataset	removed_confirmed_patients	n_patients	median	ci_low	ci_high	positive_patients	negative_patients	n_datasets
Endothelial	adjusted_full	GSE171351	4	25	-0.04245	-0.1025	-0.02991	5	20	2
Endothelial	adjusted_full	GSE238145	9	20	-0.03356	-0.04816	0.003965	7	13	2
Endothelial	adjusted_full	GSE319536	16	13	-0.02207	-0.1025	0.03371	6	7	2
Endothelial	primary	GSE171351	4	25	0.1603	0.08081	0.3332	21	4	2
Endothelial	primary	GSE238145	9	20	0.138	0.03022	0.2551	16	4	2
Endothelial	primary	GSE319536	16	13	0.1546	-0.005177	0.3235	9	4	2
Epithelial	adjusted_full	GSE171351	4	25	-0.1436	-0.2202	-0.07213	5	20	2
Epithelial	adjusted_full	GSE238145	9	20	-0.1712	-0.2741	-0.07213	3	17	2
Epithelial	adjusted_full	GSE319536	16	13	-0.04575	-0.1025	0.06586	6	7	2
Epithelial	primary	GSE171351	4	25	-0.04794	-0.1712	0.105	11	14	2
Epithelial	primary	GSE238145	9	20	-0.0898	-0.1757	0.06811	7	13	2
Epithelial	primary	GSE319536	16	13	0.105	-0.03928	0.2151	8	5	2
Myeloid	adjusted_full	GSE171351	4	25	-0.03249	-0.09463	-0.002191	6	19	2
Myeloid	adjusted_full	GSE238145	9	20	-0.03011	-0.06356	-0.004865	5	15	2
Myeloid	adjusted_full	GSE319536	16	13	-0.007539	-0.06132	0.006001	5	8	2
Myeloid	primary	GSE171351	4	25	0.1178	-0.07144	0.2389	16	9	2
Myeloid	primary	GSE238145	9	20	0.03862	-0.07636	0.1531	13	7	2
Myeloid	primary	GSE319536	16	13	0.1268	-0.0561	0.2406	9	4	2
T_NK	adjusted_full	GSE171351	4	25	-0.0253	-0.1041	0.04377	11	14	2
T_NK	adjusted_full	GSE238145	9	20	-0.0357	-0.1086	-0.023	4	16	2
T_NK	adjusted_full	GSE319536	16	13	0.03109	-0.009082	0.07432	9	4	2
T_NK	primary	GSE171351	4	25	0.2099	-0.01148	0.3777	17	8	2
T_NK	primary	GSE238145	9	20	0.09069	-0.02225	0.2343	12	8	2
T_NK	primary	GSE319536	16	13	0.1613	-0.01148	0.3558	9	4	2
```

组织标签仅3 Normal、17 Tumor、43 UNRESOLVED；Normal均来自单一来源且患者未知。本轮只保存来源×可靠组织的描述，不做泛化Normal/Tumor差异推断。

## 8. 与 CellChat 的一致性及下一步

CellChat为既有min.cells20 tumor-only结果，10/30完整保留；其样本门槛实际为群细胞数>min.cells。不同partner eligible患者分母不等；强度不能直接当作偏好。C5 vs C0在共同eligible样本计算配对差，再按患者×临床条件中位汇总，保留独立患者和患者-条件数。单细胞与Visium不做患者一一匹配。
```text
min_cells	partner	direction	paired_median_difference	patient_condition_units	independent_patients	datasets	positive_units
10	Endothelial	incoming	0.1933	8	8	2	6
10	Endothelial	outgoing	0.5329	8	8	2	8
10	Epithelial	incoming	0.1816	6	6	1	6
10	Epithelial	outgoing	0.6517	6	6	1	6
10	Myeloid	incoming	0.07503	8	8	2	7
10	Myeloid	outgoing	0.5124	8	8	2	8
10	T_NK	incoming	0.09656	8	8	2	7
10	T_NK	outgoing	0.4077	8	8	2	8
20	Endothelial	incoming	0.218	5	5	2	4
20	Endothelial	outgoing	0.4419	5	5	2	5
20	Epithelial	incoming	0.2344	4	4	1	4
20	Epithelial	outgoing	0.779	4	4	1	4
20	Myeloid	incoming	0.07861	6	6	2	5
20	Myeloid	outgoing	0.4638	6	6	2	6
20	T_NK	incoming	0.06174	6	6	2	5
20	T_NK	outgoing	0.3593	6	6	2	6
30	Endothelial	incoming	0.218	5	5	2	4
30	Endothelial	outgoing	0.4419	5	5	2	5
30	Epithelial	incoming	0.2297	3	3	1	3
30	Epithelial	outgoing	0.9514	3	3	1	3
30	Myeloid	incoming	0.1406	5	5	2	4
30	Myeloid	outgoing	0.4959	5	5	2	5
30	T_NK	incoming	0.02913	5	5	2	4
30	T_NK	outgoing	0.4588	5	5	2	5
```

本轮没有新的配体–受体推断、空间LR表达计算、CellChat、NicheNet或解卷积训练。后续具体伙伴应依据调整后患者CI、来源一致性、尺度稳定性与独立验证可行性综合选择，不能只据CellChat强度或HH面积排名。

## 统计与技术核验

- 实际输入 265002 spots；参与最近邻图统计 264497，原始孤立节点 505，缺失必需字段 0。完整63张不删除；排除spot保持spot_id并标记NA/Excluded。
- 六邻接仅允许(0,±2),(±1,±1)，有实际spot才连接；不连通块原样保留。精确Ring2/3为最短路径恰好2/3步，不是累计邻域，不换算微米。
- esda 2.7.0 / libpysal 4.13.0 / numba 0.61.2；R spdep实际版本及完整session见00_input_audit。Moran统计全部调用成熟包，没有替代自写公式。
- 全局主BH固定252个Five→目标检验；反向及每个其他方法分别252个集合。Ring每个方法×尺度分别252个。Lee每个方法分别252个。局部每切片四目标×有效spots为一个集合。BH均为探索性名义控制，不代表患者级独立显著性。
- 所有Moran permutations=999；global每次np.random.seed(123)，local函数内seed=123/n_jobs=1/keep_simulations=False。局部不保存大型置换数组；p_sim与双侧转换均保留。全组织随机置换破坏自相关，不能充分排除共同组织结构，需结合调整和平移结果。
- 原七张pilot依原先UMI接近来源中位数原则选定；技术检查通过后扩展63，不设显著性门槛。重复一张local fit逐项验证P/quadrant可复现，最大内存与耗时逐切片记录。
- 已复用CellChat的四类细胞配色。空间图Five用0–1；目标丰度及|localI|按全量共同99分位上限显示并注明截色。原始坐标上原点方向保留，无HE覆盖。另附明确标注的切片内丰度百分位图，用于观察低丰度切片内部结构，不能作为跨切片绝对丰度比较。

## 证据边界

- **空间统计事实**：完整固定输入、全量逐切片效应/局部象限与FDR、图结构及跨方法结果。
- **探索性邻近关联**：调整前后、患者CI、来源及Ring共同支持的模式；残余混杂与空间零模型限制仍存在。
- **CellChat–空间一致性**：不同数据来源两类预测/描述证据并列，不能据此证明同一细胞对发生通讯。
- **不能确认**：Five-high=真实C5、T_NK具体亚群、Epithelial=恶性上皮、因果通讯或促血管/免疫抑制机制。

## 文件入口

- GALLERY.html：全部63张四伙伴局部图与原七张代表切片。
- tables/global_results.tsv：完整2520项主/敏感性全局结果；primary严格252项。
- 03_local_bivariate_moran/*.parquet：spot_id、localI、q、原生P、双侧P、FDR与分类。
- 01_spatial_weights：节点顺序、无向边、图QC，可重建完全相同的权重。
- 06_patient_dataset_summary：患者、来源、LODO、CI与组织描述。
- 07_cellchat_integration：已有CellChat与空间对照，保留全部分母。
- scripts：复现与续跑；00_input_audit：参数、环境、来源哈希及最终检查。