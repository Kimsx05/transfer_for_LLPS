# C0 去卷积诊断、优化及外部临床分析

本轮已完成诊断、三个模型比较、signature 验证和四个队列的正式分期/KM/Cox 分析。**最终临床结果是 C0-like signature score 的探索性关联，不是 C0 细胞丰度。** 未修改 C0–C6、旧 LLPS77_AUC、原有图或配色；没有重新下载队列。

## 主要发现

原 MuSiC 的总体 Spearman 0.851 并非仅由 dataset 间差异造成：控制 dataset 后秩相关 0.813，dataset 内置换 P=0.00020。dataset 对真实/预测全局秩的解释度分别为 41.8%/52.9%，说明 dataset 贡献明显，但仍有组内排序信号。原 20 人的误差及 2000 次患者 bootstrap CI 完整保留。

比例误差涉及整个组成体系：C1 MAE=0.205、bias=−0.122；C0 MAE=0.174、bias=+0.045；C4 MAE=0.172、bias=+0.067。C0 错误与 C1/C4 的相关为 −0.307/−0.280，与 C5 为 −0.041，不能简单归结为 C0/C5 相互误分。参考 log-CPM 相关中 C0–C4=0.942，C0–C5=0.917；全表达相关受共享基因影响，不等同于真实错分证据。组成约束本身也会诱发误差负相关。

MuSiC 中 D=mean relative expression × mean UMI/cell，拟合系数在该假设下是 cell-fraction estimates；乘以 cell-size factor 后才对应模型 RNA contribution。UMI/cell 不是独立测得的细胞 RNA 总量，既有生物学差异也有技术捕获影响。显式比较真实细胞/RNA 贡献后，CNP0000460 的 RNA-fraction MAE 仍达 0.341，因此误差不能只靠把细胞比例改称 RNA 比例解释。

C6 在分析纳入细胞中有 399 个，其中 330 来自 GSE222315；21 位患者有捕获，但仅 4 位达到 20 个 C6 细胞。原失败 fold 是 **留出 GSE222315**，不是按代码排序的“第二个”；该 fold 其他数据集只有 1 位 >=20 细胞的训练患者。降低到 5 细胞增加全体覆盖至 11 人，但严格保留原 20 人后，剩余 C6/C5 训练主要依赖 GSE222315，留出该 dataset 仍无法辨识完整参考。

## 模型比较与边界

| 方法 | 原 20 人有效预测 | Spearman | MAE | RMSE | signed bias |
|---|---:|---:|---:|---:|---:|
| 历史 baseline（缓存） | 20 | 0.851 | 0.174 | 0.251 | +0.045 |
| 重拟合 baseline | 20 | 0.405 | 0.326 | 0.434 | +0.168 |
| 优化参考 MuSiC | 20 | 0.908 | 0.166 | 0.226 | −0.156 |
| 自定义两级 MuSiC | 20 | 0.886 | 0.170 | 0.236 | −0.166 |

所有新模型的训练均全局排除本次任务指定的原 20 名主评估患者，再按 dataset 留出；历史 baseline 的其他折训练池不同，因此不能把历史/新基线的差异只归因于算法。剩余 23 人是开发训练池，其中包含上一轮曾作次级低覆盖诊断预测的患者；他们不构成本轮独立验证。原 20 人已经被观察，模型比较和 signature 表现明确属于开发/评估，不能宣称新盲法验证。9 名额外跨 dataset 评估患者全部来自 GSE222315，其不支持的 fold 没有被悄悄删除，状态和分母完整记录。

优化方案：>=5 细胞/训练患者/类型、每类型至少 3 人，dataset 等权的均值/方差，限定共同芯片覆盖的 protein-coding genes，排除 MT/RPL/RPS，用训练数据的细胞类型差异相对患者变异选每类型 top 100 informative genes。实际使用约 1014–1045 基因，condition number 从重拟合基线 86–144 降至 21–33。未用任何外部分期/生存 P 值筛选。

原生 [music_prop.cluster](https://xuranw.github.io/MuSiC/reference/music_prop.cluster.html) 已检查；当前保存实现对本任务 14 个 fine types / 8 个 coarse types 的矩阵乘法维度不匹配，审查记录见 native_hierarchy_audit.tsv。采用明确的两级流程：先估计 Stroma+7 主要类型，再在同一线性 RNA 尺度扣除非 stromal 拟合贡献，正残差拟合 C0–C6，最终乘回母群比例。未把该实现称为原生 MuSiC，负残差裁剪比例已保存。两级方案没有优于优化参考。

优化参考是在本开发集合上的最佳去卷积候选，但各 dataset 的偏差并非已证明恒定尺度偏移（约 −0.315 至 −0.014），额外 dataset 无法验证，实际芯片丰度校准缺失。因此没有将其用于正式 abundance-based 临床结论，也没有再无界搜索算法；转入已授权的 expression signature 路线。

## 固定 signature 与验证

训练仅用非预留患者；原始 pseudobulk DEG 仅作为结果注释，不能参与筛选，因为它包含本轮评估患者。固定条件为 C0>=5 CPM、相对每个其他类型平均表达超过两倍、至少 70% 可用患者的配对方向为正、至少两个 dataset 方向为正；按最弱其他类型 logFC 排序取最多 20 基因。COL*/FN1 及指定 Hallmark 炎症/缺氧/TNF genes 不进入 identity 集，另存 state/ECM 候选表。没有因临床结果调整基因。

最终 20 基因：CLDN11、CACNB4、PTPRD、BCHE、STK32B、ZNF385B、ASPA、NAV3、TNXB、FGF14、ANK2、C3、CCBE1、IL17D、AGTR1、KRT222、EPHA6、NTN1、GDF10、LMO3。这是训练数据导出的 C0-like identity 优先面板，不能称为经过实验确认的纯身份、LLPS 或保护性基因集；仍含部分基质/补体功能成分。

单样本评分采用固定 13659 基因背景中的平均秩，ties=average。20 基因 score = (mean(rank)−(20+1)/2)/(13659−20)−0.5。方法与单样本 rank scoring 思路一致，但未声称调用 singscore 包。参考：[singscore 官方说明](https://davislaboratory.github.io/singscore/articles/singscore.html)。芯片先按已保存的单一 symbol 映射，重复 probes 在提供的表达尺度上取中位数；不确定是否 log2 的 GSE48075 不作线性丰度拟合，排名不依赖单调 log 转换。四队列全部 20 个 signature genes 和背景基因完整可用。

留出 patient×type profiles：19 个 C0、78 个其他 stromal profiles，C0 对其他 stromal AUC=0.919，对全部其他类型 AUC=0.970；11 个 C0/C5 配对患者的 score 差中位数 +0.165。AUC 是描述性判别量，多个类型来自同患者，未作逐细胞独立推断。额外检查了 11184 个预留患者 stromal 单细胞，按患者/cluster 汇总中位数，结果支持 C0 较高但存在重叠和掉零。

20 名 mixture score 与真实 C0/(C0–C6) 的 rho=0.815，与总 stromal 占全参考类型比例的 rho=0.729；控制 total stroma 和 dataset 后 partial rank rho=0.701。后者仅作诊断，未把 total-stroma 当作研究对象的替代分母。该 score 同时受到总体基质含量影响，因此绝不能解释成 C0 within-stroma proportion 或细胞数量。

## 外部临床结果：仅 signature-based

表达样本/原发肿瘤沿用已有队列审查：TCGA 426/407、GSE31684 93/93、GSE48075 142/142、GSE13507 256/165。正常/复发标本不纳入临床；每队列原发样本标识唯一。GEO 未公开的跨标本患者关系不凭名称猜测。GSE13507 的异常正常来源 GSM340542 不在原发肿瘤集合。未取得可靠 purity，未填补虚假数值。

分期保留各自定义，Kruskal–Wallis 为主；有序趋势仅用于可排序级别，GEO Tis 不纳入趋势。TCGA 的 stage I 只有 2 人，图中保留并标明 n。GSE48075 排除 4 个混合/范围 T 编码以及 cstage 字段中的 1 个 pT；GSE13507 只能标为 reported TNM T component，不擅自称 pT。

| 队列 | 分期定义 / n | KW P（四队列 BH FDR） | 趋势 rho / P |
|---|---|---|---|
| TCGA | pathological AJCC / 405 | 5.22e−9（2.09e−8） | +0.290 / 2.91e−9 |
| GSE31684 | cystectomy pathological T / 93 | 0.0873（0.116） | +0.150 / 0.150 |
| GSE48075 | clinical T / 137 | 0.000431（0.000863） | +0.248 / 0.00369；n=135 |
| GSE13507 | reported TNM T / 165 | 0.122（0.122） | +0.100 / 0.201 |

OS 使用此前核验字段。TCGA 原时间为天，换算 months=days/30.4375；GEO 原时间为月。TCGA OS=1 死亡/0 删失；GSE31684 DOC/DOD=死亡、NED=删失；GSE48075 使用 os censor，不用 DSS；GSE13507 使用 overall survival，不用癌症特异死亡。固定每队列有效 OS 病例的中位数，Low<=median / High>median。KM 95% CI、删失和 risk table 均输出，所有风险人数经原随访时间逐点复核。

| 队列 | OS n / 事件 | Low / High | KM log-rank P | High/Low HR（95% CI） |
|---|---:|---:|---:|---|
| TCGA | 405 / 177 | 203 / 202 | 0.00724 | 1.50（1.11–2.03） |
| GSE31684 | 93 / 65 | 47 / 46 | 0.116 | 1.49（0.90–2.47） |
| GSE48075 | 73 / 45 | 37 / 36 | 0.0439 | 1.83（1.01–3.31） |
| GSE13507 | 165 / 69 | 83 / 82 | 0.634 | 1.12（0.70–1.80） |

四队列 log-rank BH FDR 为 0.0290、0.154、0.0879、0.634。GSE48075 只有所附 73 个侵袭性病例具有有效 OS，不能把这一生存结果外推至全部 142 个肿瘤。

| 队列 | 单变量每 SD HR（95% CI）；P | 调整后每 SD HR（95% CI）；P | 调整 n / 协变量 |
|---|---|---|---|
| TCGA | 1.21（1.06–1.38）；0.00483 | 1.11（0.97–1.28）；0.120 | 403；age+stage |
| GSE31684 | 1.18（0.93–1.50）；0.170 | 1.14（0.89–1.46）；0.296 | 93；age+T |
| GSE48075 | 0.98（0.76–1.27）；0.897 | 0.91（0.66–1.24）；0.548 | 71；age+T |
| GSE13507 | 0.93（0.73–1.19）；0.567 | 0.94（0.71–1.24）；0.664 | 165；age+T+grade |

每 SD 使用同一队列全部有效 OS 病例的 SD，complete-case 调整后不重新标准化。stage/T 在主调整模型中按有序数字线性进入，其 log-hazard 线性假设需要谨慎；grade 只在覆盖、事件及共线性条件允许的 GSE13507 纳入。变量按数据可用性和事前事件阈值选择，不按 P 值挑选。所有模型完整系数、样本、cutoff、PH 检验均保存。

**PH 限制不能忽略：** TCGA 多变量 score 的 Schoenfeld PH P=0.0293，因此 1.11 不是可无条件解释的恒定 HR。补充固定 score×log(1+t/12) 模型，12/36/60 月 HR 分别为 1.06（0.91–1.23）、1.33（1.08–1.63）、1.52（1.12–2.07）。这提示可能的时间依赖风险，但属于诊断后的模型敏感性，不是重复验证的恒定独立效应。GSE31684/GSE13507 的 stage PH P=0.0080/0.00237；以 stage 分层后 score HR 为 1.09（0.85–1.39）、1.00（0.73–1.35），仍不显著。KM 二分类的 PH 检验均未发现明显违反。森林图为描述性比较，未做跨平台、不同阶段构成及非 PH 模型的强行合并。

## 生物学结论

没有得到一致、稳健的“C0 有利预后”证据；TCGA 部分分析反而指向较高表达 signal 与不良结局相关，调整及 PH 检查后解释更复杂。GSE48075 的 median 分组名义显著与连续 Cox 阴性必须同时保留，不能仅挑 KM 得出结论。没有数据支持直接声称 C0 保护性、状态转化或因果机制。

单细胞中的 C0/(C0–C6) 相对下降，与晚期 bulk 中 C0-like expression score 增加并非同一个量；bulk score 还受总基质含量、表达状态和平台影响，不能据此否定或“验证”原比例结果。PI16/SFRP2-associated、具有炎症/缺氧响应的 fibroblast 状态仍是此前单细胞证据支持的工作解释，不是由本 20 基因 score 确立的新细胞类型。

后续优先研究 C0 上游炎症/缺氧等调控及其与 LLPS 的关系，结合独立单细胞/空间、已知比例混合及正交实验，区分状态与组成变化。当前临床结果不构成启动保护性或促血管因果机制论证的依据；未启动 NicheNet。

## 文件和复现

正式临床图在 04_stage/ 与 05_survival/；诊断图、模型比较和 signature 验证分别在 01_diagnostics/、02_model_comparison/、06_signature_fallback/。所有图独立 PDF+350 dpi PNG。正式临床图为四张 Fig2E_signature_*、四张 Fig2F_signature_* 和一张 Cox forest。全部绝对路径与 SHA256 见 summary/file_index.tsv；可点击图索引见 summary/figures.md。

从 `/home/jinshengxi/bladder_metabolism` 顺序运行 scripts/01_diagnostics.R、02_models.R、03_signature.R、04_singlecell_signature.R、05_external_scores.R、06_clinical.R、07_sensitivity_and_checks.R。common.R 读取既有 32 MB 稀疏聚合及已下载数据；除单细胞 signature 描述性检查读取原 fibro 对象外，没有重新运行大型预处理。seed=123。方法决策在临床计算前写入 02_model_comparison/method_decision.md。旧 Phase 1/1B/2 输出保持原样。
