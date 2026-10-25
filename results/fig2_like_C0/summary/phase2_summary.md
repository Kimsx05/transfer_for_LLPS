# C0 Phase 2 外部验证：停在 reference validation

完成日期：2026-10-25。**状态：STOP_ABUNDANCE_CLINICAL_ANALYSIS。** 已完成四队列数据获取/临床字段整理、14 类参考构建、独立患者验证及四张独立验证图。由于比例校准和跨 dataset 稳健性未通过事先记录的验证标准，遵照用户停止规则，未计算真实 bulk 的 C0 abundance，也未生成 abundance-based Fig2E/F、KM 或 Cox 结果。没有将表达 signature 分数冒充丰度，没有生成空图或更改切点来追求显著性。

## 1. 外部数据清单

| 队列 | 平台/表达 | 下载表达样本 | 原发肿瘤样本 | 可用分期字段 | 有效 OS / 死亡事件 |
|---|---|---:|---:|---|---:|
| TCGA-BLCA | HiSeq RNASeqV2；log2(RSEM normalized counts+1) | 426 | 407 | 405 个总体 pathological stage；另保留 pT/cT | 405 / 177 |
| GSE31684 | GPL570；GCRMA，矩阵范围 2.23–16.35 | 93 | 93 | 93 个 radical-cystectomy pathological T | 93 / 65 |
| GSE48075 | GPL6947；quantile-normalized，矩阵范围 5.05–15.35 | 142 | 142 | 138 个可解析 T；其中 137 cT、1 pT，另 4 个混合/范围编码不归类 | 73 / 45 |
| GSE13507 | GPL6102；quantile-normalized log2 intensity | 256 | 165 | 165 个 TNM 的 T 分量；原 metadata 未明确临床/病理依据 | 165 / 69 |

以上是数据可用性统计，不是 C0 临床分析的样本量。TCGA 407 个原发肿瘤对应 407 个独立 participant barcode；93 人 GSE31684 由 series overall_design 明确说明。GSE48075、GSE13507 的原发肿瘤各有 142、165 个唯一标本标签，整理表保留 patient_unit；现有矩阵不足以独立核实未公开的跨标本患者关联，未凭标签猜测额外配对。GSE13507 的 23 个复发标本来自 14 人，未混入 165 个原发肿瘤集，周围组织及正常样本也未进入临床集。

TCGA OS 使用 Xena curated survival 表（OS=1 死亡、0 删失），原始时间为天；months=days/30.4375。GEO 原始时间均为月：GSE31684 的 DOC/DOD 计为 OS 死亡，NED 计为删失，保留原始字段供人工核查；GSE48075 用 OS censor 的 uncensored/censored，不用 DSS censor；GSE13507 用 overall survival 的 death/survival，不混入癌症特异死亡。缺失/非正时间不进入有效 OS。

GSE48075 的公开 processing 未明确写出 log 变换，虽数值范围符合 log-scale，未据此擅自当成线性混合表达；将来若重新开放临床丰度分析仍需确认原始处理流程。GSE13507 当前公开矩阵的 GSM340542 注释出现 MPP/mouse-mm10 RNA-seq 与人膀胱芯片来源冲突，异常原样记录于 metadata_issues.tsv；它不在原发肿瘤集合，未以猜测替换其注释。grade、T、总体 stage 分列保留，不互相转换。当前整理未获得经核验的 tumor-purity 数值。

数据来源为 [NCBI GEO FTP](https://ftp.ncbi.nlm.nih.gov/geo/series/) 和 [UCSC Xena TCGA hub](https://tcga.xenahubs.net/)。每个文件的精确 URL、大小、SHA256 在 resources/download_manifest*.tsv；初始临床 gzip 链接返回 403 后成功获取同一官方未压缩表。原始表达 gzip 各保存一份，没有再复制大型 gene-level 矩阵。

## 2. 方法、训练隔离与参考规模

当前环境未配置 CIBERSORTx 用户/token，因此采用替代方法 [MuSiC](https://xuranw.github.io/MuSiC/reference/music.basic.html)。执行的是官方仓库原样保留的 music.iter/music.basic weighted NNLS 数值函数，依赖 nnls 1.6；不是 CIBERSORTx，也不是完整安装的 MuSiC 包。公开源代码、GPL 许可证、获取时 commit 信息和 sessionInfo 均保存。方法来源：[MuSiC 原始论文](https://www.nature.com/articles/s41467-018-08023-x)、[官方代码](https://github.com/xuranw/MuSiC)。

从冻结 full atlas 的 raw counts 聚合 C0–C6 和 Epithelial、Endothelial、Myeloid、T_NK、B、Plasma、Mast，共 14 类。沿用既有 QC 排除项及同患者额外标本排除决定，最终聚合 337136 个细胞、639 个 sample×type 库；没有更改任何标签。保存一份 32 MB 稀疏聚合和 5.1 MB 最终参考统计，不复制 3.5 GB 单细胞对象。

模型假设是线性 RNA 混合、可转移的细胞类型表达、用训练患者 mean UMI/cell 近似 cell-size scaling，并按跨患者表达方差加权。它不能自动解决实际 bulk 与单细胞平台的所有偏差，也不能把未包含的细胞状态识别成新类别。既有被排除 stromal/未定谱系 Cycling 不在参考内，真实组织中若含这些成分可能被误分配。

训练按 patient×type 汇总，以 >=20 个细胞、每类至少 3 个训练患者作为覆盖条件；不以单细胞作为独立重复。均值、方差及细胞尺寸由稀疏聚合计算，在确定性小型数据上与官方 music_basis 四类统计量逐一核对，最大差均为 0。没有用全体样本 DEG 选择基因。每 fold 的基因仅按训练平均表达至少 1 CPM 筛选。

五个 dataset 逐一留出，所有同患者标本同时隔离；这也避免同一 study 内未知患者关联引起跨训练/验证泄漏。验证 mixture 为真实留出标本所有已纳入类型的 raw-count 总和。每患者事先选总保留细胞最多的一个标本；主准确度要求至少 50 个 stromal 细胞。46 个标本、34 个患者获得预测；主评估 20 人。GSE222315 fold 因 C6 仅剩 1 位合格训练患者而无效，明确保留失败记录，未删除 C6 或悄悄合并。其余所有 230 个自然混合拟合均收敛。

所有估计系数来自同一 fit 并归一化；分母为同一尺度的 C0–C6 之和。这只保证代数定义一致，不保证估计准确。

## 3. 验证结果与停止原因

在检视预测前保存的 validation_plan.md 中规定：总体 Spearman>=0.70、MAE<=0.15、RMSE<=0.20、|bias|<=0.10；足够样本的各 dataset MAE<=0.20；C0/C5 pure-type 正确分配中位数>=0.70、相互误分配<=0.20，并需要完整 dataset 覆盖。这些是本任务的操作性标准，不是领域公认绝对阈值。

| 主评估（Full RNA） | 结果 |
|---|---:|
| 独立患者数 | 20 |
| Pearson / Spearman | 0.796 / 0.851 |
| MAE | 0.1735（17.4 百分点） |
| RMSE | 0.2515（25.1 百分点） |
| 平均偏差 | +0.0454 |
| CNP0000460，n=5 | MAE 0.3117；bias +0.2953 |
| InhouseData，n=10 | MAE 0.1316；bias −0.0147 |
| GSE135337，n=2 | MAE 0.2543；样本过少，不计算相关 |
| HRA003620，n=3 | MAE 0.0291；预测饱和，相关不可估计 |

总体排序相关较好，但不代表比例数值准确。CNP0000460 存在明显高估；小样本 HRA 的低误差也不能当作成功校准，部分预测为常数，相关为 NA。错误来自跨 study/平台差异和状态可识别性等因素，当前无法分解每个因素的贡献。

C0 pure-type 控制通常可被识别：按患者平均后 C0 正确分配中位数 91.3%；C5 正确分配 65.4%，约 9.4% 分给 C0、15.9% 分给 C4。这里为各分量分别取中位数，行和不必为 100%。因此并非“C0 与 C5 完全无法区分”，但 C5 的恢复和真实混合比例校准不足。原 gate 的逐 profile 汇总 C0 为 94.0%（19 profiles/18 patients），与患者汇总方向一致；C5 8 profiles/8 patients。详细 patient-summary 表和全部 pure profiles 均保留。

只保留实际 GPL570/GPL6947/GPL6102 覆盖基因时，MAE 为 0.1543/0.1610/0.1664，RMSE 为 0.2085/0.2278/0.2320；GPL6947 再加与真实比例独立的固定基因效率扰动（log2 SD 0.5），MAE=0.1748、RMSE=0.2613。这些是覆盖/扰动敏感性，不是实际配对芯片验证；不能据此宣称已经消除跨平台效应。

**结论：未通过 abundance 临床应用的验证门槛。** 停止决定同时基于较大绝对误差、dataset 偏倚和未覆盖 fold，不只因为某个结果略低于人为阈值。未尝试更换算法或阈值寻找“显著”临床关系。

## 4. 图与交付状态

新增四张单独 PDF + 350 dpi PNG，目录为 `/home/jinshengxi/bladder_metabolism/results/fig2_like_C0/phase2_external/reference_validation/`：

- FigS_C0_deconvolution_validation：留出患者真实与估计比例。
- FigS_C0_C5_deconvolution_confusion：C0/C5 类型识别。
- FigS_C0_deconvolution_platform_sensitivity：平台覆盖/效率敏感性。
- FigS_C0_reference_patient_coverage：14 类训练患者覆盖。

四队列的 Fig2E、Fig2F、二分类/连续/多变量 Cox 均**未执行**，原因是上游 reference validation 未通过；这不是“无显著分期或生存关联”的结果。当前没有证据评判四队列临床关联是否一致，也没有可报告的 HR、95%CI 或 log-rank P。stage/ 和 survival/ 内有停止原因说明，不包含空白占位图。可选 signature-score 探索未执行，以免绕开本轮停止节点并被误读成 abundance 验证。

全量绝对路径：`/home/jinshengxi/bladder_metabolism/results/fig2_like_C0/summary/phase1B_phase2_file_index.tsv`。正式图浏览索引：`/home/jinshengxi/bladder_metabolism/results/fig2_like_C0/summary/phase1B_phase2_figures.md`。

## 5. 生物学判断和下一步

结合 Phase 1B，C0 是 PI16/SFRP2-associated fibroblast 状态，在所捕获 stromal 中相对减少，保留炎症/TNF/缺氧反应；不能由此推导保护作用、绝对细胞数减少或 tumor-specific LLPS 增强。EMT/TGF-beta 证据具有参照和样本集合依赖性，不能作为已确立机制。当前没有可靠的外部临床 abundance 关联结论。

较值得优先研究的是 **C0 的上游状态调控**：验证炎症/缺氧相关刺激是否改变 C0 状态，并用正交方法区分 LLPS 生物学变化与测序复杂度影响。现阶段 C0 对其他细胞的下游促血管效应尚缺特异性证据，不应机械沿用原论文机制。没有启动 NicheNet。

外部丰度验证下一步需要补足跨患者/跨 dataset 的 C6 及其他低覆盖类型，并用独立实测 bulk 或已知比例混合物做平台校准；若未来启用配置完整的 CIBERSORTx，需要同样的无泄漏验证，不能默认其输出准确。可以另行预先定义表达 signature 的探索性临床分析，但必须称为 expression-based score，不能解释成 C0 细胞丰度。
