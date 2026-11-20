from pathlib import Path
import pandas as pd,numpy as np,json
O=Path('results/C5_Endothelial_GSEA_nominalNicheNet_v1')
def r(p):return pd.read_csv(O/p,sep='\t')
def table(d):
 def fmt(v):
  if pd.isna(v):return 'NA'
  if isinstance(v,(float,np.floating)):return f'{v:.4g}'
  return str(v).replace('|','\\|').replace('\n',' ')
 return '| '+' | '.join(map(str,d.columns))+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in d.itertuples(index=False,name=None))
g=r('01_Endothelial_Hallmark_GSEA/all_models.tsv');wide=r('01_Endothelial_Hallmark_GSEA/model_comparison.tsv');a=r('03_nominal_NicheNet/all_activities.tsv');b=a[a.program.eq('Set_B')].sort_values('rank');rc=r('04_rank_sensitivity/rank_comparisons.tsv');lo=r('04_rank_sensitivity/LOPO_status.tsv');lr=r('04_rank_sensitivity/LOPO_rank_summary.tsv');ls=r('04_rank_sensitivity/LOPO_ligand_stability.tsv');cal=r('04_rank_sensitivity/random_program_calibration.tsv');di=r('02_nominal_response_genesets/set_diagnostics.tsv');integ=r('06_CellChat_spatial_integration/SetB_complete_evidence.tsv');cnt=r('02_nominal_response_genesets/set_counts.tsv');ident=r('02_nominal_response_genesets/identity_co_detection_patient.tsv');inf=r('04_rank_sensitivity/target_influence_summary.tsv');tchanges=r('04_rank_sensitivity/target_influence_rank_changes.tsv');gc=r('04_rank_sensitivity/SetB_gene_influence_summary.tsv');sp=r('06_CellChat_spatial_integration/new_axes_expression_audit.tsv');spat=r('06_CellChat_spatial_integration/new_axes_known_patient_description.tsv');emt=r('05_EMT_program_NicheNet/EMT_coverage.tsv');comp=r('05_EMT_program_NicheNet/program_rank_wide.tsv')
focus=['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_ANGIOGENESIS','HALLMARK_TGF_BETA_SIGNALING','HALLMARK_GLYCOLYSIS','HALLMARK_HYPOXIA','HALLMARK_TNFA_SIGNALING_VIA_NFKB']
f=g[g.pathway.isin(focus)].copy();f['NES (BH-FDR)']=f.apply(lambda x:f'{x.NES:+.3f} ({x.padj:.3g})',axis=1);f=f.pivot(index='pathway',columns='model',values='NES (BH-FDR)').reset_index()
fr=a[a.test_ligand.isin(['ADM','MDK','SEMA3C','WNT2','FN1','COL1A1','COL1A2','TGFB1','VEGFA'])&a.program.isin(['Set_A','Set_B','Set_C','EMT'])].pivot(index='test_ligand',columns='program',values='rank').reset_index()
ranktop=b.head(10).merge(ls,on='test_ligand').merge(cal[['test_ligand','empirical_exceedance','empirical_BH_across246']],on='test_ligand')
cc=integ[(integ['rank']<=10)&integ.feature.notna()].sort_values(['rank','patient_condition_median'],ascending=[True,False]).drop_duplicates('test_ligand')
spdesc=sp.groupby('ligand').agg(n_sections=('section_id','size'),median_sender_detection=('ligand_detection_sender','median'),median_joint_receptor_detection=('receptor_joint_detection_receiver','median'),median_edge_detection=('edge_joint_detection','median')).reset_index()
cross=wide[(wide.padj_primary<.1)&wide.both_major_same_as_primary];core=g[(g.model=='primary')&g.pathway.isin(cross.pathway)]
cpct=lambda c:100*ident[c].sum()/ident.n_endothelial.sum()
# Explicit top target and influential deletion evidence; no synthetic confidence intervals.
imp=tchanges[tchanges.rank_original<=10].copy();imp['rank_change']=imp.rank_drop-imp.rank_original;imp=imp.reindex(imp.rank_change.abs().sort_values(ascending=False).index).head(12)
report=f'''# C5–Endothelial GSEA 与 nominal-DEG NicheNet v1

完成日期：2026-11-20。全部分析使用现有冻结数据；随机种子 123，BLAS/OpenMP/R 单线程、CPU affinity 限于 4 核。没有修改历史结果、重新注释或重新运行 CellChat/cell2location。

## 核心结论

**没有获得主模型 FDR 支持的 Endothelial EMT、Angiogenesis 或 TGF-β 响应。** 主模型 EMT NES=1.183，BH-FDR=0.320；采样规模调整后 EMT 得到支持，但不能用敏感性模型替换主模型。标准 EMT-NicheNet 为 **REFERENCE_ONLY**，不等于观察到 EndMT。

名义 Set A/B/C 经双队列正向及固定背景筛选，分别为 **28/123/302 genes**。Set B Top10 为：{', '.join(b.head(10).test_ligand)}。这些不是 FDR 显著 DEG，也不是已验证的 C5-induced targets。完整 246 配体均输出，不只保留与 CellChat 的交集。

最值得保留的是“真实内皮转录变化下的候选排序及其不稳定性”，而非确定的 C5→内皮机制。TGFB1 在多个阈值/标准程序中排名较高，但其来源并非 C5 特异，既有空间支持仍受数据集及组成限制。

## 1. 输入、排名、资源和可重复性

- 21 名符合原标准的独立肿瘤患者：InhouseData 9、GSE222315 8、CNP0000460 3、GSE135337 1；Fibro C0–C6≥50、Endothelial≥20。同患者多个肿瘤 specimen 已聚合。
- 主模型 `~ dataset + standardized C5 fraction`，edgeR TMM → limma-voom → robust eBayes。主要 DEG 复用，不从大型单细胞对象重建。LOPO 用原小型 pseudobulk 重拟合主模型与两个主要队列模型。
- 重拟合原模型与原 t 值逐基因一致，见 `00_input_audit/model_reproduction.tsv`。
- GSEA 使用完整有限 moderated t 排名，不按 DEG 显著性截断。各模型检测背景不同，未补入缺失基因。混合 ENSG 标识原样保留在完整排名中，只对符号匹配的 Hallmark 基因计命中；未猜测别名或将其移至末尾。重复符号为 0、非有限 t 为 0；详见审查表。
- MSigDB **2026.1.Hs** Hallmark 50，与 C5 既有 GSEA 同版；fgseaMultilevel，minSize=10、maxSize=500、eps=0、nPermSimple=10000，SerialParam(RNGseed=123)。每模型全部 50 通路 BH 校正。少量 t ties 按 gene symbol 确定输入顺序，日志记录原 fgsea 警告；不存在结果导向重排。
- C5 旧 GSEA 使用 avg_log2FC、不同 detection 过滤，内皮使用患者级 moderated t，故不能将两者 NES 数值当作直接可比效应量。
- 配体/受体表达与背景沿用上一轮的全部合格肿瘤表达审查，和本轮 21 位 pseudobulk 患者的统计样本不是完全相同的集合；因此表达支持人数不能冒充本轮响应的患者分母。随机匹配也沿用原背景 pooled detection，而不是重新选择背景。
- NicheNet 固定 human v2、nichenetr 2.2.1.1、7,070 背景与 246 broad ligands。官方 `predict_ligand_activities()`；以 corrected AUPR、Pearson、gene symbol 排序。资源 SHA256 已复核，未下载或重新筛选 sender/background。
- Set A/B/C：主模型 logFC>0，P<0.01/0.05/0.10，InhouseData 与 GSE222315 logFC 均>0，且属于固定背景及 LT target universe。未在某队列通过表达过滤的基因标为未检验，不能算阴性或强行满足方向条件。
- 每套集合的精确列表见 `02_nominal_response_genesets/frozen_genesets.tsv`；原始阈值、排除清单和流程见 `00_input_audit/prespecified_parameters.json`。

## 2. Endothelial Hallmark GSEA

表中为 NES (BH-FDR)，所有 250 项检验见 `01_Endothelial_Hallmark_GSEA/all_models.tsv`。

{table(f)}

主模型 FDR<0.10 有 {len(g[(g.model=='primary')&(g.padj<.1)])} 个通路，其中 {len(cross)} 个与两个主要队列同向。以下是“主模型 FDR 支持 + 队列方向一致”，不是所有队列各自显著，也不是独立验证：

{table(core[['pathway','NES','padj']])}

E2F targets 主模型负向，而 GSE222315 正向，明确属于来源不一致。EMT 在 InhouseData BH=0.0417，但 GSE222315 BH=0.684；采样调整模型 BH=0.0341。Angiogenesis 在 GSE222315 BH=0.0941、采样调整 BH=0.0963，主模型不支持；TGF-β 仅 GSE222315 BH=0.0301。不能选择这些更有利的模型替代预定主结果。

C5 自身 EMT 在旧 5%/10%/25% 检测阈值下均为正且 FDR<0.05，是用户重点程序中最稳定的上调特征。Glycolysis 在 25% 阈值 FDR≈0.0545，稳定性较弱。C5 的 EMT leading-edge 未用于任何内皮 geneset_oi。

## 3. 名义基因集与来源稳定性

{table(cnt)}

{table(di)}

三套固定间质清单均只命中 **FAP**；保留原集，同时另存 noECM 版本 (27/122/301 genes)。这只是预定义标志排除敏感性，不能说明没有广泛 ECM 信号或污染。双队列方向条件利用了参与合并模型的患者，不构成独立复制。

21 次 LOPO 均可拟合，名义 Set B 数量 {lo.n_genes.min()}–{lo.n_genes.max()}，原集 Jaccard 中位数 {lo.jaccard.median():.3f} (范围 {lo.jaccard.min():.3f}–{lo.jaccard.max():.3f})。所有 21 次均达到预定至少 15 基因门槛；没有补齐基因。原 Set B 基因在留出中的正向比例、入选比例、最大 P 值和效应范围已逐基因保存；{int((gc.selected_fraction==1).sum())}/{len(gc)} 个在全部留出中仍入选。LOPO 中队列内方向也重新估计，未使用被留出患者。POSTN 最差留出为 InhouseData::T4295396，名次由 1 降至 114；这是单患者影响的具体例子，不能只报告中位名次仍为 1。

## 4. Set B NicheNet 与阈值/患者敏感性

{table(ranktop[['test_ligand','rank','aupr_corrected','auroc','pearson','CellChat_supported','LOPO_rank_min','LOPO_rank_max','LOPO_top10_fraction','empirical_exceedance','empirical_BH_across246']])}

ANGPTL4 在 Set B 的官方 top250-target 交集为空，未通过改阈值补造链接；activity 可计算与可展示 top250 link 是不同结果。

Top 配体 corrected AUPR 较小，且部分 AUROC 接近 0.5、Pearson 接近 0 或负值；不能只凭“第 1 名”写成配体实际活性增强。gene set 大小不同，跨 A/B/C/EMT 比较仅讨论排序，不将 AUPR 差值当作方案优劣。

{table(rc[rc.program1.isin(['Set_A','Set_B','Set_C']) & rc.program2.isin(['Set_A','Set_B','Set_C','Set_B_noECM','EMT'])])}

患者留出排名相关中位数 {lr.rank_spearman.median():.3f} (范围 {lr.rank_spearman.min():.3f}–{lr.rank_spearman.max():.3f})，Top10 与完整 Set B 交集为 {lr.top10_overlap.min()}–{lr.top10_overlap.max()}。这是影响诊断，不是独立测试性能或患者级置信区间。逐配体中位/最小/最大名次见 `LOPO_ligand_stability.tsv`。

重点配体真实名次：

{table(fr)}

尤其 ADM 从 A/B 的低排名升至 Set C 第 8，不能据此挑选宽松阈值作为唯一主分析。TGFB1 的 A/B/C 名次为 9/5/4；WNT2 为 66/12/18。MDK、SEMA3C、FN1/COLLAGEN 未一致进入名义 Top10。

50 次随机对照严格匹配 Set B 数量和内皮 pooled detection 十分位分布，每次 123 genes，无放回抽样，允许与原 Set B 重叠。`empirical_exceedance=(1+null>=observed)/(51)`；另给全部 246 配体 BH，最低经验 BH 约 0.506。分辨率最低 1/51，这些是网络/程序匹配经验校准，不是患者级显著性，更没有模拟完整名义基因筛选过程。{int((cal.empirical_BH_across246<.05).sum())}/246 配体通过此经验 BH<0.05。

预定义 FAP 删除后前三名不变。另按冻结 Set B Top10 各配体的最大先验权重 target，逐一删除唯一 target 重新评分；不是选择使结果更好的删除方案：

{table(inf)}

影响最大的名次变化 (完整表保留全部配体)：

{table(imp[['test_ligand','removed_gene','rank_original','rank_drop','rank_change']])}

这说明排名是否依赖个别 target，而不能把 LT 权重份额等同评分的精确贡献或真实靶基因激活。

## 5. 标准 EMT 程序及污染/EndMT 风险

EMT 原始 {int(emt.original_genes.iloc[0])} genes，固定背景匹配 {int(emt.matched_genes.iloc[0])}；状态 **{emt.status.iloc[0]}**。官方 NicheNet 可计算，但没有通过患者响应门槛。Top10：{', '.join(a[a.program.eq('EMT')].sort_values('rank').head(10).test_ligand)}。

EMT 与 nominal Set B 的 Top10 交集为 {int(((comp.Set_B<=10)&(comp.EMT<=10)).sum())} 个；共同候选是 {', '.join(comp.loc[(comp.Set_B<=10)&(comp.EMT<=10),'test_ligand'])}。既有 Angiogenesis/TGF-β 使用 2025.1，本轮 EMT 用 2026.1，比较图显式标注版本差异，仅比较优先级，不当作同版数值验证。

冻结内皮标注中，合格患者共有 {int(ident.n_endothelial.sum())} 个内皮细胞；至少检出 1/2 个 PECAM1/VWF/CDH5/KDR 的比例分别为 {cpct('n_endo_marker_ge1'):.1f}%/{cpct('n_endo_marker_ge2'):.1f}%。至少检出 2 个 FAP/COL1A1/COL1A2/POSTN/DCN/LUM 的比例为 {cpct('n_stromal_marker_ge2'):.1f}%，同时≥2 内皮与≥2间质标志为 {cpct('n_both_ge2'):.1f}%；≥2间质但无上述内皮标志为 {cpct('n_stromal_ge2_endo0'):.2f}%。这些细胞层比例仅为身份诊断，不作独立患者推断。共检出使用冻结 LogNormalize10000 稀疏矩阵中的非零值，与原始计数非零对应，不是新增打分。

原冻结流程已经进行 QC/doublet 过滤；本次未新聚类、未改变标签。共同检出可来自状态变化、残余 doublet 或环境 RNA，无法据此确认或排除污染；检测不足也会漏掉真正内皮标志。患者级 CPM、检测率与 specimen/患者级共同检出表均保存。适当名称是 Endothelial EMT-associated transcriptional signature，**不能宣称 EndMT 或 C5 诱导 EndMT**。

## 6. CellChat、C5 特异性与已有空间证据

Set B Top10 中有 t20 阳性 CellChat 的配体及每配体强度最高的既有轴如下；完整多轴、t10/t20/t30 和 eligible 分母均保留：

{table(cc[['test_ligand','rank','receptor','positive_patients','eligible_patients','patient_condition_median','C5_pooled_pct','Fibro_C0_median_difference','Fibro_C4_median_difference','required_subunits','Endo_min_subunit_pct']])}

POSTN 相对 C0/C4 的差值中位数约 +0.960/+0.762，患者正向比例 12/13、13/14；INHBA 约 +0.561/+0.526，正向比例 13/13、11/14。这两者相对富集较清楚，但分别存在患者/target 依赖及 Five-gene 共享基因循环性。HMGB1 也有相对表达富集，但没有本次 t20 C5→Endothelial CellChat 支持；ANGPTL4、TGFB1 则不能称为 C5 富集配体。

表达差异是共同可评估患者中 mean-log-expression 的差值中位数，不是患者级 DE FDR。只有相对 C0/C4 更高及患者方向支持的配体可称相对富集，不能推广为 C5-specific。TGFB1 并非 C5 特异；NicheNet 的 ligand-target 先验不验证某一 CellChat receptor complex。

原空间分析复用：FN1–ITGA5/ITGB1 为共享 ECM 对照，调整 Ring1 rho≈0.074 (29 位已确认患者)，TGFB1–ACVR1/TGFBR1/TGFBR2 rho≈0.046 (26 位患者)，但 TGFB1 leave-GSE319536-out CI 跨 0。原匹配平衡子集仅 1 位已确认患者，504 项平移检验无 BH<0.05；因此不能用原正相关宣称 C5 特异通讯。其他原轴结果不因本轮排序而重算或隐藏。

按冻结规则新增且仅新增 POSTN–ITGAV/ITGB5、INHBA–ACVR1B/ACVR2A、ANGPTL4–CDH5 的轻量表达审查 (63 切片)。以下为逐切片中位数，无新的空间显著性推断：

{table(spdesc)}

发送区仍为原 Five-high (>切片 90% 分位、排除零分) 且 Fibro-rich (≥75% 分位)，接收区为 Endothelial-rich (≥75% 分位)，仅实际六邻接边。受体复杂体以全部必要亚基共同检出审查，不用单个亚基代替；表达用 raw counts→log1p(CP10k)。有边的 joint detection 只是表达相容性，没有匹配或空间零模型，**不构成新轴的邻接富集证据**。INHBA 本身是 Five-gene 成员，其发送区高表达部分是定义造成的；不能视为独立空间验证。INHBA 的必要亚基单细胞最小检测比例低于 10%，完整复合体支持需谨慎。

本轮不重新运行 63 切片的完整空间统计。已确认患者按照原身份审查的两个合法状态聚合 (30 切片/29 患者)，其余切片仅描述；未将带有 UNVERIFIED 的标签误作 verified。新增轻量审查的患者描述有 {spat.patient_uid.nunique()} 位，不能与单细胞患者擅自配对。

## 7. 小文章中的证据分层与十二问回答

1. **C5 自身最稳健的重点上调程序是 EMT**，旧多检测阈值均支持；内皮主模型没有复制该显著性。
2. 内皮 EMT/Angiogenesis/TGF-β 主模型均未获 BH<0.10 支持。队列/采样调整模型出现支持，应完整呈现而不替换主分析。
3. 两主要队列对上述三程序均正向，但显著性不一致；Hypoxia、TNFA 等方向不一致。全 50 通路比较已保存。
4. A/B/C 为 28/123/302，方向筛选前后的全部数量见上表；都是名义候选。
5. Set B Top10 见第 4 节；不强制任何旧候选进入。
6. 阈值与 LOPO 稳定性为连续诊断，不能概括为全部稳定；具体相关、交集及逐配体范围已给出。
7. 重点候选完整名次见上表。ADM 的宽松阈值提升及 MDK/SEMA3C 的低/变化排名必须保留。
8. Top10 中 CellChat 支持见第 6 节；预测强度与 positive/eligible 分母不能混为普遍发生率。
9. C5 相对 C0/C1/C4/其余 Fibro 的表达差值和方向比例均在整合表，shared CAF ligand 不能称 C5-specific。
10. EMT/Set B 只有部分共同候选，主要体现 TGF 家族先验一致性；EMT 是 reference-only。
11. 间质共同检出和模型敏感性构成 EndMT 解释风险；没有充分证据定义 EndMT，也不能仅凭共同检出判定污染。
12. 正文可以展示透明的患者级 GSEA、固定 Set B 排名及关键敏感性，结论限于探索性候选；完整 A/C、LOPO、随机对照、参考 EMT 和新增空间表达审查应作为补充，避免用这些弱证据形成确定机制主线。

### FDR 支持且跨队列方向一致的观察

以第 2 节列出的主模型通路为依据，尤其 OXPHOS/MYC/DNA repair 负向及 apical surface/mitotic spindle 正向。它们是观察性患者关联，未证明由 C5 因果驱动。不同通路的表达相关结构仍会影响 GSEA，通路 FDR 不等于 DEG FDR。

### 名义 DEG 支持的探索性结果

固定 Set B 的候选排序、A/C 敏感性、21 次 LOPO、target 删除及随机程序校准；无 FDR 显著主模型 DEG。TGFB1/POSTN/INHBA/ANGPTL4 等仅按各自完整稳定性和来源支持描述，不统一宣称阳性。

### 标准功能程序 NicheNet 先验预测

完整 Hallmark EMT 及旧 Angiogenesis/TGF-β 的网络排序。EMT 不通过主门槛，参考程序不能替代观察到的内皮响应。

### CellChat 与空间表达补充证据

保留 t20 主分析、t10/t30、受体亚基、患者分母、共享 CAF 表达；旧空间混杂限制不变。新三轴只证明技术表达可评估程度，不证明空间富集、真实通讯或细胞内靶反应。

### 不支持或无法确认

不支持主模型显著 C5-associated EndMT/血管生成/TGF-β 响应，不支持真实 C5 空间身份、C5-specific ligand 或配体因果效应。正式 C5-specific 空间定位 NO-GO 不变。最有价值后续是独立患者、内皮纯化/多重原位及针对候选的扰动验证，而不是继续放宽 P 阈值。

## 文件与图说明

- `01_Endothelial_Hallmark_GSEA/`：5 模型、250 通路、完整排名、leading edges、方向比较和客观选曲线。
- `02_nominal_response_genesets/`：完整基因审查、固定列表、去间质版本、身份/污染诊断。
- `03_nominal_NicheNet/`：所有 activity、全部 Top10 links、checkpoint；`04_rank_sensitivity/`：21 LOPO、50 random、target influence。
- `05_EMT_program_NicheNet/`：REFERENCE_ONLY、原始/匹配数量、标准程序排名对比。
- `06_CellChat_spatial_integration/`：完整证据、多阈值 CellChat、冻结三轴与表达审查。
- Figures A: all Hallmark/focus，星号 * BH<0.10、** BH<0.05；B: 主模型 FDR<0.10 中最大的两个正 NES 与两个负 NES；C: 名义筛选/队列效应；D: Set B、阈值、LOPO、随机；E: Top10 prior links 的最高 max-weight 35 targets (完整 top10 links 表不作展示截断；链接定义仍为每配体官方全 target top250 与 geneset 交集，空白不等于零调控潜力)；F: 标准程序排名 (版本明确)。
- 最终自动审查见 `00_input_audit/final_validation.tsv`，核对原输入哈希、250项GSEA的BH、固定集合/配体数量、21次LOPO、50次检测率匹配随机集合、三轴63切片、患者身份、PDF矢量和PNG分辨率。
- 所有图带原始绘图数据；PDF 为矢量、PNG 300 dpi、Arial 12pt、白底；有图例时导出主体与独立图例 PDF。LOPO range 不是 CI，随机 95% 区间不是患者 CI。
'''
(O/'COMPREHENSIVE_REPORT.md').write_text(report)
(O/'README.md').write_text('''# C5–Endothelial GSEA / nominal NicheNet v1

详见 [COMPREHENSIVE_REPORT.md](COMPREHENSIVE_REPORT.md)。固定 21 肿瘤患者、7,070 背景、246 配体；主要名义 Set B，A/C 敏感性。主模型无 FDR 显著 DEG，EMT 为 REFERENCE_ONLY。

从项目根目录执行 `bash results/C5_Endothelial_GSEA_nominalNicheNet_v1/scripts/run_all.sh`。CPU affinity 最多4核，各计算进程单线程，seed=123。历史输入只读，按阶段/程序 checkpoint 续跑；不复制大矩阵。复现依赖本机既有 R 库及 Python 依赖路径，见 run_R.sh/run_python.sh 和 sessionInfo。输入哈希见 00_input_audit。

NicheNet 活性分数是先验网络的程序匹配，名义候选不称 significant DEG；队列方向与LOPO不构成外部验证。空间新增仅3轴表达审查，没有新的通讯或空间显著性推断。
''')
print('REPORT_COMPLETE')
