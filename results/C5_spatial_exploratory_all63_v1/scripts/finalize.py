import pyarrow as pa
pa.set_cpu_count(1);pa.set_io_thread_count(1)
from pathlib import Path
import json,hashlib,html,h5py
import numpy as np,pandas as pd
O=Path('results/C5_spatial_exploratory_all63_v1')
def read(p):return pd.read_csv(O/p,sep='\t')
def fmt(x):return f'{x:.3f}' if pd.notna(x) else 'NA'
def md(df):
 d=df.copy()
 for c in d.select_dtypes('number'):d[c]=d[c].map(lambda x:str(int(x)) if pd.notna(x) else 'NA') if d[c].dtype.kind in 'iu' else d[c].map(fmt)
 return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in d.fillna('NA').itertuples(index=False,name=None))
s=read('tables/section_summary.tsv');c=read('04_fibroblast_and_depth_controls/correlations.tsv');m=read('03_spatial_clustering/moran.tsv');b=read('04_fibroblast_and_depth_controls/match_balance.tsv');g=read('00_inputs_and_QC/gene_coverage_detection.tsv');p=read('07_dataset_heterogeneity/patient_effect_intervals.tsv');co=read('03_spatial_clustering/components.tsv');gene=read('04_fibroblast_and_depth_controls/A26_gene_rank_contributions.tsv');r=read('06_spatial_relationships/section_effects.tsv');sig=json.loads((O/'00_inputs_and_QC/frozen_signatures.json').read_text())
original=json.loads(Path('results/C5_spatial_identifiability_v2_relaxed_marker/01_signature_candidates/fixed_signatures.json').read_text());assert sig['A26']==original['A26'] and len(sig['A26'])==26
qc=[];allids=[];coverage_rows=[];mapping=read('00_inputs_and_QC/section_mapping.tsv').set_index('section_id')
for row in s.itertuples():
 f=O/'01_spatial_scores'/f'{row.section_id}.parquet';d=pd.read_parquet(f);allids.extend(d.spot_id);assert len(d)==row.n_spots
 for name,genes in sig.items():coverage_rows.append(dict(section_id=row.section_id,dataset=row.dataset,signature=name,requested_genes=len(genes),covered_genes=int(d[name+'_covered_genes'].iloc[0]),median_detected_genes=float(d[name+'_detected_genes'].median()),mean_detected_genes=float(d[name+'_detected_genes'].mean()),any_detected_fraction=float((d[name+'_detected_genes']>0).mean())))
 assert np.isfinite(d[[k+'_'+m for k in sig for m in ['UCell','AUCell']]].to_numpy()).all()
 for suffix in ['_overview.pdf','_overview.png','_overview_panel_only.pdf','_overview_legend_only.pdf','_HE_provisional.pdf','_HE_provisional.png']:assert (O/'02_spatial_maps'/f'{row.section_id}{suffix}').exists()
 abcols=[c for c in d if c.startswith(('major_mean__','major_q05__','sub_mean__','sub_q05__'))];assert len(abcols)==78 and np.isfinite(d[abcols].to_numpy()).all() and (d[abcols].to_numpy()>=0).all()
 with h5py.File(Path('/home/jinshengxi/LLPS/spatial_data/LLPS_spatial_transfer_manual')/mapping.loc[row.section_id,'expression_relative_path'],'r') as f:n_universe=int(f['X'].attrs['shape'][1])
 qc.append(dict(section_id=row.section_id,n_spots=len(d),A26_zero_score_fraction=float((d.A26_UCell==0).mean()),scoring_gene_universe=n_universe,UCell_maxRank=1500,AUCell_aucMaxRank=int(np.ceil(.05*n_universe)),score_status='PASS_DESCRIPTIVE_OUTPUT',map_status='PASS_OUTPUT_EXISTS',localization_status='NO_GO_UNCHANGED'))
assert len(s)==63 and len(allids)==265002 and len(set(allids))==265002
for row in read('00_inputs_and_QC/history_checksums.tsv').itertuples():assert hashlib.sha256(Path(row.path).read_bytes()).hexdigest()==row.sha256
coverage_summary=pd.DataFrame(coverage_rows);coverage_summary.to_csv(O/'00_inputs_and_QC/signature_section_detection_summary.tsv',sep='\t',index=False)
pd.DataFrame(qc).to_csv(O/'00_inputs_and_QC/final_QC.tsv',sep='\t',index=False)
pd.DataFrame(qc)[['section_id','A26_zero_score_fraction']].to_csv(O/'00_inputs_and_QC/A26_zero_score_fractions.tsv',sep='\t',index=False)
A=c[c.score=='A26_UCell'];corr=A.groupby('variable').spearman.agg(['median','min','max']).reset_index();corr.to_csv(O/'tables/A26_correlation_summary.tsv',sep='\t',index=False)
clusters=s[s.exploratory_clustering_flag];non=s[~s.exploratory_clustering_flag];res=m[m.score=='A26_adjusted_rank_residual'];known=s.dropna(subset=['patient_uid']);main=p[(p.dataset=='ALL_KNOWN_PATIENTS')];ep=r[(r.analysis=='epithelial_distance')&(r.variable.astype(str)=='0.3')]
# Per-dataset primary diagnostics, including controls, so pooled summaries never hide source heterogeneity.
dataset=s.groupby('dataset').agg(n_sections=('section_id','size'),median_moran=('moran_I','median'),A26_high_in_fib_high=('A26_high_in_fib_top_quartile','median'),A26_high_in_fib_low=('A26_high_in_fib_bottom_quartile','median'),C2C4_overlap=('C2C4_high_overlap','median')).reset_index()
for cov,label in [('major_mean__Fibroblast','rho_Fib'),('nCount','rho_UMI'),('nFeature','rho_features'),('core15_UCell','rho_core'),('five_UCell','rho_five'),('C2C4_program_UCell','rho_C2C4')]:dataset=dataset.merge(A[A.variable==cov].groupby('dataset').spearman.median().rename(label),on='dataset',how='left')
dataset.to_csv(O/'07_dataset_heterogeneity/dataset_diagnostics.tsv',sep='\t',index=False)
shares=gene.groupby('gene').high_rank_weight_share.median().sort_values(ascending=False);dominance=gene.groupby('section_id').high_rank_weight_share.apply(lambda x:x.nlargest(3).sum());missing=g[~g.present]
ringdata=read('06_spatial_relationships/rings.tsv');ring_summary=ringdata[ringdata.analysis=='matched_seed_neighborhood'].groupby(['variable','ring']).agg(n_sections=('effect','count'),median_effect=('effect','median'),n_positive=('effect',lambda x:(x>0).sum())).reset_index();ring_summary.to_csv(O/'tables/neighborhood_ring_summary.tsv',sep='\t',index=False)
sens=read('07_dataset_heterogeneity/epithelial_matching_sensitivity_summary.tsv');partial=read('06_spatial_relationships/epithelial_distance_partial.tsv');weak=s[s.moran_I<.3];fibneg=A[(A.variable=='major_mean__Fibroblast')&(A.spearman<0)]
report=f'''# C5-associated 全量空间探索报告

已完成 **63/63 张切片、265,002 spots**。全量结果用于描述性探索；既有正式 C5-specific 定位 **NO-GO 保持不变**。本轮没有重训 cell2location、没有改变 A26、没有用真实空间图反向选基因。

## 交付与数据身份

- [63张缩略图总览](figures/all63_percentile_atlas.png)、[可浏览图册](GALLERY.html)、[逐切片汇总](tables/section_summary.tsv)。每张切片10个连续值panel及辅助HE视图，PDF矢量点图、PNG 300 dpi。HE本身是原始栅格图像。
- 所有spot级原始评分、百分位、检测基因数、nCount/nFeature、30亚型和9大类mean/q05均保存在`01_spatial_scores/*.parquet`，保留原spot_id及barcode，未重复保存表达矩阵。
- 原7张pilot的五基因/三基因两种评分与历史结果逐spot完全一致，最大绝对误差0。技术通过后扩展全部63张，没有模型AUROC门槛。
- 已知患者身份切片{len(known)}张、{known.patient_uid.nunique()}位患者，来自{known.dataset.nunique()}个dataset；其余切片不虚构患者身份。可靠Normal 3、Tumor 17、UNRESOLVED 43。CanDis/Inhouse/NPJ来源标签沿用UNVERIFIED限制。
- 所有63张HE只有来源/尺寸/边界几何检查，未独立组织学核验。主图不叠加HE；辅助图明确标为provisional，并使用原始上方原点像素方向，不旋转、不镜像、不转换微米。

## 主要观察

A26呈现可复现的空间程序结构，并常见于Fibroblast较高、上皮较低的区域；但微环境关联在进一步控制上皮组成后明显减弱。当前最有支持的是**C5-associated stromal/ECM程序的描述性空间结构**，而不是C5特异细胞定位或促血管/免疫机制。不能把高Moran I或漂亮图像当成新的辨识性验证。

## 九个问题

1. **A26呈现什么分布？** 全部切片连续分数均已输出。切片A26中位数范围{fmt(s.A26_median.min())}–{fmt(s.A26_median.max())}，这是未校准的描述值，不构成跨dataset绝对丰度差异。原始score统一[0,1]显示；百分位panel与总览只比较切片内部位置；并列零分也会得到中间排名，不能将有颜色理解为已检出强表达。CanDis_2026_Sample10约59.5%的spot A26 UCell为0，完整零分比例表已保存。MMP11、丰度和UMI采用全63共享99%上限，剪裁数逐图保存。
2. **哪些有聚集？** 按预先记录的操作性定义，{len(clusters)}/63张呈正空间自相关且BH<0.05。A26 Moran I中位数{fmt(s.moran_I.median())}，范围{fmt(s.moran_I.min())}–{fmt(s.moran_I.max())}。未满足定义的切片：{', '.join(non.section_id) or '无'}。其中较弱的{len(weak)}张(描述性I<0.3，非预设生物学门槛)为：{', '.join(weak.section_id)}。统计可检测不等于视觉上形成强聚集。199次置换导致P值分辨率为0.005，不能区分所有达到该下限的切片谁更显著，应看I与连通性。

这表明混合spot程序有空间结构，不能证明离散C5 niche。严格top10%集合最大连通分量占选中spot比例的切片中位数为{fmt(co[co.top_fraction==.1].largest_fraction_high.median())}，最小{fmt(co[co.top_fraction==.1].largest_fraction_high.min())}，因此不少高分集合仍是多个分散区域。逐切片5%/10%/20%高分连通性见[组件表](03_spatial_clustering/components.tsv)。
3. **signature是否一致？** A26与core15、五基因、AUCell的逐切片相关分布见下表，包含差异而非只选最好切片。共享基因会产生部分一致性，不能当作独立验证。
4. **是否主要由Fibroblast驱动？** A26-high落在Fibroblast最高四分位的比例，切片中位数{fmt(s.A26_high_in_fib_top_quartile.median())}；落在最低四分位为{fmt(s.A26_high_in_fib_bottom_quartile.median())}。基线四分位各约25%，但丰度并列会改变比例。Fibroblast相关及调整UMI/上皮后的偏相关均保存；仍有{len(fibneg)}张呈负相关：{', '.join(fibneg.section_id)}，说明不能把所有切片解释为同一CAF模式。匹配成功率切片中位数{fmt(b[b.variable=='Fibroblast'].match_rate.median())}，Fibroblast匹配后绝对SMD中位数{fmt(b[b.variable=='Fibroblast'].standardized_mean_difference.abs().median())}、UMI为{fmt(b[b.variable=='nCount'].standardized_mean_difference.abs().median())}；但Fibroblast绝对SMD最大仍为{fmt(b[b.variable=='Fibroblast'].standardized_mean_difference.abs().max())}，不能将所有切片视为充分匹配。匹配不是细胞类型特异表达去卷积；不能凭残余结构认定C5真实存在。
5. **深度伪影？** UMI/nFeature关系见下表和逐dataset表。调整Fibroblast、UMI、nFeature和上皮后的rank残差，{int(((res.moran_I>0)&(res.FDR_within_score_63<.05)).sum())}/63张仍呈正自相关，残差Moran I中位数{fmt(res.moran_I.median())}。这排除了“全部空间结构只是这些变量线性rank关系”的简单解释，不能排除非线性、遗漏组织因素或共同程序影响。
6. **C2/C4竞争重叠？** 单细胞预定义联合C2/C4程序与A26高分集合重叠，切片中位数{fmt(s.C2C4_high_overlap.median())}，范围{fmt(s.C2C4_high_overlap.min())}–{fmt(s.C2C4_high_overlap.max())}。这是组合程序诊断，不能分别识别真实C2/C4，也不能用低重叠推翻既有C4困难负对照结果；对比程序仅覆盖所选20个基因。
7. **上皮/内皮/髓系空间关系？** 上皮mean占9大类总量≥0.3作为主定义，0.2/0.4敏感性；距离为观察spot图步数，非微米。主阈值可估计的切片{int(ep.effect.notna().sum())}张，A26与距离相关中位数{fmt(ep.effect.median())}。距离相关为正表示远离上皮富集区的score较高，为负则相反。内皮/髓系提供同spot匹配对比及ring0–3邻域曲线；下表效应正值表示高分区域周围高于匹配对照，负值表示低于对照。Fibroblast/UMI匹配的匹配率和SMD公开；相邻ring与两条件可能重叠，不能当独立重复，更不能解释为通讯/因果。
8. **跨来源/患者重复？** 下列逐dataset诊断和患者聚合CI展示方向是否一致。所有患者内多切片先聚合；未知患者只作切片描述。患者CI是探索性bootstrap，没有使用spot独立性；只重采患者，不包含dataset级抽样不确定性，已知患者仅3个来源，不能替代可靠的dataset层级推断。Normal切片患者ID未确认，无法完成可靠的患者配对Normal/Tumor检验。Normal仅来自GSE285715，不能把合并组差异解释为C5丰度变化；同来源GSE285715的3张Normal切片score中位数的中位数为0.099、3张Tumor为0.022，属于未配对描述，不能据此得出真实C5丰度在Normal更高，也不支持预设的Tumor必然更高。各组表见`07_dataset_heterogeneity/dataset_tissue_summary.tsv`。不比较固定top10%比例。留一dataset结果见[患者效应区间](07_dataset_heterogeneity/patient_effect_intervals.tsv)。
9. **值得独立验证的pattern？** 更值得验证的是多来源重复、调整后仍保留的A26-associated基质程序结构及其C5/C4细胞来源；与上皮成分分隔可作为描述性取样线索，而不是预设的血管/髓系富集机制。所有7个来源组的原始上皮距离相关中位数为正，但条件分析明显衰减。 可优先从同时具有A26连续空间结构、调整后结构仍保留、不同signature一致且多个来源重复的区域开展外部验证，但不能据此挑选/重训signature。本轮提供全部切片供人工评估，不重新给“正式C5定位PASS”。后续价值在于独立C5/C4/C2参考与多标记原位验证，检验这些区域是C5、其他CAF或共享ECM状态。

## 相关性与dataset异质性

以下均为逐切片相关的汇总，未将265,002 spots作为独立患者。

{md(corr)}

{md(dataset)}

## 匹配邻域的逐切片描述汇总

{md(ring_summary)}

## 已知患者聚合的探索性效应

同spot丰度使用log1p差，距离使用Spearman系数；不能混用单位。每位患者先汇总重复切片，bootstrap 1000次。CI没有校正跨终点多重比较，因此不作为确认性显著发现；本轮仅对空间置换P值按score内63切片BH。小于3位患者不计算CI。

{md(main[['analysis','variable','n_patients','median_effect','ci_low','ci_high','status']])}

## 上皮组成敏感性 (诊断后补充)

由于主匹配仍有上皮组成不平衡，追加了Fibroblast+UMI+Epithelial三变量匹配，未改signature或挑选空间图。下面既给全部可匹配切片，也给match rate≥0.5且三项绝对SMD≤0.1的平衡子集；后者是探索性敏感性，不能选取其中更有利的CI替代主分析。

{md(sens[sens.dataset=='ALL'])}

上皮距离的逐切片相关中位数，控制Fibroblast/UMI/nFeature后为{fmt(partial[partial.controls=='Fib_UMI_features'].partial_spearman.median())}，再控制同spot上皮丰度后为{fmt(partial[partial.controls=='Fib_UMI_features_Epi'].partial_spearman.median())}。这提示原距离关联很大程度可由局部组成解释；不能声称发现独立的上皮边界C5 niche。内皮/髓系负关联在额外上皮匹配后也明显减弱，全部匹配已知患者CI跨0。远端ring的正差异还受两组邻域重叠及残余组成影响，不支持C5特异促血管/免疫因果机制。

## 基因覆盖与少数基因贡献

签名-基因-切片记录中缺失特征{len(missing)}条，完整表记录每个基因的实际检测率。评分仅使用present genes；如缺失会改变分母，必须结合覆盖表解释。A26在63张均为26/26特征覆盖，但每切片每spot检出基因数的中位数范围为{fmt(coverage_summary[coverage_summary.signature=='A26'].median_detected_genes.min())}–{fmt(coverage_summary[coverage_summary.signature=='A26'].median_detected_genes.max())}；特征覆盖不等于实际表达检出。具体例外：GSE238145_B5竞争程序缺ADH1B/MT1M (18/20)，GSE285715_N2及T3的marker3缺COL10A1 (2/3)，GSE285715_T3 core15缺SFRP4 (14/15)。单基因MMP11只是诊断，不和多基因证据等价。

完整的UCell排名权重重建还用于逐切片核对官方评分实现，误差表见`00_inputs_and_QC/UCell_rank_reconstruction.tsv`。

在A26-high区域，按UCell正向rank权重描述，最主要的中位贡献基因为：{', '.join(f'{g} ({v:.1%})' for g,v in shares.head(6).items())}。每切片前三个基因贡献合计中位数{dominance.median():.1%}。这是排名权重分解，不是基因因果效应，也不等于RNA贡献；未因此修改A26。详见[A26贡献表](04_fibroblast_and_depth_controls/A26_gene_rank_contributions.tsv)。

## 方法边界与结论层级

**可靠的描述性事实**：输入完整、分数可复现、所显示的连续程序和观察图上的聚集/相关/重叠可直接由保存的数据复算。颜色共用尺度，但不同dataset绝对score仍未校准。

**探索性关联**：Fibroblast/测序深度条件下残余结构、上皮距离及内皮/髓系邻域差异、组织标签分组与患者聚合方向。受空间自相关、匹配不完全、组织区域未知、非线性组成混杂、来源标签及患者数量限制。

**无法确认**：真实C5细胞丰度、C5独有空间niche、恶性侵袭前沿、细胞通讯、C5导致血管生成/免疫抑制。已有低比例混合spot及跨dataset泛化限制不因本轮63张图而消失。正式定位NO-GO保持；描述性程序分析已完成。

复现参数、原始版本来源、图像坐标、检验集合见[README](README.md)、冻结基因列表及脚本。输入与历史文件哈希复核通过。
'''
(O/'COMPREHENSIVE_REPORT.md').write_text(report)
for directory,title,body in [
 ('00_inputs_and_QC','Inputs and QC','63 sections, 265002 unique spots; 7 pilot scores exactly reproduce previous five/three-gene UCell and AUCell. Source histories unchanged. See final_QC.tsv and per-gene coverage.'),
 ('03_spatial_clustering','Spatial pattern',f'{len(clusters)}/63 sections meet the exploratory positive Moran/BH flag. This is not C5-specific evidence. See moran.tsv and components.tsv for all sections and thresholds.'),
 ('04_fibroblast_and_depth_controls','Composition and depth diagnostics','Continuous correlations, rank partial correlations, matched-control balance and individual gene rank contributions. No cell-specific expression deconvolution.'),
 ('05_C4_C2_competition','CAF competitor program','Combined C2/C4 signature was defined solely from the frozen single-cell reference using the existing v2 formula. No spatial gene selection; scores do not estimate C2/C4 abundance.'),
 ('06_spatial_relationships','Exploratory associations','Epithelial graph distances and matched Endothelial/Myeloid ring summaries are descriptive. Same-spot cooccurrence, proximity and true interaction are distinct. No causal inference.'),
 ('07_dataset_heterogeneity','Patient and dataset summaries','Known patient repeated sections aggregated before bootstrap; unknown patients excluded from patient inference. Normal/Tumor labels strongly confounded with dataset. Uncertainty is exploratory.')]:
 (O/directory/'REPORT.md').write_text('# '+title+'\n\n'+body+'\n\nSee ../COMPREHENSIVE_REPORT.md for quantitative interpretation and limitations.\n')
# Portable local HTML gallery; no remote assets or tracking.
parts=['<!doctype html><meta charset="utf-8"><title>C5-associated all63 gallery</title><style>body{font-family:Arial;margin:24px;background:white}article{margin:25px 0}img{max-width:100%;height:auto}a{margin-right:15px}</style><h1>C5-associated signatures: all 63 sections</h1><p>Exploratory only. Formal C5 localization remains NO-GO. Raw score panels share scales; percentile panel is within-section only. HE overlays are geometry-checked, histology-unverified.</p><a href="COMPREHENSIVE_REPORT.md">Comprehensive report</a><a href="figures/all63_percentile_atlas.png">All-section percentile atlas</a>']
for (ds,tissue),group in s.sort_values(['dataset','tissue','section_id']).groupby(['dataset','tissue'],sort=False):
 parts.append(f'<h2>{html.escape(ds)} | {html.escape(tissue)}</h2>')
 for sid in group.section_id:
  q=html.escape(sid);parts.append(f'<article><h3>{q}</h3><a href="02_spatial_maps/{q}_overview.pdf">Vector PDF</a><a href="02_spatial_maps/{q}_overview.png">300 dpi PNG</a><a href="02_spatial_maps/{q}_HE_provisional.png">Provisional HE</a><a href="01_spatial_scores/{q}.parquet">Spot scores</a><img loading="lazy" src="02_spatial_maps/{q}_overview.png"></article>')
(O/'GALLERY.html').write_text('\n'.join(parts))
(O/'00_inputs_and_QC/FINAL_PASS.json').write_text(json.dumps(dict(sections=63,spots=265002,all_outputs_validated=True,formal_localization='NO_GO_UNCHANGED'),indent=2))
print('FINAL QC PASS',flush=True)
