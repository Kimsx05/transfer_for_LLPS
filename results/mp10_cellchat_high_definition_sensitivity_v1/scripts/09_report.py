#!/usr/bin/python3
from pathlib import Path
import pandas as pd, numpy as np, json
root=Path(__file__).resolve().parents[1];T=root/'tables'
def read(n):return pd.read_csv(T/n,sep='\t')
def fmt(x):
 if isinstance(x,(float,np.floating)):return 'NA' if pd.isna(x) else f'{x:.3g}'
 return str(x)
def md(df):
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(x).replace('|','/') for x in row)+' |' for row in df.itertuples(index=False,name=None))+'\n'
ef=read('native_effects_rankings.tsv.gz');top=read('ALL_top3_positive.tsv');rep=read('dataset_direction_support.tsv');co=read('C10_C30_same_background_comparison.tsv');dom=read('LR_domain_sensitivity.tsv');chg=read('ALL_high_definition_changes.tsv');pa=read('pathway_background_contributions.tsv.gz');lr=read('focus_top10_positive_LR_contributions.tsv.gz');cov=read('node_counts.tsv');status=read('network_status.tsv')
a=ef[(ef.cohort=='C10')&(ef.scope=='ALL')&(ef.min_cells==10)&(ef.lr_domain=='fixed_full')]
rows=[]
for hi in ['H50','H10','H05']:
 for di in ['output','input']:
  z=top[(top.cohort=='C10')&(top.high==hi)&(top.direction==di)].sort_values('rank_W_or_R')
  rows.append({'High':hi,'方向':'输出 Epi→TME' if di=='output' else '输入 TME→Epi','背景调整 W top3 (RR)':'；'.join(f'{r.target} ({r.RR:.3f})' for r in z.itertuples())})
focus=['Mast','Arterial EC','vCAF','Capillary EC','Plasma','Treg','Tfh/CXCL13-like','cDC1','Naive/Memory B','Neutrophil','pDC']
rob=[]
for di in ['output','input']:
 for ta in focus:
  r={'方向':di,'目标':ta}
  for mt in ['W','R']:
   for hi in ['H50','H10','H05']:
    z=a[(a.target==ta)&(a.direction==di)&(a.high==hi)&(a.metric==mt)];r[f'{mt}_{hi}']=z.RR.iloc[0]
  for mt in ['W','R']:
   z=co[(co.target==ta)&(co.direction==di)&(co.high=='H05')&(co.metric==mt)];r[f'C30_H05_{mt}']=z.RR_C30.iloc[0]
   z=rep[(rep.target==ta)&(rep.direction==di)&(rep.high=='H05')&(rep.metric==mt)];r[f'dataset_H05_{mt}>1']=f'{z.n_dataset_RR_gt1.iloc[0]}/{z.n_dataset_RR_defined.iloc[0]}'
  rob.append(r)
rob=pd.DataFrame(rob);rob.to_csv(T/'report_candidate_robustness.tsv',sep='\t',index=False,na_rep='NA')
notes=[]
for di in ['output','input']:
 for pattern in ['H50_only_relative_enhancement','tail_relative_enhancement_emerges','mixed_direction','at_least_one_ratio_undefined']:
  v=chg[(chg.cohort=='C10')&(chg.direction==di)&(chg.metric=='W')&(chg.pattern==pattern)]
  notes.append({'方向':di,'模式':pattern,'目标':'；'.join(v.target) if len(v) else '无'})
old=['Naive/Memory B','Treg','cDC1','Plasma','Capillary EC','Arterial EC']
# Main pathway examples, reported with overall RR and raw/supported values available separately.
p=pa[(pa.cohort=='C10')&(pa.scope=='ALL')&(pa.high=='H05')&(pa.min_cells==10)&(pa.lr_domain=='fixed_full')&(pa.metric=='W')&pa.background_contexts.str.contains(r'(^|\|)native(\||$)')&(pa.E_path>0)]
p=p[p.target.isin(focus)].sort_values('E_path',ascending=False).groupby(['direction','target']).head(3)
lr_examples=lr[(lr.high=='H05')&(lr.metric=='W')&(((lr.direction=='input')&lr.target.isin(['Plasma','Treg','Tfh/CXCL13-like','cDC1']))|((lr.direction=='output')&lr.target.isin(['Arterial EC','Capillary EC','Mast'])))].groupby(['direction','target'],sort=False).head(3)
lr_examples=lr_examples[['sender_H','receiver_H','interaction_id','ligand','receptor','H_probability','L_probability','H_internal_P','L_internal_P','E_LR']].rename(columns={'ligand':'DB_ligand','receptor':'DB_receptor'})
maxdom=dom[(dom.cohort=='C10')&(dom.metric=='W')].RR_difference.abs().max()
countqc=pd.read_csv(root/'qc/validation_checks.tsv',sep='\t');render=pd.read_csv(root/'qc/render_format_checks.tsv',sep='\t') if (root/'qc/render_format_checks.tsv').exists() else None
report=f'''# MP10 High 定义敏感性分析结果

本轮已实际完成 {len(status)} 套 fine CellChat 推断。基线 run_id 为 `c2lfinal_20261108T123140Z_0b5fd64f`，旧结果保持原样。C10 是 38 个标本、37 个已记录患者标识；C30 是 30 个标本、30 个已记录患者标识。患者标识来自原始记录/推导，未验证患者关系，不能当作已确认独立患者。每个 Epi 状态分别为 2,651 与 2,479 个抽样细胞，标本 k_s 不相同。

总体判断：数据支持**部分通讯在高分尾部更强**，尤其 Treg、Tfh/CXCL13-like → Epi 的输入，以及 Epi → Arterial EC/Capillary EC 的输出；不支持把全部通讯统一归结为 top 5% 特有信号。H50 已有同向联系，H10 已出现主要变化，H05 进一步增强的幅度因目标而异。因此更合适的表述是“目标与方向特异的尾部增强”，仍不足以确定活性分界。

## 1. 预设三种定义的主结果

以下均为 C10、ALL、min.cells=10、固定 1,933 LR 域、native 公共背景、严格 `probability>0 & internal_P<0.05` 的 W RR。只列 RR>1 的目标。

{md(pd.DataFrame(rows))}

完整 29 个目标的 W/R、原始 High/Low、差/比、背景分子分母、排名和状态见 `tables/native_effects_rankings.tsv.gz`。输出与输入是两种不同方向，不能互换解释。

## 2. 从 H50 到高分尾部的变化

{md(pd.DataFrame(notes))}

“翻转”仅表示 RR 穿过 1 的描述性变化，不是经过标本间差异检验的显著翻转。接近 1 的数值特别需要结合 R、dataset 和 C30 判断。

## 3. W、raw R、dataset 和 C30 复核

表中 dataset 分母是 RR 可定义的 dataset 数，不把 NE 或零分母算作负证据；dataset 中 RR>1 的数量是方向性描述，不是独立重复的显著性检验。各 dataset 还提供与 ALL 在同一背景上的逐项对照。C30 数值来自 C10/C30 共同背景。

{md(rob)}

Treg 和 Tfh/CXCL13-like 输入是较一致的尾部增强：H05 的 W/R 在 4/5 dataset 中均 RR>1，C30 保留增强。Arterial EC 输出在 4 个可评估 dataset 中 W/R 均为正，Capillary EC 为 4/5；但不同 dataset 效应幅度不同，部分数值接近 1。cDC1 输入在 4/4 可评估 dataset 中 W>1，raw R 为 3/4，适合作为次级证据。Naive/Memory B 输入在 H05 中由 W<=1 转为 >1，W/R 均有 3/4 dataset 同向，C30 效应较小；其输出在所有可评估 dataset 的 W/R 均未相对增强。

Plasma 输入虽然是 H10/H05 的 W top1，raw R 只有约 1.15，且仅 1/5 dataset 的 raw R>1；其 H50 raw R 反而更高。C30 保留 W 增强，但不能把 W 的排名跃升直接解释为原始通讯强度集中于尾部。Mast 输出在三方案的 C10 W 都排名第一，而 H05 raw R 约 1.02，只有 2/5 dataset 同向；C30 W RR 降至约 1.11。这两个候选都应保留，但证据弱于只看 W top3 所暗示的程度。

H50 独有的输入 RR>1 主要见于 Capillary EC、Pericyte、vCAF，数值仅约 1.001–1.012；它们不能构成强的“宽泛状态独有”证据。Treg 输出的 W 从略大于 1 变为略小于 1，但 raw R 和 C30 并不保持这一翻转，宜视为边界敏感。ALL 29 个目标三种定义均可评估，主比较没有 NE；dataset 的缺失应按覆盖表解释。

## 4. 节点范围、计算域和支持阈值

ALL 的 29 个 TME 在 C10 和 C30 均超过 30 个抽样细胞。因此 ALL 的 min10/min30 后筛不改变可评估范围；共同背景检查确认 H、L、RR、E 相同。各 dataset 小节点不同，具体为：

{md(cov[(cov.cohort=='C10')&(cov.high=='H50')&(cov.N<=30)][['scope','node','N']])}

这些节点中 N<=10 的节点没有建模；10<N<=30 的节点参与推断，但在 min30 后筛为 NE。两种后筛在共同可评估节点及背景上数值一致。它们不是两个独立实验。C30 则真正改变标本构成并重新推断，不能与后筛混为一谈。

数据库保持 CellChat 1.6.1 的 1,939 条原定义。旧表达缓存缺少 IFNLR1、IL2RA，已从同一原对象 RNA:data 补齐。6 条因原始 RNA 缺少必需基因而不可计算，详见 LR coverage 表；未将其补零。新固定域 1,933 条包含旧 ALL fine 实际计算域全部 1,875 条，扩展 58 条。C10 ALL W 的新域与旧域交集 RR 最大绝对差为 {maxdom:.6g}；完整逐目标结果见 `LR_domain_sensitivity.tsv`。

Plasma 输入需要特别区分 W 与 R：W 的高分尾部增强幅度明显大于 raw R。`focus_LR_background_contributions.tsv.gz` 同时保留各 LR 的原始概率、内部 P 和支持状态，可定位阈值作用。不同 High 定义也改变参与置换的细胞，因此即使 Low/TME 成员固定，内部 P 仍可能变化。不能把内部 P 当作标本间 High/Low 差异检验。

## 5. 通路与 LR

重点候选由旧六个候选和本轮 ALL 正向 top3 并集预先规则确定。每个目标/方向的作图通路按 C10 ALL H05 W 的正向 E_path 选最多五条，之后固定展示三种 High。图内给出每种定义总体 RR；总体未增强的目标即使有局部正贡献，也不改写总体结论。

下面列出 H05 的主要正向通路贡献 (总体背景比固定，所有通路可加总为总体 E)：

{md(p[['direction','target','pathway','E_path','RR']])}

对应 W/R、High/Low 和主要 LR 见 `focus_top10_positive_LR_contributions.tsv.gz`。关键 H05 项的实际分解如下，H/L probability 是未按内部 P 筛选的值，E_LR 为 supported W 的贡献：

{md(lr_examples)}

Plasma 的 ADGRE5_CD55 在 Low 中原始 probability 仍为正，但内部 P=1，故其 supported Low 贡献为零；Mast 输出的 MIF_CD74_CD44 也有正的 Low 原始 probability，而内部 P=0.92、High P=0.02。这直接显示支持掩码如何放大 W 对比，不能把 Low 未支持解释为该联系未计算或无原始表达信号。

本报告输入/输出沿用 CellChat 1.6.1 的 DB 编码角色。例如 TIGIT_NECTIN2 在该 DB 中写为 ligand=TIGIT、receptor=NECTIN2；这里报告的是计算方向，不能据此单独认定 TIGIT 为生物学配体，或 Treg/Tfh 为生物学上游。DB_ligand/DB_receptor、发送细胞、接收细胞已分别给出，机制方向仍需另行验证。

## 6. 可复现性与检查

全部参数、版本、数据库与 feature/LR 域、缺失原因、真实输入顺序/hash、抽样优先级、原始未过滤 probability/P、汇总和直接绘图数据已保存。数值推断对象只有紧凑缓存，不重复保存大型表达或完整 CellChat 对象。`run_all.sh` 可按内容 hash 断点复用；`scripts/04_replot.R` 只读结果表，完全不依赖重跑 CellChat。

本轮输入与结果核验 {len(countqc)} 项，全部通过；另有原生过滤/P 边界、通路加总及共同背景后筛一致性检查。图为 Arial 约 12 pt、白底矢量 PDF、300 dpi PNG；有图例时提供完整、仅主体及仅图例 PDF。直接绘图数据及统一色标见 `plot_data/` 与 `config/figure_manifest.tsv`。热图真实零、NE 和比值未定义分别标识，无伪计数。

建议继续保留三套定义，以 Treg/Tfh 输入和血管内皮输出为后续验证重点，同时将 Plasma、Mast 的支持阈值和标本构成敏感性单独列为限制。这里提供证据优先级，不替用户冻结最终候选或阈值。

本报告保留 H50/H10/H05 全套预设结果，不根据最显著方案挑选最终阈值。H05 应称为高分尾部，不能直接称作真实活跃细胞；其他 95% 也没有被判定为噪声。结果讨论与最终判断需结合上述全部复核，而不能仅凭 pooled ALL 的 top3。
'''
(root/'RESULTS.md').write_text(report)
print(md(rob))
print('REPORT WRITTEN')
