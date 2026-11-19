from pathlib import Path
import pandas as pd,html
O=Path('results/C5_Endothelial_NicheNet_spatial_v1')
# Keep plotting-table metadata synchronized after the explicit SMD constant-boundary correction.
s=pd.read_csv(O/'tables/spatial_effects.tsv',sep='\t');s[(s.top_fraction==.1)&(s.ring==1)&(s.role!='curated_program')].to_csv(O/'figures/E_sections_data.tsv',sep='\t',index=False)
r=O/'COMPREHENSIVE_REPORT.md';text=r.read_text();text+='''

## 补充：功能程序、来源稳定性与解释限制

在排除Five及入选ligand后，TGFB1与邻近内皮富集spots的标准TGF-beta程序仍有弱条件关联：患者中位rho约0.067 (95%CI 0.022–0.094)。这不是内皮细胞内SMAD响应的测量，也不是方案A失败后的替代实测DEG。
追加标准GO ECM-organization/collagen-fibril-organization去除敏感性后，TGFB1–TGF-beta程序rho约0.067 (CI 0.025–0.094)，FN1–Angiogenesis约0.228 (CI 0.161–0.304)。但本去除规则实际仅从已去重叠程序中再去除2个Angiogenesis、1个TGF-beta、6个vascular-remodeling基因，不能称已彻底排除ECM状态。完整定义、每切片结果与患者LODO见ECM_depleted*文件。这是预先要求的混杂诊断扩展，不是独立验证，也未反向改变候选排序。

WNT4在移除InhouseData后不再通过原表达候选门槛，进一步提示来源依赖；其余所列主候选在两主要来源的逐一留出表达门槛下保留。表达门槛保留不等于NicheNet响应的跨队列验证，因为本轮没有患者关联response。

单细胞患者与Visium患者没有进行未经证实的一一配对。只有三个空间来源有确认患者身份，其余来源的LODO不会提供新的患者层面独立证据。全部空间组织标签原样保留，UNRESOLVED没有改为Tumor。

### 证据分级
- 可靠表达与计算：固定资源/程序的全部配体排名、表达检测、63切片覆盖、受体亚基审查。
- 有限患者支持的探索关联：FN1共享ECM对照最稳定；TGFB1弱且依赖GSE319536。原位细胞身份和空间零模型支持尚不足。
- 先验网络和功能程序候选：HSPG2、TGFB1等配体–target关联；不能称C5诱导。
- 空间支持不足：COL3A1–ADGRG1有网络边缺口且邻接非稳定正向；LAMA1、WNT4检出/方向/来源不稳。ADM/MDK/SEMA3C/WNT2保留真实排序与覆盖，并未因不是主候选而记作阴性。

可用于小文章的是审慎的表达与描述性空间证据组合。若后续只推进一个机制，建议先补充分辨C5与其他CAF的原位证据和独立内皮响应数据，而非直接将任一当前轴写为已确立的C5特异通讯。
''';r.write_text(text)
figs=sorted((O/'figures').glob('*.png'));items=[]
for f in figs:items.append(f'<section><h2>{html.escape(f.stem)}</h2><a href="figures/{f.stem}.pdf">Vector PDF</a><br><img loading="lazy" src="figures/{f.name}" alt="{html.escape(f.stem)}"></section>')
(O/'FIGURE_GALLERY.html').write_text('<!doctype html><meta charset="utf-8"><title>C5–Endothelial exploratory evidence</title><style>body{font-family:Arial;margin:32px auto;max-width:1400px;background:white;color:#222}img{width:100%}section{margin:40px 0}h2{font-size:18px}</style><h1>C5–Endothelial NicheNet / Visium</h1><p>CURATED_EXPLORATORY. No C5-specific communication established. Spatial maps use within-section percentiles; raw values are saved. Matching balance is inadequate in most sections.</p>'+''.join(items))
(O/'README.md').write_text('''# C5–Endothelial NicheNet spatial v1

Status: analysis complete; RESPONSE_MODE=CURATED_EXPLORATORY. Formal C5-specific spatial localization NO-GO retained.

Start with COMPREHENSIVE_REPORT.md, FIGURE_GALLERY.html, tables/final_axis_summary.tsv and 06_CellChat_integration/complete_expression_activity_LR_target_evidence.tsv.

21 tumor patients meet frozen sampling criteria, but no qualifying positive endothelial response genes. Three standard programs evaluated; activation term too small. 246 broad ligands; five main axes and three ECM controls; all63 sections and265002 spots. Only2 sections achieve full matching balance (one known patient); no spatial-translation test survives global FDR.

Reproduction: run from /home/jinshengxi/bladder_metabolism using scripts/run_R.sh or scripts/run_python.sh. See numerical scripts and frozen protocol; 03d_GO extends the two Hallmark programs once, before spatial selection. 04_spatial.py resumes completed sections using DONE files. 09_balance_boundary_fix.py documents the identical-constant-covariate SMD correction (current04 already contains the corrected definition). 10/11 are standard ECM-membership program sensitivities. 07_report.py generates report; 12_finalize.py appends interpretations and builds gallery; do not repeatedly append without regenerating07 first.

Source checksums, resource versions and validation are in00_input_audit. Raw objects and historical results are read-only. No CellChat, deconvolution, Five scoring or global spatial-statistic rerun. Seed123; serial BLAS/R/Python with affinity capped at4 CPUs. Local existing human v2 NicheNet resources only.

All figures use project plotting skill styles, vector PDF and300dpi PNG. Full/panel/legend PDFs are supplied when there is a legend. Tables provide plotting values. Raw spatial expressions are log1p(CP10k); maps explicitly show within-section percentiles. No causal or cell-intrinsic interpretation.
''')
print('FINALIZED',len(figs))
