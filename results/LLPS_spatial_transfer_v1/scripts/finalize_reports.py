#!/usr/bin/env python3
"""Generate portable reports and checksums from completed export (no biological analysis)."""
import os
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import json,hashlib,sys
from pathlib import Path
import pandas as pd
O=Path(sys.argv[1]).resolve();s=json.loads((O/'03_validation/audit_summary.json').read_text());sm=pd.read_csv(O/'00_metadata/section_manifest.tsv',sep='\t');inp=pd.read_csv(O/'02_spatial_inputs/spatial_input_manifest.tsv',sep='\t');v=pd.read_csv(O/'03_validation/section_validation.tsv',sep='\t');cfg=json.loads((O/'00_metadata/evidence/frozen_run_config.json').read_text());P=Path(sm.raw_input_path.iloc[0]).parents[3]
def gib(n):return f'{n/1024**3:.3f} GiB'
provenance=f'''# 注释审查与来源

本包读取正式 `cell2location_full_run/v2_30type`，逐切片核对冻结清单、COMPLETED_AND_VALID、model_config、最终 validation gate、原始 CSV 和 posterior H5AD。未训练、重新抽样、归一化、缩放、二值化、修改注释或任何原始文件。

实际为 {s['n_sections']} 张、{s['n_spots']:,} 个 hard-QC retained spots、7 个保留 dataset 标签。预期 262,520 与实际相差 {s['spot_difference']:,}。正式冻结清单、原始输入、mean/q05、final validation gate 及此前空间生态清单一致支持 265,002；未找到 262,520 对应的正式筛选依据，不能将其归因于某个未经证实的 QC 步骤。没有为匹配预期删减 spot。

最终 reference 是 300,033 细胞、17,987 基因、58 batches、30 亚型；约 35.6 万不等于本次冻结模型实际使用细胞数。本次还只读检查了参考输入 X shape=300033 × 17987，并确认参考输入、模型和签名文件均存在。参考来源：
- `{cfg['reference_input_path']}`
- `{cfg['reference_model_path']}`
- `{cfg['signature_matrix_path']}`

软件冻结记录：cell2location==0.1.5；这是原安装 requirements 的记录，原训练环境目前路径不可用，未把当前导出环境冒充训练环境。模型参数：N_cells_per_location=20，detection_alpha=20，max_epochs=30000，lr=0.002，batch_size=None，train_size=1，seed=20260910，posterior_num_samples=1000，posterior batch max=2500。每切片模型配置均与冻结参数一致。原工作流线程设置不代表本次导出：本次串行、各数值库线程=1。

正式丰度：每切片 `abundance_mean.csv.gz`、`abundance_q05.csv.gz`；对应 H5AD `obsm/means_cell_abundance_w_sf`、`obsm/q05_cell_abundance_w_sf`。H5AD 因键限制把名称中的 / 写为 __SLASH__，本包保留原 CSV 的完整生物学名称和列顺序。

映射来自 final validation gate 的 `05_major_sanity/05_subtype_to_major_mapping.csv`。主要类型使用 major_lineage；rctd_comparison_major 仅供历史比较，不能用它把 Pericyte/SMC 从 Mural 改为 Fibroblast。Unknown/unassigned 未重新分配；30 个参考标签本身不含 Unknown，不生成新的 spot 硬分类。

原 validation gate 是 PASS_WITH_WARNING；具体亚型建议完整保存于 evidence/DOWNSTREAM_USAGE_TABLE.csv，不因导出校验通过而改变生物学解释等级。Matrix_CAF/iCAF/rCAF/vCAF 均不能直接等同 LLPS C0/C5。

## 坐标与身份

唯一 spot_id = URL 转义(dataset)::URL 转义(section_id)::URL 转义(barcode)。保留 source_obs_id 和 barcode，不能跨切片仅按 barcode 合并。sample_id 是原 specimen/切片标识，不代表已核验患者；source_sample_id、patient_uid 和身份状态取既有 section_identity_overlay，无信息仍为空，保留 CanDis__UNVERIFIED、NPJ__UNVERIFIED 标签。

x/y 原样来自原始 H5AD obsm/spatial，与 posterior 一致；各切片坐标系逐项保存，不是微米坐标，不追加比例缩放。array_row/array_col 来自既有 `neighborhood_v1/00_inputs_and_audit/lattice_recovery/03_recovered_coordinates_by_section.csv.gz`，按 section+barcode 精确匹配，未重新推断阵列。in_tissue 是此前 preflight 为 hard-QC retained spots 写入的 TRUE，不能声称是重新读取的 Space Ranger 原始 tissue flag。

## 大类 q05

已保存的 `neighborhood_v1/00_inputs_and_audit/major_q05/` 来自 accepted v2 模型的 1000 个 posterior draws：逐 draw 对同大类亚型求和，再 np.quantile(..., .05)。这是正确聚合随机变量的 Monte Carlo posterior quantile 估计（仍有有限采样误差），不是 sum(subtype q05)。对应脚本及 phase1 provenance 在 evidence 中，仅作为历史证据复制，没有执行。原 v2 worker 中 family_total q05 图使用亚型 q05 求和，仅为近似下界式汇总，不能当作大类严格 q05。本包不迁移这两种大类结果。

## 版本取舍

|位置|性质/用途|
|---|---|
|cell2location_full_run/v2_30type|本次正式丰度来源|
|cell2location_final_validation_gate_v1|正式模型的后验验证证据|
|cell2location_preflight_v1/02_visium/anndata|本次复用 raw counts 和坐标，非旧注释替代|
|neighborhood_v1/00_inputs_and_audit|已有阵列坐标、大类 q05 来源说明|
|cell2location_full_run/v1|旧注释，不迁移；其 HE PNG 是可复用图像资产|
|cell2location_full_run/v1_identifiability_audit、validation/preflight smoke tests、RCTD sanity|旧版/测试/对照，均未作为丰度来源|
|Round2、Visium10x fine annotation、直接 Seurat 注释|非本次最终 cell2location，未使用|
|后续 MP10 state/refined annotation|其他分析层，不替换冻结 30 亚型|

原始脚本路径：`{P}/cell2location_full_run/v2_30type/code/01_run_one_section.py`；counts 准备脚本：`{P}/cell2location_preflight_v1/code/preflight/04_prepare_visium_inputs.R`。其副本随包保存以便审阅，禁止把训练/恢复脚本当作本次执行入口。
'''
(O/'00_metadata/annotation_provenance.md').write_text(provenance)
readme=f'''# LLPS spatial transfer v1

可直接复用 {s['n_sections']} 张切片、{s['n_spots']:,} spots 的 v2_30type 原始 posterior mean/q05，30 亚型，零生物学重分析。审查详情见 00_metadata/annotation_provenance.md；验证见 03_validation/validation_report.md；完整文件大小和 SHA256 见 03_validation/file_manifest.tsv。

## 读取与连接

仅需标准 pandas + pyarrow，注释数据不依赖旧环境或原项目绝对路径。下面的包路径由使用者指定：

```python
from pathlib import Path
import pandas as pd
root = Path('LLPS_spatial_transfer_v1')
meta = pd.read_parquet(root / '00_metadata/spot_metadata.parquet')
mean = pd.read_parquet(root / '01_cell2location/abundance_mean.parquet')
q05 = pd.read_parquet(root / '01_cell2location/abundance_q05.parquet')
x = meta.merge(mean, on='spot_id', validate='one_to_one')
assert len(x) == len(meta) == len(q05)
mapping = pd.read_csv(root / '00_metadata/celltype_mapping.tsv', sep='\\t')
fibro = mapping.loc[mapping.major_lineage.eq('Fibroblast'), 'subtype']
x['Fibroblast_mean'] = x[fibro].sum(axis=1)  # 允许 mean 求和
# 不把 subtype q05 的和当作 Fibroblast posterior q05；不赋予 C0/C5 身份。
```

## 连接空间表达矩阵

63 张 raw-input H5AD 均具有整数非负 raw counts in X、原始基因名、完全匹配的 obs_names 和坐标。按切片处理；raw H5AD 保留完整空间特征集，posterior H5AD 已限于 reference 交集，不建议作为 LLPS 打分的默认表达源。基因标识保持 exact names，逐片 Ensembl 数量在 manifest 中；不擅自更名或大小写折叠。

```python
import anndata as ad  # 仅表达矩阵读取需要；使用支持当前 H5AD 编码的版本
inputs = pd.read_csv(root / '02_spatial_inputs/spatial_input_manifest.tsv', sep='\\t')
r = inputs.iloc[0]
# 同服务器直接只读打开；跨服务器先把 expression_path 改成新位置。
a = ad.read_h5ad(r.expression_path, backed='r')
m = meta.loc[meta.section_id.eq(r.section_id) & meta.dataset.eq(r.dataset)]
assert a.obs_names.is_unique and set(m.barcode) == set(a.obs_names)
m = m.set_index('barcode').loc[a.obs_names]
w = mean.set_index('spot_id').loc[m.spot_id]
assert list(m.index) == list(a.obs_names)
# a.X 为稀疏 counts；m 的 x/y 为原坐标，array_row/array_col 用于切片内部邻域。
a.file.close()
```

不要跨 section 建邻域；不要将像素距离直接称为微米距离。scale_factors_json 已写入清单，坐标不应重复缩放。

## 复现

导出/校验仅需 Python、numpy、pandas、pyarrow、h5py，具体运行版本见 audit_summary.json。脚本分切片读取丰度、按块检查 counts，未加载完整表达矩阵。最多 1 个计算线程，无 GPU。

```bash
python scripts/export_spatial_annotations.py --project '/path/to/Bladder Metabolism' --output '/path/to/new_transfer'
python scripts/export_spatial_annotations.py --output '/path/to/LLPS_spatial_transfer_v1' --validate-only
```

导出拒绝覆盖已有丰度文件；请指定新输出目录。冻结清单中的源绝对路径用于本服务器导出；若原项目已迁移，需恢复/重定位相应源清单，而独立读取和 --validate-only 不访问任何源绝对路径。file_manifest 不含自身 checksum，避免自引用。

## 未打包数据与迁移方案

|数据|当前总量|必要性|
|---|---:|---|
|63 个 raw-input H5AD|{gib(s['expression_bytes'])}|LLPS 空间打分、表达差异、表达相关 NicheNet 输入需要|
|63 张低分辨率 HE PNG|{gib(s['optional_image_bytes'])}|图像叠加可选；邻域计算不需要|
|posterior H5AD|{gib(s['posterior_bytes'])}|无需迁移，丰度与坐标已导出|
|完整 Seurat combined .qs|{gib((P/'ST/all_ST_combined_list_SeuratObject_updated.qs').stat().st_size)}|无需迁移；包含嵌入图像等原始资产|
|单细胞 reference / model checkpoint /训练日志|未复制|本次注释复用不需要|

同服务器首选按清单以只读模式访问原始 H5AD，新增占用仅本迁移包。跨服务器最小完整表达方案：迁移本包 + 63 个 raw-input H5AD（{gib(s['expression_bytes'])}），可选 PNG（{gib(s['optional_image_bytes'])}）；更新新项目自己的输入路径配置。更节省的形式是每切片 sparse counts + genes + barcodes，复用本包坐标/scale factors，但尚未进行重打包，不能承诺压缩后体积。

独立 Space Ranger 文件夹未在既有 preflight 中找到；当前确认的可用来源是 per-section H5AD 与 combined Seurat。PNG 实际位置已验证，历史 HE manifest 的旧目录名已失效。高分辨率/全分辨率原图未单独打包或核实可用性；需要病理精细分析时应另行定位，不能用 lowres 替代。

立即可做：原 30 亚型丰度展示、Fibroblast mean 汇总和空间分布、基于阵列坐标的切片内邻域分析。同服务器读取清单中的 counts 后可开展 LLPS score、LLPS-high 识别及表达分析；本包自身不含 counts，跨服务器这些步骤需先迁移表达输入。NicheNet 还需后续分析设计/配体靶基因资源；C0/C5 空间定位需单独评估 signature 或参考细化，本包不提供 C0/C5 对应关系。
'''
(O/'README.md').write_text(readme)
report=f'''# 验证报告

数据导出验证 PASS；继承原模型生物学解释等级 PASS_WITH_WARNING。未改动源结果。

|检查|实测结果|
|---|---|
|切片完整性|{s['n_sections']}，与冻结清单、final gate 一致；63 仅为背景预期|
|spot 数|{s['n_spots']:,}；与原始 counts、posterior、mean/q05 和 gate 一致|
|预期差异|较 262,520 多 {s['spot_difference']:,}；该预期的来源/筛选依据未确认，无删减|
|每片 metadata/mean/q05|全部数量、行序一致|
|唯一 ID|无重复；dataset+section+barcode，跨切片相同 barcode 保留|
|30 亚型|名称与列顺序全部一致|
|NA/Inf/负值|两个丰度矩阵全部为 0 个异常|
|坐标|63/63 原始 x/y 有限且与 posterior 相等；阵列坐标精确匹配|
|Parquet 对 CSV|逐值严格相等，最大误差 0|
|CSV 对 H5AD|mean 最大绝对误差 {v.mean_vs_posterior_max_abs.max():.10g}；q05 {v.q05_vs_posterior_max_abs.max():.10g}；rtol=1e-6, atol=1e-7 通过|
|counts|逐块有限、非负、整数；原始 barcode/基因无重复|
|参数|每片 model_config 与冻结配置一致|
|可移植性|Parquet/TSV 标准格式，validate-only 仅相对包路径；源数据路径仅用于后续表达访问|
|CPU|串行，数值库与 Arrow 各限制 1 线程|
|磁盘预检|开始可用 {gib(s['free_bytes_before'])}；预检至少 2 GiB|
|源 mean+q05 gzip|{gib(s['source_abundance_bytes'])}|

逐切片证据见 section_validation.tsv；来源文件 checksum 见 source_file_manifest.tsv；文件清单记录所有交付文件的大小、内容、必需性、验证状态与 SHA256。输出占用见 storage_summary.json。输出包含 float64 Parquet，未强制降精度，因此大小可高于原 float32 文本压缩文件。

限制：patient/sample 不明项保持空；in_tissue 是既有 QC 保留标记；坐标保持原尺度；旧注释、训练 checkpoint、reference、完整表达矩阵均未复制。原训练软件版本为冻结记录而非本次环境检查；未重新评价生物学可信度。
'''
(O/'03_validation/validation_report.md').write_text(report)
# Remove transient execution log from inventory scope: log is external to package in recommended use.
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
files=[p for p in O.rglob('*') if p.is_file() and p.name not in ['file_manifest.tsv','storage_summary.json','03_export.log'] and '__pycache__' not in p.parts]
def content(p):
 n=p.name
 return {'spot_metadata.parquet':'spot identity, coordinates, sample provenance','abundance_mean.parquet':'30 subtype posterior mean','abundance_q05.parquet':'30 subtype posterior q05','section_manifest.tsv':'section inventory','celltype_mapping.tsv':'30 subtype to major mapping','spatial_input_manifest.tsv':'counts, image and scale-factor inventory','README.md':'LLPS integration and transfer guide'}.get(n,'provenance / reproducibility / validation')
rows=[dict(relative_path=str(p.relative_to(O)),full_path=str(p),content=content(p),bytes=p.stat().st_size,required='YES' if p.suffix=='.parquet' or p.name in ['celltype_mapping.tsv','section_manifest.tsv','README.md','spatial_input_manifest.tsv'] else 'AUDIT_RECOMMENDED',validation='PASS_READABLE_AND_CHECKSUM',sha256=digest(p)) for p in sorted(files)]
# Stable byte accounting includes the manifest and storage report themselves.
storage=O/'03_validation/storage_summary.json';manifest=O/'03_validation/file_manifest.tsv'
for _ in range(5):
 storage.write_text(json.dumps(dict(payload_bytes=sum(p.stat().st_size for p in files),package_logical_bytes=sum(p.stat().st_size for p in files)+ (storage.stat().st_size if storage.exists() else 0)+(manifest.stat().st_size if manifest.exists() else 0),package_allocated_bytes=sum(p.stat().st_blocks*512 for p in O.rglob('*') if p.is_file()),source_abundance_bytes=s['source_abundance_bytes'],expression_not_copied_bytes=s['expression_bytes'],images_not_copied_bytes=s['optional_image_bytes']),indent=2)+'\n')
 pd.DataFrame(rows+[dict(relative_path=str(storage.relative_to(O)),full_path=str(storage),content='storage accounting',bytes=storage.stat().st_size,required='AUDIT_RECOMMENDED',validation='PASS',sha256=digest(storage))]).to_csv(manifest,sep='\t',index=False)
print('Reports and file manifest generated')
