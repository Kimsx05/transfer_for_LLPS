# LLPS spatial transfer v1

可直接复用 63 张切片、265,002 spots 的 v2_30type 原始 posterior mean/q05，30 亚型，零生物学重分析。审查详情见 00_metadata/annotation_provenance.md；验证见 03_validation/validation_report.md；完整文件大小和 SHA256 见 03_validation/file_manifest.tsv。

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
mapping = pd.read_csv(root / '00_metadata/celltype_mapping.tsv', sep='\t')
fibro = mapping.loc[mapping.major_lineage.eq('Fibroblast'), 'subtype']
x['Fibroblast_mean'] = x[fibro].sum(axis=1)  # 允许 mean 求和
# 不把 subtype q05 的和当作 Fibroblast posterior q05；不赋予 C0/C5 身份。
```

## 连接空间表达矩阵

63 张 raw-input H5AD 均具有整数非负 raw counts in X、原始基因名、完全匹配的 obs_names 和坐标。按切片处理；raw H5AD 保留完整空间特征集，posterior H5AD 已限于 reference 交集，不建议作为 LLPS 打分的默认表达源。基因标识保持 exact names，逐片 Ensembl 数量在 manifest 中；不擅自更名或大小写折叠。

```python
import anndata as ad  # 仅表达矩阵读取需要；使用支持当前 H5AD 编码的版本
inputs = pd.read_csv(root / '02_spatial_inputs/spatial_input_manifest.tsv', sep='\t')
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
|63 个 raw-input H5AD|5.475 GiB|LLPS 空间打分、表达差异、表达相关 NicheNet 输入需要|
|63 张低分辨率 HE PNG|0.185 GiB|图像叠加可选；邻域计算不需要|
|posterior H5AD|2.968 GiB|无需迁移，丰度与坐标已导出|
|完整 Seurat combined .qs|4.256 GiB|无需迁移；包含嵌入图像等原始资产|
|单细胞 reference / model checkpoint /训练日志|未复制|本次注释复用不需要|

同服务器首选按清单以只读模式访问原始 H5AD，新增占用仅本迁移包。跨服务器最小完整表达方案：迁移本包 + 63 个 raw-input H5AD（5.475 GiB），可选 PNG（0.185 GiB）；更新新项目自己的输入路径配置。更节省的形式是每切片 sparse counts + genes + barcodes，复用本包坐标/scale factors，但尚未进行重打包，不能承诺压缩后体积。

独立 Space Ranger 文件夹未在既有 preflight 中找到；当前确认的可用来源是 per-section H5AD 与 combined Seurat。PNG 实际位置已验证，历史 HE manifest 的旧目录名已失效。高分辨率/全分辨率原图未单独打包或核实可用性；需要病理精细分析时应另行定位，不能用 lowres 替代。

立即可做：原 30 亚型丰度展示、Fibroblast mean 汇总和空间分布、基于阵列坐标的切片内邻域分析。同服务器读取清单中的 counts 后可开展 LLPS score、LLPS-high 识别及表达分析；本包自身不含 counts，跨服务器这些步骤需先迁移表达输入。NicheNet 还需后续分析设计/配体靶基因资源；C0/C5 空间定位需单独评估 signature 或参考细化，本包不提供 C0/C5 对应关系。
