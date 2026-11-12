# 验证报告

数据导出验证 PASS；继承原模型生物学解释等级 PASS_WITH_WARNING。未改动源结果。

|检查|实测结果|
|---|---|
|切片完整性|63，与冻结清单、final gate 一致；63 仅为背景预期|
|spot 数|265,002；与原始 counts、posterior、mean/q05 和 gate 一致|
|预期差异|较 262,520 多 2,482；该预期的来源/筛选依据未确认，无删减|
|每片 metadata/mean/q05|全部数量、行序一致|
|唯一 ID|无重复；dataset+section+barcode，跨切片相同 barcode 保留|
|30 亚型|名称与列顺序全部一致|
|NA/Inf/负值|两个丰度矩阵全部为 0 个异常|
|坐标|63/63 原始 x/y 有限且与 posterior 相等；阵列坐标精确匹配|
|Parquet 对 CSV|逐值严格相等，最大误差 0|
|CSV 对 H5AD|mean 最大绝对误差 2.578124997e-05；q05 2.246093754e-05；rtol=1e-6, atol=1e-7 通过|
|counts|逐块有限、非负、整数；原始 barcode/基因无重复|
|参数|每片 model_config 与冻结配置一致|
|可移植性|Parquet/TSV 标准格式，validate-only 仅相对包路径；源数据路径仅用于后续表达访问|
|CPU|串行，数值库与 Arrow 各限制 1 线程|
|磁盘预检|开始可用 11343.070 GiB；预检至少 2 GiB|
|源 mean+q05 gzip|0.074 GiB|

逐切片证据见 section_validation.tsv；来源文件 checksum 见 source_file_manifest.tsv；文件清单记录所有交付文件的大小、内容、必需性、验证状态与 SHA256。输出占用见 storage_summary.json。输出包含 float64 Parquet，未强制降精度，因此大小可高于原 float32 文本压缩文件。

限制：patient/sample 不明项保持空；in_tissue 是既有 QC 保留标记；坐标保持原尺度；旧注释、训练 checkpoint、reference、完整表达矩阵均未复制。原训练软件版本为冻结记录而非本次环境检查；未重新评价生物学可信度。
