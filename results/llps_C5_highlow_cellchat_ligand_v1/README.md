# LLPS C5 High/Low CellChat 与配体关联 v1

入口：[中文综合报告](COMPREHENSIVE_REPORT.md)、[图形浏览](FIGURE_GALLERY.html)、[全文件清单](output_manifest.csv)。

复现入口：scripts/run_all.sh。原生推断由supervise_native.py调度3个单线程R worker及1个轻量监控进程，合计不超过4个进程；后续分析串行运行。run_limited.sh限制BLAS/OpenMP等线程，seed=123，不嵌套并行。运行监测中三个R进程合计RSS观测峰值6.41 GiB (30秒采样，非瞬时峰值保证)。原始输入与历史结果只读。旧launch.sh等启动器保留作执行历史，以run_all.sh为准。

主展示：specimen内中位数、min.cells=10；5为低细胞量探索，20/30仅描述。所有其他分组单列，未根据结果选择阈值。患者为生物学独立单位；NA不是0；CellChat内置P不是患者比较P。未运行NicheNet、cell2location或空间共定位。

断点：每job unfiltered.rds先原子保存，四阈值完成后COMPLETE.csv。ERROR.csv保留原因，修复后复用概率恢复导出。完全相同barcode-state分组可复用本轮新推断，另保留reuse来源。不同输入/参数须新目录；不得复用已过期缓存。脚本中所有路径为本项目实际路径，迁移时统一修改common.R的BASE与输入路径。

LR bubble和pathway图优先按阳性队列数与效应量选择，正向与反向各最多10条，未按P值筛图。全部零/阴性结果保留于统计表。不同通路类别按实际数据库注释处理，Non-protein类别中的基因表达代表数据库相关分子/酶，不应都解释为分泌蛋白配体。

图注与统计定义见综合报告；数据/模型状态不足时明确NA，不以降低门槛制造显著性。完整统计比正文展示集更广，阴性和反向结果均保留。

上传状态：GitHub连接器写入返回403，尚未上传。轻量结果包位于delivery/llps_C5_highlow_light_results.zip，范围与校验清单见delivery/UPLOAD_STATUS.md及upload_manifest.csv。
