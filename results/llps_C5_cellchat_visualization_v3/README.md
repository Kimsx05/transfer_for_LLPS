# LLPS C5 CellChat visualization v3

正式口径：完整C5，tumor-only，min.cells=20；10/30为共享上游推断的覆盖敏感性。没有重跑CellChat、NicheNet、空间分析或high/low分组。

- [中文综合报告](COMPREHENSIVE_REPORT.md)
- [图形画廊](FIGURE_GALLERY.html)
- [图形索引](FIGURE_INDEX.csv)
- [旧图逐项映射](legacy_figure_mapping.csv)
- [最终QC](FINAL_QC_REPORT.md)
- [交付状态](delivery/UPLOAD_STATUS.md)

01–09目录为完整C5及补充证据；10目录为旧图原样归档。每图都有_data.csv，统计单位和限制见索引。完整PDF为基本文件名，有图例另存_panel_only.pdf及_legend_only.pdf。不同边患者数不同，灰色NA不能当0，IQR不是CI，prevalence不是显著性。

复现：scripts/run_all.sh，单线程串行；seed123。脚本只写本v3目录，不加载大型Seurat，不复制原始表达矩阵。原始完整网络与候选全量表留在v2，本轮保留绘图和审查所需数据。09空间图只是读取原有统计表。
