# C0 去卷积优化与外部临床分析

已完成误差诊断、重拟合基线/优化参考/两级 MuSiC 比较，以及 signature 替代路线的四队列分期与生存分析。**最终临床结果是 C0-like signature score 的探索性关联，不是 C0 abundance。**

优化参考在原 20 人开发评估中的 Spearman 提升至 0.908，但比例偏差及 C6 跨 dataset 覆盖不足尚未解决。20 基因 signature 的留出 C0-vs-other-stroma AUC=0.919；四队列共输出 4 张分期图、4 张 KM+95%CI+risk table 和 1 张描述性 Cox forest。没有一致的有利预后证据，比例风险检验及敏感性结果已报告。

- 完整报告：[summary/final_report.md](summary/final_report.md)
- 正式临床图及其他独立图：[summary/figures.md](summary/figures.md)
- 全部绝对路径及 SHA256：[summary/file_index.tsv](summary/file_index.tsv)
- 方法比较：[02_model_comparison/model_comparison.tsv](02_model_comparison/model_comparison.tsv)
- 分期、KM、Cox、PH 汇总：summary/signature_*_summary.tsv
- 数值校验：[summary/numerical_checks.txt](summary/numerical_checks.txt)

运行目录为 `/home/jinshengxi/bladder_metabolism`，代码按 scripts/01_diagnostics.R 至 07_sensitivity_and_checks.R 顺序组织。输入均复用已有数据；原 20 人不进入新模型或 signature 的训练/校准，开发评估不冒称新盲法验证；seed=123。旧图、注释、LLPS AUCell、cluster 配色及已下载队列未修改。所有新图独立保存为单页矢量 PDF + 350 dpi PNG。
