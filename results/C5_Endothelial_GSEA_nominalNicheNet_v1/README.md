# C5–Endothelial GSEA / nominal NicheNet v1

详见 [COMPREHENSIVE_REPORT.md](COMPREHENSIVE_REPORT.md)。固定 21 肿瘤患者、7,070 背景、246 配体；主要名义 Set B，A/C 敏感性。主模型无 FDR 显著 DEG，EMT 为 REFERENCE_ONLY。

从项目根目录执行 `bash results/C5_Endothelial_GSEA_nominalNicheNet_v1/scripts/run_all.sh`。CPU affinity 最多4核，各计算进程单线程，seed=123。历史输入只读，按阶段/程序 checkpoint 续跑；不复制大矩阵。复现依赖本机既有 R 库及 Python 依赖路径，见 run_R.sh/run_python.sh 和 sessionInfo。输入哈希见 00_input_audit。

NicheNet 活性分数是先验网络的程序匹配，名义候选不称 significant DEG；队列方向与LOPO不构成外部验证。空间新增仅3轴表达审查，没有新的通讯或空间显著性推断。
