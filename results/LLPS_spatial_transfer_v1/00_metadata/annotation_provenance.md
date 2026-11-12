# 注释审查与来源

本包读取正式 `cell2location_full_run/v2_30type`，逐切片核对冻结清单、COMPLETED_AND_VALID、model_config、最终 validation gate、原始 CSV 和 posterior H5AD。未训练、重新抽样、归一化、缩放、二值化、修改注释或任何原始文件。

实际为 63 张、265,002 个 hard-QC retained spots、7 个保留 dataset 标签。预期 262,520 与实际相差 2,482。正式冻结清单、原始输入、mean/q05、final validation gate 及此前空间生态清单一致支持 265,002；未找到 262,520 对应的正式筛选依据，不能将其归因于某个未经证实的 QC 步骤。没有为匹配预期删减 spot。

最终 reference 是 300,033 细胞、17,987 基因、58 batches、30 亚型；约 35.6 万不等于本次冻结模型实际使用细胞数。本次还只读检查了参考输入 X shape=300033 × 17987，并确认参考输入、模型和签名文件均存在。参考来源：
- `/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/reference_training_v2/11_final_spatial_reference/reference_v2_final_30type_regression_input.h5ad`
- `/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/reference_training_v2/11_final_spatial_reference/30type_regression_model/model/model.pt`
- `/home/data/t070721/codex_workspace/Bladder Metabolism/subtype/v2/reference_training_v2/11_final_spatial_reference/30type_regression_model/learned_signatures_30type.csv.gz`

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

原始脚本路径：`/home/data/t070721/codex_workspace/Bladder Metabolism/cell2location_full_run/v2_30type/code/01_run_one_section.py`；counts 准备脚本：`/home/data/t070721/codex_workspace/Bladder Metabolism/cell2location_preflight_v1/code/preflight/04_prepare_visium_inputs.R`。其副本随包保存以便审阅，禁止把训练/恢复脚本当作本次执行入口。
