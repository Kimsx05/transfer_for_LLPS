# Exact CellChat method

CellChat 2.1.2; release v2.1.2; commit a0d3b2d231d46c8787177fffeac908270c253747. R 4.4.3.
- Human database restricted to Secreted Signaling, ECM-Receptor and Cell-Cell Contact; no category is omitted.
- One object per orig.ident biological sample. All included cells of present groups enter inference; group levels absent from that sample are dropped inside CellChat and restored to the 14x14 output afterward.
- createCellChat(expression, meta, group.by="cellchat_group", datatype="RNA"); subsetData; identifyOverExpressedGenes(do.fast=FALSE); identifyOverExpressedInteractions; computeCommunProb(type="truncatedMean", trim=0.1, raw.use=TRUE, population.size=FALSE, nboot=100, seed.use=123); filterCommunication(min.cells=30); computeCommunProbPathway(thresh=0.05); aggregateNet(thresh=0.05, remove.isolate=FALSE).
- do.fast=FALSE uses the package Wilcoxon implementation without an optional presto dependency. Other gene-selection defaults remain unchanged (thresh.pc=0, thresh.fc=0, thresh.p=0.05).
- aggregateNet and LR retained weights use bootstrap p<0.05, matching installed API. Native pathway output uses p<=0.05 in computeCommunProbPathway; this API boundary difference is retained and disclosed. Raw candidate LR probabilities/p-values are also saved.
- No projection, downsampling, pooling, stage comparison, score input, or upstream mutation. Global and per-sample seed reset to 123; sequential future plan within each sample.
- Full 14x14 outputs: source/target <30 -> NA; successful eligible undetected -> 0; failed sample -> NA with FAILED status. Failed samples never supply zeros.
- The only adaptive failure handling is catch-and-record per sample; no parameter changes or invented empty objects.

Executable source: scripts/03_cellchat_per_sample.R
