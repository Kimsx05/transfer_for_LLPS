# QC sensitivity audit summary

- Audit ID: `20260914_185530`
- Immutable source checkpoint: `/home/jinshengxi/LLPS/1.2 基本处理/scRNA_standard_pipeline_v1_20260914_151017/checkpoints/01_merged_full_20260914_153425.qs2`
- Cells audited: 410,466
- No subset, filtered object, checkpoint modification, DoubletFinder, PCA, Harmony, clustering, or UMAP was performed.

## Direct answers

1. **Is current QC obviously too strict?** 当前 Rule A 对 upper tail 有中度敏感性。 Rule A retention = 87.55%; Rule B adds 11,765 cells (2.87% of input).
2. **Lower versus upper tails:** summed flag counts (overlap possible) are lower = 14,136 and upper = 23,165. See tail tables for non-overlapping interpretation.
3. **Rule B versus A:** +11,765 retained cells.
4. **Rule C versus A:** +14,659 retained cells (3.57% of input).
5. **Most upper-tail-sensitive groups:**
   - HRA003620 / HRA003620_NC-30: B+A gain 1021
   - InhouseData / T4398343: B+A gain 903
   - InhouseData / T4471438: B+A gain 745
   - GSE222315 / GSE222315_p3_BCa: B+A gain 699
   - GSE222315 / GSE222315_p4_NAT: B+A gain 627
   - InhouseData / T4463209b: B+A gain 498
   - HRA003620 / HRA003620_NC-11: B+A gain 441
   - CNP0000460 / CNP0000460_P03T: B+A gain 433
   - GSE135337 / GSM5329919_BCN: B+A gain 422
   - GSE222315 / GSE222315_p5_NAT: B+A gain 407
6. **Why InhouseData MT is zero/no-failure:** 是：InhouseData counts layers 中没有任何 ^MT- feature，MT=0 主要是 feature naming/mapping/feature omission 问题，不能解释为真实线粒体转录为零。
7. **Does HB 5% have little impact?** Yes. It flags 189 cells (0.046%).
8. **Recommendation for future consideration:** Rule B. This audit does not apply the recommendation. Rule C is not preferred by default because it removes all upper-tail safeguards before downstream doublet assessment.

## Important implementation note

The formal v1 QC code explicitly used `stats::mad(..., constant = 1)`. This audit reproduces that exact implementation. It is stricter than R's default consistency-scaled MAD (`constant = 1.4826`) and should be considered when interpreting sensitivity.

Audit elapsed seconds: 34.68
