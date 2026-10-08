# NicheNet expression audit

- All 59 biological samples are read independently from the frozen RNA/data inputs; no pooling or resampling.
- Sender groups: Fibro_C0–C6; receiver groups: all seven major populations. CellChat will determine follow-up priority.
- Ligand/receptor universe is the verified official human NicheNet v2 LR resource. Missing RNA feature -> NA, not assumed zero. Absent/<30-cell groups -> NA. Evaluable but unexpressed genes -> true zero fraction and mean.
- Detection means RNA/data>0; mean is arithmetic mean of existing log-normalized RNA/data values. 10% expression is an explicit diagnostic flag and retains sample prevalence, not a hard biological significance cutoff.
- Cohort expression medians include eligible observed zeros. Sender/receptor summary n_eligible and n_evaluable are distinct for missing genes.
- receiver_expressed_background_candidates.tsv.gz contains per-sample expressed target-universe genes for a future background. It is NOT geneset_oi, and does not freeze the eventual program/background replication design.
- No ligand activity, DE test, receiver program selection, score input, or annotation/reduction mutation was performed.
