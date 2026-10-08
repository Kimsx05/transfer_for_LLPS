# Fibro7 CellChat + NicheNet analysis

The active methods are CellChat and official NicheNet/nichenetr. The directory name is retained from the initial request so frozen-input paths and existing links remain valid. The former NicheChat name audit is archived under `00_audit/superseded_name_correction_20261008`; it is not a computation requirement. NICHECHAT_STATUS is NOT_APPLICABLE.

CellChat uses one independent object per biological sample, a fixed 14-group output frame, min.cells=30, truncatedMean/trim=0.1, population.size=FALSE, nboot=100 and seed=123. Primary cohort strength is the median sample net$weight including eligible zeros, displayed only when n_eligible>=3 and prevalence>=0.30.

NicheNet scope in this round: official package/resource/gene-symbol audit, sample-aware sender ligand and receiver receptor expression audit, and a receiver-program template. No geneset_oi is defined, and no ligand activity or ligand-target mechanism analysis is executed. READY denotes a prepared environment and template; it does not claim completed activity inference.

See `FINAL_REPORT.md`, `00_audit/QC_REPORT.md`, `05_nichenet_preflight/NicheNet_preflight_audit.md` and `07_nichenet_templates/README.md` for results and scope. The final artifact manifest provides absolute paths to existing tables and figures.
