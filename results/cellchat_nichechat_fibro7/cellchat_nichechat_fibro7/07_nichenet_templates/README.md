# NicheNet receiver-program gate

This round prepares official NicheNet resources, expression audits and a template. It does not run ligand activity or define geneset_oi.

`NicheNet_analysis_template.R` only defines a function. The function requires an explicit receiver gene set, sample-specific expressed background, expressed candidate ligands, frozen program provenance, biological rationale, gene-selection rule and replication design. There is no default marker set, automatic top-DEG set, stage contrast or LLPS-based program.

`receiver_program_design_template.tsv` is an empty schema; no gene list exists. After CellChat, `receiver_program_design_after_CellChat.tsv` supplies supported pairs while keeping program fields undefined.

The future decision must state what receiver response is being explained, why the gene program measures that response, and how its evidence repeats across biological samples. A curated externally documented response program or a clearly justified sample-aware receiver transcriptional program can be considered at that time. Do not infer direction of causality from a sample-level association alone.

Candidate ligand activity is continuous (`pearson`, `auroc`, `aupr`, `aupr_corrected`); high activity and a reporting limit will require an explicit analysis definition. Official model regulatory-potential weights are distinct from activity scores. Retain ligand–receptor and ligand–target outputs separately, with CellChat support and expression evidence. Never classify an unrun ligand activity analysis as negative support.

Receiver background candidates from the expression audit are expressed genes in the target-matrix universe, not geneset_oi. CellChat ligand/receptor complex identifiers require an explicit subunit-aware mapping; an exact single-gene join cannot establish a whole complex interaction.

The original result-directory name is retained to preserve links to frozen inputs. The active analysis is CellChat + NicheNet; NICHECHAT_STATUS is NOT_APPLICABLE.
