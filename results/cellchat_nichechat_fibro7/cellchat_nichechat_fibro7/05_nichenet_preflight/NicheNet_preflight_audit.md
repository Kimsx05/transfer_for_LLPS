# Official NicheNet pre-flight audit

- Installed nichenetr: TRUE; version: 2.2.1.1; fixed commit: 66f90d5eeafef280b2b2f339b3fd70ffec1781dd.
- Source: [official saeyslab/nichenetr](https://github.com/saeyslab/nichenetr). Human v2 model: [Zenodo 7074291](https://zenodo.org/records/7074291), as linked by the official ligand-activity vignette.
- All three resources are locally available, loaded and checksum verified: ligand-target matrix, human ligand-receptor network, weighted networks. Exact paths, official MD5, local SHA256 and byte sizes are in verified_resource_manifest.tsv.
- Existing local proxy failed CERN TLS; direct verified TLS to Zenodo succeeded. No model replacement or disabled certificate verification was used.
- Ligand-target matrix: 33354 target genes x 1226 ligands. LR unique pairs: 4986. Weighted components: lr_sig, gr.
- RNA symbols are read directly from a frozen sample RNA/data matrix, with the same full gene universe as the authoritative atlas. Matching is exact and case-sensitive; unmatched symbols remain explicit. Expression identifiers and values are not renamed/recomputed.

      resource_component n_resource_symbols n_exact_RNA_matches n_unmatched
                  <char>              <int>               <int>       <int>
1: ligand_target_targets              33354               23248       10106
2: ligand_target_ligands               1226                1181          45
3:            LR_ligands               1226                1181          45
4:          LR_receptors               1067                1043          24
5: weighted_lr_sig_nodes              24825               19237        5588
6:     weighted_gr_nodes              30066               22780        7286
   resource_match_fraction n_RNA_genes RNA_match_fraction
                     <num>       <int>              <num>
1:               0.6970079       71683         0.32431678
2:               0.9632953       71683         0.01647531
3:               0.9632953       71683         0.01647531
4:               0.9775070       71683         0.01455017
5:               0.7749043       71683         0.26836209
6:               0.7576665       71683         0.31778804
                        matching_rule
                               <char>
1: exact case-sensitive; no remapping
2: exact case-sensitive; no remapping
3: exact case-sensitive; no remapping
4: exact case-sensitive; no remapping
5: exact case-sensitive; no remapping
6: exact case-sensitive; no remapping

## Scope and receiver-program gate

- Sender: frozen Fibro_C0–C6. Candidate receivers: Epithelial, Myeloid, Endothelial, T_NK; B, Plasma, Mast remain available when supported by CellChat. No pooled network or generic 14x14 NicheNet network.
- NicheNet links a potential ligand to a receiver target program through a prior model; ligand-target regulatory potential is distinct from inferred ligand activity. Expression gates alone do not establish receiver activity or causality.
- Per-sample sender ligand / receiver receptor expression is audited at group >=30 cells. Detection means RNA/data>0; a 10% expressed flag is included as an explicit diagnostic convention from the official workflow, alongside continuous fractions and means. It does not define geneset_oi.
- geneset_oi is UNDEFINED. No arbitrary marker list, top-DEG list, stage comparison, or CellChat pathway genes are silently used as the receiver program. Ligand-activity inference is not executed this round.
- NicheNet functions accept per-sample expressed background/potential ligands, but ordinary nichenetr does not itself implement a sample-replicated cohort test. The template records sample identity and leaves receiver program/sample evidence design explicit.
- Prioritization after CellChat will require: pair support; sender ligand expression; receiver receptor expression; receiver-program provenance; activity metrics; interpretable ligand-target links. CellChat receptor complexes remain recorded with subunits and are not naively equated to a single NicheNet receptor symbol.
- READY means environment/resources/expression audits/template are prepared; it does not mean ligand activity has run. PASS is reserved for a later completed NicheNet analysis with a frozen receiver program.
