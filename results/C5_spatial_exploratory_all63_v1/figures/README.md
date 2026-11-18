# Figure index

- `all63_percentile_atlas`: all sections ordered by dataset, audited tissue label and section ID. Colors are within-section percentiles, not absolute cross-section scores.
- `*_atlas`: one atlas per source group.
- `moran_by_dataset`, `high_score_connectivity`: spatial structure and high-score component size.
- `score_correlations`: method agreement and composition/depth diagnostics.
- `tissue_descriptive`: reliable Normal/Tumor labels only; no patient-level tissue inference.
- `epithelial_rings`, `matched_neighborhood_rings`: exploratory graph-distance profiles.
- `epithelial_matching_sensitivity`, `epithelial_distance_sensitivity`: attenuation after additional epithelial composition control.

Statistical figures have matching `_data.tsv` tables. Full plots use PDF and 300 dpi PNG; figures with legends additionally have panel-only and legend-only PDFs. The composite maps and HE views for all 63 sections are in `../02_spatial_maps/`.
