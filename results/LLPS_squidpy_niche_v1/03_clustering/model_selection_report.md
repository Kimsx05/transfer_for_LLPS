# Composition-only model selection

Balanced pilot: up to 200 spots per section, seed 123. All 12 radius/annotation/transformation combinations evaluated at k=30 and resolutions 0.2–1.0. k=15/50 refined for each annotation level at R2 using its highest-silhouette transform. Selection score = silhouette -0.10 if any niche <1% of pilot -0.03 per niche above12. 3–20 niches allowed; target5–10 not enforced. This prespecified numerical rule does not inspect signatures, clinical variables, coordinates or communications. Spatial diagnostics and reproducibility must still inform reliability ratings.

Unified training labels mapped via 15-neighbor distance-weighted classifier in the same composition feature space; maximum class probability saved per spot. Full pooled Leiden validates both annotation winners. UMAP is only visualization and is fitted on composition KNN, never used as clustering input. Igraph Leiden uses 2 iterations consistently for tractability; seed sensitivity is saved.

No Z-score or batch correction. All candidate pilot labels saved; every chosen model labels all spots. K-means K4–12 is an independent comparison in recommended composition feature space.
