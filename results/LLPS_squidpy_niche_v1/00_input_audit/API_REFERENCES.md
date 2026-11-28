# API and implementation audit

Official documentation consulted:
- https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.spatial_neighbors.html
- https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.calculate_niche_neighborhood.html

Installed Python3.10 resolves Squidpy1.6.5, whose local source uses calculate_niche(flavor='neighborhood'). The installed source was inspected directly; newer API names are not assumed available. Source: .pydeps/squidpy/gr/_niche.py. For cumulative radius composition, pass a binary, center-inclusive graph at distance1. This avoids the native distance>1 behavior that combines path-power profiles and can count repeated walks.

Native1.6.5 neighborhood implementation does not forward its random_state argument into internal Scanpy calls. Explicit adapters set both actual Scanpy random_state values to123. These do not alter source packages and are in03_native.py.

Compatibility fixes: Scanpy1.11.5 cannot consume PyNNDescent0.7 csr_array output; pin pynndescent0.5.13. Recent UMAP's typed contiguous-array requirement fails on Scanpy neighbor arrays; pin umap-learn0.5.9.post2. Parameters and biological inputs unchanged. Logs preserve errors and successful installations (primary clustering log was restarted; compatibility rationale retained here).

Exact array adjacency validates/corrects Squidpy grid edges; missing spots are not interpolated. Graph QC records both extra and missing native edges per section.
