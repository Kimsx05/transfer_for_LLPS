# LLPS Squidpy niche clustering v1

Start with COMPREHENSIVE_REPORT.md, 03_clustering/recommended_models.tsv and 04_niche_characterization/spot_metadata.parquet. Browse GALLERY.html locally; the lightweight package includes DELIVERY_GALLERY.html with its selected figures. All63-section formal maps remain available locally. Public GitHub upload is pending explicit public-payload approval; see delivery/UPLOAD_STATUS.md.

Reproduce with scripts/reproduce.sh; it invokes audit, graph/composition, clustering, native control, inference, spatial statistics, clinical source audit, method sensitivity, figures and reporting in dependency order. The supplemental script also invokes12_convergence_diagnostic.py. Python3.10 dependencies are frozen in00_input_audit/requirements.txt and installed in local.pydeps; original read-only inputs listed in input_manifest.tsv and the plotting skill must remain accessible. Existing complete R figure outputs are reused; remove only this run’s figure files if intentionally regenerating modified graphics. Input paths are frozen in input_manifest.tsv and config.json. No original data changed.

Resources:4 CPU threads,100GB limit; seed123. Main scripts use per-section spatial graphs and continuous posterior-mean composition.
