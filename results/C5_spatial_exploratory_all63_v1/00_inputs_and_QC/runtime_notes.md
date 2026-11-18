# Runtime notes

An early long-running remaining-section process exited with signal 15 (status 143); its cause was not established. No incomplete section Parquet was accepted. A subsequent R invocation encountered a parse error while the scoring script was being updated; it was restarted after the file was frozen. Completed sections were reused.

Set options(mc.cores=1L) before package loading to prevent BiocParallel from attempting to construct its default socket backend in the restricted environment. Explicit SerialParam remains in scoring calls. This does not change scoring parameters. Remaining work uses resumable batches, atomic Parquet publication and per-section logs. All final outputs undergo separate integrity checks.

HE registration is limited to source identity, dimensions and coordinate bounds. Auxiliary overlays are provisional; primary maps are coordinate-only.

The initial asyncio batch scheduler did not advance after three workers had written their final section outputs. It was stopped, and the remaining disjoint batches were launched directly. The resumable scheduler was replaced with explicit subprocess polling; finished score files were not recalculated. This scheduling issue did not alter scoring algorithms.

Figure inspection identified PostScript text-measurement warnings and default sans-serif annotations in some composite figures. A Cairo null device, global Arial theme and explicit Arial legend labels were configured; composite maps, atlases and summaries were regenerated. Main map PDFs contain vector points/color guides; HE PDFs intentionally contain the source raster image. Final validation checks PNG 300 dpi, parseable PDFs and embedded Arial fonts.

Some distance-profile points/rings are absent because a graph component has no epithelial-region seed, or no spots in a ring. These remain NA and are omitted from the relevant plotted geometry, rather than imputed as zero. Corresponding data tables retain those missing values and denominators.
