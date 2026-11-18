# Spatial maps

Every section has a 10-panel overview (PDF/300 dpi PNG), a panel-only PDF, a legend-only PDF and an auxiliary HE view with corresponding exports. Primary maps use the same cropped observed-coordinate window across panels; HE views use the complete low-resolution image window. Coordinates retain upper origin, without rotation or mirroring.

Raw signature panels share [0,1]. The A26 percentile panel is within-section only. Abundance, MMP11 log1p(CP10k) and log10(1+UMI) share global 99th-percentile upper limits; exact limits and clipping counts are in `../tables/figure_color_scales.tsv`. Common displayed limits do not calibrate cross-dataset biological scores.

HE registration has identity, dimension and bounds checks only, not independent histological validation. No image appearance was used to modify signatures or select sections. Compact plotting tables retain the data behind each panel. See `../GALLERY.html` and the all-section atlas.
