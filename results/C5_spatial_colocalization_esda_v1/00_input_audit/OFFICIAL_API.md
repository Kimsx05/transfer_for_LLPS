# Official APIs checked 2026-11-18

- https://pysal.org/esda/stable/generated/esda.Moran_BV.html
- https://pysal.org/esda/stable/generated/esda.Moran_Local_BV.html
- https://github.com/pysal/esda/blob/main/esda/moran.py
- https://r-spatial.github.io/spdep/reference/lee.mc.html

Installed source and version records are authoritative for this run (esda 2.7 compatible with Python 3.10; current development docs may include an `alternative` argument absent in the installed release). Global RNG: numpy.random.seed(123) immediately before each fit. Local conditional RNG: seed=123 inside the function, n_jobs=1. geoda_quads=False maps 1=HH,2=LH,3=LL,4=HL; each local fit is checked against standardized X and W-lagged standardized Y signs. No custom Moran or Lee statistic is implemented.

For esda 2.7 minimum-tail pseudo-P, preserve p_sim and use min(1,2*p_sim) as a conservative two-tail conversion before BH. Local BH family is all valid spots across four targets within one section. Native-p BH is supplementary globally. Lee uses alternative='two.sided', nsim=999, set.seed(123), same row-standardized weights, no diagonal/self-neighbors. Lee's package null permutes paired x/y rows jointly, different from Moran_BV's Y-only permutation and Local Moran conditional randomization; significance is not interchangeable.

No conda executable discovered in standard user/opt installation roots. System venv lacked ensurepip. Python packages were installed into isolated /tmp/c5_esda_deps; R packages into this result's .Rlib. Existing environments were not modified. Initial sandbox network failures and Python-version resolution failures are retained in installation logs.

Installed spdep::lee.mc rank-based two-sided p.value can be exactly 0 at an extreme observed rank; this is retained as the package output, never interpreted as zero true probability. Its 999-permutation resolution remains limited. All inference here is exploratory spatial randomization, not independent-patient sampling inference.
