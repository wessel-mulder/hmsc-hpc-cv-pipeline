# 3 GWR Fitting

This folder is intended for geographically weighted regression scripts used to analyse spatial variation in richness, turnover, or environmental-driver summaries.

It is currently empty. At the moment, the active GWR figure workflow lives inside `../4-data-visualisation/fig10-gwr-richness-turnover-materials/`, where the figure script can either reuse cached GWR fits or refit them when requested.

During cleanup, move reusable GWR fitting code here if it should be part of the formal analysis pipeline rather than a figure-specific helper. Cached GWR model objects and derived summaries should go under `intermediate-outputs/gwr-output/`, not inside the scripts folder.
