# 4 Data Visualisation

This folder contains scripts for generating publication and supplementary figures from HMSC summaries, spatial products, and derived intermediate outputs.

## Current Figure Scripts

- `fig1-vp-effect-sizes.R`: variance-partitioning effect-size figure.
- `fig2-trait-influence-gamma-atlas-support.R`: trait moderation / Gamma support tile figure.
- `fig3-cmw-spatial.R`: community-weighted mean spatial figure.
- `sfig-full-vs-cv-accuracy.R`: supplementary full-model vs cross-validated accuracy figure.
- `sfig2-prevalence-vs-cv-auc-tjur.R`: supplementary species prevalence vs cross-validated AUC/Tjur R2 figure.
- `sfig-raptor-species-full-tjur-prevalence-trends.R`: supplementary full-model Tjur R2 and prevalence trajectories for prevalence-increasing species, with diurnal raptors highlighted.
- `sfig-component-influence-sankey-bump.R`: supplementary sankey-bump figure aggregating variance partitioning into spatial, land-use, and climate components.
- `sfig3-diurnal-raptors-cv-accuracy.R`: supplementary diurnal-raptor vs other-species cross-validated accuracy figure.
- `sfig4-guild-cv-accuracy-temporal-change.R`: supplementary temporal-change comparison of cross-validated accuracy for diurnal raptors and four common increasing-prevalence guilds.
- `sfig4-tjur-r2-guild-prevalence-tests.R`: focused Tjur R2 temporal-change tests and prevalence-adjusted diagnostics for the diurnal-raptor comparison.
- `fig10-gwr-richness-turnover-materials/`: self-contained materials for the GWR richness and beta-turnover figure.

## Outputs

Publication-facing figure scripts write final PNGs to `manuscript-outputs/figures/`. Reusable caches and compact summaries live under `intermediate-outputs/`. Scripts that need fitted posteriors still document `HmscOutputs/` as a local-only requirement instead of tracking those bulky objects in GitHub.

All generated plots should be PNG files unless a journal explicitly requests another format.
