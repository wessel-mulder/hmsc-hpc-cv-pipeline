# Fig10 GWR Richness/Turnover Materials

This folder bundles the materials needed to reproduce the GWR figure with the panels "A Predicted richness" and "B All-site beta turnover".

Run from the repository root:

```sh
Rscript scripts/4-data-visualisation/fig10-gwr-richness-turnover-materials/fig10-gwr-richness-turnover.R
```

The script uses cached GWR objects from `intermediate-outputs/gwr-output/fig10-gwr-richness-turnover/gwr/models/` by default.

It writes regenerated CSV summaries to `intermediate-outputs/gwr-output/fig10-gwr-richness-turnover/gwr/` and the final PNG to `manuscript-outputs/figures/main/`. Use `--force-refit` only if you want to refit the GWR models from local `HmscOutputs`.

Original recovered source: git commit `9d99ffa2f7f7fa27e8d1c60029ee6de455c1fabf`.
