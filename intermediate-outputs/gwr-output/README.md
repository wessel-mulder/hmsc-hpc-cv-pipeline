# GWR Outputs

This folder stores cached GWR products used by `scripts/4-data-visualisation/fig10-gwr-richness-turnover-materials/fig10-gwr-richness-turnover.R`.

Current products live in `fig10-gwr-richness-turnover/gwr/`:

- `models/2026-03-13-gwr-fits-richness-beta-100km.rds`;
- `models/2026-03-13-gwr-analysis-inputs-100km.rds`;
- `2026-03-13-gwr-100km-summary.csv`;
- `2026-03-13-gwr-local-coefficients-100km.csv`;
- `2026-03-13-gwr-dominant-drivers-100km.csv`.

Use the figure script without `--force-refit` to rebuild the manuscript PNG from these caches. `--force-refit` requires local `HmscOutputs/` fitted models.
