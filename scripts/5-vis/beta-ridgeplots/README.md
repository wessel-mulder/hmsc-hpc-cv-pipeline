# Atlas-scaled Beta ridgeplots

Run from the project root:

```sh
Rscript scripts/5-vis/beta-ridgeplots/beta-ridgeplots-environment-only.R
```

The script uses the November 6 environment-only fits: ten environmental
predictors, 155 matched species, and all four chains × 250 retained draws.
Every variable has three ridges, labelled 1970s, 1990s, and 2010s on the left,
with increasing opacity. Variable colours match the Sankey charts. There are
no ridge outlines.

- `beta-ridgeplots-environment-only-species-means.png`: the distribution across
  species' posterior mean Betas, with equal weight per species. This is a
  descriptive distribution of species estimates, not posterior uncertainty.
- `beta-ridgeplots-environment-only-posterior-mixture.png`: the equal-weight
  mixture of every species' full posterior (155 × 1,000 values per ridge).
  It combines between-species differences with coefficient uncertainty and is
  not the posterior distribution of a community mean.

Betas are deliberately retained on the fitted, atlas-specific environmental
scales: a one-SD predictor increase **within that atlas's fitted cells**. They
are changes in the probit linear predictor, not percentage-point changes in
occurrence probability. Betas are not divided by their own posterior SD or
converted to original predictor units. The intercept is excluded.

The script verifies that the model X columns have mean zero and SD one and
reconstructs them from raw predictors using each atlas's full training sample.
The within-atlas physical SDs are retained in `predictor_scale_audit` in memory.
Consequently, temporal differences can reflect changed predictor variability
as well as changed raw-unit species responses. These plots do not isolate
response change or provide a paired species-level test of it.

All atlas ridges use a common horizontal display axis from -1 to +1 in both figures. A common KDE
bandwidth is selected from the pooled values for each variable, then used for
all three atlases. The same density-to-height multiplier is used within each
variable, with heights normalised across variables for readability. Densities
are calculated from all values on a grid extending at least four bandwidths
beyond their extrema, then the display is cropped to [-1, +1]. Some posterior
tails are outside the visible range; no values are discarded from the density
calculation. The plots remain 2D to avoid occlusion and perspective distortion.
R-hat summaries are printed and retained in
`beta_diagnostic_summary`; no coefficients are filtered based on significance
or diagnostics. Poorly mixed estimates remain descriptive.

All outputs are PNGs. Original models, VP plots, and reports are unchanged.
