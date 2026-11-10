Run from the repository root:

```sh
Rscript scripts/5-vis/env-viz/environmental-panel.R
```

`environmental-panel.png` shows all ten environmental predictors across Atlas
1 (1970s), Atlas 2 (1990s), and Atlas 3 (2010s). Each predictor also has its own
larger three-period PNG. All figures use the Denmark circular-point layout with
Bornholm insets and blue–yellow–red scales.

Reads `X` and `Design` from the portable S0 training bundle in `data-shared`.
Uses only the 1,697 training cells shared by all three atlases, without changing
the source bundle. Shows original predictor units, not atlas-specific z-scores.
Each predictor's colour limits span its observed range across all three periods
on this common sample; different predictors have different scales.

Land cover is percentage of the whole cell, including water in the denominator;
it is not rescaled to percentage of land. Habitat heterogeneity is the mean
annual Simpson diversity over all 15 HILDA classes, including ocean/freshwater.
The forest variables preserve S0's deciduous (43+44), evergreen (42+41), and
mixed/other (40+45) groups; cropland is class 22 only.

The bundle already applies S0's coverage, water and complete-predictor filters.
Atlas 1–2 coverage is Good/Average; Atlas 3 coverage is assumed Good. Atlas labels
denote periods, not single-year measurements. See `scripts/2-hmsc/S0_Data_Definitions.R`
and `data-new/2-hmsc/derived/forest_split_description.txt` for input definitions.

For a compact comparison of distributions across the same cells, run:

```sh
Rscript scripts/5-vis/env-viz/environmental-distributions.R
```

This creates `environmental-distributions.png`: ten panels showing period medians,
25th–75th percentile boxes, and full minimum–maximum whiskers (not the standard
1.5-IQR whiskers). Vertical scales differ between predictors. These describe
variation among cells, not uncertainty in the median. Tied quartiles, especially
for forest covers with many zeros, can collapse a box to a line.
`environmental-distribution-summary.md` contains the plotted values and the
percentage of cells with zero values.

For all predictors on one compact, common standard-deviation axis, run:

```sh
Rscript scripts/5-vis/env-viz/environmental-standardised-distributions.R
```

Creates `environmental-standardised-distributions.png`, with three horizontal
boxes per variable and the same atlas colours. Each predictor is transformed as
`(value - pooled mean) / pooled SD`, pooling the 5,091 cell–period observations
(1,697 cells × 3 periods). The SD includes spatial and temporal variation;
it is not the SD of just the three period means. Original-unit pooled means
and SDs appear beside the variable labels and in
`environmental-standardisation-reference.md`. Boxes retain the middle 50%,
medians, and full minimum–maximum whiskers.

For the dominant modelled land-use category in each common cell, run:

```sh
Rscript scripts/5-vis/env-viz/dominant-land-use.R
```

Creates `dominant-land-use-atlas1-2-3.png` and `dominant-land-use-summary.md`.
Uses the highest original percentage among the seven modelled land-cover
predictors, preserving the forest groups. This is the largest modelled category,
not necessarily more than 50% of the cell or the largest of all raw land-cover
classes. Water and other unmodelled classes are excluded. Tied maxima are shown
separately (tolerance 1e-8 percentage points).
Land-use colours match `variable_colours` in
`scripts/5-vis/vp-distributions/vp-geom-sankey-bump.R`.
