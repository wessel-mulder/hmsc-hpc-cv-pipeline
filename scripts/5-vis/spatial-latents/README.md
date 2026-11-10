# Spatial latent-factor maps

Run from the repository root:

```sh
Rscript scripts/5-vis/spatial-latents/plot-spatial-latents.R
```

The script reads `fitted_model$posteriors` from the three
`2026-11-02_06-55-21` atlas models under `data-new/2-hmsc/derived/models`.
It uses every retained draw (four chains, 250 draws each), all five spatial
latent factors, and all fitted cells: 1,747 in Atlas 1, 1,821 in Atlas 2,
and 1,875 in Atlas 3. No species or VP filtering is applied.

Six PNG composites are written here: posterior means and posterior SDs separately
for each atlas. Each composite contains five factor maps with Bornholm insets.
Mean maps share one zero-centred colour scale across all atlases/factors; SD
maps share one sequential scale. Full value ranges are included without clipping.

Eta rows are matched through the fitted model's `Pi`/`dfPi` site indexing.
Stored spatial coordinates are transformed from EPSG:23032 (as set in
`S0_Data_Definitions.R`) to longitude/latitude. Cell IDs and coordinate locations
are checked against the local atlas grid. Selected site/factor summaries are
checked against direct extraction from the original chains.

These maps show **Eta**, not species-specific spatial contributions
(`Eta %*% Lambda`). They retain the saved within-model alignment from posterior
import. Factor signs, orientation, and ordering are arbitrary: factor 1 in one
atlas need not correspond to factor 1 in another. No cross-atlas alignment or
factor differencing is performed. Posterior SD describes uncertainty conditional
on the saved alignment and is not an assessment of convergence. Pooled means
and SDs can conceal chain differences; these are descriptive maps, not evidence
that all latent quantities or VP estimates have converged.

Dependencies: sf, ggplot2, cowplot. Models and source data are read only.

## Climate correlations and report

Render `spatial-latents.Rmd` for the maps, climate–factor correlations, tables,
and diagnostic discussion in a self-contained HTML report:

```sh
Rscript -e 'rmarkdown::render("scripts/5-vis/spatial-latents/spatial-latents.Rmd")'
```

The separate `correlate-climate-spatial-latents.R` script computes Pearson
correlations between temperature/precipitation and each of the five factors.
It reports correlations with posterior mean maps and correlations for each of
the 1,000 draws, on both the fitted cells and the 1,697-cell intersection.
Three additional PNGs show fitted-cell intervals, shared-cell intervals, and
chain-specific mean correlations. Chain identities are retained for
rank-normalised R-hat, bulk/tail ESS, and MCSE diagnostics (`posterior` package).
No ordinary independent-cell correlation tests or p-values are calculated.

Many derived correlations have poor chain agreement. The report flags this:
pooled-draw intervals are descriptive until mixing/alignment concerns are
resolved. Factor correlations depend on the saved basis, and their signs and
factor numbers cannot be compared directly between independently fitted atlases.
Knit recomputes correlations and embeds existing map PNGs, regenerating maps
only if any are absent. Run the map script explicitly to refresh all maps.

## All predictors against the combined spatial random effect

`correlate-environment-spatial-effect.R` correlates every original environmental
predictor with each species' combined spatial effect (`Eta %*% Lambda`), using
the same 1,697 cells and all 1,000 posterior draws. This effect is invariant to
compensating factor/loadings rotations and is inferred rather than measured.

`environment-spatial-effect-correlations.png` displays the 155 species' median
correlations for all ten predictors and three atlases. The companion
`environment-spatial-effect-correlation-strength.png` summarises median absolute
correlation across species within each draw and its pooled-draw interval.
Both plots and diagnostic tables are included in `spatial-latents.Rmd`/HTML.
Cross-cell overlap is not a direct measure of coefficient bias or the
posterior compensation between fixed and spatial contributions.

## Climate spatial residualisation

`climate-spatial-residualisation.Rmd` and its matching R script inspect the
climate variation left after subtracting a fitted smooth geographical trend.
They use the same 1,697 cells in all three atlases and the original temperature
(°C) and precipitation (mm/season) values from the November preprocessing input.

```sh
Rscript -e 'rmarkdown::render("scripts/5-vis/spatial-latents/climate-spatial-residualisation.Rmd")'
```

Separate REML thin-plate GAMs are fitted per atlas and climate variable at
k = 50, 100, and 200. Five PNGs show original/trend/residual maps (main k = 100),
physical and variance-normalised semivariograms with pair counts, and
smoothing sensitivity. Tables report retained SD/variance, ranges, descriptive
Moran's I on a common symmetrised 4-neighbour graph, EDF, and basis diagnostics.
Objects remain in memory; source inputs and bird models are not modified.

The atlas-specific intercept is included in the subtracted trend, so residuals
remove atlas means as well as broad spatial gradients. Autocorrelation can
remain. This first-stage exploration does not refit Hmsc or establish that
spatial confounding is resolved, and it does not propagate smoothing uncertainty.
