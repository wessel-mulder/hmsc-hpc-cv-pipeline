# Posterior changes in Tjur R²-scaled VP

Completed 2026-11-11 11:29:55.765483. Local runtime: 76.6 seconds (1.28 minutes).

## What was calculated

Model family: spatial (2026-11-02_06-55-21).

All 155 species, ten predictors, four chains and every saved draw are retained. Each atlas uses its own fitted cells; there is no common-cell restriction or model refitting.

For each draw: compute Hmsc environmental VP; compute expected occurrence probabilities from that same draw; calculate species Tjur R²; multiply VP by Tjur R². Values are percentages and changes are percentage points. Negative draw-level fit values are retained.

Predictions include the fitted site effect. Environmental shares include the site contribution in the partition denominator. Within each draw, the environmental and site shares sum to the full Tjur R² after scaling. The legacy spatial VP export averages the fixed share and its allocation across variables separately before multiplying them.

Cells by atlas: 1747, 1821, 1875.
Negative Tjur draw/species entries by atlas: 0, 0, 6.

Independently combine whole atlas posterior draws using 10000 balanced draw triplets, seed 369. Reuse triplets across all contrasts; retain dependence among species and variables. These triplets do not create new observations or increase the original MCMC effective sample size.

## Community comparisons

Differences are later minus earlier. CI95 is a marginal credible interval. Joint95 is a joint posterior envelope across all 30 community contrasts, not a Holm adjustment or frequentist familywise-error guarantee.

Joint-envelope empirical coverage: 0.9500. 9/30 marginal intervals and 4/30 joint intervals exclude zero.
These counts must be interpreted alongside the diagnostic flags; flagged source posteriors undermine interval reliability.

| Pair | Variable | Mean_pp | CI95 | Joint95 | Pr_increase | Diagnostics |
| --- | --- | --- | --- | --- | --- | --- |
| 1970s to 1990s | tmean_breeding | 0.161 | [-0.285, 0.615] | [-0.558, 0.880] | 0.756 | source diagnostics flagged/unassessable |
| 1970s to 1990s | prec_breeding | -0.095 | [-0.445, 0.268] | [-0.664, 0.475] | 0.296 | source diagnostics flagged/unassessable |
| 1970s to 1990s | hh | 0.104 | [-0.294, 0.513] | [-0.543, 0.751] | 0.689 | source diagnostics flagged/unassessable |
| 1970s to 1990s | perc_cropland | 0.052 | [-0.348, 0.478] | [-0.615, 0.720] | 0.588 | source diagnostics flagged/unassessable |
| 1970s to 1990s | perc_urban | -0.019 | [-0.261, 0.234] | [-0.419, 0.382] | 0.431 | source diagnostics flagged/unassessable |
| 1970s to 1990s | perc_pasture | 0.099 | [0.005, 0.199] | [-0.057, 0.256] | 0.981 | source diagnostics flagged/unassessable |
| 1970s to 1990s | perc_grass_shrub | 0.147 | [-0.012, 0.305] | [-0.109, 0.403] | 0.965 | source diagnostics flagged/unassessable |
| 1970s to 1990s | perc_forest_deciduous | 0.030 | [-0.156, 0.214] | [-0.265, 0.325] | 0.631 | source diagnostics flagged/unassessable |
| 1970s to 1990s | perc_forest_evergreen | 0.324 | [0.125, 0.541] | [-0.009, 0.658] | 0.999 | source diagnostics flagged/unassessable |
| 1970s to 1990s | perc_forest_mixed_other | 0.181 | [-0.212, 0.590] | [-0.470, 0.832] | 0.808 | source diagnostics flagged/unassessable |
| 1990s to 2010s | tmean_breeding | -0.674 | [-1.072, -0.268] | [-1.315, -0.034] | 0.001 | source diagnostics flagged/unassessable |
| 1990s to 2010s | prec_breeding | 1.116 | [0.624, 1.658] | [0.269, 1.962] | 1.000 | source diagnostics flagged/unassessable |
| 1990s to 2010s | hh | -0.118 | [-0.499, 0.256] | [-0.718, 0.483] | 0.269 | source diagnostics flagged/unassessable |
| 1990s to 2010s | perc_cropland | 0.060 | [-0.368, 0.474] | [-0.612, 0.732] | 0.616 | source diagnostics flagged/unassessable |
| 1990s to 2010s | perc_urban | 0.215 | [-0.045, 0.475] | [-0.197, 0.627] | 0.949 | source diagnostics flagged/unassessable |
| 1990s to 2010s | perc_pasture | 0.034 | [-0.083, 0.157] | [-0.158, 0.226] | 0.710 | source diagnostics flagged/unassessable |
| 1990s to 2010s | perc_grass_shrub | 0.119 | [-0.086, 0.325] | [-0.208, 0.447] | 0.868 | source diagnostics flagged/unassessable |
| 1990s to 2010s | perc_forest_deciduous | -0.099 | [-0.266, 0.060] | [-0.360, 0.161] | 0.109 | source diagnostics flagged/unassessable |
| 1990s to 2010s | perc_forest_evergreen | 0.144 | [-0.131, 0.436] | [-0.304, 0.592] | 0.850 | source diagnostics flagged/unassessable |
| 1990s to 2010s | perc_forest_mixed_other | 0.084 | [-0.360, 0.510] | [-0.611, 0.778] | 0.651 | source diagnostics flagged/unassessable |
| 1970s to 2010s | tmean_breeding | -0.513 | [-0.964, -0.048] | [-1.247, 0.221] | 0.017 | source diagnostics flagged/unassessable |
| 1970s to 2010s | prec_breeding | 1.021 | [0.517, 1.584] | [0.147, 1.895] | 1.000 | source diagnostics flagged/unassessable |
| 1970s to 2010s | hh | -0.014 | [-0.389, 0.359] | [-0.618, 0.591] | 0.470 | source diagnostics flagged/unassessable |
| 1970s to 2010s | perc_cropland | 0.112 | [-0.315, 0.550] | [-0.580, 0.805] | 0.694 | source diagnostics flagged/unassessable |
| 1970s to 2010s | perc_urban | 0.196 | [-0.098, 0.499] | [-0.277, 0.669] | 0.910 | source diagnostics flagged/unassessable |
| 1970s to 2010s | perc_pasture | 0.133 | [0.018, 0.258] | [-0.059, 0.325] | 0.990 | source diagnostics flagged/unassessable |
| 1970s to 2010s | perc_grass_shrub | 0.266 | [0.072, 0.466] | [-0.049, 0.581] | 0.997 | source diagnostics flagged/unassessable |
| 1970s to 2010s | perc_forest_deciduous | -0.069 | [-0.235, 0.091] | [-0.330, 0.191] | 0.202 | source diagnostics flagged/unassessable |
| 1970s to 2010s | perc_forest_evergreen | 0.468 | [0.216, 0.745] | [0.050, 0.887] | 1.000 | source diagnostics flagged/unassessable |
| 1970s to 2010s | perc_forest_mixed_other | 0.264 | [-0.169, 0.706] | [-0.443, 0.972] | 0.882 | source diagnostics flagged/unassessable |

## Convergence and numerical uncertainty

Diagnostics use the original chain arrays. The checks require rank-normalised R-hat < 1.1, bulk and tail ESS of at least 100 per chain, and MCSE of the mean at most 5% of posterior SD. Community and distribution contrasts carry the diagnostic status of the contributing atlas summaries.

| Atlas | Quantity | Pass | Flagged | Unassessable |
| --- | --- | --- | --- | --- |
| 1 | scaled | 1273 | 277 | 0 |
| 1 | tjur |  152 |   3 | 0 |
| 1 | community |    0 |  10 | 0 |
| 1 | distribution |    3 |  57 | 0 |
| 2 | scaled | 1232 | 318 | 0 |
| 2 | tjur |  106 |  49 | 0 |
| 2 | community |    2 |   8 | 0 |
| 2 | distribution |   35 |  25 | 0 |
| 3 | scaled | 1018 | 532 | 0 |
| 3 | tjur |   85 |  70 | 0 |
| 3 | community |    0 |  10 | 0 |
| 3 | distribution |    7 |  53 | 0 |

Directional probabilities have finite MCMC and draw-pairing error. They are not p-values; 1 - Pr(increase) is not sent to Holm. Native MCMC MCSE is recorded in atlas summaries. Repeated contrast triplets do not remedy insufficient effective sample sizes.

## Product-of-means versus mean-of-products

The new posterior mean of VP × Tjur can differ from mean(VP) × mean(Tjur), because VP and fit depend on the same parameters. The legacy figure value is included separately; differences can also reflect the legacy fit estimate. This script does not overwrite that figure.

| Atlas | variable | posterior_mean_scaled_percent | product_posterior_means_percent | legacy_figure_percent |
| --- | --- | --- | --- | --- |
| 1 | tmean_breeding | 1.408 | 1.413 | 1.401 |
| 1 | prec_breeding | 0.855 | 0.858 | 0.852 |
| 1 | hh | 0.891 | 0.899 | 0.900 |
| 1 | perc_cropland | 1.307 | 1.313 | 1.326 |
| 1 | perc_urban | 0.740 | 0.744 | 0.739 |
| 1 | perc_pasture | 0.223 | 0.225 | 0.227 |
| 1 | perc_grass_shrub | 0.298 | 0.300 | 0.305 |
| 1 | perc_forest_deciduous | 0.411 | 0.413 | 0.413 |
| 1 | perc_forest_evergreen | 0.372 | 0.375 | 0.376 |
| 1 | perc_forest_mixed_other | 1.163 | 1.166 | 1.166 |
| 2 | tmean_breeding | 1.569 | 1.573 | 1.571 |
| 2 | prec_breeding | 0.760 | 0.764 | 0.758 |
| 2 | hh | 0.995 | 1.001 | 1.006 |
| 2 | perc_cropland | 1.359 | 1.364 | 1.376 |
| 2 | perc_urban | 0.721 | 0.724 | 0.721 |
| 2 | perc_pasture | 0.323 | 0.324 | 0.327 |
| 2 | perc_grass_shrub | 0.445 | 0.448 | 0.451 |
| 2 | perc_forest_deciduous | 0.441 | 0.444 | 0.441 |
| 2 | perc_forest_evergreen | 0.697 | 0.700 | 0.699 |
| 2 | perc_forest_mixed_other | 1.343 | 1.349 | 1.343 |
| 3 | tmean_breeding | 0.895 | 0.900 | 0.907 |
| 3 | prec_breeding | 1.876 | 1.883 | 1.854 |
| 3 | hh | 0.877 | 0.884 | 0.891 |
| 3 | perc_cropland | 1.419 | 1.429 | 1.436 |
| 3 | perc_urban | 0.936 | 0.943 | 0.936 |
| 3 | perc_pasture | 0.357 | 0.358 | 0.363 |
| 3 | perc_grass_shrub | 0.565 | 0.569 | 0.575 |
| 3 | perc_forest_deciduous | 0.342 | 0.345 | 0.345 |
| 3 | perc_forest_evergreen | 0.841 | 0.846 | 0.846 |
| 3 | perc_forest_mixed_other | 1.427 | 1.435 | 1.438 |

## Species and distribution results

comparison.rds contains species_summary (4,650 paired species contrasts), distribution_summary (quantile and IQR changes), cdf_difference_summary (pointwise bands for cumulative-distribution differences), and wasserstein_summary (descriptive distributional distance). Atlas RDS files contain cumulative-distribution summaries and full scaled VP draws.

Quantile changes show whether importance shifted broadly or mainly in a tail; IQR changes describe spread. Cumulative curves show the proportion of species at or below a given VP value. CDF bands are pointwise, not simultaneous over the grid or variables, and CDF-specific convergence has not been assessed. A positive Wasserstein distance is not evidence by itself against a no-change null; continuous uncertain posteriors can produce positive distances even for indistinguishable distributions.

## What this does and does not establish

This quantifies uncertainty in changes in explanatory importance conditional on the fitted models and observed data, including dependence represented by the Hmsc posterior. It does not prove causality, isolate changing response slopes, correct spatial confounding, account for every source of survey error, or establish performance at new sites. Environmental variation, fitted cell coverage, species distributions, response estimates and overall fit can all affect scaled VP.

## Reproducibility

```r
x <- readRDS("scripts/4-post/compare-scaled-vp-changes/spatial/comparison.rds")
x$community_summary
x$species_summary
x$provenance
```

Run `Rscript scripts/4-post/compare-scaled-vp-changes.R --self-test` for synthetic checks. Use `--spatial` for the 2 November spatial fits or run without arguments for the 6 November environment-only fits. Options include --model-root, --model-prefix, --model-suffix, --samples, --thin, --chains, --seed, --contrast-draws and --output-dir. Model-root defaults to CH1_MODEL_OUTPUT_DIR when set. Contrast draws must be a multiple of the retained draw count.
