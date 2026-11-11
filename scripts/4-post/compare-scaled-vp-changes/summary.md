# Posterior changes in Tjur R²-scaled VP

Completed 2026-11-11 11:29:49.995492. Local runtime: 70.9 seconds (1.18 minutes).

## What was calculated

Model family: environment-only (2026-11-06_05-11-30).

All 155 species, ten predictors, four chains and every saved draw are retained. Each atlas uses its own fitted cells; there is no common-cell restriction or model refitting.

For each draw: compute Hmsc environmental VP; compute expected occurrence probabilities from that same draw; calculate species Tjur R²; multiply VP by Tjur R². Values are percentages and changes are percentage points. Negative draw-level fit values are retained.

Cells by atlas: 1747, 1821, 1875.
Negative Tjur draw/species entries by atlas: 12, 101, 301.

Independently combine whole atlas posterior draws using 10000 balanced draw triplets, seed 369. Reuse triplets across all contrasts; retain dependence among species and variables. These triplets do not create new observations or increase the original MCMC effective sample size.

## Community comparisons

Differences are later minus earlier. CI95 is a marginal credible interval. Joint95 is a joint posterior envelope across all 30 community contrasts, not a Holm adjustment or frequentist familywise-error guarantee.

Joint-envelope empirical coverage: 0.9500. 19/30 marginal intervals and 18/30 joint intervals exclude zero.
These counts must be interpreted alongside the diagnostic flags; flagged source posteriors undermine interval reliability.

| Pair | Variable | Mean_pp | CI95 | Joint95 | Pr_increase | Diagnostics |
| --- | --- | --- | --- | --- | --- | --- |
| 1970s to 1990s | tmean_breeding | 0.326 | [0.167, 0.481] | [0.075, 0.578] | 1.000 | pass |
| 1970s to 1990s | prec_breeding | 0.221 | [0.101, 0.337] | [0.033, 0.410] | 1.000 | pass |
| 1970s to 1990s | hh | -0.061 | [-0.231, 0.103] | [-0.331, 0.210] | 0.241 | pass |
| 1970s to 1990s | perc_cropland | 0.473 | [0.283, 0.664] | [0.169, 0.777] | 1.000 | pass |
| 1970s to 1990s | perc_urban | 0.191 | [0.058, 0.330] | [-0.026, 0.408] | 0.997 | pass |
| 1970s to 1990s | perc_pasture | 0.005 | [-0.079, 0.089] | [-0.129, 0.139] | 0.554 | pass |
| 1970s to 1990s | perc_grass_shrub | 0.202 | [0.118, 0.287] | [0.067, 0.336] | 1.000 | pass |
| 1970s to 1990s | perc_forest_deciduous | -0.074 | [-0.184, 0.043] | [-0.254, 0.106] | 0.099 | pass |
| 1970s to 1990s | perc_forest_evergreen | 0.110 | [-0.010, 0.228] | [-0.078, 0.299] | 0.963 | pass |
| 1970s to 1990s | perc_forest_mixed_other | 0.379 | [0.203, 0.568] | [0.082, 0.677] | 1.000 | pass |
| 1990s to 2010s | tmean_breeding | -0.724 | [-0.855, -0.590] | [-0.935, -0.514] | 0.000 | pass |
| 1990s to 2010s | prec_breeding | 1.435 | [1.284, 1.593] | [1.186, 1.683] | 1.000 | pass |
| 1990s to 2010s | hh | 0.049 | [-0.122, 0.219] | [-0.222, 0.319] | 0.712 | pass |
| 1990s to 2010s | perc_cropland | 0.176 | [-0.032, 0.388] | [-0.158, 0.511] | 0.948 | pass |
| 1990s to 2010s | perc_urban | 0.351 | [0.195, 0.508] | [0.099, 0.604] | 1.000 | pass |
| 1990s to 2010s | perc_pasture | -0.026 | [-0.107, 0.056] | [-0.155, 0.102] | 0.263 | pass |
| 1990s to 2010s | perc_grass_shrub | 0.060 | [-0.035, 0.156] | [-0.093, 0.212] | 0.889 | pass |
| 1990s to 2010s | perc_forest_deciduous | -0.222 | [-0.329, -0.117] | [-0.392, -0.051] | 0.000 | pass |
| 1990s to 2010s | perc_forest_evergreen | 0.445 | [0.288, 0.596] | [0.198, 0.691] | 1.000 | pass |
| 1990s to 2010s | perc_forest_mixed_other | 0.090 | [-0.106, 0.279] | [-0.221, 0.400] | 0.811 | pass |
| 1970s to 2010s | tmean_breeding | -0.398 | [-0.538, -0.259] | [-0.619, -0.178] | 0.000 | pass |
| 1970s to 2010s | prec_breeding | 1.656 | [1.504, 1.812] | [1.407, 1.905] | 1.000 | pass |
| 1970s to 2010s | hh | -0.012 | [-0.183, 0.162] | [-0.288, 0.264] | 0.446 | pass |
| 1970s to 2010s | perc_cropland | 0.649 | [0.442, 0.849] | [0.323, 0.976] | 1.000 | pass |
| 1970s to 2010s | perc_urban | 0.542 | [0.393, 0.695] | [0.303, 0.781] | 1.000 | pass |
| 1970s to 2010s | perc_pasture | -0.021 | [-0.103, 0.062] | [-0.152, 0.111] | 0.309 | pass |
| 1970s to 2010s | perc_grass_shrub | 0.261 | [0.177, 0.349] | [0.124, 0.398] | 1.000 | pass |
| 1970s to 2010s | perc_forest_deciduous | -0.296 | [-0.395, -0.188] | [-0.460, -0.131] | 0.000 | pass |
| 1970s to 2010s | perc_forest_evergreen | 0.555 | [0.407, 0.695] | [0.323, 0.787] | 1.000 | pass |
| 1970s to 2010s | perc_forest_mixed_other | 0.469 | [0.279, 0.660] | [0.167, 0.771] | 1.000 | pass |

## Convergence and numerical uncertainty

Diagnostics use the original chain arrays. The checks require rank-normalised R-hat < 1.1, bulk and tail ESS of at least 100 per chain, and MCSE of the mean at most 5% of posterior SD. Community and distribution contrasts carry the diagnostic status of the contributing atlas summaries.

| Atlas | Quantity | Pass | Flagged | Unassessable |
| --- | --- | --- | --- | --- |
| 1 | scaled | 1550 | 0 | 0 |
| 1 | tjur |  155 | 0 | 0 |
| 1 | community |   10 | 0 | 0 |
| 1 | distribution |   60 | 0 | 0 |
| 2 | scaled | 1549 | 1 | 0 |
| 2 | tjur |  155 | 0 | 0 |
| 2 | community |   10 | 0 | 0 |
| 2 | distribution |   60 | 0 | 0 |
| 3 | scaled | 1543 | 7 | 0 |
| 3 | tjur |  155 | 0 | 0 |
| 3 | community |   10 | 0 | 0 |
| 3 | distribution |   60 | 0 | 0 |

Directional probabilities have finite MCMC and draw-pairing error. They are not p-values; 1 - Pr(increase) is not sent to Holm. Native MCMC MCSE is recorded in atlas summaries. Repeated contrast triplets do not remedy insufficient effective sample sizes.

## Product-of-means versus mean-of-products

The new posterior mean of VP × Tjur can differ from mean(VP) × mean(Tjur), because VP and fit depend on the same parameters. The legacy figure value is included separately; differences can also reflect the legacy fit estimate. This script does not overwrite that figure.

| Atlas | variable | posterior_mean_scaled_percent | product_posterior_means_percent | legacy_figure_percent |
| --- | --- | --- | --- | --- |
| 1 | tmean_breeding | 1.068 | 1.065 | 1.065 |
| 1 | prec_breeding | 0.748 | 0.747 | 0.747 |
| 1 | hh | 0.914 | 0.918 | 0.918 |
| 1 | perc_cropland | 1.055 | 1.054 | 1.054 |
| 1 | perc_urban | 0.686 | 0.687 | 0.687 |
| 1 | perc_pasture | 0.439 | 0.438 | 0.438 |
| 1 | perc_grass_shrub | 0.338 | 0.338 | 0.338 |
| 1 | perc_forest_deciduous | 0.744 | 0.745 | 0.745 |
| 1 | perc_forest_evergreen | 0.543 | 0.544 | 0.544 |
| 1 | perc_forest_mixed_other | 1.647 | 1.646 | 1.646 |
| 2 | tmean_breeding | 1.394 | 1.393 | 1.393 |
| 2 | prec_breeding | 0.969 | 0.969 | 0.969 |
| 2 | hh | 0.854 | 0.857 | 0.857 |
| 2 | perc_cropland | 1.528 | 1.526 | 1.526 |
| 2 | perc_urban | 0.877 | 0.877 | 0.877 |
| 2 | perc_pasture | 0.444 | 0.443 | 0.443 |
| 2 | perc_grass_shrub | 0.540 | 0.540 | 0.540 |
| 2 | perc_forest_deciduous | 0.670 | 0.671 | 0.671 |
| 2 | perc_forest_evergreen | 0.653 | 0.655 | 0.655 |
| 2 | perc_forest_mixed_other | 2.027 | 2.025 | 2.025 |
| 3 | tmean_breeding | 0.670 | 0.668 | 0.668 |
| 3 | prec_breeding | 2.403 | 2.403 | 2.403 |
| 3 | hh | 0.902 | 0.906 | 0.906 |
| 3 | perc_cropland | 1.704 | 1.702 | 1.702 |
| 3 | perc_urban | 1.229 | 1.231 | 1.231 |
| 3 | perc_pasture | 0.418 | 0.416 | 0.416 |
| 3 | perc_grass_shrub | 0.599 | 0.600 | 0.600 |
| 3 | perc_forest_deciduous | 0.449 | 0.449 | 0.449 |
| 3 | perc_forest_evergreen | 1.098 | 1.099 | 1.099 |
| 3 | perc_forest_mixed_other | 2.116 | 2.114 | 2.114 |

## Species and distribution results

comparison.rds contains species_summary (4,650 paired species contrasts), distribution_summary (quantile and IQR changes), cdf_difference_summary (pointwise bands for cumulative-distribution differences), and wasserstein_summary (descriptive distributional distance). Atlas RDS files contain cumulative-distribution summaries and full scaled VP draws.

Quantile changes show whether importance shifted broadly or mainly in a tail; IQR changes describe spread. Cumulative curves show the proportion of species at or below a given VP value. CDF bands are pointwise, not simultaneous over the grid or variables, and CDF-specific convergence has not been assessed. A positive Wasserstein distance is not evidence by itself against a no-change null; continuous uncertain posteriors can produce positive distances even for indistinguishable distributions.

## What this does and does not establish

This quantifies uncertainty in changes in explanatory importance conditional on the fitted models and observed data, including dependence represented by the Hmsc posterior. It does not prove causality, isolate changing response slopes, correct spatial confounding, account for every source of survey error, or establish performance at new sites. Environmental variation, fitted cell coverage, species distributions, response estimates and overall fit can all affect scaled VP.

## Reproducibility

```r
x <- readRDS("scripts/4-post/compare-scaled-vp-changes/comparison.rds")
x$community_summary
x$species_summary
x$provenance
```

Run `Rscript scripts/4-post/compare-scaled-vp-changes.R --self-test` for synthetic checks. Use `--spatial` for the 2 November spatial fits or run without arguments for the 6 November environment-only fits. Options include --model-root, --model-prefix, --model-suffix, --samples, --thin, --chains, --seed, --contrast-draws and --output-dir. Model-root defaults to CH1_MODEL_OUTPUT_DIR when set. Contrast draws must be a multiple of the retained draw count.
