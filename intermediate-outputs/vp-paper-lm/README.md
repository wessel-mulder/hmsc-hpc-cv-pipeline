# VP trends: direct paper-method linear regressions

Source: Antão et al. (2022), Climate change reshuffles northern species within their niches; final paragraph of Joint species distribution modelling, Methods, and Figure 4. https://doi.org/10.1038/s41558-022-01381-x

For each of the eight environmental variables, fit lm(weighted_vp ~ decade_since_1970s), where weighted_vp is the species variance-partition fraction multiplied by its own atlas-specific full-model Tjur R². This is multiplication of the response before regression, not the weights argument to lm. Each fit pools 157 species across three atlases (471 observations, residual df = 469). Model outputs are the 2026-07-27 MinOccs5 CoverageGoodAverage fits; exported VP and Tjur values were checked against VP.rds and MF.rds.

Time is coded 0, 2, 4 for the 1970s, 1990s and 2010s. Slopes are therefore changes in weighted VP per calendar decade based on these representative decades. Coding atlas as 1, 2, 3 instead doubles the slope and standard error but leaves the t statistic, p-value and R² unchanged.

Outputs are ordinary lm slopes, standard errors, 95% confidence intervals, t tests, unadjusted p-values and regression R². No species term, bootstrap, paired test or multiple-testing adjustment is added. This project has one modelled taxonomic assemblage and region, so there is one trend per environmental variable; the paper’s comparisons among independently fitted taxonomic groups and zones are not applicable to these model outputs.

This reproduces the requested point-estimate regression. Its ordinary errors treat species-atlas rows as independent; HMSC posterior uncertainty and repeated-species dependence are not incorporated.

Run Rscript scripts/4-data-visualisation/atlas-vp-paper-lm.R from the project root. CSVs and fitted lm objects are stored here; PNG figures are in manuscript-outputs/figures/exploratory/vp-paper-lm/.

| Variable | Slope / decade | 95% CI | p |
|---|---:|---:|---:|
| Temperature | -0.000370 | [-0.001662, +0.000922] | 0.5739 |
| Precipitation | +0.000401 | [-0.000376, +0.001178] | 0.31142 |
| Land-use heterogeneity | +0.000225 | [-0.000267, +0.000716] | 0.36936 |
| Urban cover | +0.001168 | [+0.000353, +0.001983] | 0.0050461 |
| Cropland cover | +0.000733 | [-0.000361, +0.001826] | 0.18848 |
| Pasture cover | +0.000928 | [+0.000577, +0.001279] | 3.007e-07 |
| Forest cover | +0.002347 | [+0.000730, +0.003964] | 0.004539 |
| Grass / shrub cover | +0.000366 | [+0.000125, +0.000607] | 0.0030361 |
