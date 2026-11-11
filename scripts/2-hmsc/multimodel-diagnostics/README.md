# S2c: model-comparison figures

Edit the settings at the top of `../S2c_multimodel_diagnostics.R`: run-date
patterns, model root, fitted-file sampling settings, diagnostic thresholds and
atlas labels. Run from the repository root:

```sh
CH1_MODEL_OUTPUT_DIR=data-new/2-hmsc/derived/models \
  Rscript scripts/2-hmsc/S2c_multimodel_diagnostics.R
```

The output is **two comparison PNGs**, each showing all selected models, under
`<model root>/multimodel-diagnostics/`. Filenames include the actual selected run
dates: for November 2 and 6, `2026-11-02-06-multimodel-convergence.png` and
`2026-11-02-06-spatial-factor-correlations.png`. Dates in different months retain
their full dates, separated by underscores.

The convergence figure follows the tile layout of
`5-vis/convergence-diagnostics/community-vp-diagnostic-flags.png`.
All models are columns in each of four panels: R-hat, bulk ESS, tail ESS and
MCSE/posterior SD. Rows show Beta overall, each predictor including the
intercept, Gamma, V, rho, spatial ranges and residual species correlations.

Each cell prints the **median individual-parameter diagnostic**, followed by
the **percentage of assessed parameters flagged**. Colour follows the displayed
median: green on the better side, white exactly at the target, red on the worse
side. The percentage failing is separate from the colour. Missing
counts are labelled `NA`; inapplicable quantities have a dash. Availability is
assessed separately for each diagnostic: undefined tail ESS does not hide
finite R-hat, bulk ESS or mean MCSE. Full counts and
benchmark deviations are retained in `comparison-summary.csv`. These are not
diagnostics of a community-average posterior or of variance partitioning.

The colour scale saturates at R-hat 1 (green) and 1 + twice the target excess
above 1 (red); ESS twice its target (green) and zero (red); relative MCSE zero
(green) and twice its target (red). Raw values and targets remain unchanged.
Unavailable/inapplicable cells are grey.

Dates match literal folder-name substrings and are deduplicated. Every selected
fit must match the requested samples, fitting thin and chains. All retained
draws are used with chain identity preserved. There is no fitting, predictive
calculation, realignment, chain exclusion or extra thinning. S3/S3a are unchanged.

`comparison.rds` stores scalar diagnostics, spatial correlation/compensation
summaries, input fingerprints, targets and model metadata. `model-comparability.csv`
audits cell/species/predictor/scaling differences. Existing scalar diagnostics
are reused only when the complete selected input set has unchanged fingerprints;
thresholds and the two figures are regenerated on every run. No per-model,
trace, distribution or supplementary PNGs are generated.

Spatial calculations use all fitted cells separately per model. Individual
factor identities are conditional on the saved basis; combined effects use
`Eta %*% Lambda`. Alpha grid indices are converted to spatial distances.
Structurally fixed parameters are excluded; constant sampled parameters remain
unavailable. Better convergence is not itself evidence of ecological accuracy.

Tests:

```sh
Rscript scripts/2-hmsc/multimodel-diagnostics/test-diagnostics.R --integration
Rscript scripts/2-hmsc/multimodel-diagnostics/test-diagnostics.R --models data-new/2-hmsc/derived/models/multimodel-diagnostics
```

The convergence integration test checks that all selected models appear in
one figure and changing thresholds reuses diagnostics and overwrites that
same date-labelled figure. The real-model check verifies actual saved-draw extraction,
spatial calculations, matching and unchanged input fingerprints.


## Spatial-factor correlation figure

Each model has its own block of columns, one per saved spatial factor. Models
with no spatial factors get a grey placeholder; additional random levels get
separate blocks. The number of columns is read from the actual posterior, not
assumed from the model name. For factors absent in some draws, the column
label reports draw coverage and intervals condition on presence.

For every retained draw, Pearson correlation is calculated across that model's
fitted cells between each environmental design column and each saved Eta
factor (using Pi to map observations to latent-score rows). The tile shows the
posterior median; `*` marks a 95% posterior interval excluding zero. `†` marks
correlation-draw diagnostics that miss the configured targets or are undefined;
interpret those summaries cautiously. These are not spatial-independence
p-value tests. The intercept is constant and its correlations are undefined.

Factor identities/signs are local to the saved fit. No new rotations or sign
alignment are applied; LF1 in different models is not necessarily the same
spatial pattern. Colours encode signed correlation, blue negative through
white zero to orange positive. Numerical intervals, directional support,
draw coverage and diagnostics are in `spatial-factor-correlations.csv` and its
RDS cache. `spatial_interval` at the top of S2c changes the interval level.

Test flexible factor counts and spatial-row mapping:

```sh
Rscript scripts/2-hmsc/multimodel-diagnostics/test-spatial-factor-correlations.R
```
