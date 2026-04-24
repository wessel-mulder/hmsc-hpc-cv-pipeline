# Danish Bird Community Assembly HMSC Pipeline

This repository contains the modelling workflow for Chapter 1 of a PhD project on
community assembly in Danish breeding birds. The analysis uses Hierarchical
Modelling of Species Communities (HMSC) to relate species occurrences across
Danish bird atlas periods to climate, land-use, landscape heterogeneity,
species traits, phylogeny, and spatial structure.

## Repository Map

- `S0_Data_Definitions.R` prepares the modelling objects: environmental data
  (`X`), occurrence data (`Y`), traits (`Tr`), and study design (`Design`).
- `S1_Model_Definitions.R` creates unfitted HMSC model objects and test-atlas
  data bundles.
- `S2*_*.R` and `S2a_Fit_HMSC_Model.sh` prepare, fit, and merge TensorFlow/HPC
  model runs.
- `S3_check_model_fit_TensorFlow.R` checks MCMC convergence and model fit.
- `S4*_*.R` and `S4a_Fit_HMSC_Model_KFold.sh` run k-fold cross-validation on
  the HPC.
- `S5_show_model_fit.R` plots explanatory and predictive model fit.
- `S6_Making_Predictions.R` makes environmental-gradient and cross-atlas
  predictions.
- `S6b_show_test_fit.R` plots temporal-transfer test performance.
- `S7_Show_parameter_estimates.R` exports variance partitioning, beta/gamma
  estimates, and residual species associations.
- `support_scripts/` contains shared helper functions used by the staged
  scripts.
- `misc-figures/` contains publication and exploratory figure scripts and
  rendered figure outputs.
- `notebooks/` contains exploratory R Markdown work that supports interpretation
  but is not part of the core reproducible pipeline.
- `docs/` contains project notes and workflow documentation.

For the complete run order, see `docs/RUN_ORDER.md`.

## Main Data Products

The primary preprocessed data objects are:

- `Data/preprocessed_data_min_occs_is_5_coverage_is_good.RData`
- `Data/preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData`

These store `X`, `Y`, `Tr`, and `Design` after applying land-cover, survey
coverage, species occurrence, trait, and taxonomy filters.

## Main Model Outputs

HMSC output folders are stored under `HmscOutputs/`. Each model folder follows:

- `Models/Unfitted/`
- `Models/INITS/`
- `Models/Sampled/`
- `Models/Fitted/`
- `Models/Temp/`
- `Models/Logs/`
- `Models/Logs_CV/`
- `Results/`
- `Tests/`

The current publication-facing model family is the March 2026 atlas comparison:

- `2026-03-13_06-58-56_Atlas1_MinOccs5_CoverageGoodAverage`
- `2026-03-13_06-58-56_Atlas2_MinOccs5_CoverageGoodAverage`
- `2026-03-13_06-58-56_Atlas3_MinOccs5_CoverageGoodAverage`

## Notes On Reproducibility

Several scripts still use a `pattern2match` variable to select model-output
folders. This keeps interactive work flexible, but before publication runs it is
worth recording the exact model folder pattern in the manuscript or supplement.

The raw and intermediate data are included in this working repository. If this
repository is prepared for public release, check data permissions for atlas,
trait, climate, and land-cover sources before pushing the data files.
