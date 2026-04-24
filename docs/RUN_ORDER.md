# Run Order

This document describes the intended order of the main analysis scripts. It is a
workflow map rather than a fully automated pipeline.

## 1. Prepare Modelling Data

Run:

```r
source("S0_Data_Definitions.R")
```

Outputs:

- `Data/preprocessed_data_min_occs_is_5_coverage_is_good.RData`
- `Data/preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData`

Main filters:

- atlas grid cells with at least 25% land
- survey coverage of `Good` or `Good + Average`
- species with at least `min_occs` occurrences in each atlas
- species names reconciled between DOF/taxonomy, trait, phylogeny, and STI data

## 2. Define HMSC Models

Run:

```r
source("S1_Model_Definitions.R")
```

Outputs per model folder:

- `Models/Unfitted/unfitted_models.RData`
- `Tests/test_atlases.RData`

The current model uses probit occurrence responses, environmental predictors,
species traits, a phylogenetic tree, and a spatial random level for site.

## 3. Prepare HPC Initial States

Run locally or on the HPC login node:

```r
source("S2_Model_Fitting_HPC_Version.R")
```

Output:

- `Models/INITS/HPC_INIT_samples_0250_thin_100_chains_4.rds`

## 4. Fit Main Models On HPC

Submit the SLURM batch script with exported model settings, for example:

```bash
sbatch --job-name="<MODEL_FOLDER>" \
  --export=NAME="<MODEL_FOLDER>",SAMP=250,THIN=100 \
  S2a_Fit_HMSC_Model.sh
```

Outputs:

- `Models/Sampled/HPC_samples_0250_thin_100_chain_<CHAIN>.rds`
- `Models/Logs/`

## 5. Merge HPC Chains

Run:

```r
source("S2b_HPC_output_merger.R")
```

Output:

- `Models/Fitted/HPC_samples_0250_thin_100_chains_4.Rdata`

## 6. Check Main Fit

Run:

```r
source("S3_check_model_fit_TensorFlow.R")
```

Typical outputs:

- MCMC convergence plots
- model-fit summary files

## 7. Cross-Validation

Create partition initial states:

```r
source("S4_Partition_Creation.R")
```

Submit k-fold sampling:

```bash
bash S4a_Fit_HMSC_Model_KFold_Launch.sh
```

Merge and evaluate k-fold predictions:

```bash
bash S4b_HPC_Post_Processing_Launch.sh
```

Outputs:

- `Models/Fitted/MF_samples_0250_thin_100_chains_4_nfolds_10.rdata`
- cross-validation logs and temporary files

## 8. Plot Fit, Predictions, And Parameters

Run:

```r
source("S5_show_model_fit.R")
source("S6_Making_Predictions.R")
source("S6b_show_test_fit.R")
source("S7_Show_parameter_estimates.R")
```

Outputs:

- model fit PDFs
- environmental-gradient predictions
- cross-atlas prediction tests
- beta, gamma, omega, and variance-partitioning exports

## 9. Publication Figures

Figure scripts are in `misc-figures/scripts/`. Rendered outputs are organized
under `misc-figures/outputs/`.
