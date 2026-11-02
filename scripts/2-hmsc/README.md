# 2 HMSC Fitting

This folder contains scripts and Slurm launchers for the HMSC model-fitting workflow.

## Main Workflow

- `S2_Model_Fitting_HPC_Version.R`: prepares HMSC model objects for TensorFlow/HPC fitting by exporting JSON/RDS objects.
- `S2a_Fit_HMSC_Model.sh`: Slurm launcher for fitting full HMSC models on the HPC.
- `S2b_HPC_output_merger.R`: merges sampled HPC output back into fitted HMSC model objects.
- `S3_check_model_fit_TensorFlow.R`: checks fitted model convergence and diagnostics.
- `S4_Partition_Creation.R`: creates cross-validation partitions and exports fold-specific model objects.
- `S4a_Fit_HMSC_Model_KFold.sh`: Slurm launcher for cross-validation fitting.
- `S4a_Fit_HMSC_Model_KFold_Launch.sh`: convenience launcher for submitting cross-validation jobs.
- `S4b_HPC_Post_Processing.R`: post-processes HPC and cross-validation output, including model-fit metrics.
- `S4b_HPC_Post_Processing_Launch.sh`: Slurm launcher for post-processing jobs.
- `S5_show_model_fit.R`: summarizes model fit.
- `S6_Making_Predictions.R`: creates model predictions and prediction-based tests.
- `S6b_show_test_fit.R`: summarizes held-out/test prediction performance.
- `S7_Show_parameter_estimates.R`: extracts variance partitioning, model fit, and parameter summaries.

## Helper Folder

- `S7_Show_parameter_estimates_scripts/`: helper scripts used by `S7_Show_parameter_estimates.R`.

## Current Path Caveat

These scripts still assume the original `HmscOutputs/` directory convention. In the publication structure, bulky fitted models and prediction objects should remain outside GitHub or in an external archive. Compact summaries needed for figures currently belong under `intermediate-outputs/hmsc-output/`.

Run these scripts from the repository root unless submitting a Slurm script from the HPC-side working directory.

## Current S0 dataset and readiness

S0 now reads the new `data-new/0-pre` inputs and writes the ten-predictor
forest-split RData plus audits to `data-new/2-hmsc/derived/`. Reproduce and
validate it with `scripts/reproduce-preprocessing.R`; see `scripts/0-pre/README.md`.
S1 and subsequent scripts retain historical paths and assumptions. They have
not been run on this revised dataset and must be adapted before use. In
particular, S1's input RData/tree paths, variable selection and predictor scaling
need review. The current audit establishes reproducibility through S0, not
fitted model results or downstream predictions.
