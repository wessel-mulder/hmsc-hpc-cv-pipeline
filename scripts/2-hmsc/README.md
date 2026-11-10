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
S1 loads the portable shared RData and its bundled tree. Model construction
has been tested for all three atlases: 10 environmental variables, 155 species
and 155 tree tips per model. No MCMC fitting or prediction has been run for the
revised dataset; subsequent scripts and scaling conventions require review.

## Portable HPC inputs

S0 also exports its final RData with bundled taxonomy tree into tracked `data-shared/`.
Local `Data/` and `data-new/` are no longer tracked. S1 loads exclusively from
`data-shared/` and creates its output root if necessary. Run S1 from the
repository root on the HPC; set `CH1_MODEL_OUTPUT_DIR` for another output path.
This supersedes the older path caveats above. Fitting and prediction scaling
still require their own checks; no MCMC sampling is performed by S1.

## Combined S3 diagnostics

Both `S3_check_model_fit_TensorFlow.R` and the compatibility wrapper
`S3a_show_model_fit.R` run the same workflow. A literal folder-name substring is
required; the other positional arguments retain their legacy defaults:

```sh
CH1_MODEL_OUTPUT_DIR=data-new/2-hmsc/derived/models \
  Rscript scripts/2-hmsc/S3_check_model_fit_TensorFlow.R 2026-11-02 250 100 4 4
Rscript scripts/2-hmsc/model-diagnostics/test-diagnostics.R
```

Arguments are folder pattern, saved samples per chain, fitting thinning interval,
chains and prediction workers. The root defaults to `HmscOutputs`. Counts and
thinning are checked against the actual fitted object; no additional thinning,
chain alignment or chain exclusion occurs. Every matching folder must contain
the requested fit. A common matrix `X` design is required; species-specific or
reduced-rank designs fail explicitly rather than silently misindexing parameters.
Dependencies: Hmsc, posterior, ggplot2, rmarkdown, knitr, Pandoc and XeLaTeX.

Each model receives `Results/model-diagnostics/HPC_samples_.../report.html`,
`report.pdf`, standalone PNGs, CSVs, `provenance.R` and session information. The model root also
receives `model-diagnostics-comparison/<pattern>-<fit>/report.html`. Reports are
self-contained; keep their adjacent CSVs for downloadable full tables. The
comparison links to model reports on the local filesystem. All spatial levels
with positive `sDim` are assessed, first on each model's fitted cells and then on
the intersection of the selected models' `studyDesign$site` values. Unique site
identifiers are required for spatial comparisons. Saved factor counts can vary;
individual-factor summaries remain conditional on the saved basis.

Hmsc exports Beta on the `X` scale (`combineParameters` reverses internal
`XScaled` scaling). Spatial projections therefore use `X`; correlations are
invariant to positive rescaling. Raw Beta is never transformed. Full and
second-half summaries use the same thresholds, with the latter clearly marked
as sensitivity results. The upper within-predictor quartile defines the broad
uncertainty overlap screen; it is a descriptive convention, not a test.

`MF_preCV_samples_...Rdata` still contains `MF` for downstream readers. Reuse
requires an additional matching `MF_provenance` object containing the input MD5,
Hmsc version and prediction mode. Legacy or stale caches are recomputed and
replaced. Bulky predictions are transient and never exported. Fitted inputs are
MD5 checked before and after the workflow. Reruns replace generated PNGs in the
specific diagnostics directory to avoid stale coefficient selections.

The R Markdown template knits both self-contained HTML and landscape PDF.
PDFs show up to ten rows and seven relevant columns per table; HTML shows up to
fifteen rows with all columns. Complete CSVs accompany both. PDF heatmaps use
at most 32 species per page; HTML retains the full species heatmaps. All figures
are PNGs, including those embedded in the PDFs. For the full synthetic legacy
wrapper check, run the test script with `--integration` (also requires XeLaTeX).
To validate completed reports against saved fitted draws and fingerprints, use:

```sh
CH1_MODEL_OUTPUT_DIR=data-new/2-hmsc/derived/models \
  Rscript scripts/2-hmsc/model-diagnostics/test-diagnostics.R --models
```
The workers argument is forwarded to Hmsc. In the installed Hmsc version,
`computePredictedValues()`'s in-sample branch calls `predict()` without forwarding
that argument, so that branch runs serially; this workflow does not patch Hmsc.
