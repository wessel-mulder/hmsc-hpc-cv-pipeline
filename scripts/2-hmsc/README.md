# Hmsc fitting workflow

Run every command from the repository root, including Slurm submissions.
`run-config.tsv` defines the run date **2026-11-02**, 250 retained samples per
chain, sampler thinning 100, transient 100000, four chains, ten random folds,
and seed 369. Environment variables with the same names override defaults.
`CH1_MODEL_OUTPUT_DIR` defaults to ignored `HmscOutputs/`.

Downstream `CH1_POSTERIOR_THIN=1` uses **all retained post-burn-in draws**:
1000 for a production four-chain fit. S3 and S7 always summarise all retained
draws. Sampler thinning 100 is unchanged. `CH1_OVERWRITE=false` prevents
accidental replacement; set it to `true` deliberately to regenerate outputs.
`CH1_NPARALLEL` defaults to one locally and the Slurm CPU allocation on HPC,
with any explicit value capped to that allocation.

## Commands

```sh
Rscript scripts/2-hmsc/S1_Model_Definitions.R
Rscript scripts/2-hmsc/S2_Model_Fitting_HPC_Version.R
bash scripts/2-hmsc/S2a_Fit_HMSC_Model_Launch.sh
# Wait for all full-fit jobs to finish successfully.
Rscript scripts/2-hmsc/S2b_HPC_output_merger.R
Rscript scripts/2-hmsc/S3_check_model_fit_TensorFlow.R
Rscript scripts/2-hmsc/S4_Partition_Creation.R
bash scripts/2-hmsc/S4a_Fit_HMSC_Model_KFold_Launch.sh
# Wait for all CV jobs to finish successfully.
sbatch scripts/2-hmsc/S4b_HPC_Post_Processing_Launch.sh
# Wait for post-processing to finish.
Rscript scripts/2-hmsc/S5_show_model_fit.R
Rscript scripts/2-hmsc/S6_Making_Predictions.R
Rscript scripts/2-hmsc/S7_Show_parameter_estimates.R
```

Inspect convergence diagnostics before interpreting models or continuing to
costly downstream work. Submission scripts do not establish Slurm dependencies.
The GPU jobs restore `myhmscstack`, source the existing cluster Conda setup,
and activate `CH1_CONDA_ENV` (default `hmsc_fix`). CPU post-processing uses
`CH1_POSTPROCESS_CONDA_ENV` (default `my_hmsc_env`), preserving the existing
cluster environments. Check that the GPU environment
contains the compatible Python Hmsc/TensorFlow versions before production use.

S1 uses the portable S0 RData in `data-shared/`, with 5443 surveys, 155 species
and ten environmental predictors. It prunes the bundled IOC taxonomic surrogate
tree to model species, treats migration and guild as factors, and passes raw
predictors to Hmsc (`XScale=TRUE`). Spatial random effects retain the Full
method and five latent factors. S1 does not run MCMC.

S2 exports initialization files. The new full-fit launcher submits all three
atlases. S2b requires exactly chains 0–3 and the configured draw counts before
merging. GPU workers record configuration, software and input checksums, write
to temporary output names, and rename only on successful completion.

S3 exports parameter-level Rhat and ESS diagnostics plus summaries. It disables
coda's additional automatic burn-in because sampler transient was already
removed. Undefined diagnostics remain explicit. All figures are PNGs.

S4 saves fold assignments and independent thread seeds and exports 40 files per
atlas. Re-running with `CH1_OVERWRITE=true` reuses those assignments/seeds and
requires the fitted input checksum to match. The CV launcher submits threads
1–40 for every atlas. S4b imports one fold at a time, concatenates chains,
checks exact prediction dimensions, and computes MF, MFCV and WAIC. Filenames
consistently use `.rdata` for model-fit results.

Validation uses **random folds**, measuring interpolation among sampled cells,
not spatially blocked transfer. Full-run predictor scaling is retained within
folds, matching Hmsc's existing CV convention. It uses no held-out responses
to estimate those predictor means/standard deviations.

S5 exports species-level training/CV metrics with finite-value plotting checks.
S6 saves environmental response curves separately per variable, using stored
training scaling. Its default uses all retained draws. Cross-atlas tests are
removed; S6b is retired and only displays an explanatory error.
S7 exports posterior means, intervals and sign supports, existing Beta/Gamma/
Omega workbooks, variance partitions, and all five latent-factor Alpha
summaries. Latent-factor identities are not comparable between independent
atlas models. Output figures use separate numbered PNG pages.

Each R stage also accepts one configured model folder name as its sole argument.
The resolver requires exactly one folder per atlas/date and rejects ambiguous
runs. Stage provenance records effective settings, input checksums and R
versions. Production outputs remain outside Git tracking.

## Local validation

```sh
Rscript scripts/2-hmsc/tests/check-S1-model-definitions.R /path/to/S1/output
Rscript scripts/2-hmsc/tests/check-pipeline.R /tmp/new-unique-smoke-output
```

The smoke test uses tiny, deliberately short R fits and exercises HPC
serialization/import with R-generated posterior draws. It is a functionality
check, not evidence of model convergence or validation of the cluster's Python
sampler. It checks all three atlases, retained draw counts, five Alpha factors,
raw-predictor scaling equivalence, missing-chain rejection, prediction
shape checks and ambiguous-folder rejection. No HPC jobs are submitted.
Real-data S1 construction and S2 initialization exports have also been checked
for all three atlases. Production MCMC and cluster execution remain to be run.
