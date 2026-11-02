source('scripts/2-hmsc/helpers/pipeline.R')
suppressPackageStartupMessages({library(Hmsc);library(jsonify)})
for (name in model_names()) {
  p <- paths(name); check_destination(p$fitted)
  models <- read_objects(p$unfitted,'models')$models
  stopifnot(length(models)==1L); validate_model(models[[1]])
  expected <- file.path(p$models,'Sampled',sprintf('HPC_samples_%04d_thin_%02d_chain_%d.rds',cfg$samples,cfg$thin,0:(cfg$chains-1L)))
  actual <- Sys.glob(file.path(p$models,'Sampled',sprintf('HPC_samples_%04d_thin_%02d_chain_*.rds',cfg$samples,cfg$thin)))
  check_files(actual,expected)
  init <- jsonify::from_json(readRDS(p$init)[[1]])
  decoded <- lapply(expected,read_hpc)
  posteriors <- importPosteriorFromHPC(models[[1]],lapply(decoded,`[[`,1),cfg$samples,cfg$thin,cfg$transient,alignPost=TRUE)
  posteriors$adaptNf <- init$adaptNf; posteriors$verbose <- init$verbose
  validate_fit(posteriors)
  fitted_model <- list(posteriors=posteriors,fit_times=lapply(decoded,`[[`,2))
  save(fitted_model,file=p$fitted)
  provenance(dirname(p$fitted),'S2b_HPC_output_merger',c(p$unfitted,p$init,expected))
  message('Merged all ',cfg$chains,' chains: ',p$fitted)
}
