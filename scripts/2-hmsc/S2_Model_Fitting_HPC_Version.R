# Export initializations; this does not run production MCMC.
source('scripts/2-hmsc/helpers/pipeline.R')
suppressPackageStartupMessages({library(Hmsc);library(jsonify)})
for (name in model_names()) {
  p <- paths(name); models <- read_objects(p$unfitted,'models')$models
  stopifnot(length(models)==1L); validate_model(models[[1]])
  check_destination(p$init)
  init <- sampleMcmc(models[[1]],samples=cfg$samples,thin=cfg$thin,
                     transient=cfg$transient,nChains=cfg$chains,verbose=1,engine='HPC')
  saveRDS(to_json(init),p$init)
  provenance(dirname(p$init),'S2_Model_Fitting_HPC_Version',p$unfitted)
  message('Exported ',p$init)
}
