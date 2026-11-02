# Construct separate, unfitted atlas models from the portable S0 inputs.
# Run from repository root: Rscript scripts/2-hmsc/S1_Model_Definitions.R
# CH1_MODEL_OUTPUT_DIR may redirect outputs (default: HmscOutputs).
source('scripts/2-hmsc/helpers/pipeline.R')
suppressPackageStartupMessages({library(Hmsc); library(ape)})

input_file <- 'data-shared/preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData'
output_dir <- cfg$root
atlas_models <- 1:3
run_id <- cfg$run_date

# S0 owns filtering, predictor selection and alignment. Load its five objects once.
if (!file.exists(input_file)) stop('S0 output not found: ', input_file)
inputs <- new.env(parent=emptyenv())
load(input_file, envir=inputs)
stopifnot(all(c('X','Y','Tr','Design','phy') %in% ls(inputs)))
X <- inputs$X; Y <- inputs$Y; Tr <- inputs$Tr; Design <- inputs$Design
expected_predictors <- c('tmean_breeding','prec_breeding','hh','perc_cropland',
                         'perc_urban','perc_pasture','perc_grass_shrub',
                         'perc_forest_deciduous','perc_forest_evergreen','perc_forest_mixed_other')
stopifnot(identical(colnames(X),expected_predictors),
          identical(rownames(X),rownames(Y)), identical(rownames(Y),rownames(Design)),
          identical(colnames(Y),rownames(Tr)), !anyDuplicated(rownames(X)),
          all(vapply(X,is.numeric,logical(1))), all(is.finite(as.matrix(X))),
          !anyNA(Y), all(as.matrix(Y) %in% c(0,1)), !anyNA(Tr),
          all(atlas_models %in% as.integer(as.character(Design$atlas))),
          all(is.finite(Design$lon)), all(is.finite(Design$lat)))

# Migration codes are categories, not a continuous gradient. Preserve reference choices.
for (trait in c('Migration_a3_DOF','foraging_guild_consensus')) {
  Tr[[trait]] <- factor(Tr[[trait]])
  reference <- names(which.max(table(Tr[[trait]])))
  Tr[[trait]] <- relevel(Tr[[trait]], ref=reference)
}
XFormula <- reformulate(colnames(X))
TrFormula <- reformulate(colnames(Tr))

# All atlas models share S0's species set, so validate and prune the tree once.
phy <- inputs$phy
if (!inherits(phy,'phylo')) stop('S0 phy must be an ape phylo tree')
missing_tips <- setdiff(colnames(Y),phy$tip.label)
if (length(missing_tips)) stop('Species missing from tree: ',paste(missing_tips,collapse=', '))
phy <- keep.tip(phy,colnames(Y))
stopifnot(setequal(phy$tip.label,colnames(Y)),length(phy$tip.label)==ncol(Y))

# Retain provenance for both batch execution and interactive sourcing.
file_arg <- grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)
script_path <- if (length(file_arg)) sub('^--file=','',file_arg[1]) else
  if (length(sys.frames()) && !is.null(sys.frames()[[1]]$ofile)) sys.frames()[[1]]$ofile else
  if (requireNamespace('rstudioapi',quietly=TRUE) && rstudioapi::isAvailable())
    rstudioapi::getSourceEditorContext()$path else NULL
if (!is.null(script_path) && file.exists(script_path)) script_path <- normalizePath(script_path)

for (atlas in atlas_models) {
  idx <- which(as.character(Design$atlas)==as.character(atlas))
  Y_sub <- Y[idx,,drop=FALSE]
  X_sub <- X[idx,,drop=FALSE]
  Design_sub <- droplevels(Design[idx,,drop=FALSE])
  stopifnot(identical(rownames(Y_sub),rownames(X_sub)),
            identical(rownames(Y_sub),rownames(Design_sub)),
            all(vapply(X_sub,sd,numeric(1))>0))

  # Projected coordinates are supplied by S0 (legacy lon/lat names, EPSG:23032).
  coordinates <- unique(data.frame(site=as.character(Design_sub$site),
                                   X=Design_sub$lon,Y=Design_sub$lat))
  if (anyDuplicated(coordinates$site)) stop('Inconsistent coordinates for the same site')
  rownames(coordinates) <- coordinates$site
  coordinates$site <- NULL
  struc_space <- HmscRandomLevel(sData=coordinates,sMethod='Full')
  struc_space <- setPriors(struc_space,nfMin=5,nfMax=5)

  # Hmsc is the single scaling owner; retain raw XData and training scale parameters.
  model <- Hmsc(Y=Y_sub,XData=X_sub,XFormula=XFormula,XScale=TRUE,
                TrData=Tr,TrFormula=TrFormula,phyloTree=phy,
                studyDesign=Design_sub[,'site',drop=FALSE],
                ranLevels=list(site=struc_space),distr='probit')
  model_description <- paste0(run_id,'_Atlas',atlas,'_MinOccs5_CoverageGoodAverage')
  # S2 loads models[[1]]; preserve this interface and the existing directory structure.
  models <- setNames(list(model),model_description)
  model_dir <- file.path(output_dir,model_description)
  folders <- c(file.path('Models',c('Fitted','Raw_HPC','Unfitted','INITS','Temp','Logs','Logs_CV','Sampled')),
               'Results','Tests')
  for (folder in folders) dir.create(file.path(model_dir,folder),recursive=TRUE,showWarnings=FALSE)
  unfitted_dir <- file.path(model_dir,'Models','Unfitted')
  destination <- file.path(unfitted_dir,'unfitted_models.RData')
  if (file.exists(destination)) stop('Refusing to overwrite an existing model: ',destination)
  save(models,file=destination)
  if (!is.null(script_path) && file.exists(script_path)) {
    if (!file.copy(script_path,file.path(unfitted_dir,paste0(run_id,'_S1_Model_Definitions.R'))))
      stop('Could not copy S1 for provenance')
  } else warning('Could not locate S1 for provenance')
  writeLines(c(paste('S0 input:',input_file),paste('Input MD5:',unname(tools::md5sum(input_file))),
               paste('Training atlas:',atlas),'Predictor scaling: Hmsc XScale=TRUE; XData is unscaled.',
               'Spatial method: Full; latent factors fixed at five. Distribution: probit.',
               'No cross-atlas test files generated.',capture.output(sessionInfo())),
             file.path(unfitted_dir,'model_definition_notes.txt'))
  provenance(unfitted_dir,'S1_Model_Definitions',input_file)
  message('Saved Atlas ',atlas,': ',nrow(Y_sub),' surveys, ',ncol(Y_sub),' species, ',ncol(X_sub),' predictors -> ',destination)
}
