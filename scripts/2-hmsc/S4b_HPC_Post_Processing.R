# Import one fold at a time to limit memory; concatenate retained chains.
source('scripts/2-hmsc/helpers/pipeline.R')
suppressPackageStartupMessages({library(Hmsc);library(jsonify)})
for (name in model_names()) {
  p <- paths(name); check_destination(p$fit)
  info <- read_objects(p$folds,c('partition','seeds','input_md5','row_ids','idfold','threads','hM1','hM'))
  hM <- info$hM; validate_fit(hM)
  stopifnot(identical(info$input_md5,unname(tools::md5sum(p$fitted))),
            identical(info$row_ids,rownames(hM$Y)),info$threads==cfg$nfolds*cfg$chains,
            length(info$partition)==hM$ny,identical(sort(unique(info$partition)),seq_len(cfg$nfolds)))
  expected <- file.path(p$temp,sprintf('Sampled_HPC_samples_%04d_thin_%02d_thread_%d.rds',cfg$samples,cfg$thin,seq_len(info$threads)))
  actual <- Sys.glob(file.path(p$temp,sprintf('Sampled_HPC_samples_%04d_thin_%02d_thread_*.rds',cfg$samples,cfg$thin)))
  check_files(actual,expected)
  ndraws <- cfg$chains*length(seq.int(1,cfg$samples,by=cfg$posterior_thin))
  predArray <- array(NA_real_,c(hM$ny,hM$ns,ndraws))
  for (fold in seq_len(cfg$nfolds)) {
    ids <- which(info$idfold==fold); stopifnot(length(ids)==cfg$chains)
    decoded <- lapply(expected[ids],read_hpc)
    m <- importPosteriorFromHPC(info$hM1[[fold]],lapply(decoded,`[[`,1),cfg$samples,cfg$thin,cfg$transient,alignPost=TRUE)
    post <- posterior_draws(m); stopifnot(length(post)==ndraws)
    val <- info$partition==fold
    predictions <- predict(m,post=post,XData=hM$XData[val,,drop=FALSE],
                           studyDesign=droplevels(hM$studyDesign[val,,drop=FALSE]),
                           mcmcStep=1,expected=TRUE,nParallel=cfg$nparallel,useSocket=FALSE)
    predArray[val,,] <- prediction_array(predictions,sum(val),hM$ns,ndraws)
    rm(m,post,decoded,predictions);gc()
    message('Predicted fold ',fold,'/',cfg$nfolds,' using ',ndraws,' draws')
  }
  stopifnot(!anyNA(predArray))
  MFCV <- evaluateModelFit(hM,predY=predArray)
  rm(predArray);gc()
  preds <- computePredictedValues(hM,thin=cfg$posterior_thin,expected=TRUE,nParallel=cfg$nparallel)
  MF <- evaluateModelFit(hM,predY=preds)
  rm(preds);gc()
  WAIC <- computeWAIC(hM)
  validation <- list(draws=ndraws,posterior_thin=cfg$posterior_thin,nfolds=cfg$nfolds,
                     method='random folds',folds_md5=unname(tools::md5sum(p$folds)))
  save(MF,MFCV,WAIC,validation,file=p$fit)
  provenance(dirname(p$fit),'S4b_HPC_Post_Processing',c(p$fitted,p$folds,expected),paste('Prediction draws:',ndraws))
}
