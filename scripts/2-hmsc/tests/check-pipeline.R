# End-to-end local smoke test. Tiny R fits exercise downstream code; no HPC jobs.
# Rscript scripts/2-hmsc/tests/check-pipeline.R [output-directory]
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else tempfile('ch1-pipeline-smoke-')
if (dir.exists(root)) stop('Use a new output directory')
Sys.setenv(CH1_MODEL_OUTPUT_DIR=root,CH1_SAMPLES=6,CH1_THIN=1,CH1_TRANSIENT=4,
           CH1_CHAINS=2,CH1_NFOLDS=2,CH1_NGRID=4,CH1_NPARALLEL=1,CH1_POSTERIOR_THIN=1)
source('scripts/2-hmsc/helpers/pipeline.R')
suppressPackageStartupMessages({library(Hmsc);library(ape);library(jsonify)})
run_stage <- function(script) {
  log <- file.path(root,paste0(script,'.log'))
  status <- system2(file.path(R.home('bin'),'Rscript'),file.path('scripts/2-hmsc',paste0(script,'.R')),stdout=log,stderr=log)
  if (status!=0) stop('Failed stage ',script,': ',paste(readLines(log,warn=FALSE),collapse='\n'))
}
expect_error <- function(expr) {
  failed <- inherits(tryCatch({force(expr);NULL},error=function(e)e),'error')
  stopifnot(failed)
}
# Python HPC draws use scaled coefficients, precision iV, a rho index and sigma SD.
# Convert R-generated draws into that representation for serialization/import checks.
as_hpc_draw <- function(draw,m) {
  s <- draw; s$iV <- solve(s$V);s$V <- NULL
  for (k in rev(seq_len(m$ncNRRR))) {
    center <- m$XScalePar[1,k]; spread <- m$XScalePar[2,k]
    if (center!=0 || spread!=1) {
      s$Beta[m$XInterceptInd,] <- s$Beta[m$XInterceptInd,]+center*s$Beta[k,]
      s$Gamma[m$XInterceptInd,] <- s$Gamma[m$XInterceptInd,]+center*s$Gamma[k,]
      s$Beta[k,] <- s$Beta[k,]*spread;s$Gamma[k,] <- s$Gamma[k,]*spread
      s$iV[k,] <- s$iV[k,]/spread;s$iV[,k] <- s$iV[,k]/spread
    }
  }
  for (k in rev(seq_len(m$nt))) {
    center <- m$TrScalePar[1,k];spread <- m$TrScalePar[2,k]
    if (center!=0 || spread!=1) {
      s$Gamma[,m$TrInterceptInd] <- s$Gamma[,m$TrInterceptInd]+center*s$Gamma[,k]
      s$Gamma[,k] <- s$Gamma[,k]*spread
    }
  }
  s$rho <- which.min(abs(m$rhopw[,1]-draw$rho));s$sigma <- sqrt(draw$sigma)
  s$BetaSel <- list()
  s
}
# Exercise missing files, dimension protection, ambiguous dates and CPU allocation.
expect_error(check_files('missing',c('missing','second')))
expect_error(prediction_array(list(matrix(1,2,2)),2,2,2))
expect_error(model_names(character()))
set.seed(cfg$seed)
ns <- 8L; ny <- 30L
species <- paste0('Species_',seq_len(ns)); phy <- ape::rtree(ns,tip.label=species)
Tr <- data.frame(Migration_a3_DOF=factor(rep(c(1,2),4)),
                 foraging_guild_consensus=factor(rep(c('A','B'),each=4)),
                 species_thermal_index=seq_len(ns),row.names=species)
for (atlas in 1:3) {
  name <- paste0(cfg$run_date,'_Atlas',atlas,'_MinOccs5_CoverageGoodAverage');p <- paths(name)
  dir.create(dirname(p$unfitted),recursive=TRUE)
  for (folder in c('INITS','Fitted','Sampled','Temp','Logs','Logs_CV')) dir.create(file.path(p$models,folder))
  dir.create(p$results)
  X <- data.frame(tmean_breeding=rnorm(ny),prec_breeding=rnorm(ny),row.names=paste0('A',atlas,'_',seq_len(ny)))
  Y <- matrix(rbinom(ny*ns,1,.5),ny,ns,dimnames=list(rownames(X),species))
  design <- data.frame(site=factor(paste0('s',seq_len(ny))),row.names=rownames(X))
  xy <- matrix(runif(ny*2,0,100000),ny,2,dimnames=list(levels(design$site),c('X','Y')))
  rL <- setPriors(HmscRandomLevel(sData=xy,sMethod='Full'),nfMin=5,nfMax=5)
  model <- Hmsc(Y=Y,XData=X,XFormula=reformulate(colnames(X)),TrData=Tr,TrFormula=reformulate(colnames(Tr)),phyloTree=phy,
                studyDesign=design,ranLevels=list(site=rL),distr='probit',XScale=TRUE)
  models <- setNames(list(model),name);save(models,file=p$unfitted)
  fit <- sampleMcmc(model,samples=cfg$samples,thin=cfg$thin,transient=cfg$transient,
                    nChains=cfg$chains,nParallel=1,verbose=0)
  for (chain in seq_len(cfg$chains)) {
    file <- file.path(p$models,'Sampled',sprintf('HPC_samples_%04d_thin_%02d_chain_%d.rds',cfg$samples,cfg$thin,chain-1L))
    saveRDS(jsonify::to_json(list(lapply(fit$postList[[chain]],as_hpc_draw,m=fit),0)),file)
  }
}
run_stage('S2_Model_Fitting_HPC_Version')
run_stage('S2b_HPC_output_merger')
# Incomplete chain sets must be rejected before import.
p <- paths(model_names(character())[1])
chain <- file.path(p$models,'Sampled',sprintf('HPC_samples_%04d_thin_%02d_chain_1.rds',cfg$samples,cfg$thin))
file.rename(chain,paste0(chain,'.hidden'))
expected <- file.path(p$models,'Sampled',sprintf('HPC_samples_%04d_thin_%02d_chain_%d.rds',cfg$samples,cfg$thin,0:1))
expect_error(check_files(Sys.glob(file.path(p$models,'Sampled','*.rds')),expected))
file.rename(paste0(chain,'.hidden'),chain)
run_stage('S4_Partition_Creation')
for (name in model_names(character())) {
  p <- paths(name);info <- read_objects(p$folds,c('hM1','partition','idfold','seeds'))
  stopifnot(all(table(info$partition)>0),length(info$seeds)==cfg$nfolds*cfg$chains)
  for (thread in seq_along(info$idfold)) {
    set.seed(info$seeds[thread])
    fit <- sampleMcmc(info$hM1[[info$idfold[thread]]],samples=cfg$samples,thin=cfg$thin,
                      transient=cfg$transient,nChains=1,verbose=0)
    file <- file.path(p$temp,sprintf('Sampled_HPC_samples_%04d_thin_%02d_thread_%d.rds',cfg$samples,cfg$thin,thread))
    saveRDS(jsonify::to_json(list(lapply(fit$postList[[1]],as_hpc_draw,m=fit),0)),file)
  }
}
run_stage('S4b_HPC_Post_Processing')
run_stage('S3_check_model_fit_TensorFlow')
run_stage('S5_show_model_fit')
run_stage('S6_Making_Predictions')
run_stage('S7_Show_parameter_estimates')
for (name in model_names(character())) {
  p <- paths(name); result <- read_objects(p$fit,'validation')$validation
  stopifnot(result$draws==cfg$samples*cfg$chains,result$posterior_thin==1)
  alpha <- read.csv(file.path(p$results,'parameter_Alpha_all_factors.csv'))
  stopifnot(nrow(alpha)==5L)
  curves <- readRDS(file.path(p$results,'Preds','response_tmean_breeding.rds'))
  stopifnot(curves$draws==cfg$samples*cfg$chains)
  m <- read_fit(p)
  # Both raw XData and its raw model matrix must receive the stored transformation.
  set.seed(42); a <- predict(m,post=posterior_draws(m),XData=m$XData,expected=TRUE,nParallel=1)
  set.seed(42); b <- predict(m,post=posterior_draws(m),X=m$X,expected=TRUE,nParallel=1)
  stopifnot(isTRUE(all.equal(a,b,tolerance=1e-10)))
}
stopifnot(!length(list.files(root,pattern='[.]pdf$',recursive=TRUE)),length(list.files(root,pattern='[.]png$',recursive=TRUE))>0)
# More than one model for an atlas/date must never silently select a run.
dir.create(file.path(root,paste0(cfg$run_date,'_12-00-00_Atlas1_MinOccs5_CoverageGoodAverage')))
expect_error(model_names(character()))
writeLines(c('PASS: all three atlas smoke pipelines',
             'R fits only; HPC serialization/import exercised with R-generated posterior draws.',
             'All retained draws used; five Alpha factors; PNG only; raw-data prediction equivalence.',
             'Missing chains, mismatched dimensions and ambiguous model folders rejected.'),file.path(root,'validation.txt'))
message('PASS. Smoke-test outputs: ',normalizePath(root))
