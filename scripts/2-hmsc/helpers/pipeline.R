# Shared configuration and IO for S1-S7. Run entry points from repository root.
.libPaths(c('~/Rlibs', .libPaths()))
config_file <- 'scripts/2-hmsc/run-config.tsv'
settings <- read.delim(config_file, colClasses='character', check.names=FALSE)
setting <- function(key) Sys.getenv(key, settings$value[match(key,settings$key)])
positive_integer <- function(key, minimum=1L) {
  x <- suppressWarnings(as.numeric(setting(key)))
  if (length(x)!=1L || !is.finite(x) || x!=floor(x) || x<minimum || x>.Machine$integer.max)
    stop('Invalid setting: ',key)
  as.integer(x)
}
cfg <- list(run_date=setting('CH1_RUN_DATE'), root=setting('CH1_MODEL_OUTPUT_DIR'),
            samples=positive_integer('CH1_SAMPLES'), thin=positive_integer('CH1_THIN'),
            transient=positive_integer('CH1_TRANSIENT',0L), chains=positive_integer('CH1_CHAINS',2L),
            nfolds=positive_integer('CH1_NFOLDS',2L), seed=positive_integer('CH1_SEED'),
            posterior_thin=positive_integer('CH1_POSTERIOR_THIN'), ngrid=positive_integer('CH1_NGRID',2L),
            overwrite=tolower(setting('CH1_OVERWRITE'))=='true')
if (!grepl('^[0-9]{4}-[0-9]{2}-[0-9]{2}$',cfg$run_date) || is.na(as.Date(cfg$run_date)))
  stop('CH1_RUN_DATE must be YYYY-MM-DD')
# Respect Slurm allocation. Outside Slurm use one worker unless explicitly requested.
cpu_setting <- Sys.getenv('CH1_NPARALLEL',Sys.getenv('SLURM_CPUS_PER_TASK','1'))
cfg$nparallel <- suppressWarnings(as.integer(cpu_setting))
if (is.na(cfg$nparallel) || cfg$nparallel<1L) stop('Invalid CPU allocation')
slurm_cpus <- suppressWarnings(as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','')))
if (!is.na(slurm_cpus)) cfg$nparallel <- min(cfg$nparallel,slurm_cpus)
set.seed(cfg$seed)

model_names <- function(selected=commandArgs(trailingOnly=TRUE)) {
  dirs <- list.dirs(cfg$root,recursive=FALSE,full.names=FALSE)
  chosen <- vapply(1:3,function(a) {
    pattern <- paste0('^',cfg$run_date,'(_[0-9]{2}-[0-9]{2}-[0-9]{2})?_Atlas',a,'_MinOccs5_CoverageGoodAverage$')
    hits <- dirs[grepl(pattern,dirs)]
    if (length(hits)!=1L) stop('Expected exactly one ',cfg$run_date,' model for Atlas ',a,
                               '; found ',length(hits),' in ',cfg$root)
    hits
  },character(1))
  if (length(selected)) {
    if (length(selected)!=1L || !selected %in% chosen) stop('Supply one model name from the configured three-atlas run')
    chosen <- selected
  }
  chosen
}
paths <- function(name) {
  base <- file.path(cfg$root,name)
  list(base=base, models=file.path(base,'Models'), results=file.path(base,'Results'),
       temp=file.path(base,'Models','Temp'), unfitted=file.path(base,'Models','Unfitted','unfitted_models.RData'),
       init=file.path(base,'Models','INITS',sprintf('HPC_INIT_samples_%04d_thin_%02d_chains_%d.rds',cfg$samples,cfg$thin,cfg$chains)),
       fitted=file.path(base,'Models','Fitted',sprintf('HPC_samples_%04d_thin_%02d_chains_%d.Rdata',cfg$samples,cfg$thin,cfg$chains)),
       folds=file.path(base,'Models','Temp',sprintf('temp_fold_info_samples_%04d_thin_%02d.rdata',cfg$samples,cfg$thin)),
       fit=file.path(base,'Models','Fitted',sprintf('MF_samples_%04d_thin_%02d_chains_%d_nfolds_%d.rdata',cfg$samples,cfg$thin,cfg$chains,cfg$nfolds)))
}
read_objects <- function(file,required) {
  if (!file.exists(file)) stop('Missing input: ',file)
  e <- new.env(parent=emptyenv()); load(file,envir=e)
  if (!all(required %in% ls(e))) stop('Missing objects in ',file,': ',paste(setdiff(required,ls(e)),collapse=', '))
  e
}
check_destination <- function(file) {
  if (file.exists(file) && !cfg$overwrite) stop('Output already exists: ',file,'; set CH1_OVERWRITE=true to replace it')
  dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE)
}
provenance <- function(directory,stage,inputs=character(),extra=character()) {
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  writeLines(c(paste('Stage:',stage),paste('Time:',Sys.time()),capture.output(str(cfg)),
               paste('Input:',inputs),if(length(inputs)) capture.output(tools::md5sum(inputs)),
               extra,capture.output(sessionInfo())),file.path(directory,paste0(stage,'_provenance.txt')))
  sources <- c(config_file,'scripts/2-hmsc/helpers/pipeline.R',file.path('scripts/2-hmsc',paste0(stage,'.R')))
  for (f in sources[file.exists(sources)]) {
    destination <- file.path(directory,paste0(stage,'_',basename(f)))
    if (normalizePath(f)!=normalizePath(destination,mustWork=FALSE)) file.copy(f,destination,overwrite=TRUE)
  }
}
validate_model <- function(m) {
  stopifnot(inherits(m,'Hmsc'),identical(rownames(m$Y),rownames(m$XData)),
            identical(colnames(m$Y),rownames(m$TrData)),all(is.finite(m$XScaled)),
            is.factor(m$TrData$Migration_a3_DOF),is.factor(m$TrData$foraging_guild_consensus),
            setequal(m$phyloTree$tip.label,colnames(m$Y)))
  invisible(m)
}
validate_fit <- function(m) {
  validate_model(m)
  if (length(m$postList)!=cfg$chains || any(lengths(m$postList)!=cfg$samples))
    stop('Fitted chains/draw counts disagree with configuration')
  if (m$thin!=cfg$thin || m$transient!=cfg$transient) stop('Fitted sampling settings disagree with configuration')
  invisible(m)
}
read_fit <- function(p) { m <- read_objects(p$fitted,'fitted_model')$fitted_model$posteriors; validate_fit(m); m }
read_hpc <- function(file,expected_draws=cfg$samples) {
  if (!file.exists(file)) stop('Missing HPC output: ',file)
  encoded <- readRDS(file)
  decoded <- jsonify::from_json(if(is.character(encoded)) encoded[[1]] else stop('Unexpected HPC serialization: ',file))
  if (!is.list(decoded) || length(decoded)<2L || length(decoded[[1]])!=expected_draws)
    stop('Incomplete HPC draws: ',file)
  for (i in seq_along(decoded[[1]])) {
    alpha <- decoded[[1]][[i]]$Alpha
    if (is.matrix(alpha)) decoded[[1]][[i]]$Alpha <- lapply(seq_len(nrow(alpha)),function(r) alpha[r,])
  }
  decoded
}
check_files <- function(actual,expected) {
  if (!setequal(normalizePath(actual,mustWork=FALSE),normalizePath(expected,mustWork=FALSE)))
    stop('Missing or unexpected chain/thread files. Missing: ',paste(setdiff(expected,actual),collapse=', '),
         '; unexpected: ',paste(setdiff(actual,expected),collapse=', '))
}
posterior_draws <- function(m) Hmsc::poolMcmcChains(m$postList,start=1,thin=cfg$posterior_thin)
prediction_array <- function(pred,ny,ns,ndraws) {
  if (length(pred)!=ndraws || !all(vapply(pred,function(x) identical(dim(x),c(as.integer(ny),as.integer(ns))),logical(1))))
    stop('Prediction dimensions/draw count do not match the receiving array')
  out <- array(unlist(pred,use.names=FALSE),c(ny,ns,ndraws))
  if (any(!is.finite(out))) stop('Non-finite predictions')
  out
}
# Filename sequences give each base graphics page its own PNG.
open_png <- function(prefix,width=2400,height=1800) {
  grDevices::png(paste0(prefix,'_%03d.png'),width=width,height=height,res=180)
}
summarise_draws <- function(draws) {
  x <- as.matrix(draws)
  data.frame(parameter=colnames(x),mean=colMeans(x),
             lower=apply(x,2,quantile,probs=.025),upper=apply(x,2,quantile,probs=.975),
             support_positive=colMeans(x>0),row.names=NULL)
}
make_fold_model <- function(hM,train) {
  # Keep the full-run transformation, as in Hmsc's computePredictedValues CV.
  m <- Hmsc::Hmsc(Y=hM$Y[train,,drop=FALSE],XData=hM$XData[train,,drop=FALSE],
                  XFormula=hM$XFormula,TrData=hM$TrData,TrFormula=hM$TrFormula,
                  phyloTree=hM$phyloTree,distr=hM$distr,
                  studyDesign=droplevels(hM$studyDesign[train,,drop=FALSE]),ranLevels=hM$rL)
  m <- Hmsc::setPriors(m,V0=hM$V0,f0=hM$f0,mGamma=hM$mGamma,UGamma=hM$UGamma,
                       aSigma=hM$aSigma,bSigma=hM$bSigma,rhopw=hM$rhowp)
  m$YScalePar <- hM$YScalePar; m$YScaled <- scale(m$Y,m$YScalePar[1,],m$YScalePar[2,])
  m$XInterceptInd <- hM$XInterceptInd; m$XScalePar <- hM$XScalePar
  m$XScaled <- scale(m$X,m$XScalePar[1,],m$XScalePar[2,])
  m$TrInterceptInd <- hM$TrInterceptInd; m$TrScalePar <- hM$TrScalePar
  m$TrScaled <- scale(m$Tr,m$TrScalePar[1,],m$TrScalePar[2,])
  stopifnot(isTRUE(all.equal(unname(m$XScaled),unname(hM$XScaled[train,,drop=FALSE]),check.attributes=FALSE)))
  m
}
