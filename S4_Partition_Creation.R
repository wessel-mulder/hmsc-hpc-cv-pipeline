cat("Script started\n"); flush.console()

remove(list=ls())
set.seed(369)
### Loading packages####
.libPaths(c("~/Rlibs", .libPaths()))

require(jsonify)
require(RColorBrewer)
require(farver)
require(scales)
require(Hmsc)
require(cli)
require(vioplot)
require(colorspace)


### Set up directories ####
pattern2match <- "2026-01-27"
  
matching_folders <- list.dirs('HmscOutputs', recursive = FALSE, full.names = F)
matching_folders <- matching_folders[grepl(pattern2match, basename(matching_folders))]

for(folders2match in matching_folders){
models_description = folders2match

getwd()
localDir = "./HmscOutputs"
ModelDir = file.path(localDir, sprintf("%s/Models/Fitted",models_description))
TempDir = file.path(localDir,sprintf("%s/Models/Temp",models_description))

samples_list = c(250)
thin_list = c(100)
nChains = 4
nfolds = 5

#Only run for the longest model run
for(Lst in length(samples_list):1){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  filename.in = file.path(ModelDir,sprintf("HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",
                                           samples, thin, nChains))
  message("Input file:", filename.in)
  filename.out = file.path(ModelDir,sprintf("MF_samples_%.4d_thin_%.2d_chains_%.1d_nfolds_%.1d.Rdata",
                                            samples, thin, nChains,nfolds))
  
  if(file.exists(filename.out)){
    cat("---------------------------\nthin = ",as.character(thin),
        "; samples = ",as.character(samples),
        "\nmodel fit had been computed already\n")
    break
  } else {
    if(file.exists(filename.in)){
      cat("---------------------------\nthin = ",
          as.character(thin),"; samples = ",as.character(samples),
          "\n",date(),"\n")
      load(file = filename.in) #models
      hM = fitted_model$posteriors
      if(is.null(hM$dfPi)){
        cat("No study Design defined, setting to NA's due to bug in function that requires it")
        hM$dfPi = data.frame(temp = factor(rep(NA, dim(hM$Y)[1])))
      }
      
      rm(fitted_model)
      
      partition = createPartition(hM, nfolds = nfolds)
      thin = 1
      start = 1
      ## STAGE 1: Basic housekeeping
      parts <- sort(unique(partition))
      nfolds <- length(parts)
      if (thin > 1 || start > 1)
        postN <- sum(sapply(hM$postList, function(z)
          length(seq(from=start, to=length(z), by=thin))))
      else
        postN <- Reduce(sum, lapply(hM$postList, length))
      ## output array
      predArray <- array(NA, c(hM$ny, hM$ns, postN))
      ## STEP 1: define new Hmsc model for each nfolds partition
      ## only implement for the simple case first
      if (is.list(hM$X))
        stop("not yet implemented for a list of model matrices")
      ## Pack setting training model into one function
      setHmsc <- function(k, hM) {
        train <- k != partition
        ## X can be a list of matrices
        XTrain <- if(is.matrix(hM$X))
          hM$X[train, , drop=FALSE]
        else if(is.list(hM$X))
          lapply(hM$X, function(a) a[train, , drop=FALSE])
        m <- Hmsc(Y = hM$Y[train, , drop=FALSE], X = XTrain,
                  XRRR = hM$XRRR[train, , drop=FALSE],
                  ncRRR = hM$ncRRR, XSelect = hM$XSelect,
                  distr = hM$distr,
                  studyDesign = droplevels(hM$dfPi[train,, drop=FALSE]),
                  Tr = hM$Tr, C = hM$C, ranLevels = hM$rL)
        ## old code calls here setPriors, but that does nothing as its
        ## result is not saved, and it is currently skipped: CHECK
        ## THIS!
        m$YScalePar <- hM$YScalePar
        m$YScaled <- scale(m$Y, m$YScalePar[1,], m$YScalePar[2,])
        m$XInterceptInd <- hM$XInterceptInd
        m$XScalePar <- hM$XScalePar
        m$XScaled <- scale(m$X, m$XScalePar[1,], m$XScalePar[2,])
        m$TrInterceptInd <- hM$TrInterceptInd
        m$TrScalePar <- hM$TrScalePar
        m$TrScaled <- scale(m$Tr, m$TrScalePar[1,], m$TrScalePar[2,])
        m
      }
      ## parallelism would only slow-down (due to start-off time)
      hM1 <- lapply(parts, function(k) setHmsc(k, hM))
      ## STEP 2: sample Hmsc models for each nfolds * nChains case
      chains <- length(hM$postList)
      threads <- nfolds * chains
      idfold <- rep(parts, each = chains)
      seeds <- sample.int(.Machine$integer.max, threads)
      ## to be called in parallel for each chain x fold, and
      ## therefore we set nChains=1, nParallel=1 within
      
      ##Loop over all the modes to make nfolds*nchains sets of training and testing data
      message(hM$transient)
      message(hM$adaptNf)
      for(i in 1:threads){
        set.seed(seeds[i])
        k <- idfold[i]
        message("Writing temp files for thread ", i, "/", threads)
        m <- sampleMcmc(hM1[[k]], samples = hM$samples, thin = hM$thin,
                        transient = hM$transient, adaptNf = hM$adaptNf,
                        initPar = hM$initPar, nChains = 1, nParallel = 1,
                        updater = hM$updater, verbose = hM$verbose,
                        alignPost = FALSE, engine = "HPC")
        #Call the gibs sampler on each thread, each thread is 1 chain so we can simple take the output and convert it to Rdata
        filename = file.path(TempDir,sprintf("temp_samples_%.4d_thin_%.2d_thread_%.1d.rds", hM$samples, hM$thin,i))
        saveRDS(to_json(m), file = filename)
      }
      save(idfold,threads, hM1, hM, partition, file = file.path(TempDir,sprintf("temp_fold_info_samples_%.4d_thin_%.2d.rdata",hM$samples, hM$thin)))
      break
    } else {
      message("No input files found for samples ",samples, " thin ", thin)
    }
  } 
}
}