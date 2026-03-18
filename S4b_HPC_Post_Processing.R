set.seed(369)
print('loading libraries')

# This ensures any parallel workers created later inherit this path
.libPaths(c("~/Rlibs", .libPaths()))

#### ECOSPAT BOYCE FUNCTION
#### Calculating Boyce index as in Hirzel et al. 2006
# fit: A vector or a SpatRaster containing the predicted suitability values 
# obs: A vector containing the predicted suitability values or xy-coordinates (if fit is a Raster-Layer) of the validation points (presence records)
# nclass : number of classes or vector with classes threshold. If nclass=0, Boyce index is calculated with a moving window (see next parameters)
# windows.w : width of the moving window (by default 1/10 of the suitability range)
# res : resolution of the moving window (by default 100 focals)
# PEplot : if True, plot the predicted to expected ratio along the suitability class
# rm.duplicate : if TRUE, the correlation exclude successive duplicated values
# method : correlation method used to compute the boyce index


ecospat.boyce <- function(fit, obs, nclass = 0, window.w = "default", res = 100, 
                          PEplot = TRUE, rm.duplicate = TRUE, method = 'spearman') {
  
  #### internal function calculating predicted-to-expected ratio for each class-interval
  boycei <- function(interval, obs, fit) {
    pi <- sum(as.numeric(obs >= interval[1] & obs <= interval[2])) / length(obs)
    ei <- sum(as.numeric(fit >= interval[1] & fit <= interval[2])) / length(fit)
    return(round(pi/ei,10))
  }

  
  mini <- min(fit,obs)
  maxi <- max(fit,obs)
  
  if(length(nclass)==1){
    if (nclass == 0) { #moving window
      if (window.w == "default") {window.w <- (max(fit) - min(fit))/10}
      vec.mov <- seq(from = mini, to = maxi - window.w, by = (maxi - mini - window.w)/res)
      vec.mov[res + 1] <- vec.mov[res + 1] + 1  #Trick to avoid error with closed interval in R
      interval <- cbind(vec.mov, vec.mov + window.w)
    } else{ #window based on nb of class
      vec.mov <- seq(from = mini, to = maxi, by = (maxi - mini)/nclass)
      interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
  } else{ #user defined window
    vec.mov <- c(mini, sort(nclass[!nclass>maxi|nclass<mini]))
    interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
  }
  
  f <- apply(interval, 1, boycei, obs, fit)
  to.keep <- which(f != "NaN")  # index to keep no NaN data
  f <- f[to.keep]
  if (length(f) < 2) {
    b <- NA  #at least two points are necessary to draw a correlation
  } else {
    r<-1:length(f)
    if(rm.duplicate == TRUE){
      r <- c(1:length(f))[f != c( f[-1],TRUE)]  #index to remove successive duplicates
    }
    b <- cor(f[r], vec.mov[to.keep][r], method = method)  # calculation of the correlation (i.e. Boyce index) after removing successive duplicated values
  }
  HS <- apply(interval, 1, sum)/2  # mean habitat suitability in the moving window
  if(length(nclass)==1 & nclass == 0) {
    HS[length(HS)] <- HS[length(HS)] - 1  #Correction of the 'trick' to deal with closed interval
  }
  HS <- HS[to.keep]  #exclude the NaN
  if (PEplot == TRUE) {
    plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio", col = "grey", cex = 0.75)
    points(HS[r], f[r], pch = 19, cex = 0.75)
  }
  return(list(F.ratio = f, cor = round(b, 3), HS = HS))
}

# Load all required libraries
library(jsonify)
library(colorspace)
library(RColorBrewer)
library(farver)
library(scales)
library(Hmsc)
library(cli)
library(vioplot)
library(parallel)
#install.packages(c("vegan", "ecodist", "dismo"),repos = "https://cloud.r-project.org")
#install.packages("ecospat",repos = "https://cloud.r-project.org")
#library(ecospat)

# Force parallel processes to use this path
options(parallel.lib.paths = "~/Rlibs")

c.Hmsc = getS3method("c","Hmsc")

# Get arguments from command line
args <- commandArgs(trailingOnly = TRUE)
# Default value if no argument is provided, otherwise use the first argument
thin2 <- if(length(args) > 0) as.numeric(args[1]) else 50
models_description <- args[2]

# ### Set up directories ####
# pattern2match <- "2026-01-27"
  
# matching_folders <- list.dirs('HmscOutputs', recursive = FALSE, full.names = F)
# matching_folders <- matching_folders[grepl(pattern2match, basename(matching_folders))]

# for(folders2match in matching_folders){
#models_description = folders2match

getwd()
localDir = "./HmscOutputs"
ModelDir = file.path(localDir, sprintf("%s/Models/Fitted",models_description))
TempDir = file.path(localDir,sprintf("%s/Models/Temp",models_description))




samples_list = c(250)
thin_list = c(100)
transient = 100000
nParallel = detectCores() - 1
print(nParallel)
nChains = 4
nfolds = 10

Lst = 1
print('starting')
for(Lst in length(samples_list):1){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  transient = transient
  filename.in = file.path(TempDir,sprintf("temp_fold_info_samples_%.4d_thin_%.2d.rdata", samples, thin))
  filename.out = file.path(ModelDir,sprintf("MF_samples_%.4d_thin_%.2d_chains_%.1d_nfolds_%.1d.rdata",
                                            samples, thin, nChains,nfolds))
  if(file.exists(filename.in)){
    ptm = proc.time()
    message("Reading in Fold information")
    load(file = filename.in)
    postN <- Reduce(sum, lapply(hM$postList, length))
    predArray <- array(NA, c(hM$ny, hM$ns, postN))
    mods = vector("list", threads)
    # for(i in 1:threads){
    #   temp_fitted = file.path(TempDir,sprintf("Sampled_HPC_samples_%.4d_thin_%.2d_thread_%.1d.rds", samples, thin,i))
    #   temp = from_json(readRDS(file = temp_fitted)[[1]])
    #   if(is.matrix(temp[[1]][[1]]$Alpha)){
    #     cat("\tAlpha is a matrix\nFixing alpha issue\n")
    #     #pb = txtProgressBar(min = 0, max = samples, initial = 0)
    #     for(z in 1:samples){
    #       temp_Alpha_Mat = temp[[1]][[z]]$Alpha
    #       temp[[1]][[z]]$Alpha = lapply(seq_len(nrow(temp_Alpha_Mat)), function(p) temp_Alpha_Mat[p,])
    #       #setTxtProgressBar(pb,i)
    #     }
    #   } else {
    #     cat("\tAlpha is not a matrix\nNo fix required\n")
    #   }
    #   k = idfold[i]
    #   m = importPosteriorFromHPC(hM1[[k]], temp[1], samples, thin, transient, alignPost = TRUE)
    #   message("finished thread ", i, "/", threads)
    #   attr(m, "fold") <- k
    #   mods[[i]] = m
    # }
    mods <- mclapply(1:threads, function(i) {
      
      temp_fitted <- file.path(TempDir, sprintf("Sampled_HPC_samples_%.4d_thin_%.2d_thread_%.1d.rds", samples, thin, i))
      
      # Loading the RDS
      temp <- from_json(readRDS(file = temp_fitted)[[1]])
      
      # Alpha fix logic
      if(is.matrix(temp[[1]][[1]]$Alpha)){
        for(z in 1:samples){
          temp_Alpha_Mat <- temp[[1]][[z]]$Alpha
          temp[[1]][[z]]$Alpha <- lapply(seq_len(nrow(temp_Alpha_Mat)), function(p) temp_Alpha_Mat[p,])
        }
      }
      
      k <- idfold[i]
      
      # THIS IS THE HEAVY PART (CPU and RAM intensive)
      m <- importPosteriorFromHPC(hM1[[k]], temp[1], samples, thin, transient, alignPost = TRUE)
      
      attr(m, "fold") <- k
      
      # Clean up 'temp' inside the function to free RAM immediately
      rm(temp)
      
      return(m)
      
    }, mc.cores = nParallel)

    rm(hM1)
    print(mods[[1]])
    #Combine predictions: this is still a loop
    idfold <- sapply(mods, attr, which = "fold")
    parts <- sort(unique(partition))
    print(partition)
    #This might be used in the future if I package up this code better, but for now its always set to NULL
    partition.sp = NULL
    Yc = NULL
    expected=TRUE

    for (p in parts) {
      message("predictions for partition ", p)
      print(partition)
      print(p)
      str(idfold)
      str(p)
      str(partition)
      val <- partition == p
      m <- do.call(c.Hmsc, mods[which(idfold == p)])
      m <- alignPosterior(m)
      postList <- poolMcmcChains(m$postList, start=1, thin=thin2)
      dfPi <- droplevels(hM$dfPi[val,, drop=FALSE])
      Xval <- if (is.matrix(hM$X)){
        hM$X[val, , drop=FALSE]
      } else{
        lapply(hM$X, function(a) a[val, , drop=FALSE])
      }
      message('starting to predict')
      pred1 <- if (is.null(partition.sp)) {
        predict(m, post=postList, X = Xval,
                XRRR = hM$XRRR[val,, drop=FALSE],
                Yc = Yc[val,, drop=FALSE], studyDesign = dfPi,
                mcmcStep = mcmcStep, expected = expected, nParallel = nParallel, useSocket = FALSE)
      } else {
        getSpeciesFoldPrediction(hM, m, val, postList, dfPi,
                                 partition.sp = partition.sp,
                                 mcmcStep = mcmcStep,
                                 expected = expected,
                                 nParallel = nParallel,
                                 useSocket = useSocket)
      }
      cat(sprintf("current mermory usage Ncels: %.1f MB Ncels: %.1f MB\n", gc()[3], gc()[4]))
      cat("High memory use section, writing predictions to array\n")
      predArray[val,,] <- simplify2array(pred1)
      cat(sprintf("current mermory usage Ncels: %.1f MB Ncels: %.1f MB\n", gc()[3], gc()[4]))
      cat("Cleaning up memory\n")
      rm(pred1,val,m,postList,Xval,dfPi)
      cat("Memory clean complete\n")
      cat(sprintf("current mermory usage Ncels: %.1f MB Ncels: %.1f MB\n", gc()[3], gc()[4]))
    }
    ### MF FOR NORMAL 
    preds = computePredictedValues(hM)
    cat("Calculating MF\n")
    MF = evaluateModelFit(hM, predY=preds)
    # --- NEW: Add Boyce Index to MF ---
    cat("Calculating Boyce Index for MF\n")
    # We calculate the mean prediction across posterior samples
    mean_preds = apply(preds, c(1,2), mean) 
    # Loop through each species to get Boyce
    MF$Boyce = sapply(1:ncol(hM$Y), function(i) {
      ecospat.boyce(fit = mean_preds[,i], 
                    obs = mean_preds[hM$Y[,i] == 1, i], 
                    nclass = 0, plot.main = FALSE)$Cor
    })
    rm(pred)
    
    ### MF FOR CV 
    cat("Calculating MFCV\n")
    MFCV = evaluateModelFit(hM, predY=predArray)
    cat("Calculating Boyce Index for MFCV\n")
    mean_predCV = apply(predArray, c(1,2), mean)
    MFCV$Boyce = sapply(1:ncol(hM$Y), function(i) {
      # We use the actual observations from hM$Y
      ecospat.boyce(fit = mean_predCV[,i], 
                    obs = mean_predCV[hM$Y[,i] == 1, i], 
                    nclass = 0, plot.main = FALSE)$Cor
    })
    
    rm(predArray)
    cat("Calculating WAIC\n")
    WAIC = computeWAIC(hM)
    computational.time = proc.time() - ptm
    cat("Time taken:", computational.time[3],"s \n\n")
    cat(sprintf("max mermory usage\n\tNcels: %.2f MB\n\tNcels: %.2f MB\n\tTotal: %.2f MB\n", gc()[11], gc()[12], gc()[11] + gc()[12]))
    save(MF,MFCV,WAIC,file = filename.out)
    break
  } else {
    message("Could not find file:", filename.in)
  }
}
