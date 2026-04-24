remove(list=ls())
.libPaths(c("~/Rlibs", .libPaths()))

require(Hmsc)
require(ggplot2)
require(parallel)
require(cli)
source(file.path("support_scripts", "project_paths.R"))
set.seed(369)
### Set up directories #### Because I run this on two difference computers this

### MY FLAGS 
ngrid <- 30
thin2 <- 50
overwrite_gradients <- TRUE
skip_covariates <- FALSE # skip the first part of the script which constructs gradients & predicts over them 
atlas_subset <- NULL # set to NULL for all atlas rows, or a number for a subset  
overwrite_prediction <- TRUE # overwrite atlas predictions if they've already been made 

# # Get arguments from command line
# args <- commandArgs(trailingOnly = TRUE)

# # Set a thinnning factor
# # Default value if no argument is provided, otherwise use the first argument

### Set up directories #### 
### Set up directories #### 
#If you are using RStudio this will set the working directory to exactly where the file is 
### Set up directories ####
pattern2match <- "2026-02-10"
  
matching_folders <- list.dirs('HmscOutputs', recursive = FALSE, full.names = F)
matching_folders <- find_model_folders(pattern = pattern2match)
folders2match <- matching_folders[1]
for(folders2match in matching_folders){
models_description = folders2match

getwd()
localDir = "./HmscOutputs"
ModelDir = file.path(localDir, sprintf("%s/Models/Fitted",models_description))
TempDir = file.path(localDir,sprintf("%s/Models/Temp",models_description))
UnfittedDir = file.path(localDir, sprintf("%s/Models/Unfitted",models_description))
ResultDir = file.path(localDir, sprintf("%s/Results",models_description))
TestDir = file.path(localDir, sprintf("%s/Tests",models_description))

samples_list = c(250)
thin_list = c(100)
transient = 100000
nChains = 4
nfolds = 5

nParallel = detectCores() - 1
print(nParallel)

nst = length(thin_list)
#Note that changing the species and traits only effect the graphical outputs and
#don't require now predictions to be calculated
#If you want to make predictions for certain species you can pass a list here,
#where the numbers relate to the position of the species in the Y matrix.
#Alternatively if you set this to NULL no species predictions will be plotted.
species.list = NULL
trait.list = 0
#Changing this does not always require recalculation of the predictions, unless
#you add a new factor since predictions are calculate for each covariant
#separately. However, the checker currently only checks for the presence of a
#file so it is a good idea to either delete the file or manually skip the check
#step to force recalculation.
env.list = NULL

#Create a Preds directory in the results directory if one does not already exist
if(!dir.exists(file.path(ResultDir,"Preds"))){
  dir.create(file.path(ResultDir,"Preds"))
}

load(file = file.path(UnfittedDir,"unfitted_models.RData"))
Lst <- nst[1]
for (Lst in nst:1) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  #Note that I use different file names for the R fitted and HPC fitted models
  #just to keep track
  
  filename = file.path(ModelDir,sprintf("HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata", samples, thin, nChains))
  #filename = file.path(ModelDir,sprintf("Fitted/FittedR_samples_%.4d_thin_%.2d_chains_%.1d.Rdata", samples, thin, nChains))
  if(file.exists(filename)){
    cli_alert_success("File {filename} exists")
    break} else{
      cli_alert_danger("File {filename} exists")
      cli_alert_info(cli_par("Run computing model fit script for:\n Thin: {thin} \t Samples: {samples} \t Chains: {nChains}"))
      cli_rule()
      }
}

if(file.exists(filename)){
  load(filename)
  #If you are using R fitted models you don't need to run the following two lines as the model is saved differently.
  m = fitted_model$posteriors
  rm(fitted_model)
  
  modelnames = models_description
  if(is.null(species.list)){
    species.list = list()
    species.list = 0
  }
  if(is.null(trait.list)){
    trait.list = list()
    trait.list = 0
  }
  if(is.null(env.list)){
    env.list = list()
    env.list = 0
  }
  
  pdf(file= file.path(ResultDir,paste0(models_description,"predictions.pdf")))
  if(all(env.list==0)){
    if(m$XFormula=="~."){
      covariates = colnames(m$XData)
    } else {
      covariates = all.vars(m$XFormula)
    }
  } else {
    covariates = env.list
  }
  
  #Change this to save the gradients, check if the file exists and if it does skip calculating them and move straight to plotting
  if(skip_covariates){
    cli_alert_info("Skipping covariate predictions as per user request")
    covariates = c()
  }
  if(length(covariates)>0){
    #Note that I use different file names for the R fitted and HPC fitted models
    #just to keep track
    outfile = file.path(ResultDir,sprintf("Preds/Preds_%s_HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",models_description, m$samples, m$thin, nChains))

    #outfile = file.path(ResultDir,sprintf("Preds/Preds_%s_R_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",model_description, m$samples, m$thin, nChains))
    cli_h1("Making predictions")
    if(file.exists(file.path(outfile)) & !overwrite_gradients){
      cli_alert_success("Predictions already calculated")
      load(outfile)
    } else {
      Preds = vector("list", length(covariates))
      k <- 1
      
      # thin2 posterior 
      postlist <- poolMcmcChains(m$postList,thin=thin2)
      ptm_tot = proc.time()
      ptm = ptm_tot
      for(k in 1:(length(covariates))){
        covariate = covariates[[k]]

        
        cli_h2("Calculating predictions for {covariate}")

        cli_progress_step("Starting to construction Gradient:")
        Gradient = constructGradient(m,focalVariable = covariate, ngrid=ngrid)

        cli_progress_step("Starting to construction Gradient2:")
        Gradient2 = constructGradient(m,focalVariable = covariate,non.focalVariables = 1,ngrid=ngrid)

        cli_progress_step("Making predictions based on Gradient 1:")
        predY = predict(m, Gradient=Gradient, expected = TRUE,nParallel = nParallel,useSocket=F,
                        post = postlist)

        cli_progress_step("Making predictions based on Gradient2:")
        predY2 = predict(m, Gradient=Gradient2, expected = TRUE,nParallel = nParallel,useSocket=F,
                         post = postlist)
        cli_process_done()
        
        Preds[[k]]$predY = predY 
        Preds[[k]]$predY2 = predY2 
        Preds[[k]]$Gradient = Gradient
        Preds[[k]]$Gradient2 = Gradient2
      }
      names(Preds) = covariates
      save(Preds, file = outfile)
      computational.time = proc.time() - ptm_tot
      cat(sprintf("Total Time taken: %.2f s \nCurrent time: %s\n\n", computational.time[3],format(Sys.time(), "%H:%M:%S")))
    }
    str(Preds)
    cli_h1("Plotting graphs")
    for(k in 1:(length(covariates))){
      par(mfrow=c(2,1))
      pl = plotGradient(m, Preds[[k]]$Gradient, pred=Preds[[k]]$predY, yshow = 0, measure="S", showData = TRUE, 
                        main = paste0(modelnames,": summed response (total effect)"))
      if(inherits(pl, "ggplot")){
        print(pl + labs(title=paste0(modelnames,": summed response (total effect)")))
      }
      pl = plotGradient(m, Preds[[k]]$Gradient2, pred=Preds[[k]]$predY2, yshow = 0, measure="S", showData = TRUE, 
                        main = paste0(modelnames,": summed response (marginal effect)"))
      if(inherits(pl, "ggplot")){
        print(pl + labs(title=paste0(modelnames,": summed response (marginal effect)")))
      }
      # only if species are supplied to list
      if(!species.list==0){
      for(l in 1:length(species.list)){
        par(mfrow=c(2,1))
        pl = plotGradient(m, Preds[[k]]$Gradient, pred=Preds[[k]]$predY, yshow = if(m$distr[1,1]==2){c(-0.1,1.1)}else{0}, measure="Y",index=species.list[l], showData = TRUE, 
                          main = paste0(modelnames,": example species (total effect)"))
        if(inherits(pl, "ggplot")){
          print(pl + labs(title=paste0(modelnames,": example species (total effect)")))
        }
        pl = plotGradient(m, Preds[[k]]$Gradient2, pred=Preds[[k]]$predY2, yshow = if(m$distr[1,1]==2){c(-0.1,1.1)}else{0}, measure="Y",index=species.list[l], showData = TRUE, 
                          main = paste0(modelnames,": example species (marginal effect)"))
        if(inherits(pl, "ggplot")){
          print(pl + labs(title=paste0(modelnames,": example species (marginal effect)")))
        }
      }}
      if(m$nt>1){
        traitSelection = 2:m$nt
        if(!all(trait.list==0)) traitSelection = trait.list
        for(l in traitSelection){
          par(mfrow=c(2,1))
          pl = plotGradient(m, Preds[[k]]$Gradient, pred=Preds[[k]]$predY, measure="T",index=l, showData = TRUE,yshow = 0,
                            main = paste0(modelnames,": community weighted mean trait (total effect)"))
          if(inherits(pl, "ggplot")){
            print(pl + labs(title=paste0(modelnames,": community weighted mean trait (total effect)")))
          }
          pl = plotGradient(m, Preds[[k]]$Gradient2, pred=Preds[[k]]$predY2, measure="T",index=l, showData = TRUE, yshow = 0,
                            main = paste0(modelnames,": community weighted mean trait (marginal effect)"))
          if(inherits(pl, "ggplot")){
            print(pl + labs(title=paste0(modelnames,": community weighted mean trait (marginal effect)")))
          }
        }
      }
    }
  }
  
  #This strange loop ensures that all output devices are turned off at the end.
  #Without it sometimes the pdf file that is create would crash on loading and
  #require the whole script to be reran
  while(!is.null(dev.list())){
    dev.off()
  }
  
  ### NOW CHECK FOR SOME ATLAS TESTS 
  filename = file.path(TestDir,'test_atlases.RData')
  #filename = file.path(ModelDir,sprintf("Fitted/FittedR_samples_%.4d_thin_%.2d_chains_%.1d.Rdata", samples, thin, nChains))
  if(file.exists(filename)){
    ptm = proc.time()

    outfiletest = file.path(TestDir,sprintf("Preds/AtlasPreds_%s_HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",models_description, m$samples, m$thin, nChains))
    
    load(filename)
    cli_h1("Running tests")
    if(file.exists(file.path(outfiletest)) & !overwrite_prediction){
      cli_alert_success("Predictions already calculated")
    } else {
      Atlases = testing_list
      atlas_test <- Atlases[[2]]
      tests<-lapply(Atlases,function(atlas_test){

        cli_h2(sprintf('Starting test for Atlas'))
        atlas_data <- atlas_test

        if(length(atlas_subset)>0){
          subset <- atlas_subset
        }else{
          subset <- nrow(atlas_data$X)
        }
        X_sub <- atlas_data$X[1:subset,]
        Design_sub <- atlas_data$Design[1:subset,]
        Y_sub <- atlas_data$Y[1:subset,]
        
        # 1. IDENTIFY TRAINING SITES
        # Get the list of sites the model was actually trained on
        train_sites <- unique(as.character(m$studyDesign$site))
        
        # 2. FILTER ATLAS TO MATCH MODEL
        # Only keep rows in the Atlas that correspond to sites in the model
        keep_rows <- as.character(Design_sub$site) %in% train_sites
        
        if(sum(keep_rows) == 0) {
          cli_alert_danger("No matching sites found between Model and Atlas! Check your site naming.")
          stop("Stopping to prevent crash.")
        }
        table(keep_rows)
        # Apply the filter
        X_sub_filtered <- X_sub[keep_rows, , drop=FALSE]
        Design_sub_filtered <- Design_sub[keep_rows, , drop=FALSE]
        Y_sub_filtered <- Y_sub[keep_rows, , drop=FALSE] # If you need Y for testing
        
        cli_alert_info(paste("Filtered prediction set from", nrow(Design_sub), "to", nrow(Design_sub_filtered), "sites matching the model."))
        
        # 1. First, make sure you are working with clean, unique data
        proj_xycoords_unique <- unique(
          data.frame(
            site = as.character(Design_sub_filtered$site), # Convert to character immediately
            X = Design_sub_filtered$lon,
            Y = Design_sub_filtered$lat,
            stringsAsFactors = FALSE
          )
        )
        
        # 2. Filter out any NAs that might have slipped in (Critical for the 'leading minor' error)
        proj_xycoords_unique <- proj_xycoords_unique[!is.na(proj_xycoords_unique$X), ]
        
        # 3. Set rownames using the character vector
        rownames(proj_xycoords_unique) <- proj_xycoords_unique$site
        
        # 4. Remove the site column so only X and Y remain for sData
        proj_xycoords_unique$site <- NULL
        
        # 5. NOW Ensure your Design_sub matches these sites exactly
        Design_sub_filtered$site <- factor(as.character(Design_sub_filtered$site), 
                                  levels = rownames(proj_xycoords_unique))
        
        struc_space <- HmscRandomLevel(sData = proj_xycoords_unique, sMethod = "Full")
        struc_space <- setPriors(struc_space,nfMin=5,nfMax=5) # set priors to limit latent factors
        
        ### GET SPATIAL PREDICTIONS 
        cli_h3('Spatial predictions')
        preds <- predict(m,
                XData = X_sub_filtered,
                studyDesign = Design_sub_filtered,
                ranLevels = list('site'=struc_space),
                nParallel = nParallel,useSocket=F)
        
        ### GET NON-SPATIAL PREDICTIONS 
        cli_h3('Non-spatial predictions')
        Design_dummy <- data.frame(site = factor(rep("dummy", nrow(X_sub_filtered))))
        struc_nospatial <- HmscRandomLevel(
          units = "dummy"
        )
        preds_nospatial <- predict(
          m,
          XData = X_sub_filtered,
          studyDesign = Design_dummy,
          ranLevels = list(site = struc_nospatial),
          expected = TRUE,
          nParallel = nParallel,
          useSocket = FALSE
        )
        

        # reduce dimentions
        EpredY = Reduce('+', preds)/length(preds) # convert to 2d 
        predArray = abind::abind(EpredY, along=3)
        # also for spatial
        Epreds_nospatialY = Reduce('+', preds_nospatial)/length(preds_nospatial)
        predArray_nonspatial = abind::abind(Epreds_nospatialY,along=3)
        
        
        ### get fake model for the fit  
        hM_test <- Hmsc(Y = Y_sub_filtered, 
                XData = X_sub_filtered, 
                studyDesign = Design_sub_filtered[,'site',drop=F],
                ranLevels = list('site'=struc_space),
                distr = m$distr) # Ensure distribution matches
        
        ### 
        test_fit = evaluateModelFit(hM_test, predArray)
        test_fit_nonspatial = evaluateModelFit(hM_test, predArray_nonspatial)
        # spatial 
        atlas_data$Ypred = predArray
        atlas_data$fit_test = test_fit
        #nonspatial
        atlas_data$Ypred_nonspatial = predArray_nonspatial
        atlas_data$fit_test_nonspatial = test_fit_nonspatial
        atlas_data
        
        })
      names(tests) <- names(Atlases)
      if(!dir.exists(file.path(TestDir,"Preds"))){
        dir.create(file.path(TestDir,"Preds"))
      }
      save(tests, file = outfiletest)
      computational.time = proc.time() - ptm
      cat("Time taken:", computational.time[3],"s \n\n")

      }
    }
  }

}

# 
# plot(tests$atlas_1$fit_test_nonspatial$AUC~tests$atlas_1$fit_test$AUC)
# plot(tests$atlas_2$fit_test_nonspatial$AUC~tests$atlas_2$fit_test$AUC)
# 
# plot(tests$atlas_1$fit_test_nonspatial$TjurR2~tests$atlas_1$fit_test$TjurR2)
# plot(tests$atlas_2$fit_test_nonspatial$TjurR2~tests$atlas_2$fit_test$TjurR2)
# 
# plot(tests$atlas_1$fit_test_nonspatial$TjurR2~tests$atlas_1$fit_test$RMSE)
# plot(tests$atlas_2$fit_test_nonspatial$TjurR2~tests$atlas_2$fit_test$RMSE)
# 
# # 1️⃣ AUC
# plot_comparison <- function(x, y, xlab = "Spatial", ylab = "Non-spatial", 
#                             main = NULL, lim = c(0,1),
#                             col_points = "black", col_line = "red") {
#   plot(y ~ x,
#        xlab = xlab, ylab = ylab,
#        xlim = lim, ylim = lim,
#        pch = 16, col = col_points,
#        main = main)
#   abline(a = 0, b = 1, col = col_line, lty = 2, lwd = 2)
# }
# 
# plot_comparison(
#   x = tests$atlas_1$fit_test$AUC,
#   y = tests$atlas_1$fit_test_nonspatial$AUC,
#   lim = c(0.3,1),
#   xlab = "Spatial AUC",
#   ylab = "Non-spatial AUC",
#   main = "Atlas 1: AUC"
# )
# 
# plot_comparison(
#   x = tests$atlas_2$fit_test$AUC,
#   y = tests$atlas_2$fit_test_nonspatial$AUC,
#   lim = c(0.3,1),
#   xlab = "Spatial AUC",
#   ylab = "Non-spatial AUC",
#   main = "Atlas 2: AUC"
# )
# 
# # 2️⃣ TjurR2
# plot_comparison(
#   x = tests$atlas_1$fit_test$TjurR2,
#   y = tests$atlas_1$fit_test_nonspatial$TjurR2,
#   lim = c(0,0.5),
#   xlab = "Spatial Tjur R2",
#   ylab = "Non-spatial Tjur R2",
#   main = "Atlas 1: Tjur R2"
# )
# 
# plot_comparison(
#   x = tests$atlas_2$fit_test$TjurR2,
#   y = tests$atlas_2$fit_test_nonspatial$TjurR2,
#   lim=c(0,0.5),
#   xlab = "Spatial Tjur R2",
#   ylab = "Non-spatial Tjur R2",
#   main = "Atlas 2: Tjur R2"
# )
# 
# # 3️⃣ TjurR2 vs RMSE
# plot_comparison(
#   x = tests$atlas_1$fit_test$RMSE,
#   y = tests$atlas_1$fit_test_nonspatial$RMSE,
#   xlab = "Spatial RMSE",
#   ylab = "Non-spatial RMSE",
#   main = "Atlas 1: RMSE vs RMSE"
# )
# 
# plot_comparison(
#   x = tests$atlas_2$fit_test$RMSE,
#   y = tests$atlas_2$fit_test_nonspatial$RMSE,
#   xlab = "Spatial RMSE",
#   ylab = "Non-spatial RMSE",
#   main = "Atlas 2: RMSE vs RMSE"
# )
