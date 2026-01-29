remove(list=ls())
.libPaths(c("~/Rlibs", .libPaths()))

require(Hmsc)
require(ggplot2)
require(parallel)
require(cli)
set.seed(369)
### Set up directories #### Because I run this on two difference computers this



### MY FLAGS 
ngrid <- 3
skip_covariates <- TRUE # skip the first part of the script which constructs gradients & predicts over them 
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
pattern2match <- "2026-01-27"
  
matching_folders <- list.dirs('HmscOutputs', recursive = FALSE, full.names = F)
matching_folders <- matching_folders[grepl(pattern2match, basename(matching_folders))]
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
trait.list = NULL
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
  m$postList = m$postList[1:2]
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
    if(file.exists(file.path(outfile))){
      cli_alert_success("Predictions already calculated")
      load(outfile)
    } else {
      Preds = vector("list", length(covariates))
      for(k in 1:(length(covariates))){
        covariate = covariates[[k]]
        cli_h2("Calculating predictions for {covariate}")
        cli_progress_step("Starting to construction Gradient:")
        sprintf('Making predictions for %s ngrids',ngrid)
        ptm_tot = proc.time()
        ptm = ptm_tot
        Gradient = constructGradient(m,focalVariable = covariate, ngrid=ngrid)
        computational.time = proc.time() - ptm
        
        cli_progress_step("Starting to construction Gradient2:")
        Gradient2 = constructGradient(m,focalVariable = covariate,non.focalVariables = 1, ngrid=ngrid)
        computational.time = proc.time() - ptm
        
        cli_progress_step("Making predictions based on Gradient 1:")
        predY = predict(m, Gradient=Gradient, expected = TRUE,nParallel = nParallel,useSocket=F)
        computational.time = proc.time() - ptm
        cli_progress_step("Making predictions based on Gradient2:")
        ptm = proc.time()
        predY2 = predict(m, Gradient=Gradient2, expected = TRUE,nParallel = nParallel,useSocket=F)
        cli_process_done()
        computational.time = proc.time() - ptm
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
        
        #rownames(X_sub) <- gsub("_2", "_3", rownames(X_sub))
        #rownames(Design_sub) <- gsub("_2", "_3", rownames(Design_sub))
        
        # get ranlevels 
  
        proj_xycoords_unique <- unique(
          data.frame(
            site = Design_sub$site,
            X    = Design_sub$lon,
            Y    = Design_sub$lat,
            stringsAsFactors = FALSE
          )
        )
        rownames(proj_xycoords_unique) <- proj_xycoords_unique$site
        proj_xycoords_unique$site <- NULL
        
        struc_space <- HmscRandomLevel(sData = proj_xycoords_unique, sMethod = "Full")
        struc_space <- setPriors(struc_space,nfMin=5,nfMax=5) # set priors to limit latent factors
        
        ### GET PREDICTIONS 
        preds <- predict(m,
                XData = X_sub,
                studyDesign = Design_sub,
                ranLevels = list('site'=struc_space),
                nParallel = nParallel,useSocket=F)
        EpredY = Reduce('+', preds)/length(preds) # convert to 2d 
        ### get fid 
        str(Design_sub)
        predArray = abind::abind(preds, along=3)
        hM_test <- Hmsc(Y = Y_sub, 
                XData = X_sub, 
                studyDesign = Design_sub[,'site',drop=F],
                ranLevels = list('site'=struc_space),
                distr = m$distr) # Ensure distribution matches
        test_fit = evaluateModelFit(hM_test, predArray)
        atlas_data$Ypred = predArray
        atlas_data$fit_test = test_fit
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