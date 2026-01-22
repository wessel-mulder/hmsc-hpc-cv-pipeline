remove(list=ls())
.libPaths(c("~/Rlibs", .libPaths()))

require(Hmsc)
require(ggplot2)
require(cli)
set.seed(369)
### Set up directories #### Because I run this on two difference computers this

### Set up directories #### 
#If you are using RStudio this will set the working directory to exactly where the file is 

guild <- 'Woodpeckers'
env_var <- 'LandusePercs'
models_description = sprintf("2026-01-20_12-40-41_%s_%s_Atlas3",guild,env_var)

getwd()
localDir = "./HmscOutputs"
ModelDir = file.path(localDir, sprintf("%s/Models/Fitted",models_description))
TempDir = file.path(localDir,sprintf("%s/Models/Temp",models_description))
ResultDir = file.path(localDir, sprintf("%s/Results",models_description))
UnfittedDir = file.path(localDir, sprintf("%s/Models/Unfitted",models_description))

samples_list = c(250)
thin_list = c(10)
transient = 100000
nParallel = 10
nChains = 4
nfolds = 5

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
  if(length(covariates)>0){
    #Note that I use different file names for the R fitted and HPC fitted models
    #just to keep track
    outfile = file.path(ResultDir,sprintf("Preds/Preds_%s_HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",models_description, m$samples, m$thin, nChains))
    outfiletest = file.path(TestDir,sprintf("Preds/AtlasPreds_%s_HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",models_description, m$samples, m$thin, nChains))

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
        ptm_tot = proc.time()
        ptm = ptm_tot
        Gradient = constructGradient(m,focalVariable = covariate, ngrid=30)
        computational.time = proc.time() - ptm
        
        cli_progress_step("Starting to construction Gradient2:")
        Gradient2 = constructGradient(m,focalVariable = covariate,non.focalVariables = 1, ngrid=30)
        computational.time = proc.time() - ptm
        
        cli_progress_step("Making predictions based on Gradient 1:")
        predY = predict(m, Gradient=Gradient, expected = TRUE)
        computational.time = proc.time() - ptm
        cli_progress_step("Making predictions based on Gradient2:")
        
        ptm = proc.time()
        predY2 = predict(m, Gradient=Gradient2, expected = TRUE)
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
}
