################################################
remove(list=ls())
.libPaths(c("~/Rlibs", .libPaths()))
require(Hmsc)
require(cli)

#Define the plotting function
plot_CV = function(type){
  #This sets the limits and values that change based on the metric, it outputs a
  #list where item 1 is abs(lower limit), item 2 is upper limit, item 3 is vline
  #value and item 4 is hline value
  limit = switch(type,
                 "RMSE"= c(rep(max(cMF[[type]],cMFCV[[type]]),2),0,0),
                 "AUC" = c(0,1,0.5,0.5),
                 "O.AUC" = c(0,1,0.5,0.5),
                 c(rep(1,2),0,0))
  plot(cMF[[type]],cMFCV[[type]],xlim=c(-limit[1],limit[2]),ylim=c(-limit[1],limit[2]),
       xlab = "explanatory power",
       ylab = "predictive power",
       main=sprintf("%s:\n%s thin = %i, samples = %i\nMF: mean(%1$s) = %.4f MFSCV: mean(%1$s) = %.4f",
                    type,modelnames,thin,samples,mean(cMF[[type]],na.rm=TRUE),
                    mean(cMFCV[[type]],na.rm=TRUE)))
  
  # 
  # main=paste0(modelnames,", thin = ",
  #             as.character(thin),
  #             ", samples = ",as.character(samples),
  #             ": Tjur R2.\n",
  #             "mean(MF) = ",as.character(mean(cMF[[type]],na.rm=TRUE)),
  #             ", mean(MFCV) = ",as.character(mean(cMFCV[[type]],na.rm=TRUE))))
  abline(0,1)
  abline(v=limit[3])
  abline(h=limit[4])
}

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


samples_list = c(250)
thin_list = c(10)
transient = 100000
nParallel = 10
nChains = 4
nfolds = 5

nst = length(thin_list)


for (Lst in nst:1) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  filename = file.path(ModelDir,sprintf("MF_samples_%.4d_thin_%.2d_chains_%.1d_nfolds_%.1d.rdata",
                                        samples, thin, nChains,nfolds))
  if(file.exists(filename)){
    cli_alert_success("Check good\nFile: {.file {filename}} exists")
    break} else{
      cli_alert_danger("Check BAD\nFile: {.file {filename}} does not exist\nRun computing model fit script for:\nThin: {.strong  {thin}} \tSamples: {.strong  {samples}}\tChains: {.strong  {nChains}}\tFolds: {.strong  {nfolds}}")
    }
}

if(file.exists(filename)){
  cli_progress_step("Loading File")
  load(filename)
  cli_progress_done()
  modelnames = models_description
  
  pdf(file = file.path(ResultDir,paste0("/",models_description,"_model_fit_nfolds_",nfolds,".pdf")))
  
  cMF = MF
  cMFCV = MFCV
  
  for(x in names(cMF)){
    cli_progress_step("Plotting {x}")
    plot_CV(x)
    cli_progress_done()
  }
  dev.off()
}
