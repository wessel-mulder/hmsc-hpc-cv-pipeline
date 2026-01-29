################################################
remove(list=ls())
.libPaths(c("~/Rlibs", .libPaths()))
require(Hmsc)
require(cli)

#Define the plotting function
type <- 'AUC'
plot_CV = function(type){
  #This sets the limits and values that change based on the metric, it outputs a
  #list where item 1 is abs(lower limit), item 2 is upper limit, item 3 is vline
  #value and item 4 is hline value
  limit = switch(type,
                 "RMSE"= c(rep(max(cMF[[type]],cMFCV[[type]]),2),0,0),
                 "AUC" = c(0,1,0.5,0.5),
                 "O.AUC" = c(0,1,0.5,0.5),
                 c(rep(1,2),0,0))
  pch = 19
  cex = 0.5
  par(mar=c(5,5,7.5,5))
  plot(cMF[[type]],cMFCV[[type]],xlim=c(-limit[1],limit[2]),ylim=c(-limit[1],limit[2]),pch=pch,cex=cex,
       xlab = "explanatory power",
       ylab = "predictive power",
       main=sprintf("%s:\n%s\nmean(explanatory) = %.2f, mean(cross-validated) = %.2f\n mean(test atlas1) = %.2f, mean(test atlas2) = %.2f",
                    type,modelnames,
                    mean(cMF[[type]],na.rm=TRUE),
                    mean(cMFCV[[type]],na.rm=TRUE),
                    mean(tMF[[1]][[type]],na.rm=T),
                    mean(tMF[[2]][[type]],na.rm=T)))

  points(cMF[[type]],tMF[[1]][[type]],col = 'red',pch=pch,cex=cex)
  points(cMF[[type]],tMF[[2]][[type]],col = 'blue',pch=pch,cex=cex)
  
  legend('topleft',legend = c('Cross validation on atlas 3','Test atlas 1','Test atlas 2'),
         fill = c('black','red','blue'))
  
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

pattern2match <- "2026-01-27"
  
matching_folders <- list.dirs('HmscOutputs', recursive = FALSE, full.names = F)
matching_folders <- matching_folders[grepl(pattern2match, basename(matching_folders))]

print('starting predictions')
for(folders2match in matching_folders){
models_description = folders2match

getwd()
localDir = "./HmscOutputs"
ModelDir = file.path(localDir, sprintf("%s/Models/Fitted",models_description))
TempDir = file.path(localDir,sprintf("%s/Models/Temp",models_description))
ResultDir = file.path(localDir, sprintf("%s/Results",models_description))
TestDir = file.path(localDir, sprintf("%s/Tests",models_description))

samples_list = c(250)
thin_list = c(100)
transient = 100000
nParallel = 10
nChains = 4
nfolds = 5

nst = length(thin_list)

# check if models exits 
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

# check for testing info 
testfilename = file.path(TestDir,sprintf("Preds/AtlasPreds_%s_HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",models_description, samples, thin, nChains))
if(file.exists(testfilename)){
  load(testfilename)
}

if(file.exists(filename)){
  #cli_progress_step("Loading File")
  load(filename)
  #cli_progress_done()
  modelnames = models_description

  pdf(file = file.path(TestDir,paste0("/",models_description,"_model_fit_nfolds_",nfolds,".pdf")))

  if(file.exists(testfilename)){
    extract_fit <- function(tests, atlas_name) {
      tests[[atlas_name]]$fit_test
    }
    tMF <- lapply(names(tests),function(test){
      extract_fit(tests,test)
     })
  }

  
  cMF = MF
  cMFCV = MFCV
  
  for(x in names(cMF)){
    cli_progress_step("Plotting {x}")
    plot_CV(x)
    cli_progress_done()
  }
  dev.off()
}
}