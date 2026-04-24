print('loading packages')
.libPaths(c("~/Rlibs", .libPaths()))
require(Hmsc)
require(cli)
source(file.path("support_scripts", "project_paths.R"))
source(file.path("support_scripts", "plot_helpers.R"))

### Set up directories #### 
#If you are using RStudio this will set the working directory to exactly where the file is 
### Set up directories ####
pattern2match <- "2026-02-10"
  
matching_folders <- list.dirs('HmscOutputs', recursive = FALSE, full.names = F)
matching_folders <- find_model_folders(pattern = pattern2match)

for(folders2match in matching_folders){
models_description = folders2match

getwd()
localDir = "./HmscOutputs"
ModelDir = file.path(localDir, sprintf("%s/Models/Fitted",models_description))
TempDir = file.path(localDir,sprintf("%s/Models/Temp",models_description))
ResultDir = file.path(localDir, sprintf("%s/Results",models_description))


samples_list = c(250)
thin_list = c(100)
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
    plot_model_fit_cv(x, cMF, cMFCV, modelnames, thin, samples)
    cli_progress_done()
  }
  dev.off()
}
}
