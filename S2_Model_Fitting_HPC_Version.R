#### Hmsc analyses on ####
#General cleaning of the workspace
#remove(list=ls())
print('loading libraries')

# 1. SET THE LIBPATH GLOBALLY FIRST
# This ensures any parallel workers created later inherit this path
.libPaths(c("~/Rlibs", .libPaths()))
require(Hmsc)
require(jsonify)
source(file.path("support_scripts", "project_paths.R"))

### Set up directories #### 

#If you are using RStudio this will set the working directory to exactly where the file is 
#setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path)))
pattern2match <- "2026-03-13_06-58-56"
  
matching_folders <- list.dirs('HmscOutputs', recursive = FALSE, full.names = F)
matching_folders <- find_model_folders(pattern = pattern2match)
folders2match <- matching_folders[1]
for(folders2match in matching_folders){
model_description = folders2match
localDir = sprintf("./HmscOutputs/%s",model_description)
ModelDir = file.path(localDir, "Models")

### Read in the unfitted models ####
load(file = file.path(ModelDir, "Unfitted/unfitted_models.RData"))

samples_list = 250
thin_list = 100
transient = 100000
nChains = 4
Lst = 1
verbose = 1

while(Lst <= length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  transient = transient
  
  filename = file.path(ModelDir,sprintf("INITS/HPC_INIT_samples_%.4d_thin_%.2d_chains_%.1d.rds",samples,thin,nChains))
  
  # Check if there is a mismatch in the study design
  m = sampleMcmc(models[[1]], samples = samples, thin=thin,
                 #adaptNf=rep(transient,models[[1]]$nr), 
                 transient = transient,
                 nChains = nChains,
                 verbose = verbose, 
                 engine = "HPC")
  
  saveRDS(to_json(m), filename)
  Lst = Lst + 1
}
}
