remove(list=ls())
set.seed(369)
require(Hmsc)

### Set up directories #### 

#If you are using RStudio this will set the working directory to exactly where the file is 
setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path)))
model_description = "Example1_R_x1_TRphy_Site"
localDir = sprintf("./Hmsc Outputs/%s",model_description)
ModelDir = file.path(localDir, "Models")

### Read in the unfitted models ####
load(file = file.path(ModelDir, "Unfitted/unfitted_models.RData"))

#Now we can start fitting the model here, we set up the number of model runs
samples_list = c(100, 250, 500)
thin_list = c(10, 20, 20)
nChains = 4
nParallel = 4
Lst = 1
verbose = 1

while(Lst <= length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  transient = ceiling(0.5*samples*thin)
  filename = file.path(ModelDir,sprintf("Fitted/FittedR_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",samples,thin,nChains))
  if(file.exists(filename)){
    cat("model had been fitted already\n\n")
  } else {
    ptm = proc.time()
    cat("Starting model fit\nRunning:\t",nParallel,"Parallel Process", ifelse(nParallel>1,"\nParallel processing, no outputs",
                                                                              paste0("\nShowing every:\t", verbose, " itteration")),
        "\nRunning:\n\t", transient, "transient\n\t",samples, "samples\n\t", thin,"sample thinning\n\t",
        transient+samples*thin,"total\n\n")
    m = sampleMcmc(models[[1]], thin = thin, samples = samples,
                   transient = transient, nChains = nChains,
                   verbose = verbose, nParallel = nParallel)
    computational.time = proc.time() - ptm
    cat("Time taken:", computational.time[3],"s \nCurrent time:",format(Sys.time(), "%H:%M:%S"),"\n\n")
    save(m,file=filename)
  }
  Lst = Lst + 1
}
