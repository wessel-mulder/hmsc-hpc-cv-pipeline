#### Hmsc analyses on ####
#General cleaning of the workspace
remove(list=ls())
gc()

require(Hmsc)
#### Set up directories #### 
#If you are using RStudio this will set the working directory to exactly where the file is 
setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path)))
localDir = "./Hmsc Outputs"
dataDir = "./Data"
date <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

#### Model specs ####
guild_models <- c('Plovers','Woodpeckers')
variable_models <- c('Climate','LandusePercs')
atlas_models <- 3
for(guild2run in guild_models){
for(variable2run in variable_models){
# assign correct names 
guild <- guild2run
variables <- variable2run
atlas <- atlas_models
# make description
model_description = paste(as.character(date),
                          guild,
                          variables,
                          paste('Atlas',atlas,sep=''),
                          sep='_')

#### Edit the data ####
# If the data are saved as an RData file
# These files contain the preprocessed data values. Filtered for >5 occurrences in each atlas and 
# at least 25% of the atlas grid cell is land 
load(file.path(dataDir,"preprocessed_data.RData"))

# Subset environmental data 
X <- X %>%
  {
    if (variables == "LandusePercs") {
      select(., matches("^perc_"))
    } else if (variables == 'Climate') {
      select(., tmean_year,prec_year)
    } else {
      .
    }
  }
# Grab species from the guild 
Y <- Y %>%
  select(any_of(rownames(Tr)[Tr$foraging_guild_consensus == guild]))

#### Prepare for model ####
# Define model formulas for environmental and trait data
XFormula <- as.formula(paste("~", paste(colnames(X), collapse = "+"), sep = " "))
#TrFormula <- as.formula(paste("~", paste(colnames(Tr), collapse = "+"), sep = " "))

# get random effects for space
proj_xycoords_unique <- distinct(data.frame(X = Design$lon,
                                            Y = Design$lat))
rownames(proj_xycoords_unique) <- unique(Design$site)
struc_space <- HmscRandomLevel(sData = proj_xycoords_unique, sMethod = "Full")
struc_space <- setPriors(struc_space,nfMin=5,nfMax=5) # set priors to limit latent factors

#### Grab training atlas data ####
pattern <- paste0("_(", paste(atlas, collapse = "|"), ")$")

Y_sub <- Y[rownames(Y)[grep(pattern, rownames(Y))],,drop=F]
X_sub <- X[rownames(X)[grep(pattern, rownames(X))],,drop=F]
Design_sub <- Design[rownames(Design)[grep(pattern, rownames(Design))],,drop=F]
Design_sub$atlas <- droplevels(Design_sub$atlas)

model = Hmsc(Y = Y_sub,
                     XData = X_sub,
                     XFormula = XFormula,
                     # TrData = Tr,
                     # TrFormula = TrFormula,
                     # phyloTree = phy,
                     studyDesign = Design_sub[,c('site'),drop=F],
                     ranLevels = list('site'=struc_space),
                     distr='probit')

#### Save model ####
#These lists are left over from the orginal pipeline which had multiple models
#defined at the same time
models = list(model)
names(models) = c(model_description)

#Check if the model with this description has been created before, if not create
#both the model and results directories and all sub folders
ModelDir = file.path(localDir, sprintf("%s",model_description))
save.dir = file.path(ModelDir, "Models")
results.dir = file.path(ModelDir, "Results")
test.dir = file.path(ModelDir, "Tests")

if(!dir.exists(save.dir)){
  dir.create(ModelDir)
  dir.create(save.dir)
  for(x in c("Fitted","Raw_HPC","Unfitted","INITS","Temp")){dir.create(file.path(save.dir,x))}
  if(!dir.exists(results.dir)){dir.create(results.dir)}
  if(!dir.exists(test.dir)){dir.create(test.dir)}
} else {
  cat(sprintf("\n%1$s\nWARNING\n%1$s\nThis model has been ran before, unfitted model file being over written\n",strrep("-",10)))
}
save(models, file = file.path(save.dir, "/Unfitted/unfitted_models.RData"))

#### Store a copy of the script #### 
# --- Backup this script in Unfitted folder ---
args <- commandArgs(trailingOnly = FALSE)
script_path <- NULL

# Detect current script path
file_arg <- grep("--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  # Works for Rscript or batch jobs
  script_path <- normalizePath(sub("--file=", "", file_arg))
} else if (requireNamespace("rstudioapi", quietly = TRUE) &&
           rstudioapi::isAvailable()) {
  # Fallback for interactive runs in RStudio
  script_path <- normalizePath(rstudioapi::getSourceEditorContext()$path)
} else if (!is.null(sys.frames()) && !is.null(sys.frames()[[1]]$ofile)) {
  # Another fallback
  script_path <- normalizePath(sys.frames()[[1]]$ofile)
}

if (!is.null(script_path) && file.exists(script_path)) {
  backup_file <- file.path(
    save.dir, "Unfitted",
    paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S_"), basename(script_path))
  )
  file.copy(script_path, backup_file, overwrite = TRUE)
  message("Script copied to: ", backup_file)
} else {
  warning("Could not determine script path for backup")
}

#### Grab data from the other atlases for testing later ####
tests <- setdiff(c(1,2,3),atlas)
test<-1
# Create a list of testing datasets
testing_list <- lapply(tests, function(test) {
  
  # pattern for rownames ending with _1, _2, etc.
  pattern <- paste0("_(", paste(test, collapse = "|"), ")$")
  
  # subset Y, X, Design
  Y_sub <- Y[grep(pattern, rownames(Y)), , drop = FALSE]
  X_sub <- X[grep(pattern, rownames(X)), , drop = FALSE]
  Design_sub <- Design[grep(pattern, rownames(Design)), , drop = FALSE]
  
  # drop unused factor levels in atlas
  Design_sub$atlas <- droplevels(Design_sub$atlas)
  
  # return a list with your 3 datasets
  list(
    Y = Y_sub,
    X = X_sub,
    Design = Design_sub
  )
})
# Name each element by atlas
names(testing_list) <- paste0("atlas_", tests)
save(testing_list,
     file = file.path(test.dir,'test_atlases.RData'))

}
}
