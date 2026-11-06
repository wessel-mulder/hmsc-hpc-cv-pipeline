#### Hmsc analyses on ####
#General cleaning of the workspace
remove(list=ls())
print('loading libraries')

# 1. SET THE LIBPATH GLOBALLY FIRST
# This ensures any parallel workers created later inherit this path
.libPaths(c("~/Rlibs", .libPaths()))

require(Hmsc)
require(dplyr)
require(ape)
#### Set up directories #### 
#If you are using RStudio this will set the working directory to exactly where the file is 
localDir = Sys.getenv("CH1_MODEL_OUTPUT_DIR", "./HmscOutputs")
dir.create(localDir, recursive = TRUE, showWarnings = FALSE)
dataDir = "./data-shared"
date <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

#### Model specs ####
guild_models <- c('All')
variable_models <- c('All')   # 10 S0 predictors: climate, heterogeneity, seven cover variables 
coverages <- c('GoodAverage')
atlas_models <- c(1,2,3)
min_occs2run <- c(5)
# Choose c("spatial", "none") to prepare both versions for each atlas.
# Use "spatial" alone for environment + space, or "none" for environment only.
# "none" removes all site latent factors; it does not substitute non-spatial ones.
spatial_models <- c("none")
spatial_n_factors <- 5L
stopifnot(length(spatial_models) > 0L, !anyDuplicated(spatial_models),
          all(spatial_models %in% c("spatial", "none")),
          length(spatial_n_factors) == 1L, is.finite(spatial_n_factors),
          spatial_n_factors >= 1L, spatial_n_factors == as.integer(spatial_n_factors))
guild2run <- guild_models[[1]]
variable2run <- variable_models[[1]]
atlas_model <- atlas_models[[1]]
min_occs <- min_occs2run[[1]]
coverage <- coverages[[1]]

for(coverage in coverages){
for(min_occs in min_occs2run){
for(atlas_model in atlas_models){
for(guild2run in guild_models){
for(variable2run in variable_models){
for(spatial_model in spatial_models){
# assign correct names 
guild <- guild2run
variables <- variable2run
atlas <- atlas_model
# make description
model_description = paste(as.character(date),
                          #guild,
                          #variables,
                          paste('Atlas',atlas,sep=''),
                          paste('MinOccs',min_occs,sep=''),
                          paste('Coverage',coverage,sep=''),
                          if (spatial_model == "spatial") {
                            paste0("EnvironmentSpace_nf", spatial_n_factors)
                          } else {
                            "EnvironmentOnly"
                          },
                          sep='_')

#### Edit the data ####
# If the data are saved as an RData file
# S0 retains >=5 occurrences in each atlas and ocean + freshwater <=75%,
# with the configured survey-coverage rule and complete environmental predictors.
if(coverage=='Good'){coveragerename='good'}
if(coverage=='GoodAverage'){coveragerename='good_and_average'}
input_file <- file.path(dataDir, sprintf("preprocessed_data_min_occs_is_%s_coverage_is_%s_forest_split.RData",
                                         min_occs, coveragerename))
if (!file.exists(input_file)) stop("S0 output not found: ", input_file, ". Run scripts/2-hmsc/S0_Data_Definitions.R first.")
load(input_file)
stopifnot(all(c("X", "Y", "Tr", "Design", "phy") %in% ls()),
          identical(rownames(X), rownames(Y)),
          identical(rownames(Y), rownames(Design)),
          identical(colnames(Y), rownames(Tr)))

# Subset environmental data 
X <- X %>%
  {
    if (variables == "LandusePercs") {
      select(., matches("^perc_"))
    } else if (variables == 'Climate') {
      select(., tmean_breeding,prec_breeding)
    } else if (variables == 'All') {
      select(., tmean_breeding,prec_breeding, hh,matches("^perc_"))
      
    }
  }
# Grab species from the guild 
if(guild != 'All'){
Y <- Y %>%
  select(any_of(rownames(Tr)[Tr$foraging_guild_consensus == guild]))
}
#### Prepare for model ####
# Define model formulas for environmental and trait data
XFormula <- as.formula(paste("~", paste(colnames(X), collapse = "+"), sep = " "))
TrFormula <- as.formula(paste("~", paste(colnames(Tr), collapse = "+"), sep = " "))

print(XFormula)
### GET PHYLOGENY
# phy is bundled by S0; reload above restores the full tree for each model.
if (!inherits(phy, "phylo")) stop("S0 phy must be an ape phylo tree")
missing_tips <- setdiff(colnames(Y), phy$tip.label)
if (length(missing_tips)) stop("Model species missing from bundled tree: ", paste(missing_tips, collapse=", "))
phy <- keep.tip(phy,colnames(Y))
stopifnot(setequal(phy$tip.label,colnames(Y)),length(phy$tip.label)==ncol(Y))

pd_matrix <- cophenetic.phylo(phy)
pd_matrix <- pd_matrix[sort(rownames(pd_matrix)), sort(colnames(pd_matrix))]
pd_matrix

# check if species-lists are identical
setdiff(rownames(pd_matrix),names(Y))
# stunning 

#### Grab training atlas data ####
pattern <- paste0("_(", paste(atlas, collapse = "|"), ")$")

Y_sub <- Y[rownames(Y)[grep(pattern, rownames(Y))],,drop=F]
X_sub <- X[rownames(X)[grep(pattern, rownames(X))],,drop=F]
Design_sub <- Design[rownames(Design)[grep(pattern, rownames(Design))],,drop=F]
Design_sub$atlas <- droplevels(Design_sub$atlas)
Design_sub$site <- droplevels(Design_sub$site)

# Are the rownames identical (including order)?
identical(rownames(Y_sub), rownames(X_sub))
identical(rownames(Y_sub), rownames(Design_sub))

setequal(rownames(Y_sub), rownames(X_sub))
setequal(rownames(Y_sub), rownames(Design_sub))

# set identical
X_sub      <- X_sub[rownames(Y_sub), , drop = FALSE]
Design_sub <- Design_sub[rownames(Y_sub), , drop = FALSE]

stopifnot(
  identical(rownames(Y_sub), rownames(X_sub)),
  identical(rownames(Y_sub), rownames(Design_sub))
)

setdiff(colnames(Y_sub),
        rownames(Tr))

# Only construct spatial random effects for the environment + space version.
# Retain site identifiers in studyDesign for both versions, but an empty
# ranLevels list means Hmsc has zero random levels in the environment-only model.
random_levels <- list()
if (spatial_model == "spatial") {
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
struc_space <- setPriors(struc_space, nfMin=spatial_n_factors, nfMax=spatial_n_factors)
nrow(struc_space$s)
random_levels <- list(site = struc_space)
}

# set factors of traits to determine level orders 
# mig
Tr$Migration_a3_DOF <- as.factor(Tr$Migration_a3_DOF)
most_common <- names(which.max(table(Tr$Migration_a3_DOF)))
Tr$Migration_a3_DOF <- relevel(Tr$Migration_a3_DOF, ref = most_common)
# guild
Tr$foraging_guild_consensus <- as.factor(Tr$foraging_guild_consensus)
most_common <- names(which.max(table(Tr$foraging_guild_consensus)))
Tr$foraging_guild_consensus <- relevel(Tr$foraging_guild_consensus, ref = most_common)

# scale 
vars_to_scale <- colnames(X_sub)
X_sub[vars_to_scale] <- scale(X_sub[vars_to_scale])

model = Hmsc(Y = Y_sub,
                     XData = X_sub,
                     XFormula = XFormula,
                     TrData = Tr,
                     TrFormula = TrFormula,
                     phyloTree = phy,
                     studyDesign = Design_sub[,c('site'),drop=F],
                     ranLevels = random_levels,
                     distr='probit')
stopifnot(model$nr == if (spatial_model == "spatial") 1L else 0L)
print(model)
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
  for(x in c("Fitted","Raw_HPC","Unfitted","INITS","Temp","Logs","Logs_CV","Sampled")){dir.create(file.path(save.dir,x))}
  if(!dir.exists(results.dir)){dir.create(results.dir)}
  if(!dir.exists(test.dir)){dir.create(test.dir)}
} else {
  cat(sprintf("\n%1$s\nWARNING\n%1$s\nThis model has been ran before, unfitted model file being over written\n",strrep("-",10)))
}
save(models, file = file.path(save.dir, "Unfitted", "unfitted_models.RData"))

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

}
}
}
}
}
}


