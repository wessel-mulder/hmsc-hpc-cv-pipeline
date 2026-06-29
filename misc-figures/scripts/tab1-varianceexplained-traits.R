rm(list = ls())

# GETTING STARTED  --------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
source(file.path("support_scripts", "figure_data_helpers.R"))

# get folders 
pacman::p_load(tidyverse, Hmsc)
pattern <- "2026-03-13"
folders <- figure_model_folders(pattern)
# get models 
mods <- load_hmsc_posteriors(folders)
# KABLE  ------------------------------------------------------------------

for(m1 in mods){

VP = computeVariancePartitioning(m1)
print(VP$R2T)

mpost <- convertToCodaObject(m1)
print(summary(mpost$Rho)$quant)
}

lapply

plotVariancePartitioning(mods[[1]], VP,)
