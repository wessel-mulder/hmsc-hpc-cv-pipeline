rm(list = ls())

# GETTING STARTED ---------------------------------------------------------
if (interactive() && Sys.getenv("RSTUDIO") == "1") {
  message("Running in RStudio")
  library(Hmsc)
  library(jsonify)
  library(knitr)
  library(corrplot)
  library(ape)
  library(dplyr)
  library(sp)
  library(sf)
  library(terra)
  library(tidyverse)
  input <- '.'
  python <- file.path(getwd(), 'hmsc-hpc-main',"hmsc-venv", "bin", "python3.11")
  flagInit = 0
  flagFitR = 1
} else {
  message("Running from terminal or non-interactive environment")
  library(RColorBrewer,lib="~/Rlibs")
  library(farver,lib="~/Rlibs")
  library(scales,lib="~/Rlibs")
  library(jsonify,lib="~/Rlibs")
  library(ape,lib="~/Rlibs")
  library(dplyr,lib="~/Rlibs")
  library(Hmsc,lib="~/Rlibs")
  library(dplyr,lib="~/Rlibs")
  library(withr,lib="~/Rlibs")
  library(sp,lib="~/Rlibs")
  library(terra,lib="~/Rlibs")
  library(sf,lib="~/Rlibs")
  library(tidyverse,lib="~/Rlibs")
  
  
  input <- '~/home/projects/hmsc-danishbirds'
  python <- '/maps/projects/cmec/people/bhr597/projects/hmsc-danishbirds/hmsc-venv'
  flagInit = 1
  flagFitR = 0
}

# check if accurately installed
summary(TD)

# LOADING DATA -----------------------------------------------------

#### ENVIRONMENT ####
# OCEAN THRESHOLD
grids_thresholds <- st_read(file.path(input,'Data/data/1_preprocessing/atlas-grids/grids-ocean-thresholds/grids_ocean_thresholds.shp'))
thresholds <- grids_thresholds$kvdrtkd[grids_thresholds$pct_lnd>=25]
# NAMES TO MATCH LANDUSE 
lulc_lookup <- c(
  "0"  = "ocean",
  "11" = "urban",
  "22" = "cropland",
  "33" = "pasture",
  "44" = "forest",
  "55" = "grass_shrub",
  "66" = "other",
  "77" = "water"
)
# GET ENVIRONMENTAL VALUES 
X <- read.csv(file.path(input,'Data/data/1_preprocessing/X_environmental/X_Environmental.csv'),row.names=1) %>%
  rename_with(
    .cols = matches("^LULC_"),
    .fn = function(col) {
      lulc_code <- stringr::str_extract(col, "\\d+")
      paste0("perc_", lulc_lookup[lulc_code])
    }
  ) %>%
  mutate(
    perc_fresh_saltwater = perc_water + perc_ocean
  ) %>%
  select(-perc_other, -perc_ocean, -perc_water) %>%
  mutate(
    dominant = factor(
      names(select(., starts_with("perc_")))[
        max.col(select(., starts_with("perc_")), ties.method = "first")
      ]
    )
  ) %>%
  # keep only grids passing land threshold
  filter(
    sub("_[123]$", "", rownames(.)) %in% thresholds
  ) %>%
  # drop rows with any NA
  drop_na()

#### OCCURRENCES ####  
Y <- read.csv(
  file.path(input, "Data/data/1_preprocessing/Y_occurrences/Y_occurrences.csv"),
  row.names = 1
) %>%
  filter(row.names(.) %in% row.names(X))

### filter based on richness per atlas 
effort_list <- c(
  '1' = .03,
  '2' = .01,
  '3' = 0
)
number='1'
for(number in c('1','2','3')){
quantile <- effort_list[[number]]
Y_sub <- Y[rownames(Y)[grep(paste0("_",number,"$"), rownames(Y))],,drop=F]
if(quantile > 0){
  print(sprintf('In atlas %s the quantile is %s percent',
          number,quantile))
  th = quantile(rowSums(Y_sub),quantile)
  tofilter = rownames(Y_sub)[rowSums(Y_sub)<th]
  Y <- Y[!rownames(Y) %in% tofilter,]
  print(sprintf("Removed %s sites from atlas %s",
                length(tofilter),number))
}
}

min_occs <- 5
for(number in c('1','2','3')){
  Y_sub <- Y[rownames(Y)[grep(paste0("_",number,"$"), rownames(Y))],,drop=F]
  if(any(colSums(Y_sub, na.rm =T)<min_occs)){
    print(paste0('In atlas ',number,' these species: '))
    print(names(which(colSums(Y_sub, na.rm =T)<min_occs)))
    tofilter <- names(which(colSums(Y_sub, na.rm =T)<min_occs))
    print('have less than 5 occurrences')
    Y <- Y[, !(colnames(Y) %in% tofilter)]
    print('and are now filtered ')
  }
}

# now refilter X again 
X <- X %>%
  filter(rownames(.) %in% rownames(Y))

##### TRAITS #####
Tr <- read.csv(
  file.path(input, "Data/data/1_preprocessing/Tr_aits/traits-guild_migration.csv"),
  row.names = 2
) %>%
  dplyr::select(2, 3) %>%
  .[colnames(Y), , drop = FALSE]

#### DESIGN ####
Design <- read.csv(
  file.path(input, "Data/data/1_preprocessing/design/studyDesign.csv"),
  row.names = 5
) %>%
  dplyr::filter(row.names(.) %in% row.names(X)) %>%
  .[sort(row.names(.)), , drop = FALSE] %>%
  dplyr::mutate(
    site  = factor(site),
    atlas = factor(atlas),
    year  = dplyr::case_when(
      atlas == "1" ~ 1971,
      atlas == "2" ~ 1992,
      atlas == "3" ~ 2014
    )
  )

# project cordinates 
xycoords <- Design %>%
  dplyr::select(lon, lat)
v <- terra::vect(xycoords, geom = c("lon", "lat"), crs = "EPSG:4326") %>%
  terra::project("EPSG:23032")
proj_xycoords <- terra::crds(v)

# merge back 
Design <- Design %>%
  dplyr::mutate(
    lon = proj_xycoords[, 1],
    lat = proj_xycoords[, 2]
  )

save(
  X,
  Y,
  Tr,
  Design,
  file = file.path(sprintf("Data/preprocessed_data_minoccs%s_atlasrichnessfilterbyeffort.RData",min_occs))
)

# 
# 
# # TAKE A LOOK AT LANDUSE VARS -------------------------------------------------
# 
# 
# # PREPARING MODEL BUILD ---------------------------------------------------
# Define model types:
# date <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
# 
# #atlasnr <- c('1')
# # loop over different atlases
# for(atlasnr in atlases){
# 
#   ### SAME ACROSS ALL MODELS
#   # Define model formulas for environmental and trait data
#   XFormula <- as.formula(paste("~", paste(colnames(X), collapse = "+"), sep = " "))
#   TrFormula <- as.formula(paste("~", paste(colnames(Tr), collapse = "+"), sep = " "))
# 
#   # get random effects for space
#   proj_xycoords_unique <- distinct(data.frame(X = Design$lon,
#                                               Y = Design$lat))
#   rownames(proj_xycoords_unique) <- unique(Design$site)
#   struc_space <- HmscRandomLevel(sData = proj_xycoords_unique, sMethod = "Full")
#   struc_space <- setPriors(struc_space,nfMin=5,nfMax=5) # set priors to limit latent factors
# 
#   # keep only atlas 1,2,3
#   pattern <- paste0("_(", paste(atlasnr, collapse = "|"), ")$")
# 
#   Y_sub <- Y[rownames(Y)[grep(pattern, rownames(Y))],,drop=F]
#   X_sub <- X[rownames(X)[grep(pattern, rownames(X))],,drop=F]
#   Design_sub <- Design[rownames(Design)[grep(pattern, rownames(Design))],,drop=F]
#   Design_sub$atlas <- droplevels(Design_sub$atlas)
#   #Design_sub <- Design_sub[,c('site','year'),drop=F]
# 
#   # if it's just 1 atlas
#   if(length(atlasnr)==1){
# 
#     m <-Hmsc(Y = Y_sub,
#              XData = X_sub,
#              XFormula = XFormula,
#              TrData = Tr,
#              TrFormula = TrFormula,
#              phyloTree = phy,
#              studyDesign = Design_sub[,c('site'),drop=F],
#              ranLevels = list('site'=struc_space),
#              distr='probit')
# 
#   }else{
#     # and time
#     years_unique <- distinct(data.frame(Year = Design_sub$year))
#     rownames(years_unique) <- unique(Design_sub$atlas)
#     struc_time <- HmscRandomLevel(sData = years_unique, sMethod = "Full")
# 
#     m <-Hmsc(Y = Y_sub,
#              XData = X_sub,
#              XFormula = XFormula,
#              TrData = Tr,
#              TrFormula = TrFormula,
#              phyloTree = phy,
#              studyDesign = Design_sub[,c('site','atlas'),drop=F],
#              ranLevels = list('site'=struc_space,
#                               'atlas'=struc_time),
#              distr='probit')
#   }
# 
#   print(head(Y_sub[,1:5]))
#   print(tail(Y_sub[,1:5]))
# 
#   summary(m)
#   m$rLNames
#   m$ranLevels
#   m$ranLevels$site
#   m$ranLevels$atlas
#   m$studyDesign
#   ### IN RSTUDIO START SAMPLING
#   if(flagFitR){
#     print('Rstudio stuff executed')
#     init_obj <- sampleMcmc(m, samples=nSamples, thin=thin,
#                            transient=transient, nChains=nChains,
#                            verbose = verbose,
#                            engine="HPC")
#   }
#   ### IN HPC ENVIORNMENT SET UP INIT
#   if(flagInit){
#     # initiate mcmc sampling
#     init_obj <- sampleMcmc(m, samples=nSamples, thin=thin,
#                            transient=transient, nChains=nChains,
#                            verbose = verbose,
#                            engine="HPC")
# 
#     dir_name <- paste0(date,name_of_dir,paste(atlasnr, collapse = ""))
#     dir.create(file.path(input,'tmp_rds',dir_name))
# 
#     init_file_path = file.path(input,'tmp_rds',dir_name, "init_file.rds")
#     m_file_path = file.path(input,'tmp_rds',dir_name, "m_object.rds")
#     param_file_path = file.path(input,'tmp_rds',dir_name, "params.rds")
#     lines <- paste(names(params), format(unlist(params), scientific = FALSE, trim=T), sep = "=")
#     writeLines(lines, file.path(input,'tmp_rds',dir_name, "params.txt"))
# 
#     # save as json
#     saveRDS(to_json(init_obj), file=init_file_path)
#     saveRDS(m,file=m_file_path)
#     saveRDS(params,file=param_file_path)
# 
#     # operates in python, so formulate the required call
#     post_file_path = file.path(input,'tmp_rds',dir_name, "post_file.rds")
#     python_cmd_args = paste("-m hmsc.run_gibbs_sampler",
#                             "--input", shQuote(init_file_path),
#                             "--output", shQuote(post_file_path),
#                             "--samples", nSamples,
#                             "--transient", format(transient,scientific=F),
#                             "--thin", thin,
#                             "--verbose", verbose)
#     cat(paste(shQuote(python), python_cmd_args), "\n")
#     print('Init files created')
# 



