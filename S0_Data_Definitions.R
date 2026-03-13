
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
input <- '.'
coveragelevels <- 'good_and_average'
min_occs <- 5


#### support scripts  ######
get_effort <- function(dir){
  atlas1 <- read_excel(file.path(dir, 'effort', 'dækning atlas1+2.xlsx'), sheet = 'kvadrater atlas1')
  atlas2 <- read_excel(file.path(dir, 'effort', 'dækning atlas1+2.xlsx'), sheet = 'kvadrater atlas2')
  
  atlas2 <- atlas2 %>% 
    select(-UNDS_KODE_) %>% 
    rename(UNDS_KODE_ = UNDS_KOD_1)
  
  # Load design (Note: 'design' object is overwritten here in original code)
  # load('HmscOutputs/2026-02-10_07-03-53_All_All_Atlas3_MinOccs5/Models/Fitted/HPC_samples_0250_thin_100_chains_4.Rdata')
  # design <- read.csv('Data/data/1_preprocessing/design/studyDesign.csv')
  # xycoords <- design[design$atlas == 2, c('site', 'lat', 'lon')]
  
  key <- read_excel(file.path(dir, 'effort', 'dækning atlas1+2.xlsx'), sheet = 'dækning') %>% 
    mutate(coverage = case_when(
      UNDSKODE == 0 ~ 'Not understood',
      UNDSKODE == 1 ~ 'Good',
      UNDSKODE == 2 ~ 'Average',
      UNDSKODE == 3 ~ 'Bad'
    ))
  
  effort <- rbind(
    atlas1 %>% mutate(atlas = 1),
    atlas2 %>% mutate(atlas = 2)
  ) %>% 
    left_join(key %>% rename(UNDS_KODE_ = UNDSKODE), by = "UNDS_KODE_") %>% 
    rename(site = KVADRATKOD) #%>% 
  #left_join(xycoords, by = 'site') 
  
  effort$coverage <- factor(effort$coverage, levels = c('Good','Average','Bad','Not understood'))
  
  effort <- effort %>% 
    mutate(survey = paste(.$site, .$atlas, sep='_')) %>% 
    select(survey, coverage) 
  
  return(effort)
} ### PREPARE COVERAGE 

# LOADING DATA -----------------------------------------------------
#### ENVIRONMENT ####
# OCEAN THRESHOLD
grids_thresholds <- st_read(file.path(input,'Data/data/1_preprocessing/atlas-grids/grids-ocean-thresholds/grids_ocean_thresholds.shp'))
thresholds <- grids_thresholds$kvdrtkd[grids_thresholds$pct_lnd>=25] ### FILTERING GRID CELLS WITH LESS

effort <- get_effort(file.path(input,'Data')) %>% 
  ### MODIFY THIS TO EDIT THE ATLAS THRESHOLD FOR 3 WHEN I GET IT
  rbind(., data.frame(
    survey = paste(thresholds, 3, sep='_'),
    coverage = 'Good'
  ))

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
X_inter <- X

### GET LIST OF SURVEYS WITH BAD / UNKOWN SURVEYS
if(coveragelevels=='good_and_average'){surveys2keep <- effort$survey[effort$coverage%in%c('Good','Average')]}
if(coveragelevels=='good'){surveys2keep <- effort$survey[effort$coverage%in%c('Good')]}

X <- X_inter %>% 
  rownames_to_column(.,var='survey') %>% 
  filter(survey %in% surveys2keep) %>% 
  column_to_rownames(.,var = 'survey')

# get rid of it 
rm(X_inter)
  
#### OCCURRENCEcoverage#### OCCURRENCES ####  
Y <- read.csv(
  file.path(input, "Data/data/1_preprocessing/Y_occurrences/Y_occurrences.csv"),
  row.names = 1
) %>%
  filter(row.names(.) %in% row.names(X))

### filter based on coverage per atlas 
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

# get rid of it 
rm(Y_sub)

# now refilter X again 
X <- X %>%
  filter(rownames(.) %in% rownames(Y))

#### MERGE IN THE STI ####
sti <- read.csv('Data/sti/STI_Devictor.csv')

my_list <- colnames(Y)
sti$SPECIES <- gsub(" ", "_", sti$SPECIES)

setdiff(my_list,sti$SPECIES)
# astur_gentilis missing 
grep('gentilis',sti$SPECIES,value=T)  # accipiter gentilis --> astur gentilis
# Anarhynchus_alexandrinus missing 
grep('alexandrinus',sti$SPECIES,value=T)  # charadrius alexandrinus --> Anarhynchus_alexandrinus
# Chroicocephalus_ridibundus
grep("ridibundus", sti$SPECIES, value = TRUE) # Larus ridibundus --> Chroicocephalus_ridibundus
# Coloeus_monedula
grep("monedula", sti$SPECIES, value = TRUE)  # corvus monedula --> coloeus_monedula
# corvus cornix
grep("cornix", sti$SPECIES, value = TRUE)    # nada
grep("corone", sti$SPECIES, value = TRUE)    # corvus corone --> subspecies designation now corvus cornix
# Curruca_communis
grep("communis", sti$SPECIES, value = TRUE)  #  Sylvia communis
# Curruca_curruca
grep("curruca", sti$SPECIES, value = TRUE)   # often Sylvia curruca

# 1. Create the translation mapping
# Format: "Their Name" = "Your Name"
name_map <- c(
  "Accipiter_gentilis"       = "Astur_gentilis",
  "Charadrius_alexandrinus"  = "Anarhynchus_alexandrinus",
  "Larus_ridibundus"         = "Chroicocephalus_ridibundus",
  "Corvus_monedula"          = "Coloeus_monedula",
  "Sylvia_communis"          = "Curruca_communis",
  "Sylvia_curruca"           = "Curruca_curruca"
)

# 2. Perform the replacement
# We use 'match' to find which indices in their_data need updating
sti$SPECIES <- ifelse(sti$SPECIES %in% names(name_map), 
                      name_map[sti$SPECIES], 
                      sti$SPECIES)

sti_inter <- sti

# 3. Add Corvus_cornix using the value from Corvus corone
# This keeps 'Corvus corone' intact and creates a new entry for 'Corvus_cornix'
if ("Corvus_corone" %in% sti_inter$SPECIES) {
  corone_data <- sti_inter[sti_inter$SPECIES == "Corvus_corone", ]
  corone_data$SPECIES <- "Corvus_cornix"
  # Append this new row to the dataframe
  sti <- rbind(sti_inter, corone_data)
}

### CHECK IF THIS IS RIHT 
sti[sti$SPECIES == 'Corvus_cornix',]
sti[sti$SPECIES == 'Corvus_corone',]


their_final_list <- unique(sti$SPECIES)
missing <- setdiff(my_list, their_final_list)
print(paste("Still missing:", paste(missing, collapse=", ")))

#### MERGE WITH Tr
##### TRAITS #####
Tr <- read.csv(
  file.path(input, "Data/data/1_preprocessing/Tr_aits/traits-guild_migration.csv"),
  row.names = 2
) %>%
  dplyr::select(2, 3) %>%
  .[colnames(Y), , drop = FALSE] %>% 
  rownames_to_column(.,var = 'SPECIES') %>% 
  left_join(sti %>% 
              select(SPECIES,STI) %>% 
              rename(species_thermal_index = STI),
            by = 'SPECIES') %>% 
  column_to_rownames(.,var = 'SPECIES')

sti$STI[sti$SPECIES == 'Accipiter_nisus']
Tr$species_thermal_index[Tr$SPECIES == 'Accipiter_nisus']

sti$STI[sti$SPECIES == 'Picus_viridis']
Tr$species_thermal_index[Tr$SPECIES == 'Picus_viridis']

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
  file = file.path(sprintf("Data/preprocessed_data_min_occs_is_%s_coverage_is_%s.RData",min_occs,coveragelevels))
)
