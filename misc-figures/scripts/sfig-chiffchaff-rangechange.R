rm(list=ls())
#### GETTING STARTED ####
if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,Hmsc,RColorBrewer,ggplot2,
               rnaturalearth,rnaturalearthdata,
               gridExtra,patchwork,sf,cowplot,
               terra)

#### LOAD MODELS #### 
dir <- './HmscOutputs'
pattern <- '2026-03-13'

matching_folders <- list.dirs(dir, recursive = FALSE, full.names = F)
matching_folders <- matching_folders[grepl(pattern,
                                           basename(matching_folders))]
model <- matching_folders[1]
models_nums <- as.numeric(str_extract(matching_folders, "(?<=Atlas)\\d+"))

# load models
mods <- lapply(matching_folders,function(model){
  # read objects
  load(file.path(dir,model,'Models','Fitted','HPC_samples_0250_thin_100_chains_4.Rdata'))
  m <- fitted_model$posteriors
})
names(mods) <- models_nums
#mod <- mods[[1]]
# load designs
designs <- lapply(mods,function(mod){
  mod$studyDesign %>% 
    rownames_to_column(.,var='survey') %>% 
    left_join(.,
              mod$ranLevels$site$s %>% rownames_to_column(.,var='site'),
              by = 'site')
})

m <- mods[[1]]
species <- 'Phylloscopus_collybita'
eng_name <- 'CommonChiffchaff'

m$Y
m[['Y']]
# prepare for plotting 
dfs <- map2(mods, designs, ~ {
  # Calculate row sums and keep it as a data frame
  .x[['Y']] %>% 
    as.data.frame() %>% 
    select(all_of(species)) %>% 
    rownames_to_column(var = "survey") %>% 
    left_join(.y, by = 'survey')
})

denmark <- ne_countries(scale = "large", country = "denmark", returnclass = "sf")
year_lookup <- c("1" = "T1", "2" = "T2", "3" = "T3")
shape <- vect('~/box/PhD/logistics/data/distributions/DK5km_ED50grid_approx_kvadrkod_DOF/DK5km_ED50grid_approx_kvadrkod_DOF.shp')
shape_sf <- st_as_sf(shape)

df <- dfs[[1]]
mainland_bbox <- st_bbox(c(xmin = 400000,  # Moved further West
                           xmax = 750000,  # Moved further East
                           ymin = 6000000, # Moved further South
                           ymax = 6450000), # Moved further North
                         crs = st_crs(25832))

library(jpeg)
library(grid)
img <- readJPEG(paste0('../../../../box/PhD/old-projects/data/other/images_birds/',eng_name,'.jpg'))  # Use your own image path
img_grob <- rasterGrob(img, interpolate = TRUE)

name <- names(dfs)[[1]]
plots <- imap(dfs, function(df, name) {
  
  year_label <- year_lookup[[name]]
  
  # Ensure species column is a factor for manual scaling
  plot_data <- shape_sf %>%
    left_join(df, by = c("kvadratkod" = "site")) %>%
    mutate(Phylloscopus_collybita = as.factor(Phylloscopus_collybita))
  
  ggplot() +
    # 1. Background: The full country outline
    geom_sf(data = denmark, fill = 'grey90', color = "grey80", size = 0.3) +
    
    # 2. The Grid: Only '1' will be visible
    geom_sf(data = plot_data, aes(fill = Phylloscopus_collybita), color = NA) +
    
    # 3. Scale: Set 0 and NA to transparent
    scale_fill_manual(
      values = c("1" = "red", "0" = "transparent"),
      na.value = "transparent"
    ) +
    
    # 4. Zoom: Crop to mainland (removes Bornholm from view)
    # Adjust xlim/ylim based on your CRS (these are approx for WGS84)
    coord_sf(xlim = c(8, 13), ylim = c(54.5, 58), expand = FALSE) +
    
    labs(x = NULL, y = NULL, title = year_label) +
    theme_minimal() +
    theme(
      legend.position = 'none',
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      #axis.text = element_blank()
    ) +
  
  # 1. Add the Text (replacing text())
  # x and y are map coordinates; 'label' is the year/index
  annotate("text", x = 13.5, y = 57.7, label = eng_name, 
           size = 6, fontface = "bold", hjust = 1) +
    
    # 2. Add the Image (replacing grid.raster)
    # Adjust these map coordinates to place the image where you want it
    annotation_custom(img_grob, 
                      xmin = 12.5, xmax = 13.5, 
                      ymin = 57.0, ymax = 58.0) 
})

plots[[1]]
plots[[2]]
plots[[3]]


'../../../../box/'


plots[[1]]
