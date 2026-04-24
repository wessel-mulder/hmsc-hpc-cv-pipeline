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
# load MFs
MFs <- lapply(matching_folders,function(model){
  # read objects
  load(file.path(dir,
                    model,
                 'Models',
                    'Fitted',
                    'MF_samples_0250_thin_100_chains_4_nfolds_10.rdata'))
  list <- list(MF,MFCV)
  names(list) <- c('MF','MFCV')
  list
})
names(MFs) <- models_nums
str(MFs)


#### GET PREDS #### 
preds <- lapply(mods,pcomputePredictedValues)
predsY <- lapply(preds,rowMeans,dims=2)
rm(preds)

# prepare for plotting 
dfs <- map2(predsY, designs, ~ {
  # Calculate row sums and keep it as a data frame
  .x %>% 
    as.data.frame() %>% 
    mutate(richness = rowSums(.)) %>% 
    rownames_to_column(var = 'survey') %>% 
    select(survey,richness) %>% 
    # Now join with .y (the design list element)
    left_join(.y, by = 'survey')
})

# get richness limits & aesthetics 
lims <- dfs %>%
  imap_dfr(~ data.frame(
    name = .y,
    min = min(.x$richness, na.rm = TRUE),
    max = max(.x$richness, na.rm = TRUE)
  ))
lim <- c(min(lims$min),max(lims$max))
my_breaks <- seq(lim[1], lim[2], length.out = 5)

# colors 
pal <- colorRampPalette(rev(brewer.pal(11,'RdYlBu')))

#background 
denmark <- ne_countries(scale = "large", country = "denmark", returnclass = "sf")
# positions for annotation 
inset <- 12
df <- dfs[[3]]
x_pos <- max(df$X) - (diff(range(df$X)) / inset)
y_pos <- max(df$Y) - (diff(range(df$Y)) / inset)
# years 
year_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")

#### READ ATLAS POLYGON ####
shape <- vect('~/box/PhD/logistics/data/distributions/DK5km_ED50grid_approx_kvadrkod_DOF/DK5km_ED50grid_approx_kvadrkod_DOF.shp')
shape_sf <- st_as_sf(shape)
str(dfs[[1]])

# Define a bounding box for Bornholm (approximate UTM 32N)
# You might need to tweak these numbers based on your specific grid extent
# Increase the range between min and max to "zoom out"
mainland_bbox <- st_bbox(c(xmin = 400000,
                           xmax = 750000,
                           ymin = 6000000,
                           ymax = 6450000),
                         crs = st_crs(25832))
bornholm_bbox <- st_bbox(c(xmin = 855000, xmax = 905000, 
                           ymin = 6100000, ymax = 6160000), 
                         crs = st_crs(25832))

# Compute inset size as fraction of main plot dimensions
mainland_width  <- as.numeric(mainland_bbox["xmax"] - mainland_bbox["xmin"])
mainland_height <- as.numeric(mainland_bbox["ymax"] - mainland_bbox["ymin"])
bornholm_width  <- as.numeric(bornholm_bbox["xmax"] - bornholm_bbox["xmin"])
bornholm_height <- as.numeric(bornholm_bbox["ymax"] - bornholm_bbox["ymin"])

inset_w <- bornholm_width  / mainland_width
inset_h <- bornholm_height / mainland_height

#### MAPS ####
plots <- imap(dfs, function(df, name) {
  
  year_label <- year_lookup[[name]]
  
  plot_data <- shape_sf %>%
    left_join(df, by = c("kvadratkod" = "site"))
  
  # --- shared geom builder to keep fill scale consistent ---
  make_geom <- function(data) {
    list(
      geom_sf(data = data, aes(fill = richness), color = 'grey30', size = 0.1),
      scale_fill_gradientn(
        colors = pal(10),
        limits = lim,
        breaks = my_breaks,
        labels = round(my_breaks, 0),
        na.value = "transparent"
      )
    )
  }
  
  # --- Main plot (mainland) ---
  p_main <- ggplot() +
    make_geom(plot_data) +
    coord_sf(
      xlim = c(mainland_bbox["xmin"], mainland_bbox["xmax"]),
      ylim = c(mainland_bbox["ymin"], mainland_bbox["ymax"]),
      expand = FALSE
    ) +
    labs(x = NULL, y = NULL, title = year_label, fill = 'Richness') +
    theme_minimal() +
    theme(
      legend.position = 'none',
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_blank(),
      panel.grid = element_blank()
    )
  
  # --- Inset plot (Bornholm) ---
  p_inset <- ggplot() +
    make_geom(plot_data) +
    coord_sf(
      xlim = c(bornholm_bbox["xmin"], bornholm_bbox["xmax"]),
      ylim = c(bornholm_bbox["ymin"], bornholm_bbox["ymax"]),
      expand = FALSE
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.5)
    )
  
  # --- Combine with cowplot ---
  p_final <- ggdraw(p_main) +
    draw_plot(
      p_inset,
      x = 1 - inset_w - 0.2,
      y = 1 - inset_h - 0.2,
      width  = inset_w,
      height = inset_h
    )
  
  p_final <- as_grob(p_final)  # convert the complete object, not mid-pipe  
  plot(p_final)
  list(
    grob = as_grob(p_final),
    ggplot = p_main  # keep original for legend extraction
  )
  
})



# ACCURACY BOXPLOTS  ------------------------------------------------------
#### PREPARE FOR PLOTTING ####
df_plot_prep <- imap_dfr(MFs, function(atlas_data, atlas_name) {
  
  # Loop through MF and MFCV
  imap_dfr(atlas_data, function(metric_list, model_type) {
    
    # Convert the 3 metric vectors into a data frame
    as_tibble(metric_list) %>%
      mutate(
        atlas = atlas_name,
        evaluation = model_type,
        obs_id = row_number() # Keep track of the 157 observations
      )
  })
}) %>%
  # Pivot metrics into one column for easy faceting
  pivot_longer(cols = c(RMSE, AUC, TjurR2), 
               names_to = "metric", 
               values_to = "value")



# rename
df_plot <- df_plot_prep %>%
  mutate(
    # # 1. Rename and Order Metrics
    # metric = case_when(
    #   metric == "AUC" ~ "Area Under the Curve",
    #   metric == "TjurR2" ~ "Tjur R²",
    #   metric == "RMSE" ~ "Root Mean Squared Error"
    # ),
    # metric = factor(metric, levels = c("Area Under the Curve", "Tjur R²", "Root Mean Squared Error")),
    # 
    # 1. Rename and Order Metrics
    metric = case_when(
      metric == "AUC" ~ "AUC",
      metric == "TjurR2" ~ "Tjur R²",
      metric == "RMSE" ~ "RMSE"
    ),
    metric = factor(metric, levels = c("AUC", "Tjur R²", "RMSE")),

    # 2. Rename Evaluation types for X-axis
    evaluation = case_when(
      evaluation == "MF" ~ "MF",
      evaluation == "MFCV" ~ "CV"
    ),
    evaluation = factor(evaluation, levels = c("MF",'CV')),
    
    
    # 3. Rename Atlas for the legend
    atlas = case_when(
      atlas == "1" ~ "1970s",
      atlas == "2" ~ "1990s",
      atlas == "3" ~ "2010s"
    )
  )

#### PLOTTING ####
plots[[4]] <- ggplot(df_plot, aes(y = value, alpha = atlas, fill = metric, x = atlas)) +
  geom_boxplot(fill = 'grey40',
               position = position_dodge(width = 0.8),
               color = "grey30",
               outlier.shape = NA) +
  facet_grid(~metric) +
  scale_alpha_manual(
    values = c("1970s" = 0.2, "1990s" = 0.5, "2010s" = 0.9),
    name = "Atlas Period"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none",   # <-- changed
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_blank(),
    legend.title = element_text(face = "bold")
  )

plots[[4]]$grob <- as_grob(plots[[4]])
plots[[4]]$ggplot <- plots[[4]]  # so get_legend() can access it downstream

# ASSEMBLE ----------------------------------------------------------------
#### GET LEGENDS ####
# Extract the Richness legend from one of your map plots
# (Make sure plots[[1]] has legend.position = "right" temporarily to catch it)

layout <- "
AB
CD
EG
"

final_layout <- wrap_plots(
  A = plots[[1]]$grob,
  B = plots[[2]]$grob,
  C = plots[[3]]$grob,
  D = plots[[4]]$grob,
  E = as_grob(get_legend(plots[[1]]$ggplot + theme(legend.position = "bottom"))),
  G = as_grob(get_legend(plots[[4]]$ggplot + theme(legend.position = "bottom")))
  
) +
  plot_layout(design = layout, heights = c(1, 1,0.05)) +
  plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D','',''))) &
  theme(plot.tag = element_text(face = "bold", size = 14))

final_layout

plot(plots[[1]]$grob)
plot(plots[[2]]$grob)
plot(plots[[4]]$grob)



ggsave(sprintf('misc-figures/%s-fig1-richness-accuracies.png',pattern),width = 7.5,height= 7.5)

