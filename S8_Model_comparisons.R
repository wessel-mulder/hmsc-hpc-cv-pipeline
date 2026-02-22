remove(list=ls())
.libPaths(c("~/Rlibs", .libPaths()))
require(Hmsc)
require(colorspace)
require(corrplot)
require(writexl)
require(cli)

set.seed(369)

### Set up directories #### 
#If you are using RStudio this will set the working directory to exactly where the file is 
### Set up directories ####

### LOAD VP DATA 
pattern2match <- "2026-02-10"

matching_folders <- list.dirs('HmscOutputs', recursive = FALSE, full.names = F)
matching_folders <- matching_folders[grepl(pattern2match, basename(matching_folders))]

model<-matching_folders[1]
models_nums <- as.numeric(str_extract(matching_folders, "(?<=Atlas)\\d+"))
VPs_scaled <- lapply(matching_folders,function(model){
  atlas_num <- as.numeric(str_extract(model, "(?<=Atlas)\\d+"))
  df <- read.csv(file.path('HmscOutputs',
                     model,
                     'Results',
                     sprintf('%sparameter_estimates_VP_.csv',model)
           ) # load the info
  ) %>% column_to_rownames(var='X') # set rownames
  
  
  # # 2. Extract the TjurR2 row as a named vector
  # r2_values <- df["TjurR2", ] %>% as.numeric()
  # 
  # # 3. Filter out the TjurR2 row and multiply the rest
  # # 3. Clean and Multiply
  # df_scaled <- df %>%
  #   filter(rownames(.) != "TjurR2") %>%
  #   # Use mapply or map2 style logic to multiply columns by the vector
  #   map2_dfc(r2_values, ~ .x * .y) %>%
  #   # Put the row names back because map2_dfc drops them
  #   mutate(Factor = rownames(df)[rownames(df) != "TjurR2"]) %>%
  #   column_to_rownames(var = "Factor")
  
  return(df)
})
names(VPs_scaled) <- models_nums

### GET AVERAGES
VPs_means <- lapply(VPs_scaled,function(VP){
  VP %>%
    mutate(average = rowMeans(pick(everything()),na.rm=T)) %>%
    select(average)
})

### PREPARE DATA 
plot_data <- bind_rows(VPs_means, .id = "Atlas") %>%
  rownames_to_column("variable") %>%
  mutate(variable = str_remove(variable, "\\.\\.\\.\\d+$"))

### SET COLORS 
base_colors <- c(
  "Temperature"     = "red",
  "Precipitation"   = "blue",
  "Land-use heterogeneity"   = "pink",
  "Urban (% coverage)"           = "grey",
  "Cropland (% coverage)"        = "yellow3", # Slightly darker so it doesn't vanish
  "Pasture (% coverage)"         = "orange",
  "Forest (% coverage)"          = "darkgreen",
  "Grass/Shrubland (% coverage)" = "lightgreen"
)

### FILTER SITE FOR THIS PLOT
plot_data_clean <- plot_data %>%
  filter(!(variable %in% c("Random: site", "TjurR2"))) %>%  
  mutate(
    variable_clean = case_match(variable,
                                "tmean_year"       ~ "Temperature",
                                "prec_year"        ~ "Precipitation",
                                "hh"               ~ "Land-use heterogeneity",
                                "perc_urban"       ~ "Urban (% coverage)",
                                "perc_cropland"    ~ "Cropland (% coverage)",
                                "perc_pasture"     ~ "Pasture (% coverage)",
                                "perc_forest"      ~ "Forest (% coverage)",
                                "perc_grass_shrub" ~ "Grass/Shrubland (% coverage)",
                                .default = variable
    ),
    # Set the manual order for the X-axis
    variable_clean = factor(variable_clean, levels = c(
      "Temperature", "Precipitation", "Land-use heterogeneity", 
      "Urban (% coverage)", "Cropland (% coverage)", 
      "Pasture (% coverage)", "Forest (% coverage)", "Grass/Shrubland (% coverage)"
    )),
    atlas = factor(Atlas)
  )

# 2. Create the beautiful plot
ggplot(plot_data_clean, aes(x = variable_clean, y = average, fill = variable_clean,alpha=atlas)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  # Clean, publication-ready theme
  theme_minimal(base_size = 14) + 
  scale_y_continuous(labels = scales::percent) +
  # Set the specific colors for the categories
  scale_fill_manual(values = base_colors, guide = "none") +
  # Set the "shades" using transparency (1.0 = dark/solid, 0.4 = light)
  scale_alpha_manual(values = c("1" = 0.4, "2" = 0.7, "3" = 1.0), 
                     name = "Atlas Period",
                     labels = c("1970s", "1990s", "2010s")) +
  
  labs(
    title = "% of explained variance by environmental variables",
    subtitle = "Average across all species (Spatial effects excluded)",
    x = NULL, # Variable names are self-explanatory
    y = "% of explained variance (R²)"
  ) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )




# CHECK OUT STACKED -------------------------------------------------------

library(tidyverse)

# 1. Prepare data (Including Site, excluding TjurR2)
plot_data_stacked <- plot_data %>%
  filter(variable != "TjurR2") %>%
  mutate(
    variable_clean = case_match(variable,
                                "tmean_year"       ~ "Temperature",
                                "prec_year"        ~ "Precipitation",
                                "hh"               ~ "Land-use heterogeneity",
                                "perc_urban"       ~ "Urban",
                                "perc_cropland"    ~ "Cropland",
                                "perc_pasture"     ~ "Pasture",
                                "perc_forest"      ~ "Forest",
                                "perc_grass_shrub" ~ "Grass/Shrubland",
                                "Random: site"     ~ "Random: Site Effect",
                                .default = variable
    ),
    # Ensure the order is correct (Environmental first, then Random at the top)
    variable_clean = factor(variable_clean, levels = rev(c(
      "Temperature", "Precipitation", "Land-use heterogeneity", 
      "Urban", "Cropland", "Pasture", "Forest", "Grass/Shrubland",
      "Random: Site Effect"
    )))
  )

# 2. Define colors (Adding a neutral grey/tan for the Random effect)
stack_colors <- c(
  "Temperature"            = "red",
  "Precipitation"          = "blue",
  "Land-use heterogeneity" = "pink",
  "Urban"                  = "grey40",
  "Cropland"               = "yellow3",
  "Pasture"                = "orange",
  "Forest"                 = "darkgreen",
  "Grass/Shrubland"        = "lightgreen",
  "Random: Site Effect"    = "grey80" # Neutral tone for the "background" noise
)

# 3. Plot
ggplot(plot_data_stacked, aes(x = Atlas, y = average, fill = variable_clean)) +
  # position = "fill" scales everything to 100% (0 to 1)
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.1) +
  
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = stack_colors) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Relative Contribution to Explained Variance",
    subtitle = "Stacked to 100% (Includes Site Random Effect)",
    x = "Atlas Period",
    y = "Proportion of Explained Variance",
    fill = "Variable"
  ) +
  theme(panel.grid.minor = element_blank())

