remove(list=ls()) 
.libPaths(c("~/Rlibs", .libPaths()))
require(Hmsc)
require(colorspace)
require(corrplot)
require(writexl)
require(cli)
require(tidyverse)


### Set up directories #### 
#If you are using RStudio this will set the working directory to exactly where the file is 
### Set up directories ####

### LOAD VP DATA 
patterns <- c('2026-02-10','2026-03-13','2026-03-24')
names <- c('all','goodaverage','good')
VPs <- lapply(patterns,function(pattern){
pattern2match <- pattern
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

names(VPs_scaled)

# Use imap_dfr to iterate over the list and its names/indices
combined_df <- imap_dfr(VPs_scaled, function(df, name) {
  df %>%
    rownames_to_column(var = 'variable') %>%
    pivot_longer(
      cols = -variable,
      names_to = 'species',
      values_to = 'VP'
    ) %>%
    mutate(atlas = name) # Add the list name/index to the 'atlas' column
}) 
})
names(VPs) <- names

## SET COLORS 
base_colors <- c(
  "Random: site" = "ivory2",
  "Temperature"     = "firebrick3",
  "Precipitation"   = "dodgerblue3",
  "Land-use heterogeneity"   = "orchid3",
  "Urban (% coverage)"           = "snow3",
  "Cropland (% coverage)"        = "goldenrod1", # Slightly darker so it doesn't vanish
  "Pasture (% coverage)"         = "darkorange",
  "Forest (% coverage)"          = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2"
)

str(VPs)
### COMBINE 
# This creates a new column called 'source' as the very first column
combined_df <- bind_rows(VPs, .id = "source")
unique(combined_df$variable)
### FILTER SITE FOR THIS PLOT
plot_data_clean <- combined_df %>%
  filter(!(variable %in% c("TjurR2","Random: site"))) %>%  
  mutate(
    variable_clean = case_match(variable,
                                "tmean_breeding"       ~ "Temperature",
                                "prec_breeding"        ~ "Precipitation",
                                "hh"               ~ "Land-use heterogeneity",
                                "perc_urban"       ~ "Urban (% coverage)",
                                "perc_cropland"    ~ "Cropland (% coverage)",
                                "perc_pasture"     ~ "Pasture (% coverage)",
                                "perc_forest"      ~ "Forest (% coverage)",
                                "perc_grass_shrub" ~ "Grass/Shrubland (% coverage)",
                                "Random: site" ~ "Site-level random effect",
                                .default = variable
    ),
    # Set the manual order for the X-axis
    variable_clean = factor(variable_clean, levels = rev(c(
      "Site-level random effect",
      "Temperature", "Precipitation", "Land-use heterogeneity", 
      "Urban (% coverage)", "Cropland (% coverage)", 
      "Pasture (% coverage)", "Forest (% coverage)", "Grass/Shrubland (% coverage)"
    ))),
    atlas = factor(atlas,levels = c('1','2','3')),
    source = factor(source,levels = c('all','goodaverage','good'))
  )


unique(plot_data_clean$variable)

### FILTER OUT 3 
plot_data_clean_sub <- plot_data_clean[plot_data_clean$atlas%in%c(1,2),]
# 2. Create the beautiful plot
ggplot(plot_data_clean_sub, 
       aes(x = variable_clean, y = VP, fill = variable_clean,alpha=source)) +
  geom_boxplot() +
  # Clean, publication-ready theme
  theme_minimal(base_size = 14) + 
  #scale_y_continuous(labels = scales::percent) +
  # Set the specific colors for the categories
  scale_fill_manual(values = base_colors, guide = "none") +
  # Set the "shades" using transparency (1.0 = dark/solid, 0.4 = light)
  scale_alpha_manual(values = c("all" = 0.4, "goodaverage" = 0.7, "good" = 1.0),
                     name = "Coverage",
                     labels = c("All", "Good/Average", "Good")) +
  facet_wrap(~atlas) +
  
  labs(
    title = "% of explained variance by environmental variables",
    subtitle = "",
    x = "% of explained variance (R²)",
    y = NULL,
  ) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = 'top'
  )   +
  guides(alpha = guide_legend(override.aes = list(alpha = c(0.4,0.7,1),
                                                   fill='firebrick3'))) 

plot_data_clean$atlas[plot_data_clean$VP > 0.6]
ggsave("misc-figures/env-vp-atlas-boxplots.png")

