rm(list=ls())
#### GETTING STARTED ####
if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,Hmsc,RColorBrewer,ggplot2,
               rnaturalearth,rnaturalearthdata,
               gridExtra,patchwork,
               readxl,cowplot)

devtools::install_github("davidsjoberg/ggbump")
devtools::install_github("davidsjoberg/ggsankey")
library(ggbump)
library(ggsankey)

scaled = F # set to NULL to not scale 

#### LOAD VP DATA ####
pattern <- "2026-03-13"
# get folders
matching_folders <- list.dirs('HmscOutputs', recursive = FALSE, full.names = F)
matching_folders <- matching_folders[grepl(pattern, basename(matching_folders))]
# get model numbers
model<-matching_folders[1]
models_nums <- as.numeric(str_extract(matching_folders, "(?<=Atlas)\\d+"))

# load VPs for each model
VPs <- lapply(matching_folders,function(model){
  atlas_num <- as.numeric(str_extract(model, "(?<=Atlas)\\d+"))
  df <- read.csv(file.path('HmscOutputs',
                           model,
                           'Results',
                           sprintf('%sparameter_estimates_VP_.csv',model)
  ) # load the info
  ) %>% column_to_rownames(var='X') # set rownames
  if(scaled==T){
    
  # 2. Extract the TjurR2 row as a named vector
  r2_values <- df["TjurR2", ] %>% as.numeric()
  # multiple with tjur value 
  df <- df %>%
    filter(rownames(.) != "TjurR2") %>%
    # Use mapply or map2 style logic to multiply columns by the vector
    map2_dfc(r2_values, ~ .x * .y) %>%
    # Put the row names back because map2_dfc drops them
    mutate(Factor = rownames(df)[rownames(df) != "TjurR2"]) %>%
    column_to_rownames(var = "Factor")
  }
  return(df)
})
# add appropriate atlas name
names(VPs) <- models_nums

#### PREPARE FOR PLOTTING ####
# Calculate rowwise medians for each dataframe in the list
median_list <- map(VPs, function(df) {
  # apply(df, 1, ...) calculates the median across every row
  row_medians <- apply(df, 1, median, na.rm = TRUE)
  row_means <- apply(df, 1, mean, na.rm = TRUE)
  
  
  # Return as a tibble for easy plotting/table making
  # Replace 'Variable_1:10' with your actual variable names if you have them
  tibble(
    variable = rownames(df), 
    median_VP = row_medians,
    mean_VP = row_means
  )
})

# Combine them into one dataframe for your final table
median_df <- bind_rows(median_list, .id = "atlas") %>%
  mutate(atlas = case_when(
    atlas == "1" ~ "1970s",
    atlas == "2" ~ "1990s",
    atlas == "3" ~ "2010s", # In case you added the fake 3rd one
    TRUE ~ atlas
  ))

# set colors 
base_colors <- c(
  "Site-level random effect" = "ivory2",
  "Temperature"     = "firebrick3",
  "Precipitation"   = "dodgerblue3",
  "Land-use heterogeneity"   = "orchid3",
  "Urban (% coverage)"           = "snow3",
  "Cropland (% coverage)"        = "goldenrod1", # Slightly darker so it doesn't vanish
  "Pasture (% coverage)"         = "darkorange",
  "Forest (% coverage)"          = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2"
)

# rename appropriately
plot_data_clean <- median_df %>%
  filter(!(variable %in% c("TjurR2",'Random: site'))) %>%  
  mutate(
    variable = case_match(variable,
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
    variable = factor(variable, levels = c(
      "Site-level random effect",
      "Temperature", "Precipitation", "Land-use heterogeneity", 
      "Urban (% coverage)", "Cropland (% coverage)", 
      "Pasture (% coverage)", "Forest (% coverage)", "Grass/Shrubland (% coverage)"
    )),
    atlas = factor(atlas,levels = c('1970s','1990s','2010s'))
  )

#### AND PLOTS ####
# prepare labels
y_label <- if (scaled==T) "Absolute variance explained" else "Relative variance explained"

plots <- list()
plots[[1]] <- ggplot(plot_data_clean, aes(x = atlas, 
                            node = variable, 
                            fill = variable, 
                            value = mean_VP)) +
  # The bump chart lines
  geom_sankey_bump(space = 0,type='alluvial',alpha=0.8) +
  

  # Fix the axis: Sankey bumps usually don't need Y-axis numbers
  # Removing them prevents the "negative" confusion
  #scale_y_continuous() + 
  scale_fill_manual(values = base_colors) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  
  # 
  labs(
    title = "Variable importance through time",
    #subtitle = "Relative contribution to Species Richness (VP)",
    x = NULL,
    fill = "Variable",
    y = y_label
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "none", # Hide legend here if you're putting it in Patchwork
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.background = element_rect(fill = "grey99", color = NA),
    axis.text.x = element_text(face = "bold", size = 12)
  )
plots[[1]]


#### LOAD EFFECT SIZES ####  
# load effects 
make_effect_sizes <- function(dir, pattern, effect = "Beta") {
  # 1. Input Validation and Setup
  # Standardize effect name to Title case for file matching
  effect_type <- ifelse(tolower(effect) == "gamma", "Gamma", "Beta")
  id_col_name <- ifelse(effect_type == "Gamma", "Traits", "Species")
  rename_to   <- ifelse(effect_type == "Gamma", "traits", "species")
  
  matching_folders <- list.dirs(dir, recursive = FALSE, full.names = FALSE)
  matching_folders <- matching_folders[grepl(pattern, basename(matching_folders))]
  
  # 2. Iterate through folders
  outputs <- lapply(matching_folders, function(model) {
    
    # Construct file path
    file_path <- file.path(dir, model, 'Results', 
                           sprintf('%sparameter_estimates_%s_.xlsx', model, effect_type))
    
    # Read Posterior Mean
    df_mean <- read_excel(file_path, sheet = 'Posterior mean')
    
    # Read Significance values [Pr(x>0)]
    df_sigs <- read_excel(file_path, sheet = "Pr(x>0)")
    
    # Optional: Load model posteriors if needed elsewhere in your script
    # load(file.path(dir, model, 'Models/Fitted', 'HPC_samples_0250_thin_100_chains_4.Rdata'))
    
    # 3. Pivot and Join
    long_effects <- df_mean %>%
      pivot_longer(cols = -all_of(id_col_name), 
                   names_to = "variable", 
                   values_to = "effect_size")
    
    long_sigs <- df_sigs %>%
      pivot_longer(cols = -all_of(id_col_name), 
                   names_to = "variable", 
                   values_to = "sig_val")
    
    # 4. Merge and Mask
    final_df <- long_effects %>%
      left_join(long_sigs, by = c(id_col_name, "variable")) %>%
      # Mask: sig <= 0.05 OR sig >= 0.95
      mutate(effect_size = ifelse(sig_val > 0.05 & sig_val < 0.95, NA, effect_size)) %>%
      # Standardize the ID column name (traits or species)
      rename(!!rename_to := all_of(id_col_name)) %>%
      select(-sig_val)
    
    # 5. Add Metadata
    atlas_num <- as.numeric(sub(".*Atlas([0-9]+).*", "\\1", model))
    final_df$atlas <- atlas_num
    
    return(final_df)
  })
  
  # 6. Combine all models
  merged <- bind_rows(outputs)
  return(merged)
}

beta <- make_effect_sizes('HmscOutputs','2026-03-13',effect='Beta')
gamma <- make_effect_sizes('HmscOutputs','2026-03-13',effect='Gamma')

n_species_total <- length(unique(beta$species))

# prepare betas
beta_processed <- beta %>%
  filter(variable != "(Intercept)", !is.na(effect_size)) %>%
  mutate(
    # Ensure Positive is the first level so it plots to the left of Negative
    direction = factor(ifelse(effect_size > 0, "Positive", "Negative"), 
                       levels = c("Positive", "Negative")),
    atlas = case_when(
      atlas == 1 ~ "1970s",
      atlas == 2 ~ "1990s",
      atlas == 3 ~ "2010s"
    )
  ) %>%
  group_by(atlas, variable, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  # Use the dynamic count here
  mutate(perc = (n / n_species_total) * 100)

# add placeholder info and prepare for plotting
beta_clean <- beta_processed %>%
  mutate(
    variable = case_match(variable,
                          "tmean_breeding"     ~ "Temperature",
                          "prec_breeding"      ~ "Precipitation",
                          "hh"                 ~ "Land-use heterogeneity",
                          "perc_urban"         ~ "Urban (% coverage)",
                          "perc_cropland"      ~ "Cropland (% coverage)",
                          "perc_pasture"       ~ "Pasture (% coverage)",
                          "perc_forest"        ~ "Forest (% coverage)",
                          "perc_grass_shrub"   ~ "Grass/Shrubland (% coverage)",
                          "Random: site"       ~ "Site-level random effect",
                          .default = variable
    ),
    plot_val = ifelse(direction == "Negative", -perc, perc),
    variable = factor(variable, levels = names(base_colors)),
    atlas = factor(atlas, levels = c('1970s', '1990s', '2010s'))
  )

# prepare aesthetics 
# set colors 
base_colors <- c(
  "Site-level random effect" = "ivory2",
  "Temperature"     = "firebrick3",
  "Precipitation"   = "dodgerblue3",
  "Land-use heterogeneity"   = "orchid3",
  "Urban (% coverage)"           = "snow3",
  "Cropland (% coverage)"        = "goldenrod1", # Slightly darker so it doesn't vanish
  "Pasture (% coverage)"         = "darkorange",
  "Forest (% coverage)"          = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2"
)

# get limits 
# 1. Get the highest absolute percentage in your data
current_max <- max(max(beta_clean$plot_val), abs(min(beta_clean$plot_val)))

# 2. Round up to the nearest multiple of 5
# Examples: 61.1 -> 65 | 66.0 -> 70 | 32.0 -> 35
plot_limit <- ceiling(current_max / 5) * 5

# and plot
plots[[2]] <- ggplot(beta_clean, aes(x = variable, y = plot_val, fill = variable,group = atlas,alpha = atlas)) +
  # Draw the bars
  geom_col(position = position_dodge(width = 0.9), color = "white", linewidth = 0.1) +
  
  # Add a horizontal line at 0
  geom_hline(yintercept = 0, color = "black") +
  
  # Use your base color scheme
  scale_fill_manual(values = base_colors) +
  
  scale_alpha_manual(values=c(0.4,0.7,1))+
  
  # Fix the Y-axis to show absolute numbers (no negative signs)
  # Dynamic Y-axis limits
  #scale_y_continuous(limits = c(-plot_limit, plot_limit)) +
  
  # Corner Annotations
  annotate("text", x = 4.5, y = 55, 
           label = "Positive effect", 
           hjust = 0.5, vjust = 0.5,
           fontface = "bold", size = 5, color = "black") +
  annotate("text", x = 4.5, y = -55, 
           label = "Negative effect", 
           hjust = 0.5, vjust = -0.5,
           fontface = "bold", size = 5, color = "black") +
  
  
  labs(
    title = 'Significant variable effects across modelled species',
    x = NULL,
    y = "Percentage of species",
    alpha = 'Atlas',
    fill = 'Variable'
  ) +
  
  theme_minimal() +
  theme(
    legend.direction = "vertical",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank()
  )

plots[[2]]

# ASSEMBLE ----------------------------------------------------------------
# get legend
# Shared Variable Legend (for Maps A & B and the Bar/Sankey plots)

# Alpha Legend (from plot 2 - Significance/Effect)
leg <- get_legend(
  plots[[2]] + 
    theme(legend.position = "right", 
          legend.direction = "vertical")+
    labs(fill = 'Variable',alpha = "Atlas")
  
)

layout <- "
AL
BL
"

final_layout <- wrap_plots(
  A = plots[[1]] + theme(legend.position = "none"), 
  B = plots[[2]] + theme(legend.position = "none"),
  L = leg,
  design = layout
) +
  # Adjust heights: Maps and Bars get most space, S (legend) gets very little
  plot_layout(heights = c(1, 1, 0.1), widths = c(1, 0.4)) + 
  
  # Tag only the actual data plots, not the legends
  plot_annotation(tag_levels = list(c('A', 'B', '', ''))) & 
  theme(plot.tag = element_text(face = "bold", size = 14))

final_layout

# Save with a slightly wider aspect ratio to accommodate the side legend
name <- if (scaled==T) sprintf('misc-figures/%s-fig2-scaled-vp-effect-sizes.png',pattern) else sprintf('misc-figures/%s-fig2-unscaled-vp-effect-sizes.png',pattern)

ggsave(name, final_layout, width = 10, height = 6)

