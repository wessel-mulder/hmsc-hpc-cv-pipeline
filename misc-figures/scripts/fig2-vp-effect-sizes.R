rm(list=ls())
#### GETTING STARTED ####
if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,Hmsc,RColorBrewer,ggplot2,
               rnaturalearth,rnaturalearthdata,
               gridExtra,patchwork,
               readxl,cowplot,scales)
source(file.path("support_scripts", "figure_data_helpers.R"))

if (!requireNamespace("ggbump", quietly = TRUE) || !requireNamespace("ggsankey", quietly = TRUE)) {
  stop(
    "Install the `ggbump` and `ggsankey` packages before running this script. ",
    "They are used for the variance-partitioning alluvial panel.",
    call. = FALSE
  )
}
library(ggbump)
library(ggsankey)

scaled = T # set to NULL to not scale 

predictor_sds_by_atlas <- function(models, variables) {
  imap_dfr(models, function(model, atlas) {
    tibble(
      atlas = as.numeric(atlas),
      variable = variables,
      predictor_sd = map_dbl(variables, ~ sd(model$XData[[.x]], na.rm = TRUE))
    )
  })
}

#### LOAD VP DATA ####
pattern <- "2026-03-13"
matching_folders <- figure_model_folders(pattern = pattern)
model <- matching_folders[1]
models_nums <- atlas_numbers(matching_folders)
VPs <- load_vp_estimates(matching_folders, scaled = scaled)

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
beta <- read_parameter_effects(pattern, effect = "Beta")
gamma <- read_parameter_effects(pattern, effect = "Gamma")

n_species_total <- length(unique(beta$species))

# These SDs convert beta coefficients into effects for a one-standard-deviation
# change in the environmental variable. The SD is calculated separately within
# each atlas, matching the atlas-specific model used to estimate each beta.
mods <- load_hmsc_posteriors(matching_folders)
beta_variables <- intersect(setdiff(unique(beta$variable), "(Intercept)"), colnames(mods[[1]]$XData))
beta_predictor_sds <- predictor_sds_by_atlas(mods, beta_variables)

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

# Prepare a second version of the beta bar panel. Instead of counting the
# percentage of species with supported effects, this sums the magnitude of the
# atlas-specific SD-scaled beta effects across species. Unsupported beta
# coefficients have already been removed by `read_parameter_effects()`, which
# keeps only effects with Pr(x > 0) >= 0.95 or Pr(x > 0) <= 0.05.
beta_total_effect <- beta %>%
  filter(variable != "(Intercept)", !is.na(effect_size)) %>%
  inner_join(beta_predictor_sds, by = c("atlas", "variable")) %>%
  mutate(
    sd_scaled_effect_size = effect_size * predictor_sd,
    direction = factor(
      ifelse(sd_scaled_effect_size > 0, "Positive", "Negative"),
      levels = c("Positive", "Negative")
    ),
    atlas = case_when(
      atlas == 1 ~ "1970s",
      atlas == 2 ~ "1990s",
      atlas == 3 ~ "2010s"
    )
  ) %>%
  group_by(atlas, variable, direction) %>%
  summarise(
    n_species = n_distinct(species),
    total_abs_sd_scaled_effect = sum(abs(sd_scaled_effect_size), na.rm = TRUE),
    mean_abs_sd_scaled_effect = mean(abs(sd_scaled_effect_size), na.rm = TRUE),
    .groups = "drop"
  )

beta_total_effect_clean <- beta_total_effect %>%
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
    plot_val = ifelse(direction == "Negative", -total_abs_sd_scaled_effect, total_abs_sd_scaled_effect),
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

# The total-effect version keeps the same signed layout as the percentage plot,
# but the bar height is now the summed SD-scaled beta magnitude among species
# with at least 95% posterior support for the effect direction.
total_effect_max <- max(abs(beta_total_effect_clean$plot_val), na.rm = TRUE)
total_effect_limit <- pretty(c(-total_effect_max, total_effect_max), n = 6)
total_effect_limit <- max(abs(total_effect_limit))

plots[[3]] <- ggplot(beta_total_effect_clean, aes(x = variable, y = plot_val, fill = variable, group = atlas, alpha = atlas)) +
  geom_col(position = position_dodge(width = 0.9), color = "white", linewidth = 0.1) +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(values = base_colors) +
  scale_alpha_manual(values=c(0.4,0.7,1))+
  scale_y_continuous(
    limits = c(-total_effect_limit, total_effect_limit),
    labels = label_number(accuracy = 0.1)
  ) +
  annotate("text", x = 4.5, y = total_effect_limit * 0.82,
           label = "Positive effect",
           hjust = 0.5, vjust = 0.5,
           fontface = "bold", size = 5, color = "black") +
  annotate("text", x = 4.5, y = -total_effect_limit * 0.82,
           label = "Negative effect",
           hjust = 0.5, vjust = 0.5,
           fontface = "bold", size = 5, color = "black") +
  labs(
    title = 'Total supported variable effects across modelled species',
    x = NULL,
    y = "Total SD-scaled effect size",
    alpha = 'Atlas',
    fill = 'Variable'
  ) +
  theme_minimal() +
  theme(
    legend.direction = "vertical",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank()
  )

plots[[3]]

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

total_effect_leg <- get_legend(
  plots[[3]] +
    theme(legend.position = "right",
          legend.direction = "vertical")+
    labs(fill = 'Variable',alpha = "Atlas")
)

total_effect_layout <- wrap_plots(
  A = plots[[1]] + theme(legend.position = "none"),
  B = plots[[3]] + theme(legend.position = "none"),
  L = total_effect_leg,
  design = layout
) +
  plot_layout(heights = c(1, 1, 0.1), widths = c(1, 0.4)) +
  plot_annotation(tag_levels = list(c('A', 'B', '', ''))) &
  theme(plot.tag = element_text(face = "bold", size = 14))

total_effect_layout

# Save with a slightly wider aspect ratio to accommodate the side legend
name <- if (scaled==T) sprintf('misc-figures/%s-fig2-scaled-vp-effect-sizes.png',pattern) else sprintf('misc-figures/%s-fig2-unscaled-vp-effect-sizes.png',pattern)

ggsave(name, final_layout, width = 10, height = 6)

total_effect_name <- if (scaled==T) {
  sprintf('misc-figures/%s-fig2-scaled-vp-total-supported-effect-sizes.png', pattern)
} else {
  sprintf('misc-figures/%s-fig2-unscaled-vp-total-supported-effect-sizes.png', pattern)
}

ggsave(total_effect_name, total_effect_layout, width = 10, height = 6)

write_csv(
  beta_total_effect_clean,
  sprintf('misc-figures/%s-fig2-total-supported-sd-scaled-beta-effects.csv', pattern)
)
