remove(list=ls()) 
.libPaths(c("~/Rlibs", .libPaths()))
require(Hmsc)
require(colorspace)
require(corrplot)
require(writexl)
require(cli)
require(tidyverse)
require(readxl)

source(file.path("support_scripts", "figure_data_helpers.R"))


### Set up directories #### 
#If you are using RStudio this will set the working directory to exactly where the file is 
### Set up directories ####

### LOAD VP DATA 
pattern2match <- "2026-03-13"

matching_folders <- list.dirs('HmscOutputs', recursive = FALSE, full.names = F)
matching_folders <- matching_folders[grepl(pattern2match, basename(matching_folders))]


model<-matching_folders[1]
models_nums <- as.numeric(str_extract(matching_folders, "(?<=Atlas)\\d+"))

#### GET VPs ####
VPs_scaled <- lapply(matching_folders,function(model){
  atlas_num <- as.numeric(str_extract(model, "(?<=Atlas)\\d+"))
  df <- read.csv(file.path('HmscOutputs',
                           model,
                           'Results',
                           sprintf('%sparameter_estimates_VP_.csv',model)
  ) # load the info
  ) %>% column_to_rownames(var='X') # set rownames
  
  # 2. Extract the TjurR2 row as a named vector
  r2_values <- df["TjurR2", ] %>% as.numeric()

  # 3. Filter out the TjurR2 row and multiply the rest
  # 3. Clean and Multiply
  df_scaled <- df %>%
    filter(rownames(.) != "TjurR2") %>%
    # Use mapply or map2 style logic to multiply columns by the vector
    map2_dfc(r2_values, ~ .x * .y) %>%
    # Put the row names back because map2_dfc drops them
    mutate(Factor = rownames(df)[rownames(df) != "TjurR2"]) %>%
    column_to_rownames(var = "Factor")
  
  return(df_scaled)
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

#### GET TRAITS ####
load(file.path('HmscOutputs',
               model,
               'Models/Fitted',
               sprintf('HPC_samples_0250_thin_100_chains_4.Rdata')))
m <- fitted_model$posteriors
Tr <- m$TrData

#### PLOTTING PARAMS ####
## SET COLORS 
base_colors <- c(
  "Site-level random effects" = "ivory2",
  "Temperature"     = "firebrick3",
  "Precipitation"   = "dodgerblue3",
  "Land-use heterogeneity"   = "orchid3",
  "Urban (% coverage)"           = "snow3",
  "Cropland (% coverage)"        = "goldenrod1", # Slightly darker so it doesn't vanish
  "Pasture (% coverage)"         = "darkorange",
  "Forest (% coverage)"          = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2"
)




### FILTER SITE FOR THIS PLOT
plot_data_clean <- combined_df %>%
  filter(!(variable %in% c("TjurR2"))) %>%  
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
                                "Random: site" ~ "Site-level random effects",
                                .default = variable
    ),
    # Set the manual order for the X-axis
    variable_clean = factor(variable_clean, levels = rev(c(
      "Temperature", "Precipitation", "Land-use heterogeneity", 
      "Urban (% coverage)", "Cropland (% coverage)", 
      "Pasture (% coverage)", "Forest (% coverage)", "Grass/Shrubland (% coverage)",
      "Site-level random effects"
    ))),
    atlas = factor(atlas,levels = rev(c('1','2','3')))
  )

#### PREPARE PLOTTING ####
atlas <- 1
subset <- plot_data_clean[plot_data_clean$atlas==atlas,] %>% 
  left_join(.,
            Tr %>% 
              rownames_to_column(var = 'species'),
            by = 'species')

sorted_guilds <- c(
  # 1. Terrestrial Vertivores
  "Diurnal raptors", "Owls", "Shrikes",
  
  # 2. Terrestrial Invertivores
  "Aerial insectivores", "Nigthjars", "Cuckoos", "Tit-like birds", 
  "Open-land insectivores", "Low flycatching feeders", "Foliage gleaners", 
  "Flycatchers", "Woodpeckers",
  
  # 3. Terrestrial Granivores/Frugivores
  "Columbids", "Gallinaceous birds", "Passerine seedeaters", 
  "Buntings", "Orioles",
  
  # Omnivores
  "Gulls", "Omnivorous corvids", "Starlings", "Thrushes",
  
  # 4. Aquatic/Wetland/Coastal Vertivores
  "Diving ducks", "Aquatic pursuers", "Terns", "Kingfishers", "Wading birds",
  
  # 5. Aquatic/Wetland/Coastal Invertivores
  "Dabbling ducks", "Plovers", "Scolopacids", "Rails", 
  "Marsh warblers", "Dippers", "Reedlings",
  
  # 6. Grazers
  "Grazing waterfowl"
)

subset <- subset %>%
  mutate(species = str_replace(species, "^([A-Z])[a-z]+_", "\\1. "))



subset_ordered <- subset %>%
  # Calculate total VP per species (so all rows for one species stay together)
  filter(!variable_clean == 'Site-level random effects') %>% 
  group_by(species) %>%
  mutate(total_species_VP = sum(VP, na.rm = TRUE)) %>%
  ungroup() %>%
  # Convert guild to factor using our custom levels
  mutate(foraging_guild_consensus = factor(foraging_guild_consensus, levels = sorted_guilds)) %>%
  # Arrange: 1st by the Guild levels, 2nd by Total VP (descending)
  # Change to total_species_VP (no desc) if you want lowest VP first
  arrange(foraging_guild_consensus, desc(total_species_VP)) %>%
  # Optional: Remove the helper column if you don't need it anymore
  select(-total_species_VP)

subset$species <- factor(subset$species,
                         levels = unique(subset_ordered$species))

#### AND PLOT ####
# ggplot(subset_ordered, 
#        aes(x = species, y = VP, fill = variable_clean)) +
#   geom_point(pch=21,size=1.5) +
#   # Clean, publication-ready theme
#   theme_minimal(base_size = 14) + 
#   #scale_y_continuous(labels = scales::percent) +
#   # Set the specific colors for the categories
#   scale_color_manual(values = base_colors) +
# 
#   labs(
#     title = "% of explained variance by environmental variables",
#     subtitle = "",
#     x = "% of explained variance (R²)",
#     y = NULL,
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 40, hjust = 1),
#     panel.grid.minor = element_blank(),
#     legend.position = 'top'
#   ) 
# 
# plot_data_clean$atlas[plot_data_clean$VP > 0.6]
# ggsave("misc-figures/env-vp-species-points.png")
if(atlas==1){yearlabel = '1970s'}
if(atlas==2){yearlabel = '1990s'}
if(atlas==3){yearlabel = '2010s'}


ggplot(subset, 
       aes(x = species, y = VP, fill = variable_clean)) +
  geom_bar(stat='identity',position='stack',width=1,color='white',
           linewidth = 0.1) +
  # Clean, publication-ready theme
  theme_minimal(base_size = 14) + 
  #scale_y_continuous(labels = scales::percent) +
  # Set the specific colors for the categories
  scale_fill_manual(values = base_colors) +

  #facet_grid(~foraging_guild_consensus,scales = 'free_x')+

  labs(
    title = "Proportions of variance explained",
    subtitle = yearlabel,
    y = "Absolute variance explained",
    x = NULL,
    fill = 'Variable'
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust =1,size = 7),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = 'bottom',
    guide_legend(reverse=T)
  ) 
ggsave(sprintf("misc-figures/%s-env-vp-species-stacked-atlas-%s.png",pattern2match,atlas))


# guild-seperated --------------------------------------------------------------------
# 
# # 1. Ensure levels are set correctly
# subset_ordered$foraging_guild_consensus <- factor(subset_ordered$foraging_guild_consensus, levels = sorted_guilds)
# subset_ordered$species <- factor(subset_ordered$species, levels = rev(unique(subset_ordered$species))) 
# # Note: rev() keeps the highest VP at the top of the vertical list
# 
# # 2. Horizontal Guild Plot
# ggplot(subset_ordered, 
#        aes(x = foraging_guild_consensus, y = VP, fill = variable_clean)) +
#   # Using a boxplot to summarize the species within each guild
#   geom_boxplot(outlier.size = 1, alpha = 0.7) + 
#   
#   theme_minimal(base_size = 14) + 
#   scale_fill_manual(values = base_colors) +
#   
#   labs(
#     title = "% Explained Variance by Foraging Guild",
#     subtitle = "Summarized across species within each guild",
#     x = NULL, # Guild names are on the axis, so we don't need a title
#     y = "Explained variance (R²)",
#     fill = "Environmental Variable"
#   ) +
#   
#   theme(
#     # Rotate guild names so they don't overlap
#     axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
#     
#     # Legend and Grid cleanup
#     panel.grid.minor = element_blank(),
#     legend.position = 'top',
#     legend.justification = "left"
#   )
# 
# ggsave("misc-figures/env-vp-species-byguild-horizontal-box.png",
#        height = 10,width = 25)
# 

# boxplots, edian per guild -----------------------------------------------

# 1. Ensure levels are still set correctly
subset_ordered$foraging_guild_consensus <- factor(
  subset_ordered$foraging_guild_consensus, 
  levels = sorted_guilds
)

# 2. Horizontal Stacked Median Barplot
ggplot(subset_ordered, 
       aes(x = foraging_guild_consensus, y = VP, fill = variable_clean)) +
  
  # Calculate the median for each variable per guild and stack them
  stat_summary(fun = median, geom = "bar", position = "stack", color = "white", linewidth = 0.2) +
  
  theme_minimal(base_size = 14) + 
  scale_fill_manual(values = base_colors) +
  
  labs(
    title = "Median Explained Variance by Foraging Guild",
    subtitle = "Stacked median VP values across species",
    x = NULL,
    y = "Cumulative Median VP",
    fill = "Environmental Variable"
  ) +
  
  theme(
    # Rotate guild names
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
    
    # Legend and Grid cleanup
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(), # Remove vertical lines for cleaner bars
    legend.position = 'top',
    legend.justification = "left"
  )

ggsave(sprintf("misc-figures/%s-xfig-env-vp-median-byguild-vertical-stack-atlas-%s.png",pattern2match,atlas),
       height = 10,width = 25)



# MEDIAN GUILD TIME DOTPLOT -----------------------------------------
subset_ordered <- plot_data_clean %>%
  left_join(.,
            Tr %>% 
              rownames_to_column(var = 'species'),
            by = 'species') %>% 
  # Calculate total VP per species (so all rows for one species stay together)
  filter(!variable_clean == 'Site-level random effects') %>% 
  group_by(species) %>%
  mutate(total_species_VP = sum(VP, na.rm = TRUE)) %>%
  ungroup() %>%
  # Convert guild to factor using our custom levels
  mutate(foraging_guild_consensus = factor(foraging_guild_consensus, levels = sorted_guilds)) %>%
  # Arrange: 1st by the Guild levels, 2nd by Total VP (descending)
  # Change to total_species_VP (no desc) if you want lowest VP first
  arrange(foraging_guild_consensus, desc(total_species_VP)) %>%
  # Optional: Remove the helper column if you don't need it anymore
  select(-total_species_VP)

### SUMMARISE TRAIT FUNCTION 
summarise_trait_group <- function(
    df,
    trait = c("foraging_guild_consensus", "Migration_a3_DOF", "species_thermal_index"),
    stat = c("median", "mean"),
    min_species = 5
) {
  
  trait <- match.arg(trait)
  stat  <- match.arg(stat)
  
  # Add year label
  df <- df %>%
    mutate(
      yearlabel = factor(atlas, levels = 1:3, labels = c("1970s", "1990s", "2010s"))
    )
  
  # For STI, bin into groups first (it's continuous)
  if (trait == "species_thermal_index") {
    df <- df %>%
      mutate(trait_group = cut(species_thermal_index, breaks = 5, dig.lab = 2))
  } else {
    df <- df %>%
      mutate(trait_group = .data[[trait]])
  }
  
  sum_fn <- if (stat == "median") median else mean
  
  # Step 1: per species x atlas x variable — one VP value each (already the case)
  # Step 2: per trait_group x variable x atlas — average across species
  per_atlas <- df %>%
    group_by(trait_group, variable_clean, atlas, yearlabel) %>%
    summarise(
      n_species   = n_distinct(species),
      vp_per_atlas = sum_fn(VP, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_species >= min_species)
  
  # Step 3: across atlases — central tendency + variability
  result <- per_atlas %>%
    group_by(trait_group, variable_clean) %>%
    summarise(
      n_species_min    = min(n_species),        # most conservative species count
      n_species_max    = max(n_species),
      vp_central       = sum_fn(vp_per_atlas, na.rm = TRUE),   # median/mean across atlases
      sd_across_atlas  = sd(vp_per_atlas, na.rm = TRUE),       # temporal variability
      .groups = "drop"
    )
  
  # Step 4: variability across species (within trait group, pooling all atlases)
  species_var <- df %>%
    mutate(trait_group = if (trait == "species_thermal_index")
      cut(species_thermal_index, breaks = 5, dig.lab = 2)
      else .data[[trait]]) %>%
    group_by(trait_group, variable_clean, species) %>%
    summarise(vp_species = sum_fn(VP, na.rm = TRUE), .groups = "drop") %>%
    group_by(trait_group, variable_clean) %>%
    summarise(sd_across_species = sd(vp_species, na.rm = TRUE), .groups = "drop")
  
  # Join both variability measures
  result <- result %>%
    left_join(species_var, by = c("trait_group", "variable_clean"))
  
  return(result)
}

# Foraging guild, median, min 5 species
fg_summary <- summarise_trait_group(subset_ordered, 
                                    trait = "foraging_guild_consensus",
                                    stat = "mean",
                                    min_species = 5)

# plot_trait_variance <- function(summary_df, trait_label = "Foraging guild",
#                                 guild_order = row_order) {
#   
#   # Color palette for variables — adjust to match your base_colors
#   var_colors <- c(
#     "Site-level random effect" = "ivory2",
#     "Temperature"     = "firebrick3",
#     "Precipitation"   = "dodgerblue3",
#     "Land-use heterogeneity"   = "orchid3",
#     "Urban (% coverage)"           = "snow3",
#     "Cropland (% coverage)"        = "goldenrod1", # Slightly darker so it doesn't vanish
#     "Pasture (% coverage)"         = "darkorange",
#     "Forest (% coverage)"          = "springgreen4",
#     "Grass/Shrubland (% coverage)" = "springgreen2"
#   )
#   
#   # Force panel order via factor levels
#   if (!is.null(guild_order)) {
#     summary_df <- summary_df %>%
#       mutate(trait_group = factor(trait_group, levels = guild_order))
#   }
#   
#   ggplot(summary_df,
#          aes(x     = sd_across_atlas,
#              y     = sd_across_species,
#              color = variable_clean,
#              size  = vp_central)) +
#     geom_point(alpha = 0.85, stroke = 0.3) +
#     scale_color_manual(values = var_colors) +
#     scale_size_continuous(
#       range  = c(2, 12),
#       name   = "Median VP (R²)"
#     ) +
#     facet_wrap(~ trait_group, scales = "fixed") +
#     labs(
#       title    = paste("Variance partitioning —", trait_label),
#       subtitle = "x: temporal variability (SD across atlases) · y: species variability (SD across species)",
#       x        = "Temporal variability (SD across atlases)",
#       y        = "Species variability (SD across species)",
#       color    = "Variable"
#     ) +
#     theme_minimal(base_size = 13) +
#     theme(
#       panel.grid.minor  = element_blank(),
#       panel.grid.major  = element_line(linewidth = 0.3, color = "grey90"),
#       legend.position   = "bottom",
#       legend.box        = "vertical",
#       strip.text        = element_text(face = "bold", size = 11),
#       axis.text         = element_text(size = 9),
#       plot.subtitle     = element_text(size = 10, color = "grey50")
#     ) +
#     guides(
#       color = guide_legend(override.aes = list(size = 4), nrow = 2),
#       size  = guide_legend(nrow = 1)
#     )
# }
# 
# # Usage
# fg_summary <- summarise_trait_group(subset_ordered,
#                                     trait       = "foraging_guild_consensus",
#                                     stat        = "mean",
#                                     min_species = 5)
# unique(fg_summary$trait_group)
# row_order
# plot_trait_variance(fg_summary, trait_label = "Foraging guild",
#                     guild_order = row_order)
# 
# ggsave(sprintf('misc-figures/%s-xfig-guilds-guilds-vp-points-variability.png', pattern), 
#        width = 7.5, height = 7.5)


# TRY TILE INSTEAD --------------------------------------------------------
plot_trait_tile <- function(summary_df,
                            trait_label  = "Foraging guild",
                            guild_order  = NULL,
                            var_order    = NULL,
                            var_colors   = NULL,
                            sd_type      = c("species", "temporal")) {  # <-- new
  library(scales)
  library(purrr)
  
  sd_type <- match.arg(sd_type)
  sd_col  <- if (sd_type == "species") "sd_across_species" else "sd_across_atlas"
  sd_label <- if (sd_type == "species") "intra-guild variability" else "temporal variability"
  
  mean_mat <- summary_df %>%
    select(trait_group, variable_clean, vp_central) %>%
    pivot_wider(names_from = variable_clean, values_from = vp_central) %>%
    column_to_rownames("trait_group")
  
  sd_mat <- summary_df %>%
    select(trait_group, variable_clean, all_of(sd_col)) %>%
    pivot_wider(names_from = variable_clean, values_from = all_of(sd_col)) %>%
    column_to_rownames("trait_group")
  
  row_order <- if (!is.null(guild_order)) guild_order[guild_order %in% rownames(mean_mat)] else rownames(mean_mat)
  col_order <- if (!is.null(var_order))   var_order[var_order %in% colnames(mean_mat)]     else colnames(mean_mat)
  
  mean_long <- mean_mat[row_order, col_order] %>%
    as.data.frame() %>%
    rownames_to_column("trait_group") %>%
    pivot_longer(-trait_group, names_to = "variable_clean", values_to = "vp_central")
  
  sd_long <- sd_mat[row_order, col_order] %>%
    as.data.frame() %>%
    rownames_to_column("trait_group") %>%
    pivot_longer(-trait_group, names_to = "variable_clean", values_to = "sd")
  
  plot_df <- left_join(mean_long, sd_long, by = c("trait_group", "variable_clean")) %>%
    mutate(
      trait_group    = factor(trait_group,    levels = rev(row_order)),
      variable_clean = factor(variable_clean, levels = col_order)
    )
  
  plot_df <- plot_df %>%
    group_by(variable_clean) %>%
    mutate(vp_scaled = rescale(vp_central, to = c(0, 1))) %>%
    ungroup()
  
  if (is.null(var_colors)) {
    vars <- levels(plot_df$variable_clean)
    default_pal <- c("firebrick3", "dodgerblue3", "orchid3", "snow3",
                     "goldenrod1", "darkorange", "springgreen4", "springgreen2", "ivory2")
    var_colors <- setNames(default_pal[seq_along(vars)], vars)
  }
  
  plot_df <- plot_df %>%
    mutate(
      var_hex    = var_colors[as.character(variable_clean)],
      fill_color = map2_chr(var_hex, vp_scaled,
                            ~ colorRampPalette(c("white", .x))(101)[round(.y * 100) + 1]),
      alpha_val  = 1 - rescale(sd, to = c(0.15, 0.9))
    )
  
  ggplot(plot_df, aes(x = variable_clean, y = trait_group)) +
    geom_tile(aes(fill = fill_color, alpha = alpha_val),
              color = "white", linewidth = 0.5) +
    scale_fill_identity() +
    scale_alpha_identity(
      name  = paste0("certainty\n(1 – scaled ", sd_label, ")"),
      guide = guide_legend()
    ) +
    scale_x_discrete(position = "top") +
    labs(title    = paste("Variance partitioning —", trait_label),
         subtitle = paste("Transparency encodes", sd_label),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 0, vjust = 0),
      axis.text.y     = element_text(size = 10),
      panel.grid      = element_blank(),
      legend.position = "right",
      plot.title      = element_text(size = 13, face = "bold"),
      plot.subtitle   = element_text(size = 10, color = "grey50")
    )
}

fg_summary <- summarise_trait_group(subset_ordered,
                                    trait       = "foraging_guild_consensus",
                                    stat        = "median",
                                    min_species = 5)
fg_summary$variable_clean <- gsub(' (% coverage)','',fg_summary$variable_clean,
                                  fixed = T)
var_colors <- c(
  "Site-level random effect" = "ivory2",
  "Temperature"     = "firebrick3",
  "Precipitation"   = "dodgerblue3",
  "Land-use heterogeneity"   = "orchid3",
  "Urban"           = "snow3",
  "Cropland"        = "goldenrod1", # Slightly darker so it doesn't vanish
  "Pasture"         = "darkorange",
  "Forest"          = "springgreen4",
  "Grass/Shrubland" = "springgreen2"
)

col_order <- names(var_colors)[names(var_colors) %in% unique(fg_summary$variable_clean)]

clean_variable_names <- function(x) {
  case_match(x,
             "tmean_breeding"       ~ "Temperature",
             "prec_breeding"        ~ "Precipitation",
             "hh"                   ~ "Land-use heterogeneity",
             "perc_urban"           ~ "Urban",
             "perc_cropland"        ~ "Cropland",
             "perc_pasture"         ~ "Pasture",
             "perc_forest"          ~ "Forest",
             "perc_grass_shrub"     ~ "Grass/Shrubland",
             "Random: site"         ~ "Site-level random effects",
             .default = x)
}

summarise_trait_direction <- function(df,
                                      trait = "foraging_guild_consensus",
                                      stat = c("median", "mean"),
                                      min_species = 5,
                                      pattern = pattern2match,
                                      dominance_threshold = 0.1) {
  stat <- match.arg(stat)
  sum_fn <- if (stat == "median") median else mean
  
  beta <- read_parameter_effects(pattern, effect = "Beta") %>%
    filter(variable != "(Intercept)", !is.na(effect_size)) %>%
    mutate(
      atlas = as.character(atlas),
      variable_clean = clean_variable_names(variable),
      effect_sign = sign(effect_size)
    ) %>%
    select(species, atlas, variable_clean, effect_sign)
  
  direction_df <- df %>%
    mutate(
      atlas = as.character(atlas),
      trait_group = .data[[trait]],
      variable_clean = as.character(variable_clean)
    ) %>%
    filter(!is.na(trait_group), !variable_clean %in% c("Site-level random effects", "Site-level random effect")) %>%
    mutate(variable_clean = gsub(" (% coverage)", "", variable_clean, fixed = TRUE)) %>%
    left_join(beta, by = c("species", "atlas", "variable_clean")) %>%
    group_by(trait_group, variable_clean, atlas) %>%
    summarise(
      n_species = n_distinct(species),
      vp_per_atlas = sum_fn(VP, na.rm = TRUE),
      signed_vp = sum(VP * effect_sign, na.rm = TRUE),
      supported_vp = sum(ifelse(is.na(effect_sign), 0, VP), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_species >= min_species)
  
  direction_df %>%
    group_by(trait_group, variable_clean) %>%
    summarise(
      direction_score = ifelse(
        sum(supported_vp, na.rm = TRUE) > 0,
        sum(signed_vp, na.rm = TRUE) / sum(supported_vp, na.rm = TRUE),
        NA_real_
      ),
      supported_vp_share = ifelse(
        sum(vp_per_atlas, na.rm = TRUE) > 0,
        sum(supported_vp, na.rm = TRUE) / sum(vp_per_atlas, na.rm = TRUE),
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      effect_direction = case_when(
        is.na(direction_score) ~ "Not supported",
        direction_score >= dominance_threshold ~ "Higher at max",
        direction_score <= -dominance_threshold ~ "Lower at max",
        TRUE ~ "Mixed"
      ),
      effect_direction = factor(
        effect_direction,
        levels = c("Higher at max", "Lower at max", "Mixed", "Not supported")
      )
    )
}

fg_direction <- summarise_trait_direction(subset_ordered,
                                          trait = "foraging_guild_consensus",
                                          stat = "median",
                                          min_species = 5,
                                          pattern = pattern2match)

fg_summary_directed <- fg_summary %>%
  left_join(fg_direction, by = c("trait_group", "variable_clean")) %>%
  mutate(
    direction_score = replace_na(direction_score, 0),
    signed_vp_central = vp_central * direction_score
  )

pca_order_trait_groups <- function(summary_df, value_col = "signed_vp_central",
                                   var_order = NULL, fallback_order = NULL) {
  input <- summary_df %>%
    select(trait_group, variable_clean, value = all_of(value_col)) %>%
    mutate(value = replace_na(value, 0)) %>%
    pivot_wider(names_from = variable_clean, values_from = value, values_fill = 0) %>%
    column_to_rownames("trait_group")
  
  if (!is.null(var_order)) {
    input <- input[, var_order[var_order %in% colnames(input)], drop = FALSE]
  }
  
  nonzero_cols <- vapply(input, sd, numeric(1), na.rm = TRUE) > 0
  input <- input[, nonzero_cols, drop = FALSE]
  
  if (nrow(input) < 2 || ncol(input) < 2) {
    return(fallback_order[fallback_order %in% rownames(input)])
  }
  
  pca <- prcomp(input, center = TRUE, scale. = TRUE)
  n_pc <- min(3, ncol(pca$x))
  pc_scores <- pca$x[, seq_len(n_pc), drop = FALSE]
  clustering <- hclust(dist(pc_scores), method = "ward.D2")
  
  clustering$labels[clustering$order]
}

row_order <- pca_order_trait_groups(
  fg_summary_directed,
  value_col = "signed_vp_central",
  var_order = col_order,
  fallback_order = sorted_guilds
)
target_guilds <- row_order

# Let clustering decide row order
plot_trait_tile(fg_summary,
                trait_label = "Foraging guild",
                guild_order = row_order,
                var_order   = col_order,
                var_colors  = var_colors,
                sd_type      = 'species')


# BUBBLE ------------------------------------------------------------------
plot_trait_bubble <- function(summary_df,
                              trait_label  = "Foraging guild",
                              guild_order  = NULL,
                              var_order    = NULL,
                              var_colors   = NULL,
                              sd_type      = c("species", "temporal"),
                              show_direction = FALSE) {
  library(scales)
  library(purrr)
  
  sd_type  <- match.arg(sd_type)
  sd_col   <- if (sd_type == "species") "sd_across_species" else "sd_across_atlas"
  sd_label <- if (sd_type == "species") "intra-guild variability" else "temporal variability"
  
  # Apply row/col ordering
  row_order <- if (!is.null(guild_order)) guild_order[guild_order %in% summary_df$trait_group] else unique(summary_df$trait_group)
  col_order <- if (!is.null(var_order))   var_order[var_order %in% summary_df$variable_clean]  else unique(summary_df$variable_clean)
  
  # Default colors
  if (is.null(var_colors)) {
    vars <- col_order
    default_pal <- c("firebrick3", "dodgerblue3", "orchid3", "snow3",
                     "goldenrod1", "darkorange", "springgreen4", "springgreen2", "ivory2")
    var_colors <- setNames(default_pal[seq_along(vars)], vars)
  }
  
  plot_df <- summary_df %>%
    filter(trait_group %in% row_order, variable_clean %in% col_order) %>%
    rename(sd = all_of(sd_col)) %>%
    mutate(
      alpha_val      = 1 - rescale(sd, to = c(0.15, 0.9)),
      trait_group    = factor(trait_group,    levels = rev(row_order)),
      variable_clean = factor(variable_clean, levels = col_order),
      fill_color     = var_colors[as.character(variable_clean)],
      effect_direction = if ("effect_direction" %in% names(.)) {
        fct_na_value_to_level(effect_direction, level = "Not supported")
      } else {
        factor("Not supported", levels = c("Higher at max", "Lower at max", "Mixed", "Not supported"))
      }
    )
  
  p <- ggplot(plot_df, aes(x = variable_clean, y = trait_group))
  
  if (show_direction) {
    p <- p +
      geom_point(
        aes(
          size = vp_central,
          color = variable_clean,
          alpha = alpha_val,
          shape = effect_direction
        ),
        stroke = 1.1
      )
  } else {
    p <- p +
      geom_point(aes(size = vp_central, color = variable_clean, alpha = alpha_val),
                 shape = 16) +
      geom_point(aes(size = vp_central, color = variable_clean),
                 shape = 21, fill = NA, stroke = 0.4, alpha = 0.6)
  }
  
  p <- p +
    scale_color_manual(values = var_colors, guide = "none") +
    scale_size_continuous(
      range  = c(1, 12),
      name   = "Median VP (R²)",
      labels = label_number(accuracy = 0.01)
    ) +
    scale_alpha_identity(
      name  = paste0("certainty\n(1 – scaled\n", sd_label, ")"),
      guide = guide_legend()
    )

  if (show_direction) {
    p <- p +
    scale_shape_manual(
      values = c("Higher at max" = 16, "Lower at max" = 4, "Mixed" = 1, "Not supported" = 3),
      name = "Direction at\nhigh predictor value",
      drop = FALSE
    )
  }

  p <- p +
    scale_x_discrete(position = "top") +
    labs(
      title    = paste("Variance partitioning —", trait_label),
      subtitle = if (show_direction) {
        paste("Size: median VP  ·  Transparency:", sd_label, " ·  Symbol: supported Beta direction")
      } else {
        paste("Size: median VP  ·  Transparency:", sd_label)
      },
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 0, vjust = 0),
      axis.text.y      = element_text(size = 10),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      plot.title       = element_text(size = 13, face = "bold"),
      plot.subtitle    = element_text(size = 10, color = "grey50")
    )
  
  p
}

plot_trait_bubble(fg_summary,
                  trait_label = "Foraging guild",
                  guild_order = row_order,
                  var_order   = col_order,
                  var_colors  = var_colors,
                  sd_type     = "temporal")
ggsave(sprintf("misc-figures/%s-env-vp-guild-bubble-temporalvariabaility.png",pattern2match))

plot_trait_bubble(fg_summary_directed,
                  trait_label = "Foraging guild",
                  guild_order = row_order,
                  var_order   = col_order,
                  var_colors  = var_colors,
                  sd_type     = "temporal",
                  show_direction = TRUE)
ggsave(sprintf("misc-figures/%s-env-vp-guild-bubble-temporalvariabaility-direction.png", pattern2match),
       width = 9, height = 7.5, dpi = 300)

# GUILD - ENVIRONMENT - DOTPLOT - ATLAS SEPEARTE --------------------------

subset_ordered <- plot_data_clean %>%
  left_join(.,
            Tr %>% 
              rownames_to_column(var = 'species'),
            by = 'species') %>% 
  filter(!variable_clean == 'Site-level random effects') %>% 
  group_by(species) %>%
  mutate(total_species_VP = sum(VP, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(foraging_guild_consensus = factor(foraging_guild_consensus, levels = sorted_guilds)) %>%
  arrange(foraging_guild_consensus, desc(total_species_VP)) %>%
  select(-total_species_VP)

# Strip coverage label
subset_ordered$variable_clean <- gsub(' (% coverage)', '', subset_ordered$variable_clean, fixed = TRUE)

VP_mean <- subset_ordered %>% 
  group_by(atlas,variable_clean,foraging_guild_consensus) %>% 
  summarize(mean_vp = mean(VP)) %>% 
  filter(foraging_guild_consensus%in%target_guilds) %>% 
  mutate(foraging_guild_consensus = factor(foraging_guild_consensus,
                                              levels = rev(row_order)),
         variable_clean = factor(variable_clean,
                                    levels=col_order),
         atlas = factor(atlas,
                        levels = c(1,2,3)))


ggplot(VP_mean,
       aes(x = atlas,
           y = foraging_guild_consensus,
           group = variable_clean,
           size = mean_vp,
           fill = variable_clean)) +

    geom_point(shape = 21) +
  scale_fill_manual(values = var_colors, guide = "none") +

  facet_grid(~variable_clean) +
  theme_minimal()
ggsave(sprintf("misc-figures/%s-env-vp-guild-bubble-temporalvariabaility-atlas-seperated.png",pattern2match))
