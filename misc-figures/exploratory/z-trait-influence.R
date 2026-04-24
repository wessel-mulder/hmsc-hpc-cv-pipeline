rm(list=ls())
#### GETTING STARTED ####
if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,Hmsc,RColorBrewer,ggplot2,
               rnaturalearth,rnaturalearthdata,
               gridExtra,patchwork,sf,cowplot,
               terra,packcircles,readxl)

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


# INFLUENCE OF TRAITS  ----------------------------------------------------
m <- mods[[1]]
m
VPs <- lapply(mods,computeVariancePartitioning)

b <- VPs[[1]]$R2T$Beta
b1 <- as.data.frame(b) %>% rename('value'='b') %>% rownames_to_column(var='')

# Generate the layout. This function return a dataframe with one line per bubble. 
# It gives its center (x and y) and its radius, proportional of the value
packing <- circleProgressiveLayout(b1, sizetype='area')

# We can add these packing information to the initial data frame
data <- cbind(b1, packing)

# Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
# plot(data$radius, data$value)

# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=50)

# Make the plot
ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill=as.factor(id)), colour = "black", alpha = 0.6) +
  
  # Add text in the center of each bubble + control its size
  geom_text(data = b, aes(x, y, size=value, label = group)) +
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

##### FIGURE OUT TRAIT PLOTS  -----------------------------------------------------------
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

gamma <- make_effect_sizes('HmscOutputs','2026-03-13',effect='Gamma') %>% 
  drop_na()

print(gamma,n=100)

# ── Species counts from TrData ─────────────────────────────────────────────────

tr <- mods[[1]]$TrData

migration_keep <- names(table(tr$Migration_a3_DOF))[table(tr$Migration_a3_DOF) >= 5]
foraging_keep  <- names(table(tr$foraging_guild_consensus))[table(tr$foraging_guild_consensus) >= 5]


gamma_heatmap <- gamma %>%
  mutate(
    trait_category = case_when(
      str_starts(traits, "Migration")      ~ "Migration",
      str_starts(traits, "foraging_guild") ~ "Foraging guild",
      str_starts(traits, "species_therm")  ~ "Thermal index",
      TRUE                                 ~ NA_character_
    ),
    trait_clean = traits %>%
      str_remove("Migration_a3_DOF") %>%
      str_remove("foraging_guild_consensus") %>%
      str_trim()
  ) %>%
  # Remove intercept
  filter(traits != "(Intercept)") %>%
  filter(!is.na(trait_category)) %>%
  # Filter small categories
  filter(
    (trait_category == "Migration"      & trait_clean %in% migration_keep) |
      (trait_category == "Foraging guild" & trait_clean %in% foraging_keep)  |
      (trait_category == "Thermal index")
  ) %>%
  mutate(
    direction = case_when(
      is.na(effect_size)  ~ "non-significant",
      effect_size > 0     ~ "positive",
      effect_size < 0     ~ "negative"
    ),
    atlas = factor(atlas, levels = c(1, 2, 3))
  )

# ── 2. Cluster traits and environments ─────────────────────────────────────────

cluster_mat <- gamma_heatmap %>%
  group_by(trait_clean, variable) %>%
  summarise(mean_effect = mean(effect_size, na.rm = TRUE), .groups = "drop") %>%
  mutate(mean_effect = replace_na(mean_effect, 0)) %>%
  pivot_wider(names_from = variable, values_from = mean_effect, values_fill = 0) %>%
  tibble::column_to_rownames("trait_clean")

trait_order <- rownames(cluster_mat)[hclust(dist(cluster_mat))$order]
env_order   <- colnames(cluster_mat)[hclust(dist(t(cluster_mat)))$order]

gamma_heatmap <- gamma_heatmap %>%
  mutate(
    trait_clean = factor(trait_clean, levels = trait_order),
    variable    = factor(variable, levels = env_order)
  )

# ── 3. Triangle geometry ───────────────────────────────────────────────────────

make_triangles <- function(df) {
  df %>%
    mutate(
      x = as.numeric(variable),
      y = as.numeric(trait_clean)
    ) %>%
    group_by(trait_clean, variable, atlas) %>%
    summarise(
      x              = first(x),
      y              = first(y),
      direction      = first(direction),
      trait_category = first(trait_category),
      .groups        = "drop"
    ) %>%
    mutate(h = 0.5) %>%
    mutate(
      x1 = case_when(atlas == 1 ~ x - h, atlas == 2 ~ x,     atlas == 3 ~ x - h),
      y1 = case_when(atlas == 1 ~ y + h, atlas == 2 ~ y + h, atlas == 3 ~ y - h),
      x2 = case_when(atlas == 1 ~ x - h, atlas == 2 ~ x + h, atlas == 3 ~ x + h),
      y2 = case_when(atlas == 1 ~ y - h, atlas == 2 ~ y + h, atlas == 3 ~ y - h),
      x3 = x,
      y3 = y
    )
}

tri      <- make_triangles(gamma_heatmap)

tri_long <- tri %>%
  mutate(group = paste(trait_clean, variable, atlas, sep = "_")) %>%
  pivot_longer(cols = c(x1, x2, x3), names_to = "point", values_to = "px") %>%
  mutate(py = case_when(
    point == "x1" ~ y1,
    point == "x2" ~ y2,
    point == "x3" ~ y3
  ))

# ── 4. Category strip data ─────────────────────────────────────────────────────

n_envs <- length(levels(gamma_heatmap$variable))

category_labels <- gamma_heatmap %>%
  distinct(trait_clean, trait_category) %>%
  mutate(y = as.numeric(trait_clean)) %>%
  group_by(trait_category) %>%
  summarise(
    ymin = min(y) - 0.5,
    ymax = max(y) + 0.5,
    ymid = mean(y),
    .groups = "drop"
  )

# ── 5. Colours ─────────────────────────────────────────────────────────────────

direction_colors <- c(
  "positive"        = "#378ADD",
  "negative"        = "#D85A30",
  "non-significant" = "grey92"
)

category_colors <- c(
  "Migration"      = "#7F77DD",
  "Foraging guild" = "#D85A30",
  "Thermal index"  = "#3B6D11"
)

all_colors <- c(direction_colors, category_colors)

# ── 6. Plot ────────────────────────────────────────────────────────────────────

p <- ggplot() +
  
  # Triangles
  geom_polygon(
    data = tri_long,
    aes(x = px, y = py, group = group, fill = direction),
    color = "white", linewidth = 0.3
  ) +
  
  # Cell borders
  geom_rect(
    data = gamma_heatmap %>%
      distinct(trait_clean, variable) %>%
      mutate(x = as.numeric(variable), y = as.numeric(trait_clean)),
    aes(xmin = x - 0.5, xmax = x + 0.5, ymin = y - 0.5, ymax = y + 0.5),
    fill = NA, color = "grey80", linewidth = 0.3
  ) +
  
  # Category colour strip
  geom_rect(
    data = category_labels,
    aes(
      xmin = n_envs + 0.6, xmax = n_envs + 0.9,
      ymin = ymin, ymax = ymax,
      fill = trait_category
    ),
    color = NA
  ) +
  
  # Category text
  geom_text(
    data = category_labels,
    aes(
      x = n_envs + 1.1, y = ymid,
      label = trait_category,
      color = trait_category
    ),
    hjust = 0, size = 3, fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = all_colors,
    breaks = names(direction_colors),
    labels = c("Positive", "Negative", "Non-significant"),
    name   = "Effect direction"
  ) +
  scale_color_manual(values = category_colors, guide = "none") +
  
  scale_x_continuous(
    breaks = seq_along(levels(gamma_heatmap$variable)),
    labels = levels(gamma_heatmap$variable),
    expand = expansion(add = c(0.5, 3))
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(gamma_heatmap$trait_clean)),
    labels = levels(gamma_heatmap$trait_clean),
    expand = expansion(add = 0.5)
  ) +
  
  coord_fixed() +
  
  annotate(
    "text",
    x = 0.5,
    y = length(levels(gamma_heatmap$trait_clean)) + 1.2,
    label = "◤ Atlas 1     ◥ Atlas 2     ▼ Atlas 3",
    hjust = 0, size = 3, color = "grey40"
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 8),
    axis.title      = element_blank(),
    panel.grid      = element_blank(),
    legend.position = "bottom",
    plot.margin     = margin(10, 100, 10, 10)
  )

print(p)


# V2 ----------------------------------------------------------------------


library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

# ── 1. Variable rename + colour mapping ────────────────────────────────────────

var_rename <- c(
  "tmean_breeding"   = "Temperature",
  "prec_breeding"    = "Precipitation",
  "hh"               = "Land-use heterogeneity",
  "perc_urban"       = "Urban (% coverage)",
  "perc_cropland"    = "Cropland (% coverage)",
  "perc_pasture"     = "Pasture (% coverage)",
  "perc_forest"      = "Forest (% coverage)",
  "perc_grass_shrub" = "Grass/Shrubland (% coverage)"
)

base_colors <- c(
  "Temperature"                  = "firebrick3",
  "Precipitation"                = "dodgerblue3",
  "Land-use heterogeneity"       = "orchid3",
  "Urban (% coverage)"           = "snow3",
  "Cropland (% coverage)"        = "goldenrod1",
  "Pasture (% coverage)"         = "darkorange",
  "Forest (% coverage)"          = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2"
)

atlas_rename <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")

# ── 2. Prepare data ────────────────────────────────────────────────────────────

tr <- mods[[1]]$TrData
migration_keep <- names(table(tr$Migration_a3_DOF))[table(tr$Migration_a3_DOF) >= 5]
foraging_keep  <- names(table(tr$foraging_guild_consensus))[table(tr$foraging_guild_consensus) >= 5]

gamma_plot <- gamma %>%
  filter(traits != "(Intercept)") %>%
  filter(variable != "(Intercept)") %>% 
  mutate(
    trait_category = case_when(
      str_starts(traits, "Migration")      ~ "Migration",
      str_starts(traits, "foraging_guild") ~ "Foraging guild",
      str_starts(traits, "species_therm")  ~ "Thermal index",
      TRUE                                 ~ NA_character_
    ),
    trait_clean = traits %>%
      str_remove("Migration_a3_DOF") %>%
      str_remove("foraging_guild_consensus") %>%
      str_trim()
  ) %>%
  filter(!is.na(trait_category)) %>%
  filter(
    (trait_category == "Migration"      & trait_clean %in% migration_keep) |
    (trait_category == "Foraging guild" & trait_clean %in% foraging_keep)  |
    (trait_category == "Thermal index")
  ) %>%
  filter(!is.na(effect_size)) %>%
  mutate(
    variable    = var_rename[variable],
    variable    = factor(variable, levels = names(base_colors)),
    direction   = ifelse(effect_size > 0, "Positive", "Negative"),
    atlas       = factor(atlas_rename[as.character(atlas)],
                         levels = c("1970s", "1990s", "2010s")),
    # Trait category order: Thermal index first, then Foraging guild, then Migration
    trait_category = factor(trait_category,
                            levels = c("Thermal index", "Foraging guild", "Migration")),
    trait_clean = factor(trait_clean)
  )

# ── 3. Cluster traits within each category ────────────────────────────────────

cluster_within_category <- function(df) {
  cats <- levels(df$trait_category)
  ordered_traits <- c()
  
  for (cat in cats) {
    sub <- df %>%
      filter(trait_category == cat) %>%
      group_by(trait_clean, variable) %>%
      summarise(mean_effect = mean(effect_size, na.rm = TRUE), .groups = "drop") %>%
      mutate(mean_effect = replace_na(mean_effect, 0)) %>%
      pivot_wider(names_from = variable, values_from = mean_effect, values_fill = 0) %>%
      tibble::column_to_rownames("trait_clean")
    
    if (nrow(sub) > 1) {
      ordered_traits <- c(ordered_traits,
                          rownames(sub)[hclust(dist(sub))$order])
    } else {
      ordered_traits <- c(ordered_traits, rownames(sub))
    }
  }
  ordered_traits
}

trait_order <- cluster_within_category(gamma_plot)

# Cluster environments
env_mat <- gamma_plot %>%
  group_by(trait_clean, variable) %>%
  summarise(mean_effect = mean(effect_size, na.rm = TRUE), .groups = "drop") %>%
  mutate(mean_effect = replace_na(mean_effect, 0)) %>%
  pivot_wider(names_from = variable, values_from = mean_effect, values_fill = 0) %>%
  tibble::column_to_rownames("trait_clean")

env_order <- colnames(env_mat)[hclust(dist(t(env_mat)))$order]

gamma_plot <- gamma_plot %>%
  mutate(
    trait_clean = factor(trait_clean, levels = trait_order),
    variable    = factor(variable, levels = env_order)
  )

present_vars <- unique(gamma_plot$variable)
base_colors_sub <- base_colors[names(base_colors) %in% present_vars]

# ── 4. Category divider positions ─────────────────────────────────────────────

category_dividers <- gamma_plot %>%
  distinct(trait_clean, trait_category) %>%
  mutate(y = as.numeric(trait_clean)) %>%
  group_by(trait_category) %>%
  summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5, ymid = mean(y), .groups = "drop")

n_envs <- length(levels(gamma_plot$variable))

# ── 5. Circle sizes per atlas ─────────────────────────────────────────────────
# Atlas 1 = small filled circle
# Atlas 2 = medium ring
# Atlas 3 = large ring

atlas_sizes <- c("1970s" = 2.5, "1990s" = 5.5, "2010s" = 9)
atlas_fills <- c("1970s" = "filled", "1990s" = "ring", "2010s" = "ring")

# Split into filled (atlas 1) and rings (atlas 2, 3)
gamma_filled <- gamma_plot %>% filter(atlas == "1970s")
gamma_rings  <- gamma_plot %>% filter(atlas != "1970s")

# ── 6. Category colour strip ──────────────────────────────────────────────────

category_colors <- c(
  "Thermal index"  = "#3B6D11",
  "Foraging guild" = "#D85A30",
  "Migration"      = "#7F77DD"
)

# ── 7. Plot ───────────────────────────────────────────────────────────────────

make_panel <- function(direction_filter, panel_title) {
  
  gf <- gamma_filled %>% filter(direction == direction_filter) %>% droplevels()
  gr <- gamma_rings  %>% filter(direction == direction_filter) %>% droplevels()
  
  # Also resubset base_colors to only present variables
  present_vars <- levels(gf$variable) %>% c(levels(gr$variable)) %>% unique()
  present_vars <- present_vars[present_vars %in% names(base_colors)]
  base_colors_panel <- base_colors[present_vars]
  
  ggplot() +
    
    
    
    # Background cells
    geom_rect(
      data = gamma_plot %>%
        distinct(trait_clean, variable) %>%
        mutate(x = as.numeric(variable), y = as.numeric(trait_clean)),
      aes(xmin = x - 0.5, xmax = x + 0.5, ymin = y - 0.5, ymax = y + 0.5),
      fill = "grey97", color = "grey88", linewidth = 0.2
    ) +
    
    # Category divider lines
    geom_hline(
      data = category_dividers %>% filter(trait_category != last(levels(gamma_plot$trait_category))),
      aes(yintercept = ymax),
      color = "grey60", linewidth = 0.4, linetype = "dashed"
    ) +
    
    # Atlas 3 ring (largest, draw first so it's behind)
geom_point(
  data = gr %>% filter(atlas == "2010s"),
  aes(x = as.numeric(variable), y = as.numeric(trait_clean), color = variable),
  size = atlas_sizes["2010s"], shape = 21,
  fill = NA, stroke = 1.8
) +
  
  # Atlas 2 ring (medium)
  geom_point(
    data = gr %>% filter(atlas == "1990s"),
    aes(x = as.numeric(variable), y = as.numeric(trait_clean), color = variable),
    size = atlas_sizes["1990s"], shape = 21,
    fill = NA, stroke = 1.5
  ) +
  
  # Atlas 1 filled circle (smallest, on top)
  geom_point(
    data = gf,
    aes(x = as.numeric(variable), y = as.numeric(trait_clean),
        color = variable, fill = variable),
    size = atlas_sizes["1970s"], shape = 21, stroke = 0.5
  ) +
  
  # Category colour strip on right
  geom_rect(
    data = category_dividers,
    aes(
      xmin = n_envs + 0.6, xmax = n_envs + 0.85,
      ymin = ymin, ymax = ymax,
      fill = trait_category
    ),
    color = NA
  ) +
  
  # Category labels
  geom_text(
    data = category_dividers,
    aes(x = n_envs + 1.05, y = ymid, label = trait_category,
        color = trait_category),
    hjust = 0, size = 2.8, fontface = "bold"
  ) +
  
    # Then in make_panel, replace scale_color_manual and scale_fill_manual with:
    scale_color_manual(
      values = c(base_colors_panel, category_colors),
      breaks = names(base_colors_panel),
      name   = "Environment"
    ) +
    scale_fill_manual(
      values = c(base_colors_panel, category_colors),
      breaks = names(base_colors_panel),
      name   = "Environment"
    ) +
  
  scale_x_continuous(
    breaks = seq_along(levels(gamma_plot$variable)),
    labels = levels(gamma_plot$variable),
    expand = expansion(add = c(0.5, 3.5))
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(gamma_plot$trait_clean)),
    labels = levels(gamma_plot$trait_clean),
    expand = expansion(add = 0.5)
  ) +
  
  coord_fixed() +
  ggtitle(panel_title) +
  
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 8),
    axis.title      = element_blank(),
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", size = 11, hjust = 0.5),
    legend.position = "bottom",
    plot.margin     = margin(10, 80, 10, 10)
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      override.aes = list(size = 3, shape = 16, fill = unname(base_colors))
    ),
    fill = "none"
  )
}

p_pos <- make_panel("Positive", "Positive effects")
p_neg <- make_panel("Negative", "Negative effects")

# Combine panels
library(patchwork)

p_combined <- p_pos + p_neg +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(p_combined)


# V3 ----------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(patchwork)

# ── 1. Variable rename + colour mapping ───────────────────────────────────────

var_rename <- c(
  "tmean_breeding"   = "Temperature",
  "prec_breeding"    = "Precipitation",
  "hh"               = "Land-use heterogeneity",
  "perc_urban"       = "Urban (% coverage)",
  "perc_cropland"    = "Cropland (% coverage)",
  "perc_pasture"     = "Pasture (% coverage)",
  "perc_forest"      = "Forest (% coverage)",
  "perc_grass_shrub" = "Grass/Shrubland (% coverage)"
)

base_colors <- c(
  "Temperature"                  = "firebrick3",
  "Precipitation"                = "dodgerblue3",
  "Land-use heterogeneity"       = "orchid3",
  "Urban (% coverage)"           = "snow3",
  "Cropland (% coverage)"        = "goldenrod1",
  "Pasture (% coverage)"         = "darkorange",
  "Forest (% coverage)"          = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2"
)

category_colors <- c(
  "Thermal index"  = "#3B6D11",
  "Foraging guild" = "#D85A30",
  "Migration"      = "#7F77DD"
)

atlas_rename <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")
atlas_sizes  <- c("1970s" = 2.5, "1990s" = 5.5, "2010s" = 9)

# ── 2. Prepare data ───────────────────────────────────────────────────────────

tr             <- mods[[1]]$TrData
migration_keep <- names(table(tr$Migration_a3_DOF))[table(tr$Migration_a3_DOF) >= 5]
foraging_keep  <- names(table(tr$foraging_guild_consensus))[table(tr$foraging_guild_consensus) >= 5]

gamma_plot <- gamma %>%
  filter(traits != "(Intercept)", variable != "(Intercept)") %>%
  mutate(
    trait_category = case_when(
      str_starts(traits, "Migration")      ~ "Migration",
      str_starts(traits, "foraging_guild") ~ "Foraging guild",
      str_starts(traits, "species_therm")  ~ "Thermal index",
      TRUE                                 ~ NA_character_
    ),
    trait_clean = traits %>%
      str_remove("Migration_a3_DOF") %>%
      str_remove("foraging_guild_consensus") %>%
      str_trim(),
    trait_clean = ifelse(trait_clean == "species_thermal_index", "Thermal index", trait_clean)
  ) %>%
  filter(!is.na(trait_category)) %>%
  filter(
    (trait_category == "Migration"      & trait_clean %in% migration_keep) |
      (trait_category == "Foraging guild" & trait_clean %in% foraging_keep)  |
      (trait_category == "Thermal index")
  ) %>%
  filter(!is.na(effect_size)) %>%
  # Recode variable BEFORE making factor, then drop NAs
  mutate(variable = var_rename[variable]) %>%
  filter(!is.na(variable)) %>%
  mutate(
    variable       = factor(variable, levels = names(base_colors)),
    variable       = droplevels(variable),
    direction      = ifelse(effect_size > 0, "Positive", "Negative"),
    atlas          = factor(atlas_rename[as.character(atlas)],
                            levels = c("1970s", "1990s", "2010s")),
    trait_category = factor(trait_category,
                            levels = c("Thermal index", "Foraging guild", "Migration")),
    trait_clean    = factor(trait_clean)
  )

# ── 3. Cluster traits within each category ────────────────────────────────────

cluster_within_category <- function(df) {
  cats           <- levels(df$trait_category)
  ordered_traits <- c()
  for (cat in cats) {
    sub <- df %>%
      filter(trait_category == cat) %>%
      group_by(trait_clean, variable) %>%
      summarise(mean_effect = mean(effect_size, na.rm = TRUE), .groups = "drop") %>%
      mutate(mean_effect = replace_na(mean_effect, 0)) %>%
      pivot_wider(names_from = variable, values_from = mean_effect, values_fill = 0) %>%
      tibble::column_to_rownames("trait_clean")
    if (nrow(sub) > 1) {
      ordered_traits <- c(ordered_traits, rownames(sub)[hclust(dist(sub))$order])
    } else {
      ordered_traits <- c(ordered_traits, rownames(sub))
    }
  }
  ordered_traits
}

trait_order <- cluster_within_category(gamma_plot)

env_mat <- gamma_plot %>%
  group_by(trait_clean, variable) %>%
  summarise(mean_effect = mean(effect_size, na.rm = TRUE), .groups = "drop") %>%
  mutate(mean_effect = replace_na(mean_effect, 0)) %>%
  pivot_wider(names_from = variable, values_from = mean_effect, values_fill = 0) %>%
  tibble::column_to_rownames("trait_clean")

env_order <- colnames(env_mat)[hclust(dist(t(env_mat)))$order]

gamma_plot <- gamma_plot %>%
  mutate(
    trait_clean = factor(trait_clean, levels = trait_order),
    variable    = factor(variable,    levels = env_order)
  )

# Colours actually present across the whole dataset
present_vars   <- levels(droplevels(gamma_plot$variable))
base_colors_sub <- base_colors[names(base_colors) %in% present_vars]

# ── 4. Shared layout objects ──────────────────────────────────────────────────

category_dividers <- gamma_plot %>%
  distinct(trait_clean, trait_category) %>%
  mutate(y = as.numeric(trait_clean)) %>%
  group_by(trait_category) %>%
  summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5, ymid = mean(y), .groups = "drop")

n_envs         <- length(present_vars)
gamma_filled   <- gamma_plot %>% filter(atlas == "1970s")
gamma_rings    <- gamma_plot %>% filter(atlas != "1970s")

# ── 5. Panel function ─────────────────────────────────────────────────────────

make_panel <- function(direction_filter, panel_title) {
  
  gf <- gamma_filled %>% filter(direction == direction_filter) %>% droplevels()
  gr <- gamma_rings  %>% filter(direction == direction_filter) %>% droplevels()
  
  # Colors present in THIS panel
  pv <- unique(c(as.character(gf$variable), as.character(gr$variable)))
  pv <- pv[!is.na(pv) & pv %in% names(base_colors)]
  bc <- base_colors[pv]
  
  ggplot() +
    
    # Background cells (always use full gamma_plot for consistent grid)
    geom_rect(
      data = gamma_plot %>%
        distinct(trait_clean, variable) %>%
        mutate(x = as.numeric(variable), y = as.numeric(trait_clean)),
      aes(xmin = x - 0.5, xmax = x + 0.5, ymin = y - 0.5, ymax = y + 0.5),
      fill = "grey97", color = "grey88", linewidth = 0.2
    ) +
    
    # Category divider lines
    geom_hline(
      data = category_dividers %>%
        filter(trait_category != last(levels(gamma_plot$trait_category))),
      aes(yintercept = ymax),
      color = "grey60", linewidth = 0.4, linetype = "dashed"
    ) +
    
    # Atlas 3 ring (largest, behind)
    geom_point(
      data = gr %>% filter(atlas == "2010s"),
      aes(x = as.numeric(variable), y = as.numeric(trait_clean), color = variable),
      size = atlas_sizes["2010s"], shape = 21, fill = NA, stroke = 1.8
    ) +
    
    # Atlas 2 ring (medium)
    geom_point(
      data = gr %>% filter(atlas == "1990s"),
      aes(x = as.numeric(variable), y = as.numeric(trait_clean), color = variable),
      size = atlas_sizes["1990s"], shape = 21, fill = NA, stroke = 1.5
    ) +
    
    # Atlas 1 filled circle (smallest, front)
    geom_point(
      data = gf,
      aes(x = as.numeric(variable), y = as.numeric(trait_clean),
          color = variable, fill = variable),
      size = atlas_sizes["1970s"], shape = 21, stroke = 0.5
    ) +
    
    # Category colour strip
    geom_rect(
      data = category_dividers,
      aes(
        xmin = n_envs + 0.6, xmax = n_envs + 0.85,
        ymin = ymin, ymax = ymax,
        fill = trait_category
      ),
      color = NA
    ) +
    
    # Category labels
    geom_text(
      data = category_dividers,
      aes(x = n_envs + 1.05, y = ymid, label = trait_category,
          color = trait_category),
      hjust = 0, size = 2.8, fontface = "bold"
    ) +
    
    scale_color_manual(
      values = c(bc, category_colors),
      breaks = names(bc),
      name   = "Environment"
    ) +
    scale_fill_manual(
      values = c(bc, category_colors),
      breaks = names(bc),
      name   = "Environment"
    ) +
    
    scale_x_continuous(
      breaks = seq_along(present_vars),
      labels = present_vars,
      expand = expansion(add = c(0.5, 3.5))
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(gamma_plot$trait_clean)),
      labels = levels(gamma_plot$trait_clean),
      expand = expansion(add = 0.5)
    ) +
    
    coord_fixed() +
    ggtitle(panel_title) +
    
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x  = element_text(angle = 35, hjust = 1, size = 8),
      axis.text.y  = element_text(size = 8),
      axis.title   = element_blank(),
      panel.grid   = element_blank(),
      plot.title   = element_text(face = "bold", size = 11, hjust = 0.5),
      legend.position = "bottom",
      plot.margin  = margin(10, 80, 10, 10)
    ) +
    guides(
      color = guide_legend(
        nrow = 2,
        override.aes = list(size = 3, shape = 16, fill = unname(bc))  # use bc not base_colors
      ),
      fill = "none"
    )
}

# ── 6. Combine ────────────────────────────────────────────────────────────────

p_pos <- make_panel("Positive", "Positive effects")
p_neg <- make_panel("Negative", "Negative effects")

p_combined <- p_pos + p_neg +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(p_combined)

gamma_plot

ggsave("gamma_circles.pdf", p_combined, width = 18, height = 10)
ggsave("misc-figures/gamma_circles.png", p_combined, width = 18, height = 10, dpi = 300)

