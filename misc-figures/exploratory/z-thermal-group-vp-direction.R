remove(list = ls())
.libPaths(c("~/Rlibs", .libPaths()))

library(tidyverse)
library(readxl)
library(scales)

source(file.path("support_scripts", "figure_data_helpers.R"))

pattern2match <- "2026-03-13"
base_dir <- "HmscOutputs"
out_dir <- file.path("misc-figures", "outputs", "exploratory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

thermal_levels <- c("very cold", "cold", "medium", "warm", "very warm")
atlas_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")

var_colors <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban" = "snow3",
  "Cropland" = "goldenrod1",
  "Pasture" = "darkorange",
  "Forest" = "springgreen4",
  "Grass/Shrubland" = "springgreen2"
)
var_order <- names(var_colors)

clean_variable_names <- function(x) {
  case_match(
    x,
    "tmean_breeding" ~ "Temperature",
    "prec_breeding" ~ "Precipitation",
    "hh" ~ "Land-use heterogeneity",
    "perc_urban" ~ "Urban",
    "perc_cropland" ~ "Cropland",
    "perc_pasture" ~ "Pasture",
    "perc_forest" ~ "Forest",
    "perc_grass_shrub" ~ "Grass/Shrubland",
    "Random: site" ~ "Site-level random effects",
    .default = x
  ) |>
    str_replace(" \\(% coverage\\)$", "")
}

thermal_groups <- function(sti) {
  cut(
    sti,
    breaks = quantile(sti, seq(0, 1, 0.2), na.rm = TRUE),
    include.lowest = TRUE,
    labels = thermal_levels
  )
}

message("Loading variance partitioning estimates and model traits...")
model_folders <- figure_model_folders(pattern = pattern2match, base_dir = base_dir)
vp_scaled <- load_vp_estimates(model_folders, base_dir = base_dir, scaled = TRUE)

vp_long <- imap_dfr(vp_scaled, function(df, atlas) {
  df |>
    rownames_to_column("variable") |>
    pivot_longer(-variable, names_to = "species", values_to = "VP") |>
    mutate(atlas = as.character(atlas))
}) |>
  filter(variable != "TjurR2") |>
  mutate(variable_clean = clean_variable_names(variable)) |>
  filter(variable_clean %in% var_order)

load(file.path(
  base_dir,
  model_folders[[1]],
  "Models/Fitted",
  "HPC_samples_0250_thin_100_chains_4.Rdata"
))

traits <- fitted_model$posteriors$TrData |>
  as.data.frame() |>
  rownames_to_column("species") |>
  select(species, Migration_a3_DOF, foraging_guild_consensus, species_thermal_index) |>
  filter(!is.na(species_thermal_index)) |>
  mutate(
    thermal_group = thermal_groups(species_thermal_index),
    thermal_group = factor(thermal_group, levels = thermal_levels)
  )

thermal_group_summary <- traits |>
  group_by(thermal_group) |>
  summarise(
    n_species = n_distinct(species),
    min_sti = min(species_thermal_index, na.rm = TRUE),
    median_sti = median(species_thermal_index, na.rm = TRUE),
    max_sti = max(species_thermal_index, na.rm = TRUE),
    .groups = "drop"
  )

message("Loading supported species-level Beta effects...")
beta <- read_parameter_effects(pattern2match, effect = "Beta", base_dir = base_dir) |>
  filter(variable != "(Intercept)") |>
  mutate(
    atlas = as.character(atlas),
    variable_clean = clean_variable_names(variable),
    effect_sign = sign(effect_size)
  ) |>
  filter(variable_clean %in% var_order) |>
  select(species, atlas, variable_clean, effect_size, effect_sign)

thermal_vp_by_species <- vp_long |>
  inner_join(traits, by = "species") |>
  left_join(beta, by = c("species", "atlas", "variable_clean")) |>
  mutate(
    atlas_label = recode(atlas, !!!atlas_lookup),
    supported = !is.na(effect_sign),
    signed_vp = ifelse(supported, VP * effect_sign, NA_real_)
  )

thermal_vp_summary <- thermal_vp_by_species |>
  group_by(thermal_group, variable_clean, atlas, atlas_label) |>
  summarise(
    n_species = n_distinct(species),
    mean_vp = mean(VP, na.rm = TRUE),
    median_vp = median(VP, na.rm = TRUE),
    signed_vp_mean = mean(signed_vp, na.rm = TRUE),
    supported_species = n_distinct(species[supported]),
    supported_vp = sum(ifelse(supported, VP, 0), na.rm = TRUE),
    total_vp = sum(VP, na.rm = TRUE),
    direction_score = ifelse(supported_vp > 0, sum(VP * effect_sign, na.rm = TRUE) / supported_vp, NA_real_),
    positive_species = n_distinct(species[effect_sign > 0]),
    negative_species = n_distinct(species[effect_sign < 0]),
    .groups = "drop"
  ) |>
  mutate(
    supported_vp_share = ifelse(total_vp > 0, supported_vp / total_vp, NA_real_),
    effect_direction = case_when(
      is.na(direction_score) ~ "Not supported",
      direction_score >= 0.1 ~ "Higher at max",
      direction_score <= -0.1 ~ "Lower at max",
      TRUE ~ "Mixed"
    )
  )

thermal_vp_overview <- thermal_vp_summary |>
  group_by(thermal_group, variable_clean) |>
  summarise(
    n_species_min = min(n_species, na.rm = TRUE),
    vp_central = median(median_vp, na.rm = TRUE),
    mean_vp_central = mean(mean_vp, na.rm = TRUE),
    sd_across_atlas = sd(median_vp, na.rm = TRUE),
    supported_vp_share = mean(supported_vp_share, na.rm = TRUE),
    direction_score = ifelse(
      sum(supported_vp, na.rm = TRUE) > 0,
      sum(supported_vp * direction_score, na.rm = TRUE) / sum(supported_vp, na.rm = TRUE),
      NA_real_
    ),
    positive_species_mean = mean(positive_species, na.rm = TRUE),
    negative_species_mean = mean(negative_species, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    signed_vp_central = vp_central * direction_score,
    effect_direction = case_when(
      is.na(direction_score) ~ "Not supported",
      direction_score >= 0.1 ~ "Higher at max",
      direction_score <= -0.1 ~ "Lower at max",
      TRUE ~ "Mixed"
    ),
    thermal_group = factor(thermal_group, levels = thermal_levels),
    variable_clean = factor(variable_clean, levels = var_order),
    effect_direction = factor(
      effect_direction,
      levels = c("Higher at max", "Lower at max", "Mixed", "Not supported")
    )
  ) |>
  left_join(thermal_group_summary, by = "thermal_group")

write_csv(
  thermal_group_summary,
  file.path(out_dir, paste0(pattern2match, "-thermal-group-species-summary.csv"))
)
write_csv(
  thermal_vp_summary,
  file.path(out_dir, paste0(pattern2match, "-thermal-group-vp-direction-by-atlas.csv"))
)
write_csv(
  thermal_vp_overview,
  file.path(out_dir, paste0(pattern2match, "-thermal-group-vp-direction-overview.csv"))
)

message("Plotting thermal group VP overview...")
bubble_df <- thermal_vp_overview |>
  mutate(
    alpha_val = 1 - rescale(sd_across_atlas, to = c(0.15, 0.85), from = range(sd_across_atlas, na.rm = TRUE)),
    thermal_label = paste0(thermal_group, "\n(n=", n_species, ")"),
    thermal_label = factor(thermal_label, levels = rev(paste0(thermal_levels, "\n(n=", thermal_group_summary$n_species, ")")))
  )

p_bubble <- ggplot(bubble_df, aes(variable_clean, thermal_label)) +
  geom_point(
    aes(size = vp_central, colour = variable_clean, alpha = alpha_val, shape = effect_direction),
    stroke = 1.1
  ) +
  scale_colour_manual(values = var_colors, guide = "none") +
  scale_shape_manual(
    values = c("Higher at max" = 16, "Lower at max" = 4, "Mixed" = 1, "Not supported" = 3),
    drop = FALSE,
    name = "Direction at\nhigh predictor value"
  ) +
  scale_size_continuous(
    range = c(2, 13),
    labels = label_number(accuracy = 0.001),
    name = "Median VP (R2)"
  ) +
  scale_alpha_identity(
    name = "Certainty\n(1 - scaled\ntemporal SD)",
    guide = guide_legend()
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title = "Variance partitioning across thermal-affinity groups",
    subtitle = "Size: median species-level VP across atlas periods; symbol: VP-weighted supported Beta direction",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey40", size = 10)
  )

ggsave(
  file.path(out_dir, paste0(pattern2match, "-thermal-group-vp-direction-bubble.png")),
  p_bubble,
  width = 9,
  height = 4.8,
  dpi = 300
)

bar_df <- thermal_vp_overview |>
  mutate(
    variable_clean = fct_reorder(variable_clean, signed_vp_central, .fun = median, na.rm = TRUE),
    sign_class = ifelse(signed_vp_central >= 0, "Higher at max", "Lower at max")
  )

p_bar <- ggplot(bar_df, aes(signed_vp_central, variable_clean, fill = variable_clean)) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.35) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  facet_wrap(~ thermal_group, ncol = 5) +
  scale_fill_manual(values = var_colors, guide = "none") +
  scale_x_continuous(labels = label_number(accuracy = 0.001)) +
  labs(
    title = "Signed environmental relationships by thermal-affinity group",
    subtitle = "Positive means higher predicted occurrence at high environmental-variable values; magnitude is median VP weighted by supported Beta direction",
    x = "Signed median VP (R2)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey40", size = 10)
  )

ggsave(
  file.path(out_dir, paste0(pattern2match, "-thermal-group-signed-vp-bars.png")),
  p_bar,
  width = 12,
  height = 5.6,
  dpi = 300
)

message("Finished thermal group VP direction overview.")

message("Summarising guild and migration composition inside thermal groups...")

summarise_trait_composition <- function(df, trait_col, trait_label) {
  trait_col <- enquo(trait_col)

  df |>
    filter(!is.na(!!trait_col), !is.na(thermal_group)) |>
    count(thermal_group, trait_value = !!trait_col, name = "n_species") |>
    group_by(thermal_group) |>
    mutate(
      thermal_group_n = sum(n_species),
      prop_species = n_species / thermal_group_n
    ) |>
    ungroup() |>
    group_by(trait_value) |>
    mutate(
      trait_total_n = sum(n_species),
      weighted_thermal_position = weighted.mean(as.integer(thermal_group), n_species)
    ) |>
    ungroup() |>
    mutate(
      trait_family = trait_label,
      thermal_group = factor(thermal_group, levels = thermal_levels)
    )
}

guild_composition <- summarise_trait_composition(
  traits,
  foraging_guild_consensus,
  "Foraging guild"
)

migration_composition <- summarise_trait_composition(
  traits,
  Migration_a3_DOF,
  "Migratory strategy"
)

write_csv(
  guild_composition,
  file.path(out_dir, paste0(pattern2match, "-thermal-group-foraging-guild-composition.csv"))
)
write_csv(
  migration_composition,
  file.path(out_dir, paste0(pattern2match, "-thermal-group-migratory-strategy-composition.csv"))
)

plot_composition_heatmap <- function(composition_df, title, min_total_n = 1) {
  trait_order <- composition_df |>
    filter(trait_total_n >= min_total_n) |>
    distinct(trait_value, trait_total_n, weighted_thermal_position) |>
    arrange(weighted_thermal_position, desc(trait_total_n), trait_value) |>
    pull(trait_value)

  plot_df <- composition_df |>
    filter(trait_value %in% trait_order) |>
    mutate(
      trait_value = factor(trait_value, levels = rev(trait_order)),
      label = ifelse(n_species > 0, as.character(n_species), "")
    )

  ggplot(plot_df, aes(thermal_group, trait_value)) +
    geom_tile(aes(fill = prop_species), colour = "white", linewidth = 0.35) +
    geom_text(aes(label = label), size = 3.2, colour = "grey15") +
    scale_fill_gradient(
      low = "white",
      high = "#2166ac",
      labels = percent_format(accuracy = 1),
      name = "Share within\nthermal group"
    ) +
    labs(
      title = title,
      subtitle = "Numbers are species counts; fill is the proportion of each thermal-affinity group",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey40", size = 10)
    )
}

plot_composition_stack <- function(composition_df, title, palette = NULL) {
  trait_order <- composition_df |>
    distinct(trait_value, trait_total_n, weighted_thermal_position) |>
    arrange(weighted_thermal_position, desc(trait_total_n), trait_value) |>
    pull(trait_value)

  plot_df <- composition_df |>
    mutate(trait_value = factor(trait_value, levels = trait_order))

  if (is.null(palette)) {
    palette <- setNames(hue_pal()(length(trait_order)), trait_order)
  }

  ggplot(plot_df, aes(thermal_group, prop_species, fill = trait_value)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = palette, name = NULL) +
    labs(
      title = title,
      subtitle = "Stack height is the species composition within each thermal-affinity group",
      x = NULL,
      y = "Species share"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey40", size = 10)
    )
}

guild_palette <- setNames(hue_pal(l = 60, c = 80)(n_distinct(guild_composition$trait_value)),
                          sort(unique(guild_composition$trait_value)))
migration_palette <- setNames(
  c("#4daf4a", "#377eb8", "#984ea3", "#ff7f00", "#e41a1c")[seq_len(n_distinct(migration_composition$trait_value))],
  sort(unique(migration_composition$trait_value))
)

p_guild_heatmap <- plot_composition_heatmap(
  guild_composition,
  "Foraging guild composition across thermal-affinity groups"
)
p_migration_heatmap <- plot_composition_heatmap(
  migration_composition,
  "Migratory strategy composition across thermal-affinity groups"
)
p_guild_stack <- plot_composition_stack(
  guild_composition,
  "Foraging guild shares across thermal-affinity groups",
  guild_palette
)
p_migration_stack <- plot_composition_stack(
  migration_composition,
  "Migratory strategy shares across thermal-affinity groups",
  migration_palette
)

ggsave(
  file.path(out_dir, paste0(pattern2match, "-thermal-group-foraging-guild-composition-heatmap.png")),
  p_guild_heatmap,
  width = 8.6,
  height = 8.8,
  dpi = 300
)
ggsave(
  file.path(out_dir, paste0(pattern2match, "-thermal-group-migratory-strategy-composition-heatmap.png")),
  p_migration_heatmap,
  width = 8.2,
  height = 4.4,
  dpi = 300
)
ggsave(
  file.path(out_dir, paste0(pattern2match, "-thermal-group-foraging-guild-composition-stacked.png")),
  p_guild_stack,
  width = 10,
  height = 6.2,
  dpi = 300
)
ggsave(
  file.path(out_dir, paste0(pattern2match, "-thermal-group-migratory-strategy-composition-stacked.png")),
  p_migration_stack,
  width = 8.8,
  height = 4.8,
  dpi = 300
)

message("Finished thermal group trait-composition overview.")
