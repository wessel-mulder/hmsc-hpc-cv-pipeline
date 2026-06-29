# Beta-parameter trends by trait group
#
# This exploratory script mirrors the temporal variance-partitioning trend
# figure, but uses species-level HMSC beta parameters as the response. It
# summarizes every thermal-affinity group, migratory strategy, and foraging
# guild for each environmental predictor. It writes matched versions of the
# figures for raw beta parameters and beta effect sizes scaled by the within-
# atlas standard deviation of each environmental variable.
#
# Output policy for this project: all plots written by this script are PNGs.

remove(list = ls())
.libPaths(c("~/Rlibs", .libPaths()))

required_packages <- c("tidyverse", "readxl", "patchwork", "scales", "Hmsc")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(tidyverse)
library(readxl)
library(patchwork)
library(scales)
library(Hmsc)

source(file.path("support_scripts", "figure_data_helpers.R"))

# ---- Configuration ----
pattern <- "2026-03-13"
base_dir <- "HmscOutputs"
out_dir <- file.path("misc-figures", "outputs", "exploratory", "beta-parameter-trait-trends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Keep all trait groups, including small guilds, because the goal here is a
# complete scan rather than a filtered manuscript panel.
min_species_per_trait_level <- 1

period_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")
period_levels <- unname(period_lookup)
thermal_levels <- c("very cold", "cold", "medium", "warm", "very warm")
migration_levels <- c(
  "long-distance",
  "short-and long-distance",
  "short-distance",
  "sedentary and short-distance",
  "sedentary"
)

variable_order <- c(
  "tmean_breeding",
  "prec_breeding",
  "hh",
  "perc_urban",
  "perc_cropland",
  "perc_pasture",
  "perc_forest",
  "perc_grass_shrub"
)

variable_label_lookup <- c(
  "tmean_breeding" = "Temperature",
  "prec_breeding" = "Precipitation",
  "hh" = "Habitat\nheterogeneity",
  "perc_urban" = "Urban",
  "perc_cropland" = "Cropland",
  "perc_pasture" = "Pasture",
  "perc_forest" = "Forest",
  "perc_grass_shrub" = "Grass/\nshrubland"
)

effect_scale_order <- c("Unscaled beta", "SD-scaled effect size")
effect_scale_slugs <- c(
  "Unscaled beta" = "unscaled-beta",
  "SD-scaled effect size" = "sd-scaled-effect-size"
)
effect_scale_legend_titles <- c(
  "Unscaled beta" = "Temporal beta slope",
  "SD-scaled effect size" = "Temporal effect-size slope"
)
effect_scale_subtitles <- c(
  "Unscaled beta" = "Slope is fitted across atlas period using species-level certainty-weighted beta values.",
  "SD-scaled effect size" = "Slope is fitted across atlas period using species-level certainty-weighted beta values multiplied by each atlas-specific predictor SD."
)

clean_variable_label <- function(variable) {
  recode(
    variable,
    !!!variable_label_lookup,
    .default = str_to_sentence(str_replace_all(variable, "_", " "))
  )
}

thermal_groups <- function(sti) {
  cut(
    sti,
    breaks = quantile(sti, seq(0, 1, 0.2), na.rm = TRUE),
    include.lowest = TRUE,
    labels = thermal_levels
  )
}

support_state <- function(prob_positive) {
  case_when(
    prob_positive >= 0.95 ~ "Supported positive",
    prob_positive <= 0.05 ~ "Supported negative",
    TRUE ~ "Unsupported"
  )
}

read_beta_parameter_workbook <- function(model_folder, base_dir = "HmscOutputs") {
  beta_path <- file.path(
    base_dir,
    model_folder,
    "Results",
    paste0(model_folder, "parameter_estimates_Beta_.xlsx")
  )

  posterior_mean <- read_excel(beta_path, sheet = "Posterior mean") |>
    pivot_longer(
      cols = -Species,
      names_to = "variable",
      values_to = "beta_mean"
    )

  posterior_positive <- read_excel(beta_path, sheet = "Pr(x>0)") |>
    pivot_longer(
      cols = -Species,
      names_to = "variable",
      values_to = "prob_positive"
    )

  posterior_mean |>
    left_join(posterior_positive, by = c("Species", "variable")) |>
    transmute(
      species = Species,
      variable,
      beta_mean,
      prob_positive,
      atlas = as.character(sub(".*Atlas([0-9]+).*", "\\1", model_folder)),
      model_folder = model_folder
    )
}

predictor_sds_by_atlas <- function(models, variables) {
  imap_dfr(models, function(model, atlas) {
    tibble(
      atlas = as.character(atlas),
      variable = variables,
      predictor_sd = map_dbl(variables, ~ sd(model$XData[[.x]], na.rm = TRUE))
    )
  })
}

fit_temporal_slope <- function(df, response_col) {
  response_name <- rlang::as_string(rlang::ensym(response_col))

  model_df <- df |>
    filter(is.finite(.data[[response_name]]), is.finite(atlas_numeric))

  if (nrow(model_df) < 3 || n_distinct(model_df$atlas_numeric) < 2) {
    return(tibble(
      slope = NA_real_,
      p_value = NA_real_,
      r_squared = NA_real_
    ))
  }

  model <- lm(reformulate("atlas_numeric", response = response_name), data = model_df)
  model_summary <- summary(model)
  coef_table <- coef(model_summary)

  tibble(
    slope = coef_table["atlas_numeric", "Estimate"],
    p_value = coef_table["atlas_numeric", "Pr(>|t|)"],
    r_squared = model_summary$r.squared
  )
}

make_slope_panel <- function(plot_data,
                             title,
                             y_levels,
                             slope_limit,
                             legend_title,
                             show_x_labels = FALSE) {
  plot_data |>
    mutate(
      trait_value = factor(as.character(trait_value), levels = y_levels),
      variable_label = factor(variable_label, levels = clean_variable_label(variable_order))
    ) |>
    ggplot(aes(x = variable_label, y = trait_value, fill = signed_beta_slope)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sig_label), colour = "grey10", fontface = "bold", size = 3) +
    scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      limits = c(-slope_limit, slope_limit),
      oob = squish,
      breaks = pretty_breaks(n = 5),
      labels = label_number(accuracy = 0.01),
      name = legend_title,
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(3.2, "in"),
        barheight = grid::unit(0.18, "in")
      )
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = if (show_x_labels) {
        element_text(size = 8.4, angle = 35, hjust = 1, vjust = 1)
      } else {
        element_blank()
      },
      axis.text.y = element_text(size = 8.6),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.title = element_text(face = "bold")
    )
}

make_period_heatmap <- function(plot_data,
                                trait_family_name,
                                y_levels,
                                output_slug,
                                plot_height,
                                effect_scale_name,
                                also_save_as = character()) {
  scale_data <- plot_data |>
    filter(effect_scale == effect_scale_name, trait_family == trait_family_name)

  beta_limit <- max(abs(scale_data$median_weighted_beta_value), na.rm = TRUE)
  beta_limit <- ifelse(is.finite(beta_limit) && beta_limit > 0, beta_limit * 1.05, 1)

  p <- scale_data |>
    mutate(
      trait_value = factor(as.character(trait_value), levels = y_levels),
      period = factor(period, levels = period_levels),
      variable_label = factor(variable_label, levels = clean_variable_label(variable_order))
    ) |>
    ggplot(aes(x = period, y = trait_value, fill = median_weighted_beta_value)) +
    geom_tile(colour = "white", linewidth = 0.18) +
    facet_wrap(~ variable_label, ncol = 4) +
    scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      limits = c(-beta_limit, beta_limit),
      oob = squish,
      labels = label_number(accuracy = 0.01),
      name = if_else(effect_scale_name == "SD-scaled effect size", "Median effect size", "Median beta")
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(
      title = paste0(trait_family_name, " beta-parameter trends: ", effect_scale_name),
      subtitle = paste0(
        "Tiles show median certainty-weighted ",
        if_else(effect_scale_name == "SD-scaled effect size", "SD-scaled beta effect size", "raw beta"),
        " by atlas period."
      ),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey35"),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  output_slugs <- c(output_slug, also_save_as)
  walk(output_slugs, function(slug) {
    ggsave(
      filename = file.path(out_dir, paste0(pattern, "-", slug, ".png")),
      plot = p,
      width = 11,
      height = plot_height,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  })

  invisible(p)
}

cluster_trait_levels_by_slopes <- function(slope_df, trait_family_name, fallback_levels) {
  variable_levels <- clean_variable_label(variable_order)

  slope_matrix <- slope_df |>
    filter(trait_family == trait_family_name, trait_value %in% fallback_levels) |>
    mutate(
      trait_value = factor(as.character(trait_value), levels = fallback_levels),
      variable_label = factor(as.character(variable_label), levels = variable_levels)
    ) |>
    select(trait_value, variable_label, signed_beta_slope) |>
    complete(
      trait_value = fallback_levels,
      variable_label = variable_levels,
      fill = list(signed_beta_slope = 0)
    ) |>
    pivot_wider(names_from = variable_label, values_from = signed_beta_slope) |>
    arrange(factor(trait_value, levels = fallback_levels)) |>
    column_to_rownames("trait_value") |>
    as.matrix()

  slope_matrix[!is.finite(slope_matrix)] <- 0

  if (nrow(slope_matrix) < 3) {
    return(fallback_levels)
  }

  # Scale each environmental-variable column before clustering so guilds are
  # grouped by the shape of their trend profiles, rather than being dominated by
  # whichever predictor happens to have the widest slope range.
  scaled_matrix <- scale(slope_matrix)
  scaled_matrix[!is.finite(scaled_matrix)] <- 0

  rownames(scaled_matrix)[hclust(dist(scaled_matrix), method = "ward.D2")$order]
}

# ---- Load beta coefficients and species traits ----
message("Loading model folders and beta coefficients...")
model_folders <- figure_model_folders(pattern = pattern, base_dir = base_dir)
model_atlases <- atlas_numbers(model_folders)

mods <- load_hmsc_posteriors(model_folders, base_dir = base_dir)
names(mods) <- as.character(model_atlases)

beta_variables <- intersect(variable_order, colnames(mods[[1]]$XData))
variable_labels <- clean_variable_label(beta_variables)

beta_raw <- map_dfr(model_folders, read_beta_parameter_workbook, base_dir = base_dir) |>
  filter(variable %in% beta_variables)

beta_sds <- predictor_sds_by_atlas(mods, beta_variables)

beta_long <- beta_raw |>
  left_join(beta_sds, by = c("atlas", "variable")) |>
  mutate(
    atlas_numeric = as.numeric(atlas),
    period = factor(period_lookup[atlas], levels = period_levels),
    certainty_weight = 2 * abs(prob_positive - 0.5),
    unscaled_beta = beta_mean,
    weighted_unscaled_beta = unscaled_beta * certainty_weight,
    abs_weighted_unscaled_beta = abs(weighted_unscaled_beta),
    sd_scaled_effect_size = beta_mean * predictor_sd,
    weighted_sd_scaled_effect_size = sd_scaled_effect_size * certainty_weight,
    abs_weighted_sd_scaled_effect_size = abs(weighted_sd_scaled_effect_size),
    # Keep the older column names as aliases so previous exploratory notebooks
    # do not break when this script is rerun.
    standardized_beta = sd_scaled_effect_size,
    weighted_standardized_beta = weighted_sd_scaled_effect_size,
    abs_weighted_standardized_beta = abs_weighted_sd_scaled_effect_size,
    beta_support_state = support_state(prob_positive),
    variable_label = clean_variable_label(variable),
    variable_label = factor(variable_label, levels = variable_labels)
  )

message("Loading trait metadata...")
load(file.path("Data", "preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData"))

species_traits <- Tr |>
  as.data.frame() |>
  rownames_to_column("species") |>
  transmute(
    species,
    thermal_group = thermal_groups(species_thermal_index),
    migratory_strategy = Migration_a3_DOF,
    foraging_guild = foraging_guild_consensus,
    species_thermal_index
  )

guild_levels <- species_traits |>
  filter(!is.na(foraging_guild), foraging_guild != "") |>
  count(foraging_guild, sort = TRUE, name = "n_species") |>
  filter(n_species >= min_species_per_trait_level) |>
  arrange(desc(n_species), foraging_guild) |>
  pull(foraging_guild)

trait_membership <- bind_rows(
  species_traits |>
    transmute(
      species,
      trait_family = "All species",
      trait_value = "All species"
    ),
  species_traits |>
    transmute(
      species,
      trait_family = "Thermal affinity",
      trait_value = as.character(thermal_group)
    ),
  species_traits |>
    transmute(
      species,
      trait_family = "Migratory strategy",
      trait_value = as.character(migratory_strategy)
    ),
  species_traits |>
    transmute(
      species,
      trait_family = "Foraging guild",
      trait_value = as.character(foraging_guild)
    )
) |>
  filter(!is.na(trait_value), trait_value != "")

beta_trait_species <- beta_long |>
  inner_join(trait_membership, by = "species", relationship = "many-to-many")

beta_trait_species_effects <- bind_rows(
  beta_trait_species |>
    mutate(
      effect_scale = "Unscaled beta",
      weighted_beta_value = weighted_unscaled_beta,
      abs_weighted_beta_value = abs_weighted_unscaled_beta
    ),
  beta_trait_species |>
    mutate(
      effect_scale = "SD-scaled effect size",
      weighted_beta_value = weighted_sd_scaled_effect_size,
      abs_weighted_beta_value = abs_weighted_sd_scaled_effect_size
    )
) |>
  mutate(effect_scale = factor(effect_scale, levels = effect_scale_order))

# ---- Group-period summaries ----
message("Summarising beta trends by trait group and period...")
beta_group_period_summary <- beta_trait_species_effects |>
  group_by(effect_scale, trait_family, trait_value, variable, variable_label, atlas, atlas_numeric, period) |>
  summarise(
    n_species = n_distinct(species),
    mean_weighted_beta_value = mean(weighted_beta_value, na.rm = TRUE),
    median_weighted_beta_value = median(weighted_beta_value, na.rm = TRUE),
    mean_abs_weighted_beta_value = mean(abs_weighted_beta_value, na.rm = TRUE),
    median_abs_weighted_beta_value = median(abs_weighted_beta_value, na.rm = TRUE),
    supported_positive_share = mean(beta_support_state == "Supported positive", na.rm = TRUE),
    supported_negative_share = mean(beta_support_state == "Supported negative", na.rm = TRUE),
    unsupported_share = mean(beta_support_state == "Unsupported", na.rm = TRUE),
    .groups = "drop"
  )

# ---- Group temporal slopes ----
beta_group_temporal_slopes <- beta_trait_species_effects |>
  group_by(effect_scale, trait_family, trait_value, variable, variable_label) |>
  group_modify(~ {
    signed_fit <- fit_temporal_slope(.x, weighted_beta_value)
    strength_fit <- fit_temporal_slope(.x, abs_weighted_beta_value)

    tibble(
      n_observations = nrow(.x),
      n_species = n_distinct(.x$species),
      signed_beta_slope = signed_fit$slope,
      signed_beta_p_value = signed_fit$p_value,
      signed_beta_r_squared = signed_fit$r_squared,
      beta_strength_slope = strength_fit$slope,
      beta_strength_p_value = strength_fit$p_value,
      beta_strength_r_squared = strength_fit$r_squared,
      mean_1970s = mean(.x$weighted_beta_value[.x$atlas == "1"], na.rm = TRUE),
      mean_2010s = mean(.x$weighted_beta_value[.x$atlas == "3"], na.rm = TRUE),
      mean_2010s_minus_1970s = mean_2010s - mean_1970s
    )
  }) |>
  ungroup() |>
  mutate(
    significance = case_when(
      signed_beta_p_value < 0.001 ~ "***",
      signed_beta_p_value < 0.01 ~ "**",
      signed_beta_p_value < 0.05 ~ "*",
      signed_beta_p_value < 0.1 ~ ".",
      TRUE ~ "ns"
    ),
    sig_label = if_else(significance == "ns" | is.na(significance), "", significance),
    effect_scale = factor(effect_scale, levels = effect_scale_order)
  )

guild_levels_by_scale <- set_names(
  map(effect_scale_order, function(effect_scale_name) {
    beta_group_temporal_slopes |>
      filter(effect_scale == effect_scale_name) |>
      cluster_trait_levels_by_slopes(
        trait_family_name = "Foraging guild",
        fallback_levels = guild_levels
      )
  }),
  effect_scale_order
)

guild_cluster_order <- map_dfr(effect_scale_order, function(effect_scale_name) {
  ordered_guilds <- guild_levels_by_scale[[effect_scale_name]]

  beta_group_temporal_slopes |>
    filter(
      effect_scale == effect_scale_name,
      trait_family == "Foraging guild",
      trait_value %in% ordered_guilds
    ) |>
    distinct(effect_scale, trait_value, n_species) |>
    mutate(cluster_order = match(trait_value, ordered_guilds)) |>
    arrange(cluster_order)
})

write_csv(
  beta_long,
  file.path(out_dir, paste0(pattern, "-beta-standardized-effects.csv"))
)
write_csv(
  beta_long,
  file.path(out_dir, paste0(pattern, "-beta-species-effects.csv"))
)
write_csv(
  beta_trait_species_effects,
  file.path(out_dir, paste0(pattern, "-beta-trait-species-effects-long.csv"))
)
write_csv(
  beta_group_period_summary,
  file.path(out_dir, paste0(pattern, "-beta-trait-group-period-summary.csv"))
)
write_csv(
  beta_group_temporal_slopes,
  file.path(out_dir, paste0(pattern, "-beta-trait-group-temporal-slopes.csv"))
)
write_csv(
  guild_cluster_order,
  file.path(out_dir, paste0(pattern, "-beta-foraging-guild-trend-cluster-order.csv"))
)

# ---- Main slope figure, matching the VP trend style ----
message("Plotting temporal slope heatmap...")

make_slope_figure <- function(effect_scale_name, output_slug, also_save_as = character()) {
  slope_data <- beta_group_temporal_slopes |>
    filter(effect_scale == effect_scale_name)

  slope_limit <- max(abs(slope_data$signed_beta_slope), na.rm = TRUE)
  slope_limit <- ifelse(is.finite(slope_limit) && slope_limit > 0, slope_limit * 1.05, 1)

  all_species_panel <- slope_data |>
    filter(trait_family == "All species") |>
    make_slope_panel(
      title = "All species",
      y_levels = "All species",
      slope_limit = slope_limit,
      legend_title = effect_scale_legend_titles[[effect_scale_name]],
      show_x_labels = FALSE
    )

  thermal_panel <- slope_data |>
    filter(trait_family == "Thermal affinity") |>
    make_slope_panel(
      title = "Thermal affinity groups",
      y_levels = thermal_levels,
      slope_limit = slope_limit,
      legend_title = effect_scale_legend_titles[[effect_scale_name]],
      show_x_labels = FALSE
    )

  migration_panel <- slope_data |>
    filter(trait_family == "Migratory strategy") |>
    make_slope_panel(
      title = "Migratory strategies",
      y_levels = rev(migration_levels),
      slope_limit = slope_limit,
      legend_title = effect_scale_legend_titles[[effect_scale_name]],
      show_x_labels = FALSE
    )

  guild_panel <- slope_data |>
    filter(trait_family == "Foraging guild") |>
    make_slope_panel(
      title = "Foraging guilds",
      y_levels = rev(guild_levels_by_scale[[effect_scale_name]]),
      slope_limit = slope_limit,
      legend_title = effect_scale_legend_titles[[effect_scale_name]],
      show_x_labels = TRUE
    )

  beta_slope_panels <- all_species_panel / thermal_panel / migration_panel / guild_panel +
    plot_layout(heights = c(0.42, 0.8, 0.9, 3.8), guides = "collect") +
    plot_annotation(
      title = paste0("Trait-group trends in HMSC beta parameters: ", effect_scale_name),
      subtitle = effect_scale_subtitles[[effect_scale_name]],
      tag_levels = "A"
    ) &
    theme(
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey35"),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  output_slugs <- c(output_slug, also_save_as)
  output_paths <- file.path(out_dir, paste0(pattern, "-", output_slugs, ".png"))

  walk(output_paths, function(output_path) {
    ggsave(
      output_path,
      beta_slope_panels,
      width = 11,
      height = 16,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  })

  invisible(output_paths)
}

slope_output_paths <- c(
  make_slope_figure(
    effect_scale_name = "Unscaled beta",
    output_slug = "beta-trait-group-temporal-slope-panels-unscaled-beta"
  ),
  make_slope_figure(
    effect_scale_name = "SD-scaled effect size",
    output_slug = "beta-trait-group-temporal-slope-panels-sd-scaled-effect-size",
    also_save_as = "beta-trait-group-temporal-slope-panels"
  )
)

# ---- Period-by-period trend heatmaps ----
message("Plotting period-level trend heatmaps...")
walk(effect_scale_order, function(effect_scale_name) {
  scale_slug <- effect_scale_slugs[[effect_scale_name]]
  legacy_slug <- function(slug) {
    if (effect_scale_name == "SD-scaled effect size") slug else character()
  }

  make_period_heatmap(
    beta_group_period_summary,
    trait_family_name = "All species",
    y_levels = "All species",
    output_slug = paste0("beta-period-trends-all-species-", scale_slug),
    plot_height = 3.4,
    effect_scale_name = effect_scale_name,
    also_save_as = legacy_slug("beta-period-trends-all-species")
  )

  make_period_heatmap(
    beta_group_period_summary,
    trait_family_name = "Thermal affinity",
    y_levels = thermal_levels,
    output_slug = paste0("beta-period-trends-thermal-affinity-", scale_slug),
    plot_height = 6,
    effect_scale_name = effect_scale_name,
    also_save_as = legacy_slug("beta-period-trends-thermal-affinity")
  )

  make_period_heatmap(
    beta_group_period_summary,
    trait_family_name = "Migratory strategy",
    y_levels = rev(migration_levels),
    output_slug = paste0("beta-period-trends-migratory-strategy-", scale_slug),
    plot_height = 6.2,
    effect_scale_name = effect_scale_name,
    also_save_as = legacy_slug("beta-period-trends-migratory-strategy")
  )

  make_period_heatmap(
    beta_group_period_summary,
    trait_family_name = "Foraging guild",
    y_levels = rev(guild_levels_by_scale[[effect_scale_name]]),
    output_slug = paste0("beta-period-trends-foraging-guild-", scale_slug),
    plot_height = 13.5,
    effect_scale_name = effect_scale_name,
    also_save_as = legacy_slug("beta-period-trends-foraging-guild")
  )
})

message("Saved outputs under: ", out_dir)
