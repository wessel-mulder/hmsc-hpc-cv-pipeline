rm(list = ls())

# Figure 8 draft: species-level environmental dependency structure of 1970s
# avian bioregion clusters.
#
# This is a lightweight plotting script. It does not rerun the HMSC models,
# recompute VP, or refit the PCA. Instead, it reads the saved exploratory outputs
# from `z-avian-bioregion-cluster-vp-effects.Rmd` and rebuilds the manuscript
# candidate figure:
#
#   1. Weighted waffles of supported positive/negative environmental dependency.
#   2. Species dependency profiles projected into the Cluster 3 reference PCA.
#   3. Environmental loading plot showing which positive/negative features define
#      the PCA axes.
#
# The saved inputs are deliberately explicit CSVs, so the plotted quantities can
# be inspected without stepping through the exploratory notebook.

# ---- Packages ----
# `tidyverse` handles tabular reshaping and plotting, `patchwork` assembles the
# multi-panel figure, and `scales` formats percentages on legends/axes.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, patchwork, scales, grid, ggsankey, cowplot, ggrepel)

source(file.path("support_scripts", "figure_data_helpers.R"))

# ---- Paths ----
# The pattern and period identify the saved outputs created by the exploratory
# bioregion-cluster VP/Beta notebook.
pattern <- "2026-03-13"
period <- "1970s"
input_prefix <- paste0(pattern, "-", period)
base_dir <- "HmscOutputs"

in_dir <- file.path("notebooks", "exploratory", "outputs", "bioregion-cluster-vp-effects")
out_dir <- file.path("misc-figures", "outputs", "main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

waffle_tile_path <- file.path(in_dir, paste0(input_prefix, "-cluster-species-driver-waffle-tiles.csv"))
waffle_summary_path <- file.path(in_dir, paste0(input_prefix, "-cluster-species-driver-waffle-summary.csv"))
pca_score_path <- file.path(in_dir, paste0(input_prefix, "-cluster3-reference-dependency-pca-scores.csv"))
pca_loading_path <- file.path(in_dir, paste0(input_prefix, "-cluster3-reference-dependency-pca-loadings-annotated.csv"))

required_paths <- c(waffle_tile_path, waffle_summary_path, pca_score_path, pca_loading_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop(
    "Missing required saved output(s):\n",
    paste(missing_paths, collapse = "\n"),
    "\nRender notebooks/exploratory/z-avian-bioregion-cluster-vp-effects.Rmd first.",
    call. = FALSE
  )
}

# ---- Shared labels and palettes ----
# These colours match the environmental variable colours used in the Fig. 2
# variable-importance figures and the exploratory bioregion VP/Beta notebook.
cluster_levels <- paste("Cluster", 1:3)
environment_order <- c(
  "Temperature",
  "Precipitation",
  "Land-use heterogeneity",
  "Urban (% coverage)",
  "Cropland (% coverage)",
  "Pasture (% coverage)",
  "Forest (% coverage)",
  "Grass/Shrubland (% coverage)"
)

environment_colours <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban (% coverage)" = "snow3",
  "Cropland (% coverage)" = "goldenrod1",
  "Pasture (% coverage)" = "darkorange",
  "Forest (% coverage)" = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2"
)

var_labels <- c(
  "tmean_breeding" = "Temperature",
  "prec_breeding" = "Precipitation",
  "hh" = "Land-use heterogeneity",
  "perc_urban" = "Urban (% coverage)",
  "perc_cropland" = "Cropland (% coverage)",
  "perc_pasture" = "Pasture (% coverage)",
  "perc_forest" = "Forest (% coverage)",
  "perc_grass_shrub" = "Grass/Shrubland (% coverage)"
)

# These cluster colours match the reordered bioregion map/PCA figures:
# Cluster 1 = light blue, Cluster 2 = light yellow, Cluster 3 = dark yellow.
cluster_palette <- c(
  "Cluster 1" = "#76B7D8",
  "Cluster 2" = "#F6E58D",
  "Cluster 3" = "#C49A00"
)

# ---- Load saved plotting inputs ----
# Waffle tiles: each row is one tile. Within each Cluster x Beta-sign waffle,
# the 100 tiles represent the relative composition of weighted supported
# environmental dependency. The tile colour is the environmental variable.
waffle_tiles <- read_csv(waffle_tile_path, show_col_types = FALSE) |>
  mutate(
    bioregion = factor(.data$bioregion, levels = cluster_levels),
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative")),
    variable_label = factor(.data$variable_label, levels = environment_order)
  )

# Waffle summary: this is the same information as the waffle tiles, but already
# aggregated to Cluster x Beta sign x variable. It is the right input for the
# sankey-bump version because the value column carries the proportional weighted
# variable importance within each cluster and sign.
waffle_summary <- read_csv(waffle_summary_path, show_col_types = FALSE) |>
  mutate(
    bioregion = factor(.data$bioregion, levels = cluster_levels),
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative")),
    variable_label = factor(.data$variable_label, levels = environment_order)
  )

# ---- 1970s versus 2010s dependency waffles ----
# The original manuscript candidate above uses the saved 1970s outputs. For the
# temporal diagnostic requested here, we recompute the same dependency summaries
# for Atlas 1 and Atlas 3. Cluster labels are still the fixed bioregion-PCA
# identities, but species weights and VP/Beta effects are period specific.
cluster_recode <- c(
  "Cluster 1" = "Cluster 1",
  "Cluster 3" = "Cluster 2",
  "Cluster 2" = "Cluster 3"
)

relabel_clusters <- function(x) {
  factor(unname(cluster_recode[as.character(x)]), levels = cluster_levels)
}

allocate_waffle_tiles <- function(data, n_tiles = 100) {
  data <- data |>
    arrange(.data$variable_label) |>
    mutate(
      raw_tiles = .data$dependency_share * n_tiles,
      tile_n = floor(.data$raw_tiles),
      remainder = .data$raw_tiles - .data$tile_n
    )

  missing_tiles <- n_tiles - sum(data$tile_n, na.rm = TRUE)
  if (missing_tiles > 0 && sum(data$dependency_share, na.rm = TRUE) > 0) {
    add_index <- data |>
      arrange(desc(.data$remainder), desc(.data$dependency_share)) |>
      slice_head(n = missing_tiles) |>
      mutate(add_tile = 1) |>
      select(variable_label, add_tile)

    data <- data |>
      left_join(add_index, by = "variable_label") |>
      mutate(tile_n = .data$tile_n + replace_na(.data$add_tile, 0)) |>
      select(-add_tile)
  }

  data |>
    filter(.data$tile_n > 0) |>
    uncount(.data$tile_n, .id = "tile_within_variable") |>
    arrange(.data$variable_label, .data$tile_within_variable) |>
    mutate(
      tile_id = row_number(),
      tile_col = ((.data$tile_id - 1) %% 10) + 1,
      tile_row = 10 - ((.data$tile_id - 1) %/% 10)
    )
}

allocate_partial_waffle_tiles <- function(data, n_tiles = 100) {
  target_tiles <- round(sum(data$dependency_share, na.rm = TRUE) * n_tiles)
  if (!is.finite(target_tiles) || target_tiles <= 0) {
    return(data |> filter(FALSE))
  }

  data <- data |>
    arrange(.data$variable_label) |>
    mutate(
      raw_tiles = .data$dependency_share * n_tiles,
      tile_n = floor(.data$raw_tiles),
      remainder = .data$raw_tiles - .data$tile_n
    )

  missing_tiles <- target_tiles - sum(data$tile_n, na.rm = TRUE)
  if (missing_tiles > 0 && sum(data$dependency_share, na.rm = TRUE) > 0) {
    add_index <- data |>
      arrange(desc(.data$remainder), desc(.data$dependency_share)) |>
      slice_head(n = missing_tiles) |>
      mutate(add_tile = 1) |>
      select(variable_label, add_tile)

    data <- data |>
      left_join(add_index, by = "variable_label") |>
      mutate(tile_n = .data$tile_n + replace_na(.data$add_tile, 0)) |>
      select(-add_tile)
  }

  data |>
    filter(.data$tile_n > 0) |>
    uncount(.data$tile_n, .id = "tile_within_variable") |>
    arrange(.data$variable_label, .data$tile_within_variable) |>
    mutate(
      tile_id = row_number(),
      tile_col = ((.data$tile_id - 1) %% 10) + 1,
      tile_row = 10 - ((.data$tile_id - 1) %/% 10)
    )
}

matching_folders <- figure_model_folders(pattern = pattern, base_dir = base_dir)
model_nums <- as.character(atlas_numbers(matching_folders))
names(matching_folders) <- model_nums

bioregion_assignments_all <- read_csv(
  file.path(
    "notebooks", "exploratory", "outputs", "bioregion-pca",
    paste0(pattern, "-bioregion-assignments.csv")
  ),
  show_col_types = FALSE
) |>
  mutate(atlas = as.character(.data$atlas))

vp_scaled_all <- load_vp_estimates(matching_folders, base_dir = base_dir, scaled = TRUE)
beta_effects_all <- read_parameter_effects(pattern, effect = "Beta", base_dir = base_dir)

build_period_dependency <- function(atlas_id, period_label) {
  prediction_path <- file.path(
    base_dir,
    matching_folders[[atlas_id]],
    "Results", "Preds", "predicted-values.rdata"
  )
  if (!file.exists(prediction_path)) {
    stop("Missing fitted-site prediction file: ", prediction_path, call. = FALSE)
  }

  predictions <- readRDS(prediction_path)
  pred_surveys <- rownames(predictions)

  assignments <- bioregion_assignments_all |>
    filter(.data$atlas == .env$atlas_id) |>
    mutate(
      bioregion_original = .data$bioregion,
      bioregion = relabel_clusters(.data$bioregion)
    ) |>
    select(survey, site, bioregion_original, bioregion)

  missing_prediction_assignments <- setdiff(pred_surveys, assignments$survey)
  missing_assignment_predictions <- setdiff(assignments$survey, pred_surveys)
  if (length(missing_prediction_assignments) > 0 ||
      length(missing_assignment_predictions) > 0) {
    stop(
      "Prediction rows and bioregion assignments do not align for Atlas ",
      atlas_id,
      ". Missing assignments: ", length(missing_prediction_assignments),
      "; missing predictions: ", length(missing_assignment_predictions),
      call. = FALSE
    )
  }

  assignments_ordered <- assignments |>
    slice(match(pred_surveys, .data$survey))

  cluster_species_sums <- rowsum(
    predictions,
    group = assignments_ordered$bioregion,
    reorder = FALSE
  )
  cluster_totals <- rowSums(cluster_species_sums, na.rm = TRUE)

  species_weights <- cluster_species_sums |>
    as.data.frame() |>
    rownames_to_column("bioregion") |>
    pivot_longer(
      cols = -bioregion,
      names_to = "species",
      values_to = "expected_occurrence_sum"
    ) |>
    mutate(
      period = period_label,
      atlas = atlas_id,
      bioregion = factor(.data$bioregion, levels = cluster_levels),
      species_weight = .data$expected_occurrence_sum / cluster_totals[as.character(.data$bioregion)]
    )

  vp_long <- vp_scaled_all[[atlas_id]] |>
    as.data.frame() |>
    rownames_to_column("variable") |>
    filter(.data$variable %in% names(var_labels)) |>
    pivot_longer(
      cols = -variable,
      names_to = "species",
      values_to = "scaled_vp"
    ) |>
    mutate(
      variable_label = unname(var_labels[.data$variable]),
      variable_label = factor(.data$variable_label, levels = environment_order)
    )

  beta_period <- beta_effects_all |>
    filter(
      .data$atlas == as.numeric(.env$atlas_id),
      .data$variable %in% names(var_labels)
    ) |>
    mutate(
      supported = !is.na(.data$effect_size),
      effect_sign = case_when(
        .data$effect_size > 0 ~ 1,
        .data$effect_size < 0 ~ -1,
        TRUE ~ NA_real_
      ),
      effect_direction = case_when(
        .data$effect_size > 0 ~ "Positive",
        .data$effect_size < 0 ~ "Negative",
        TRUE ~ "Unsupported"
      )
    )

  species_weights |>
    inner_join(vp_long, by = "species", relationship = "many-to-many") |>
    left_join(
      beta_period |>
        select(species, variable, effect_size, supported, effect_sign, effect_direction),
      by = c("species", "variable")
    ) |>
    mutate(
      supported = replace_na(.data$supported, FALSE),
      effect_direction = replace_na(.data$effect_direction, "Unsupported"),
      weighted_scaled_vp = .data$species_weight * .data$scaled_vp,
      weighted_supported_vp = if_else(.data$supported, .data$weighted_scaled_vp, 0)
    )
}

temporal_dependency_data <- bind_rows(
  build_period_dependency("1", "1970s"),
  build_period_dependency("3", "2010s")
) |>
  filter(.data$supported, .data$effect_direction %in% c("Positive", "Negative")) |>
  mutate(
    period = factor(.data$period, levels = c("1970s", "2010s")),
    bioregion = factor(.data$bioregion, levels = cluster_levels),
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative")),
    variable_label = factor(.data$variable_label, levels = environment_order)
  )

temporal_waffle_summary <- temporal_dependency_data |>
  group_by(.data$period, .data$atlas, .data$bioregion, .data$effect_direction, .data$variable_label) |>
  summarise(
    n_supported_species = n_distinct(.data$species),
    weighted_dependency = sum(.data$weighted_supported_vp, na.rm = TRUE),
    mean_species_weight = mean(.data$species_weight, na.rm = TRUE),
    .groups = "drop"
  ) |>
  complete(
    period = factor(c("1970s", "2010s"), levels = c("1970s", "2010s")),
    bioregion = factor(cluster_levels, levels = cluster_levels),
    effect_direction = factor(c("Positive", "Negative"), levels = c("Positive", "Negative")),
    variable_label = factor(environment_order, levels = environment_order),
    fill = list(
      atlas = NA_character_,
      n_supported_species = 0,
      weighted_dependency = 0,
      mean_species_weight = NA_real_
    )
  ) |>
  mutate(
    atlas = case_when(
      .data$period == "1970s" ~ "1",
      .data$period == "2010s" ~ "3",
      TRUE ~ .data$atlas
    )
  ) |>
  group_by(.data$period, .data$bioregion, .data$effect_direction) |>
  mutate(
    total_weighted_dependency = sum(.data$weighted_dependency, na.rm = TRUE),
    dependency_share = if_else(
      .data$total_weighted_dependency > 0,
      .data$weighted_dependency / .data$total_weighted_dependency,
      0
    )
  ) |>
  ungroup()

temporal_waffle_tiles <- temporal_waffle_summary |>
  group_by(.data$period, .data$bioregion, .data$effect_direction) |>
  group_modify(~allocate_waffle_tiles(.x, n_tiles = 100)) |>
  ungroup() |>
  mutate(
    period = factor(.data$period, levels = c("1970s", "2010s")),
    bioregion = factor(.data$bioregion, levels = cluster_levels),
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative")),
    variable_label = factor(.data$variable_label, levels = environment_order),
    row_label = factor(
      paste(.data$period, .data$effect_direction, sep = " | "),
      levels = c(
        "1970s | Positive", "1970s | Negative",
        "2010s | Positive", "2010s | Negative"
      )
    )
  )

temporal_waffle_plot <- temporal_waffle_tiles |>
  ggplot(aes(x = .data$tile_col, y = .data$tile_row, fill = .data$variable_label)) +
  geom_tile(
    width = 0.92,
    height = 0.92,
    colour = "white",
    linewidth = 0.2
  ) +
  facet_grid(.data$row_label ~ .data$bioregion) +
  scale_fill_manual(
    values = environment_colours,
    breaks = environment_order,
    drop = FALSE,
    name = "Environmental variable"
  ) +
  coord_equal(clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0), breaks = NULL) +
  scale_y_continuous(expand = expansion(mult = 0), breaks = NULL) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.text.y = element_text(angle = 0),
    panel.spacing = unit(0.48, "lines"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.35, "cm")
  )

# PCA scores: each row is a Cluster x species profile projected into the PCA
# fitted on Cluster 3 species. Species weights are retained only for point size.
pca_scores <- read_csv(pca_score_path, show_col_types = FALSE) |>
  mutate(
    bioregion = factor(.data$bioregion, levels = cluster_levels),
    is_reference_cluster = .data$bioregion == "Cluster 3"
  )

# PCA loadings: each row is one feature, e.g. "Temperature | Positive". The
# annotated output already separates the environmental variable from the Beta
# sign, which keeps this plotting script small and checkable.
pca_loadings <- read_csv(pca_loading_path, show_col_types = FALSE) |>
  mutate(
    variable_label = factor(.data$variable_label, levels = environment_order),
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative")),
    loading_strength = sqrt(.data$PC1^2 + .data$PC2^2),
    loading_label = if_else(
      is.na(.data$loading_label),
      paste0(.data$variable_label, if_else(.data$effect_direction == "Positive", " +", " -")),
      .data$loading_label
    )
  )

# ---- PCA axis labels ----
# The PCA was fitted on Cluster 3 in the exploratory notebook. Because the score
# CSV stores all projected PCs, the explained variance can be recovered from the
# variance of Cluster 3 scores across all PC axes.
pca_cols <- names(pca_scores)[str_detect(names(pca_scores), "^PC[0-9]+$")]
cluster3_pca_variance <- pca_scores |>
  filter(.data$bioregion == "Cluster 3") |>
  select(all_of(pca_cols)) |>
  summarise(across(everything(), ~var(.x, na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "pc", values_to = "variance") |>
  mutate(proportion = .data$variance / sum(.data$variance, na.rm = TRUE))

pca_axis_labels <- c(
  paste0("PC1 (", percent(cluster3_pca_variance$proportion[cluster3_pca_variance$pc == "PC1"], accuracy = 0.1), ")"),
  paste0("PC2 (", percent(cluster3_pca_variance$proportion[cluster3_pca_variance$pc == "PC2"], accuracy = 0.1), ")")
)

# ---- Waffle panel ----
# Layout requested by the user: clusters are columns and supported Beta sign is
# the row. Each small square is one proportional tile, so tile count rather than
# point size carries the weighted dependency share.
waffle_plot <- waffle_tiles |>
  ggplot(aes(x = .data$tile_col, y = .data$tile_row, fill = .data$variable_label)) +
  geom_tile(
    width = 0.92,
    height = 0.92,
    colour = "white",
    linewidth = 0.2
  ) +
  facet_grid(.data$effect_direction ~ .data$bioregion) +
  scale_fill_manual(
    values = environment_colours,
    breaks = environment_order,
    drop = FALSE,
    name = "Environmental variable"
  ) +
  coord_equal(clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0), breaks = NULL) +
  scale_y_continuous(expand = expansion(mult = 0), breaks = NULL) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.55, "lines"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.35, "cm")
  )

# ---- Corrected 1970s dependency-balance panels ----
# The original waffle above normalizes the positive and negative rows
# separately. That is useful for composition within each sign, but it makes the
# two signs look equally important by construction. The corrected summaries
# below normalize across both signs within each cluster, so positive and
# negative effects are on one shared scale.
corrected_dependency_summary <- waffle_summary |>
  mutate(
    weighted_dependency = replace_na(.data$weighted_dependency, 0),
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative")),
    variable_label = factor(.data$variable_label, levels = environment_order)
  ) |>
  group_by(.data$bioregion) |>
  mutate(
    cluster_total_dependency = sum(.data$weighted_dependency, na.rm = TRUE),
    cluster_dependency_share = if_else(
      .data$cluster_total_dependency > 0,
      .data$weighted_dependency / .data$cluster_total_dependency,
      0
    )
  ) |>
  ungroup()

corrected_direction_labels <- corrected_dependency_summary |>
  group_by(.data$bioregion, .data$effect_direction) |>
  summarise(
    direction_share = sum(.data$cluster_dependency_share, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    bioregion = factor(.data$bioregion, levels = cluster_levels),
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative")),
    label = percent(.data$direction_share, accuracy = 1),
    tile_col = 10.8,
    tile_row = 9.2
  )

corrected_waffle_summary <- corrected_dependency_summary |>
  left_join(
    corrected_direction_labels |>
      select(.data$bioregion, .data$effect_direction, .data$direction_share),
    by = c("bioregion", "effect_direction")
  ) |>
  group_by(.data$bioregion, .data$effect_direction) |>
  mutate(
    # Reuse the existing tile allocator, but pass the cluster-level share.
    dependency_share = .data$cluster_dependency_share
  ) |>
  ungroup()

corrected_waffle_tiles <- corrected_waffle_summary |>
  group_by(.data$bioregion, .data$effect_direction) |>
  group_modify(~allocate_partial_waffle_tiles(.x, n_tiles = 100)) |>
  ungroup() |>
  mutate(
    bioregion = factor(.data$bioregion, levels = cluster_levels),
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative")),
    variable_label = factor(.data$variable_label, levels = environment_order)
  )

corrected_waffle_plot <- corrected_waffle_tiles |>
  ggplot(aes(x = .data$tile_col, y = .data$tile_row, fill = .data$variable_label)) +
  geom_tile(
    width = 0.92,
    height = 0.92,
    colour = "white",
    linewidth = 0.2
  ) +
  geom_text(
    data = corrected_direction_labels,
    aes(x = .data$tile_col, y = .data$tile_row, label = .data$label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.2,
    fontface = "bold"
  ) +
  facet_grid(.data$effect_direction ~ .data$bioregion) +
  scale_fill_manual(
    values = environment_colours,
    breaks = environment_order,
    drop = FALSE,
    name = "Environmental variable"
  ) +
  coord_equal(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)), breaks = NULL) +
  scale_y_continuous(expand = expansion(mult = 0), breaks = NULL) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.text.y = element_text(angle = 0),
    panel.spacing = unit(0.55, "lines"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.35, "cm")
  )

diverging_dependency_data <- corrected_dependency_summary |>
  mutate(
    signed_dependency_share = if_else(
      .data$effect_direction == "Positive",
      .data$cluster_dependency_share,
      -.data$cluster_dependency_share
    )
  )

diverging_dependency_plot <- diverging_dependency_data |>
  ggplot(aes(
    x = .data$bioregion,
    y = .data$signed_dependency_share,
    fill = .data$variable_label
  )) +
  geom_hline(yintercept = 0, colour = "grey25", linewidth = 0.35) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = environment_colours,
    breaks = environment_order,
    drop = FALSE,
    name = "Environmental variable"
  ) +
  scale_y_continuous(
    labels = function(x) percent(abs(x), accuracy = 1),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Share of total supported dependency"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.35, "cm")
  )

# ---- Variable-lane corrected waffle ----
# This version keeps the corrected one-budget-per-cluster interpretation, but
# makes the plot easier to audit. Every environmental variable gets a fixed
# lane within the positive and negative panels. If a variable needs more than
# ten tiles, the extra tiles wrap to a second lane directly underneath. Empty
# lanes therefore mean no supported weighted dependency for that variable.
allocate_cluster_feature_tiles <- function(data, n_tiles = 100) {
  data <- data |>
    arrange(.data$effect_direction, .data$variable_label) |>
    mutate(
      raw_tiles = .data$cluster_dependency_share * n_tiles,
      tile_n = floor(.data$raw_tiles),
      remainder = .data$raw_tiles - .data$tile_n
    )

  missing_tiles <- n_tiles - sum(data$tile_n, na.rm = TRUE)
  if (missing_tiles > 0 && sum(data$cluster_dependency_share, na.rm = TRUE) > 0) {
    add_index <- data |>
      arrange(desc(.data$remainder), desc(.data$cluster_dependency_share)) |>
      slice_head(n = missing_tiles) |>
      mutate(add_tile = 1) |>
      select(effect_direction, variable_label, add_tile)

    data <- data |>
      left_join(add_index, by = c("effect_direction", "variable_label")) |>
      mutate(tile_n = .data$tile_n + replace_na(.data$add_tile, 0)) |>
      select(-add_tile)
  }

  data
}

variable_lane_counts <- corrected_dependency_summary |>
  mutate(
    bioregion = factor(.data$bioregion, levels = cluster_levels),
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative")),
    variable_label = factor(.data$variable_label, levels = environment_order)
  ) |>
  group_by(.data$bioregion) |>
  group_modify(~allocate_cluster_feature_tiles(.x, n_tiles = 100)) |>
  ungroup()

variable_lane_layout <- variable_lane_counts |>
  group_by(.data$variable_label) |>
  summarise(
    max_tiles = max(.data$tile_n, na.rm = TRUE),
    rows_needed = pmax(1, ceiling(.data$max_tiles / 10)),
    .groups = "drop"
  ) |>
  mutate(
    variable_label = factor(.data$variable_label, levels = environment_order)
  ) |>
  arrange(.data$variable_label) |>
  mutate(
    row_start = sum(.data$rows_needed) - lag(cumsum(.data$rows_needed), default = 0),
    label_y = .data$row_start - (.data$rows_needed - 1) / 2
  )

variable_lane_background <- variable_lane_layout |>
  mutate(row_offsets = map(.data$rows_needed, ~seq.int(0, .x - 1))) |>
  unnest(.data$row_offsets) |>
  crossing(
    bioregion = factor(cluster_levels, levels = cluster_levels),
    effect_direction = factor(c("Positive", "Negative"), levels = c("Positive", "Negative")),
    tile_col = seq_len(10)
  ) |>
  mutate(tile_y = .data$row_start - .data$row_offsets)

variable_lane_tiles <- variable_lane_counts |>
  filter(.data$tile_n > 0) |>
  uncount(.data$tile_n, .id = "tile_within_feature") |>
  left_join(
    variable_lane_layout |> select(.data$variable_label, .data$row_start),
    by = "variable_label"
  ) |>
  mutate(
    tile_col = ((.data$tile_within_feature - 1) %% 10) + 1,
    row_offset = (.data$tile_within_feature - 1) %/% 10,
    tile_y = .data$row_start - .data$row_offset,
    variable_label = factor(.data$variable_label, levels = environment_order)
  )

variable_lane_waffle_plot <- ggplot() +
  geom_tile(
    data = variable_lane_background,
    aes(x = .data$tile_col, y = .data$tile_y),
    fill = NA,
    colour = "grey90",
    linewidth = 0.18,
    width = 0.88,
    height = 0.88
  ) +
  geom_tile(
    data = variable_lane_tiles,
    aes(x = .data$tile_col, y = .data$tile_y, fill = .data$variable_label),
    colour = "white",
    linewidth = 0.2,
    width = 0.88,
    height = 0.88
  ) +
  facet_grid(.data$effect_direction ~ .data$bioregion, switch = "y") +
  scale_fill_manual(
    values = environment_colours,
    breaks = environment_order,
    drop = FALSE,
    name = "Environmental variable"
  ) +
  scale_x_continuous(expand = expansion(mult = 0), breaks = NULL) +
  scale_y_continuous(
    breaks = variable_lane_layout$label_y,
    labels = as.character(variable_lane_layout$variable_label),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_equal(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8.8),
    axis.title = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.text.y.left = element_text(angle = 0),
    strip.placement = "outside",
    panel.spacing = unit(0.65, "lines"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.35, "cm")
  )

# ---- Mirrored directionality bars ----
# This draft keeps the corrected one-budget-per-cluster scale, but replaces the
# waffle grid with a central axis for each cluster. Positive supported
# dependency extends left of zero, while negative supported dependency extends
# right. The environmental colours are retained so the panel can be read beside
# the waffle versions.
mirrored_dependency_data <- corrected_dependency_summary |>
  mutate(
    variable_label = factor(.data$variable_label, levels = rev(environment_order)),
    signed_axis_share = if_else(
      .data$effect_direction == "Positive",
      -.data$cluster_dependency_share,
      .data$cluster_dependency_share
    )
  )

mirrored_dependency_limit <- max(abs(mirrored_dependency_data$signed_axis_share), na.rm = TRUE)

mirrored_dependency_plot <- mirrored_dependency_data |>
  ggplot(aes(
    x = .data$signed_axis_share,
    y = .data$variable_label,
    fill = .data$variable_label
  )) +
  geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.45) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.2) +
  facet_wrap(~bioregion, nrow = 1) +
  annotate(
    "text",
    x = -mirrored_dependency_limit * 0.86,
    y = length(environment_order) + 0.7,
    label = "Positive",
    hjust = 0.5,
    size = 3.2,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = mirrored_dependency_limit * 0.86,
    y = length(environment_order) + 0.7,
    label = "Negative",
    hjust = 0.5,
    size = 3.2,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = environment_colours,
    breaks = environment_order,
    drop = FALSE,
    guide = "none"
  ) +
  scale_x_continuous(
    labels = function(x) percent(abs(x), accuracy = 1),
    limits = c(-mirrored_dependency_limit, mirrored_dependency_limit),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Share of total supported dependency",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9),
    plot.margin = margin(18, 8, 8, 8)
  )

# ---- Sign-split pie charts ----
# Each cluster has one pie. The two main sectors are positive and negative
# supported dependency, and the environmental variables are the coloured slices
# inside those sectors. Labels at the sector midpoints show where the positive
# and negative parts of the pie sit.
pie_dependency_data <- corrected_dependency_summary |>
  mutate(
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative")),
    variable_label = factor(.data$variable_label, levels = environment_order)
  ) |>
  arrange(.data$bioregion, .data$effect_direction, .data$variable_label) |>
  group_by(.data$bioregion) |>
  mutate(
    slice_end = cumsum(.data$cluster_dependency_share),
    slice_start = lag(.data$slice_end, default = 0),
    slice_mid = (.data$slice_start + .data$slice_end) / 2
  ) |>
  ungroup()

pie_direction_labels <- pie_dependency_data |>
  group_by(.data$bioregion, .data$effect_direction) |>
  summarise(
    direction_start = min(.data$slice_start, na.rm = TRUE),
    direction_end = max(.data$slice_end, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    label_position = (.data$direction_start + .data$direction_end) / 2,
    label = as.character(.data$effect_direction)
  )

pie_dependency_plot <- pie_dependency_data |>
  ggplot(aes(
    x = 1,
    y = .data$cluster_dependency_share,
    fill = .data$variable_label
  )) +
  geom_col(
    width = 1,
    colour = "white",
    linewidth = 0.35
  ) +
  geom_hline(
    data = pie_direction_labels,
    aes(yintercept = .data$direction_start),
    inherit.aes = FALSE,
    colour = "white",
    linewidth = 0.9
  ) +
  geom_text(
    data = pie_direction_labels,
    aes(x = 1.62, y = .data$label_position, label = .data$label),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "bold"
  ) +
  facet_wrap(~bioregion, nrow = 1) +
  coord_polar(theta = "y", clip = "off") +
  scale_fill_manual(
    values = environment_colours,
    breaks = environment_order,
    drop = FALSE,
    name = "Environmental variable"
  ) +
  xlim(0.2, 1.75) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.35, "cm"),
    plot.margin = margin(8, 28, 8, 8)
  )

# ---- Nested sign-split donut charts ----
# The simple pie makes the sign sections hard to see because the environmental
# colours dominate. This nested version gives sign its own inner ring and keeps
# the environmental-variable slices on the outside. Positive/negative is thus
# visible as a structural layer rather than only as two labels.
nested_pie_direction_data <- pie_direction_labels |>
  mutate(
    direction_share = .data$direction_end - .data$direction_start,
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative"))
  )

sign_colours <- c(
  "Positive" = "grey25",
  "Negative" = "grey78"
)

nested_pie_plot <- ggplot() +
  geom_col(
    data = nested_pie_direction_data,
    aes(
      x = 0.78,
      y = .data$direction_share,
      fill = .data$effect_direction
    ),
    width = 0.25,
    colour = "white",
    linewidth = 0.45
  ) +
  geom_col(
    data = pie_dependency_data,
    aes(
      x = 1.08,
      y = .data$cluster_dependency_share,
      fill = .data$variable_label
    ),
    width = 0.34,
    colour = "white",
    linewidth = 0.28
  ) +
  geom_hline(
    data = pie_direction_labels,
    aes(yintercept = .data$direction_start),
    inherit.aes = FALSE,
    colour = "white",
    linewidth = 1.1
  ) +
  geom_text(
    data = pie_direction_labels,
    aes(x = 1.52, y = .data$label_position, label = .data$label),
    inherit.aes = FALSE,
    size = 3.1,
    fontface = "bold"
  ) +
  facet_wrap(~bioregion, nrow = 1) +
  coord_polar(theta = "y", clip = "off") +
  scale_fill_manual(
    values = c(sign_colours, environment_colours),
    breaks = c("Positive", "Negative", environment_order),
    name = "Dependency layer"
  ) +
  xlim(0.35, 1.7) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.35, "cm"),
    plot.margin = margin(8, 28, 8, 8)
  )

# ---- Scaled semicircle dependency charts ----
# Positive and negative are shown as two separate half-circles. Their area is
# proportional to total supported dependency in that direction, while the
# coloured slices within each half-circle show the environmental-variable
# composition of that direction.
semicircle_dependency_data <- corrected_dependency_summary |>
  mutate(
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative")),
    variable_label = factor(.data$variable_label, levels = environment_order)
  ) |>
  arrange(.data$bioregion, .data$effect_direction, .data$variable_label) |>
  group_by(.data$bioregion, .data$effect_direction) |>
  mutate(
    direction_share = sum(.data$cluster_dependency_share, na.rm = TRUE),
    direction_slice_share = if_else(
      .data$direction_share > 0,
      .data$cluster_dependency_share / .data$direction_share,
      0
    ),
    direction_slice_end = cumsum(.data$direction_slice_share),
    direction_slice_start = lag(.data$direction_slice_end, default = 0),
    radius = sqrt(.data$direction_share),
    centre_x = if_else(.data$effect_direction == "Positive", -0.07, 0.07),
    centre_y = 0,
    angle_start = case_when(
      .data$effect_direction == "Positive" ~ pi / 2 + .data$direction_slice_start * pi,
      TRUE ~ -pi / 2 + .data$direction_slice_start * pi
    ),
    angle_end = case_when(
      .data$effect_direction == "Positive" ~ pi / 2 + .data$direction_slice_end * pi,
      TRUE ~ -pi / 2 + .data$direction_slice_end * pi
    ),
    slice_id = paste(.data$bioregion, .data$effect_direction, .data$variable_label, sep = " | ")
  ) |>
  ungroup()

build_semicircle_slice <- function(row, n_points = 60) {
  angles <- seq(row$angle_start, row$angle_end, length.out = n_points)
  tibble(
    bioregion = row$bioregion,
    effect_direction = row$effect_direction,
    variable_label = row$variable_label,
    slice_id = row$slice_id,
    direction_share = row$direction_share,
    x = c(row$centre_x, row$centre_x + row$radius * cos(angles), row$centre_x),
    y = c(row$centre_y, row$centre_y + row$radius * sin(angles), row$centre_y)
  )
}

semicircle_slice_polygons <- semicircle_dependency_data |>
  split(seq_len(nrow(semicircle_dependency_data))) |>
  map_dfr(~build_semicircle_slice(.x))

semicircle_direction_labels <- semicircle_dependency_data |>
  distinct(.data$bioregion, .data$effect_direction, .data$direction_share, .data$centre_x, .data$radius) |>
  mutate(
    label = as.character(.data$effect_direction),
    label_x = .data$centre_x + if_else(.data$effect_direction == "Positive", -.data$radius * 0.55, .data$radius * 0.55),
    label_y = -0.98
  )

semicircle_dependency_plot <- ggplot(semicircle_slice_polygons) +
  geom_polygon(
    aes(
      x = .data$x,
      y = .data$y,
      group = .data$slice_id,
      fill = .data$variable_label
    ),
    colour = "white",
    linewidth = 0.35
  ) +
  geom_segment(
    data = semicircle_direction_labels,
    aes(
      x = .data$centre_x,
      xend = .data$centre_x,
      y = -.data$radius,
      yend = .data$radius
    ),
    inherit.aes = FALSE,
    colour = "white",
    linewidth = 0.7
  ) +
  geom_text(
    data = semicircle_direction_labels,
    aes(x = .data$label_x, y = .data$label_y, label = .data$label),
    inherit.aes = FALSE,
    size = 3.1,
    fontface = "bold"
  ) +
  facet_wrap(~bioregion, nrow = 1) +
  scale_fill_manual(
    values = environment_colours,
    breaks = environment_order,
    drop = FALSE,
    name = "Environmental variable"
  ) +
  coord_equal(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.35, "cm"),
    plot.margin = margin(8, 18, 8, 8)
  )

# ---- 1970s-reference species dependency PCA and Atlas 1 cluster map ----
# To pair the dependency chart with a trajectory plot, we build one
# environmental-dependency profile per period x cluster x species. Each profile
# is a 16-feature vector: 8 variables x positive/negative supported effect
# direction. The PCA is fitted on all 1970s species profiles, and the 2010s
# species profiles are projected into that fixed 1970s dependency space. The
# arrows show species-weighted cluster centroids, so common/characteristic
# species have more influence on the plotted cluster trajectory.
dependency_feature_order <- as.vector(outer(
  c("Positive", "Negative"),
  environment_order,
  paste,
  sep = " | "
))

species_dependency_profile_long <- temporal_dependency_data |>
  mutate(
    period = factor(.data$period, levels = c("1970s", "2010s")),
    bioregion = factor(.data$bioregion, levels = cluster_levels),
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative")),
    variable_label = factor(.data$variable_label, levels = environment_order),
    dependency_feature = paste(.data$effect_direction, .data$variable_label, sep = " | "),
    dependency_value = .data$weighted_supported_vp
  ) |>
  group_by(.data$period, .data$bioregion, .data$species, .data$species_weight, .data$dependency_feature) |>
  summarise(
    dependency_value = sum(.data$dependency_value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  complete(
    nesting(period, bioregion, species, species_weight),
    dependency_feature = dependency_feature_order,
    fill = list(dependency_value = 0)
  )

species_dependency_matrix <- species_dependency_profile_long |>
  mutate(
    dependency_feature = factor(.data$dependency_feature, levels = dependency_feature_order)
  ) |>
  pivot_wider(
    id_cols = c(.data$period, .data$bioregion, .data$species, .data$species_weight),
    names_from = .data$dependency_feature,
    values_from = .data$dependency_value,
    values_fill = 0
  ) |>
  arrange(.data$period, .data$bioregion, .data$species)

dependency_feature_cols <- intersect(dependency_feature_order, names(species_dependency_matrix))
reference_dependency_matrix <- species_dependency_matrix |>
  filter(.data$period == "1970s") |>
  select(all_of(dependency_feature_cols)) |>
  as.matrix()

dependency_feature_sds <- apply(reference_dependency_matrix, 2, sd, na.rm = TRUE)
dependency_pca_cols <- names(dependency_feature_sds)[
  is.finite(dependency_feature_sds) & dependency_feature_sds > 0
]

dependency_pca <- prcomp(
  reference_dependency_matrix[, dependency_pca_cols, drop = FALSE],
  center = TRUE,
  scale. = TRUE
)

dependency_pca_scores <- predict(
  dependency_pca,
  newdata = species_dependency_matrix |>
    select(all_of(dependency_pca_cols)) |>
    as.matrix()
) |>
  as.data.frame() |>
  bind_cols(species_dependency_matrix |> select(.data$period, .data$bioregion, .data$species, .data$species_weight)) |>
  mutate(
    period = factor(.data$period, levels = c("1970s", "2010s")),
    bioregion = factor(.data$bioregion, levels = cluster_levels)
  )

dependency_pca_variance <- dependency_pca$sdev^2 / sum(dependency_pca$sdev^2)
dependency_pca_axis_labels <- c(
  paste0("PC1 (", percent(dependency_pca_variance[[1]], accuracy = 0.1), ")"),
  paste0("PC2 (", percent(dependency_pca_variance[[2]], accuracy = 0.1), ")")
)

dependency_pca_centroids <- dependency_pca_scores |>
  group_by(.data$period, .data$bioregion) |>
  summarise(
    PC1 = weighted.mean(.data$PC1, .data$species_weight, na.rm = TRUE),
    PC2 = weighted.mean(.data$PC2, .data$species_weight, na.rm = TRUE),
    n_species = n_distinct(.data$species),
    .groups = "drop"
  ) |>
  mutate(
    period = factor(.data$period, levels = c("1970s", "2010s")),
    bioregion = factor(.data$bioregion, levels = cluster_levels)
  )

dependency_pca_trajectory <- dependency_pca_centroids |>
  arrange(.data$bioregion, .data$period)

dependency_pca_x_limits <- range(
  c(
    quantile(dependency_pca_scores$PC1, c(0.01, 0.975), na.rm = TRUE),
    dependency_pca_centroids$PC1
  ),
  na.rm = TRUE
)
dependency_pca_y_limits <- range(
  c(
    quantile(dependency_pca_scores$PC2, c(0.01, 0.985), na.rm = TRUE),
    dependency_pca_centroids$PC2
  ),
  na.rm = TRUE
)
dependency_pca_x_pad <- diff(dependency_pca_x_limits) * 0.08
dependency_pca_y_pad <- diff(dependency_pca_y_limits) * 0.10
dependency_pca_x_limits <- dependency_pca_x_limits + c(-dependency_pca_x_pad, dependency_pca_x_pad)
dependency_pca_y_limits <- dependency_pca_y_limits + c(-dependency_pca_y_pad, dependency_pca_y_pad)

dependency_pca_plot <- ggplot(
  dependency_pca_scores,
  aes(x = .data$PC1, y = .data$PC2, colour = .data$bioregion)
) +
  geom_hline(yintercept = 0, colour = "grey86", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey86", linewidth = 0.3) +
  geom_point(alpha = 0.10, size = 0.75, show.legend = FALSE) +
  geom_path(
    data = dependency_pca_trajectory,
    aes(group = .data$bioregion),
    linewidth = 0.85,
    arrow = arrow(length = unit(0.13, "cm")),
    alpha = 0.86
  ) +
  geom_point(
    data = dependency_pca_centroids,
    aes(shape = .data$period),
    size = 3.8,
    stroke = 0.85
  ) +
  ggrepel::geom_text_repel(
    data = dependency_pca_centroids,
    aes(label = paste0(str_replace(as.character(.data$bioregion), "Cluster ", "C"), " ", .data$period)),
    size = 3,
    show.legend = FALSE,
    min.segment.length = 0,
    segment.colour = "grey70",
    segment.linewidth = 0.25,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf
  ) +
  scale_colour_manual(values = cluster_palette, drop = FALSE, name = "Cluster") +
  scale_shape_manual(values = c("1970s" = 16, "2010s" = 17), name = "Period") +
  labs(
    x = dependency_pca_axis_labels[[1]],
    y = dependency_pca_axis_labels[[2]]
  ) +
  coord_cartesian(
    xlim = dependency_pca_x_limits,
    ylim = dependency_pca_y_limits,
    clip = "off"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(7, 18, 7, 7)
  )

mainland_bbox <- c(xmin = 400000, xmax = 750000, ymin = 6000000, ymax = 6450000)
bornholm_bbox <- c(xmin = 855000, xmax = 905000, ymin = 6100000, ymax = 6160000)
mainland_width <- mainland_bbox[["xmax"]] - mainland_bbox[["xmin"]]
mainland_height <- mainland_bbox[["ymax"]] - mainland_bbox[["ymin"]]
bornholm_width <- bornholm_bbox[["xmax"]] - bornholm_bbox[["xmin"]]
bornholm_height <- bornholm_bbox[["ymax"]] - bornholm_bbox[["ymin"]]
bornholm_inset_width <- bornholm_width / mainland_width
bornholm_inset_height <- bornholm_height / mainland_height

map_theme <- function() {
  theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(0, 0, 0, 0)
    )
}

plot_cluster_map_base <- function(df, bbox, border = FALSE) {
  ggplot(df) +
    geom_point(aes(x = .data$X, y = .data$Y, colour = .data$bioregion), size = 1.25, alpha = 0.95) +
    scale_colour_manual(values = cluster_palette, drop = FALSE, name = "Cluster") +
    coord_fixed(
      xlim = c(bbox[["xmin"]], bbox[["xmax"]]),
      ylim = c(bbox[["ymin"]], bbox[["ymax"]]),
      expand = FALSE
    ) +
    map_theme() +
    theme(
      panel.border = if (border) {
        element_rect(colour = "grey35", fill = NA, linewidth = 0.55)
      } else {
        element_blank()
      }
    )
}

atlas1_cluster_map_data <- bioregion_assignments_all |>
  filter(.data$atlas == "1") |>
  mutate(bioregion = relabel_clusters(.data$bioregion))

atlas1_cluster_map <- {
  p_main <- plot_cluster_map_base(atlas1_cluster_map_data, mainland_bbox)
  p_inset <- plot_cluster_map_base(atlas1_cluster_map_data, bornholm_bbox, border = TRUE) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.7)
    )

  ggdraw(p_main) +
    draw_plot(
      p_inset,
      x = 1 - bornholm_inset_width - 0.2,
      y = 1 - bornholm_inset_height - 0.2,
      width = bornholm_inset_width,
      height = bornholm_inset_height
    )
}

dependency_semicircle_pca_map_composite <- plot_grid(
  semicircle_dependency_plot +
    theme(legend.position = "right"),
  plot_grid(
    dependency_pca_plot + theme(legend.position = "bottom"),
    atlas1_cluster_map,
    labels = c("B", "C"),
    label_size = 15,
    label_fontface = "bold",
    nrow = 1,
    rel_widths = c(1.15, 0.85)
  ),
  labels = c("A", ""),
  label_size = 15,
  label_fontface = "bold",
  ncol = 1,
  rel_heights = c(1, 1.05)
)

# ---- Sankey-bump dependency-flow panel ----
# This is an alternate display for the same positive/negative proportional
# dependency data used in the waffles. Instead of showing each Cluster x sign as
# a 100-tile block, the sankey-bump chart connects the same variable through
# Cluster 1 -> Cluster 2 -> Cluster 3. Wider ribbons mean that a larger share of
# the weighted supported dependency is assigned to that variable. If two
# clusters are similar, the variable ribbons maintain comparable thickness and
# rank; if they differ, the ribbons expand, contract, or cross.
sankey_dependency_data <- waffle_summary |>
  mutate(
    dependency_share = replace_na(.data$dependency_share, 0),
    # Keep all variables visible even when a dependency share is exactly zero.
    # A tiny plotting value avoids dropped ribbons without altering the CSV
    # values or the labels used for interpretation.
    plot_dependency_share = pmax(.data$dependency_share, 0.0001)
  )

sankey_bump_plot <- ggplot(
  sankey_dependency_data,
  aes(
    x = .data$bioregion,
    node = .data$variable_label,
    group = .data$variable_label,
    fill = .data$variable_label,
    value = .data$plot_dependency_share
  )
) +
  geom_sankey_bump(
    space = 0.01,
    type = "alluvial",
    alpha = 0.86,
    colour = "white",
    linewidth = 0.2
  ) +
  facet_grid(.data$effect_direction ~ .) +
  scale_fill_manual(
    values = environment_colours,
    breaks = environment_order,
    drop = FALSE,
    name = "Environmental variable"
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.04, 0.04))) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Proportional weighted dependency"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(face = "bold", angle = 90),
    axis.text.x = element_text(face = "bold"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.35, "cm")
  )

# ---- PCA score panel ----
# The large centroid circles are intentionally unweighted means. Since Cluster 3
# defines the PCA, its unweighted centroid sits at the reference centre. The
# small points remain sized by expected occurrence share so common species are
# visually emphasized without moving the mathematical centroid.
pca_centroids <- pca_scores |>
  group_by(.data$bioregion) |>
  summarise(
    PC1 = mean(.data$PC1, na.rm = TRUE),
    PC2 = mean(.data$PC2, na.rm = TRUE),
    .groups = "drop"
  )

# A few projected species can be far from the main cloud. The plotted panel uses
# a quantile zoom for readability; the score CSV still contains full values.
pca_x_limits <- quantile(pca_scores$PC1, c(0.01, 0.985), na.rm = TRUE)
pca_y_limits <- quantile(pca_scores$PC2, c(0.01, 0.99), na.rm = TRUE)

pca_score_plot <- pca_scores |>
  ggplot(aes(x = .data$PC1, y = .data$PC2, colour = .data$bioregion)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(
    aes(size = .data$species_weight, alpha = .data$is_reference_cluster),
    shape = 16
  ) +
  geom_point(
    data = pca_centroids,
    aes(x = .data$PC1, y = .data$PC2, fill = .data$bioregion),
    inherit.aes = FALSE,
    shape = 21,
    colour = "grey15",
    size = 5,
    stroke = 0.55
  ) +
  geom_text(
    data = pca_centroids,
    aes(x = .data$PC1, y = .data$PC2, label = .data$bioregion),
    inherit.aes = FALSE,
    position = position_nudge(
      x = c(0.45, -0.4, -0.45),
      y = c(0.28, -0.32, 0.32)
    ),
    size = 3.3,
    fontface = "bold"
  ) +
  scale_colour_manual(values = cluster_palette, name = "Cluster") +
  scale_fill_manual(values = cluster_palette, guide = "none") +
  scale_alpha_manual(values = c("TRUE" = 0.75, "FALSE" = 0.32), guide = "none") +
  scale_size_continuous(
    range = c(0.8, 4.5),
    labels = percent_format(accuracy = 0.1),
    name = "Expected occurrence\nshare"
  ) +
  coord_cartesian(xlim = pca_x_limits, ylim = pca_y_limits, clip = "off") +
  labs(x = pca_axis_labels[[1]], y = pca_axis_labels[[2]]) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# ---- PCA loading panel ----
# Loadings are shown in the same PC1/PC2 coordinate system as the score panel,
# but as a separate small biplot so the score cloud is not overprinted. Arrows
# from the origin show feature direction. Solid/filled points are positive Beta
# features, while dotted/open points are negative Beta features.
pca_loading_limit <- max(abs(c(pca_loadings$PC1, pca_loadings$PC2)), na.rm = TRUE) * 1.18

pca_loading_plot <- ggplot(pca_loadings, aes(x = .data$PC1, y = .data$PC2)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_segment(
    aes(
      x = 0,
      y = 0,
      xend = .data$PC1,
      yend = .data$PC2,
      colour = .data$variable_label,
      linetype = .data$effect_direction
    ),
    linewidth = 0.45,
    arrow = arrow(length = unit(0.08, "in"))
  ) +
  geom_point(
    data = filter(pca_loadings, .data$effect_direction == "Positive"),
    aes(size = .data$loading_strength, colour = .data$variable_label, fill = .data$variable_label),
    shape = 21,
    stroke = 0.5
  ) +
  geom_point(
    data = filter(pca_loadings, .data$effect_direction == "Negative"),
    aes(size = .data$loading_strength, colour = .data$variable_label),
    shape = 21,
    fill = "white",
    stroke = 0.8
  ) +
  geom_text(
    aes(label = .data$loading_label),
    size = 2.65,
    hjust = if_else(pca_loadings$PC1 >= 0, -0.08, 1.08),
    vjust = if_else(pca_loadings$PC2 >= 0, -0.15, 1.15),
    check_overlap = TRUE
  ) +
  scale_colour_manual(values = environment_colours, guide = "none") +
  scale_fill_manual(values = environment_colours, guide = "none") +
  scale_linetype_manual(values = c("Positive" = "solid", "Negative" = "dotted"), guide = "none") +
  scale_size_continuous(range = c(2, 7), guide = "none") +
  coord_equal(
    xlim = c(-pca_loading_limit, pca_loading_limit),
    ylim = c(-pca_loading_limit, pca_loading_limit),
    clip = "off"
  ) +
  labs(
    title = "Environmental loadings",
    subtitle = "Solid/filled = positive Beta feature; dotted/open = negative Beta feature.",
    x = pca_axis_labels[[1]],
    y = pca_axis_labels[[2]]
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(5.5, 18, 5.5, 5.5)
  )

# ---- Composite figure ----
# The requested layout puts the waffle section on top. The PCA score panel and
# PCA loading panel are placed underneath, side-by-side, so the reader can see
# where clusters fall and which environmental dependency features define the
# axes.
pca_row <- pca_score_plot | pca_loading_plot + plot_layout(widths = c(1, 0.72))
composite_plot <- waffle_plot / pca_row + plot_layout(heights = c(1.2, 1), guides = "keep")

# The alternate composite swaps the waffle blocks for the sankey-bump panel.
# This version emphasizes how the same environmental dependency categories
# change their proportional importance across the three 1970s communities.
sankey_composite_plot <- sankey_bump_plot / pca_row +
  plot_layout(heights = c(1.05, 1), guides = "keep")

# ---- Alternative Cluster 3 reference displays ----
# The sankey-bump chart shows the full flow of dependency proportions across
# clusters, but the manuscript argument may be sharper if Cluster 3 is treated
# explicitly as the reference. The next three plots all ask the same question in
# different visual forms:
#
#   "How far are Cluster 1 and Cluster 2 from Cluster 3 in species-level
#    environmental dependency structure?"
#
# These are intentionally exploratory alternatives. They all use the same
# `dependency_share` values as the waffles and sankey-bump chart.

comparison_palette <- c(
  "Cluster 1 - Cluster 3" = "#76B7D8",
  "Cluster 2 - Cluster 3" = "#C49A00"
)

reference_delta_data <- waffle_summary |>
  select(bioregion, effect_direction, variable_label, dependency_share) |>
  pivot_wider(
    names_from = bioregion,
    values_from = dependency_share
  ) |>
  mutate(
    `Cluster 1 - Cluster 3` = .data$`Cluster 1` - .data$`Cluster 3`,
    `Cluster 2 - Cluster 3` = .data$`Cluster 2` - .data$`Cluster 3`
  ) |>
  pivot_longer(
    cols = c("Cluster 1 - Cluster 3", "Cluster 2 - Cluster 3"),
    names_to = "comparison",
    values_to = "signed_delta"
  ) |>
  mutate(
    comparison = factor(.data$comparison, levels = names(comparison_palette)),
    abs_delta = abs(.data$signed_delta),
    variable_label = factor(.data$variable_label, levels = rev(environment_order))
  )

# Alternative 1: signed difference bars.
# Values to the right mean the focal cluster has more proportional dependency
# than Cluster 3 for that variable/sign; values to the left mean less.
signed_delta_limit <- max(abs(reference_delta_data$signed_delta), na.rm = TRUE)

reference_signed_difference_plot <- reference_delta_data |>
  ggplot(aes(
    x = .data$signed_delta,
    y = .data$variable_label,
    fill = .data$variable_label
  )) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.15) +
  facet_grid(.data$effect_direction ~ .data$comparison) +
  scale_fill_manual(values = environment_colours, guide = "none") +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(-signed_delta_limit, signed_delta_limit)
  ) +
  labs(
    title = "A. Signed difference from Cluster 3",
    subtitle = "Right = higher proportional dependency than Cluster 3; left = lower.",
    x = "Difference in proportional weighted dependency",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# Alternative 2: absolute distance lollipops.
# This removes the direction of the difference and makes the size of the
# mismatch easier to compare between Cluster 1 and Cluster 2.
reference_distance_plot <- reference_delta_data |>
  group_by(.data$effect_direction, .data$variable_label) |>
  mutate(sort_distance = max(.data$abs_delta, na.rm = TRUE)) |>
  ungroup() |>
  mutate(variable_label = fct_reorder(.data$variable_label, .data$sort_distance)) |>
  ggplot(aes(
    x = .data$abs_delta,
    y = .data$variable_label,
    colour = .data$comparison
  )) +
  geom_line(aes(group = interaction(.data$effect_direction, .data$variable_label)), colour = "grey80") +
  geom_point(size = 3.2) +
  facet_wrap(~effect_direction, ncol = 1) +
  scale_colour_manual(values = comparison_palette, name = "Comparison") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "B. Absolute distance from Cluster 3",
    subtitle = "Shorter distances indicate dependency profiles more similar to Cluster 3.",
    x = "Absolute proportional difference",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

# Alternative 3: profile-level similarity to Cluster 3.
# Bray-Curtis similarity is 1 - Bray-Curtis dissimilarity, so higher values mean
# a profile is compositionally closer to Cluster 3. Positive, negative, and
# combined dependency profiles are shown separately.
bray_similarity <- function(x, y) {
  denominator <- sum(x + y, na.rm = TRUE)
  if (!is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }
  1 - sum(abs(x - y), na.rm = TRUE) / denominator
}

similarity_profiles <- waffle_summary |>
  mutate(
    feature = paste(.data$effect_direction, .data$variable_label, sep = " | "),
    profile = "Combined"
  ) |>
  select(profile, bioregion, feature, dependency_share) |>
  bind_rows(
    waffle_summary |>
      transmute(
        profile = as.character(.data$effect_direction),
        bioregion,
        feature = as.character(.data$variable_label),
        dependency_share
      )
  )

similarity_to_cluster3 <- similarity_profiles |>
  pivot_wider(
    names_from = bioregion,
    values_from = dependency_share,
    values_fill = 0
  ) |>
  group_by(.data$profile) |>
  summarise(
    `Cluster 1 vs Cluster 3` = bray_similarity(.data$`Cluster 1`, .data$`Cluster 3`),
    `Cluster 2 vs Cluster 3` = bray_similarity(.data$`Cluster 2`, .data$`Cluster 3`),
    .groups = "drop"
  ) |>
  pivot_longer(
    cols = -profile,
    names_to = "comparison",
    values_to = "similarity"
  ) |>
  mutate(
    profile = factor(.data$profile, levels = c("Positive", "Negative", "Combined")),
    comparison = factor(.data$comparison, levels = c("Cluster 1 vs Cluster 3", "Cluster 2 vs Cluster 3"))
  )

similarity_palette <- c(
  "Cluster 1 vs Cluster 3" = "#76B7D8",
  "Cluster 2 vs Cluster 3" = "#C49A00"
)

cluster3_similarity_plot <- similarity_to_cluster3 |>
  ggplot(aes(
    x = .data$similarity,
    y = .data$profile,
    colour = .data$comparison
  )) +
  geom_line(aes(group = .data$profile), colour = "grey78", linewidth = 0.5) +
  geom_point(size = 4.3) +
  geom_text(
    aes(label = percent(.data$similarity, accuracy = 1)),
    nudge_x = 0.035,
    size = 3.2,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = similarity_palette, name = "Comparison") +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(
    title = "C. Similarity to Cluster 3",
    subtitle = "Bray-Curtis similarity of proportional dependency profiles.",
    x = "Similarity to Cluster 3",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

cluster3_comparison_concepts <- (
  reference_signed_difference_plot |
    (reference_distance_plot / cluster3_similarity_plot + plot_layout(heights = c(1.35, 1)))
) +
  plot_layout(widths = c(1.25, 1), guides = "collect") &
  theme(legend.position = "bottom")

# ---- Export ----
# Save both raster and vector versions. The PNG is convenient for quick review;
# the PDF is useful for manuscript assembly or manual editing.
png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-waffle-pca.png"))
pdf_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-waffle-pca.pdf"))
sankey_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-sankey-bump.png"))
sankey_pdf_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-sankey-bump.pdf"))
sankey_composite_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-sankey-bump-pca.png"))
sankey_composite_pdf_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-sankey-bump-pca.pdf"))
reference_comparison_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-cluster3-comparison-concepts.png"))
reference_comparison_pdf_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-cluster3-comparison-concepts.pdf"))
reference_signed_difference_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-cluster3-signed-differences.png"))
reference_distance_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-cluster3-distances.png"))
reference_similarity_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-cluster3-similarity.png"))
temporal_waffle_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-2010s-waffles.png"))
temporal_waffle_pdf_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-2010s-waffles.pdf"))
temporal_waffle_summary_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-2010s-waffle-summary.csv"))
temporal_waffle_tile_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-2010s-waffle-tiles.csv"))
corrected_waffle_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-corrected-waffles.png"))
corrected_waffle_pdf_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-corrected-waffles.pdf"))
variable_lane_waffle_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-variable-lane-waffles.png"))
variable_lane_waffle_pdf_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-variable-lane-waffles.pdf"))
mirrored_dependency_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-mirrored-direction-bars.png"))
mirrored_dependency_pdf_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-mirrored-direction-bars.pdf"))
pie_dependency_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-sign-split-pies.png"))
pie_dependency_pdf_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-sign-split-pies.pdf"))
nested_pie_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-nested-sign-variable-donuts.png"))
nested_pie_pdf_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-nested-sign-variable-donuts.pdf"))
semicircle_dependency_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-scaled-semicircles.png"))
semicircle_dependency_pdf_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-scaled-semicircles.pdf"))
semicircle_pca_map_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-semicircles-pca-map.png"))
semicircle_pca_map_pdf_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-semicircles-pca-map.pdf"))
dependency_pca_scores_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-pca-scores.csv"))
dependency_pca_centroids_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-pca-centroids.csv"))
diverging_dependency_png_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-diverging-dependency.png"))
diverging_dependency_pdf_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-diverging-dependency.pdf"))
corrected_dependency_summary_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-corrected-summary.csv"))
variable_lane_tiles_path <- file.path(out_dir, paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-variable-lane-waffle-tiles.csv"))

write_csv(temporal_waffle_summary, temporal_waffle_summary_path)
write_csv(temporal_waffle_tiles, temporal_waffle_tile_path)
write_csv(corrected_dependency_summary, corrected_dependency_summary_path)
write_csv(variable_lane_tiles, variable_lane_tiles_path)
write_csv(dependency_pca_scores, dependency_pca_scores_path)
write_csv(dependency_pca_centroids, dependency_pca_centroids_path)

ggsave(
  temporal_waffle_png_path,
  temporal_waffle_plot,
  width = 11.5,
  height = 9.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  temporal_waffle_pdf_path,
  temporal_waffle_plot,
  width = 11.5,
  height = 9.2,
  units = "in",
  bg = "white"
)

ggsave(
  corrected_waffle_png_path,
  corrected_waffle_plot,
  width = 10.2,
  height = 6.7,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  corrected_waffle_pdf_path,
  corrected_waffle_plot,
  width = 10.2,
  height = 6.7,
  units = "in",
  bg = "white"
)

ggsave(
  variable_lane_waffle_png_path,
  variable_lane_waffle_plot,
  width = 11.4,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  variable_lane_waffle_pdf_path,
  variable_lane_waffle_plot,
  width = 11.4,
  height = 9,
  units = "in",
  bg = "white"
)

ggsave(
  mirrored_dependency_png_path,
  mirrored_dependency_plot,
  width = 10.8,
  height = 4.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  mirrored_dependency_pdf_path,
  mirrored_dependency_plot,
  width = 10.8,
  height = 4.8,
  units = "in",
  bg = "white"
)

ggsave(
  pie_dependency_png_path,
  pie_dependency_plot,
  width = 9.8,
  height = 4.4,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  pie_dependency_pdf_path,
  pie_dependency_plot,
  width = 9.8,
  height = 4.4,
  units = "in",
  bg = "white"
)

ggsave(
  nested_pie_png_path,
  nested_pie_plot,
  width = 10.2,
  height = 4.6,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  nested_pie_pdf_path,
  nested_pie_plot,
  width = 10.2,
  height = 4.6,
  units = "in",
  bg = "white"
)

ggsave(
  semicircle_dependency_png_path,
  semicircle_dependency_plot,
  width = 10.2,
  height = 4.3,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  semicircle_dependency_pdf_path,
  semicircle_dependency_plot,
  width = 10.2,
  height = 4.3,
  units = "in",
  bg = "white"
)

ggsave(
  semicircle_pca_map_png_path,
  dependency_semicircle_pca_map_composite,
  width = 12.8,
  height = 8.4,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  semicircle_pca_map_pdf_path,
  dependency_semicircle_pca_map_composite,
  width = 12.8,
  height = 8.4,
  units = "in",
  bg = "white"
)

ggsave(
  diverging_dependency_png_path,
  diverging_dependency_plot,
  width = 8.4,
  height = 5.4,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  diverging_dependency_pdf_path,
  diverging_dependency_plot,
  width = 8.4,
  height = 5.4,
  units = "in",
  bg = "white"
)

ggsave(
  png_path,
  composite_plot,
  width = 13,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  pdf_path,
  composite_plot,
  width = 13,
  height = 12,
  units = "in",
  bg = "white"
)

suppressWarnings(
  ggsave(
    sankey_png_path,
    sankey_bump_plot,
    width = 10.5,
    height = 7,
    units = "in",
    dpi = 300,
    bg = "white"
  )
)

suppressWarnings(
  ggsave(
    sankey_pdf_path,
    sankey_bump_plot,
    width = 10.5,
    height = 7,
    units = "in",
    bg = "white"
  )
)

suppressWarnings(
  ggsave(
    sankey_composite_png_path,
    sankey_composite_plot,
    width = 13,
    height = 11.5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
)

suppressWarnings(
  ggsave(
    sankey_composite_pdf_path,
    sankey_composite_plot,
    width = 13,
    height = 11.5,
    units = "in",
    bg = "white"
  )
)

ggsave(
  reference_comparison_png_path,
  cluster3_comparison_concepts,
  width = 14,
  height = 9.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  reference_comparison_pdf_path,
  cluster3_comparison_concepts,
  width = 14,
  height = 9.5,
  units = "in",
  bg = "white"
)

ggsave(
  reference_signed_difference_png_path,
  reference_signed_difference_plot,
  width = 9.5,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  reference_distance_png_path,
  reference_distance_plot,
  width = 6.8,
  height = 7.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  reference_similarity_png_path,
  cluster3_similarity_plot,
  width = 6.8,
  height = 4.3,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Saved: ", png_path)
message("Saved: ", pdf_path)
message("Saved: ", sankey_png_path)
message("Saved: ", sankey_pdf_path)
message("Saved: ", sankey_composite_png_path)
message("Saved: ", sankey_composite_pdf_path)
message("Saved: ", reference_comparison_png_path)
message("Saved: ", reference_comparison_pdf_path)
message("Saved: ", reference_signed_difference_png_path)
message("Saved: ", reference_distance_png_path)
message("Saved: ", reference_similarity_png_path)
message("Saved: ", temporal_waffle_png_path)
message("Saved: ", temporal_waffle_pdf_path)
message("Saved: ", corrected_waffle_png_path)
message("Saved: ", corrected_waffle_pdf_path)
message("Saved: ", variable_lane_waffle_png_path)
message("Saved: ", variable_lane_waffle_pdf_path)
message("Saved: ", mirrored_dependency_png_path)
message("Saved: ", mirrored_dependency_pdf_path)
message("Saved: ", pie_dependency_png_path)
message("Saved: ", pie_dependency_pdf_path)
message("Saved: ", nested_pie_png_path)
message("Saved: ", nested_pie_pdf_path)
message("Saved: ", semicircle_dependency_png_path)
message("Saved: ", semicircle_dependency_pdf_path)
message("Saved: ", semicircle_pca_map_png_path)
message("Saved: ", semicircle_pca_map_pdf_path)
message("Saved: ", diverging_dependency_png_path)
message("Saved: ", diverging_dependency_pdf_path)
