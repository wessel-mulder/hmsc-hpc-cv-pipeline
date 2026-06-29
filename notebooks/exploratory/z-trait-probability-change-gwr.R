rm(list = ls())

.libPaths(c("~/Rlibs", .libPaths()))

library(tidyverse)
library(Hmsc)
library(GWmodel)
library(sp)
library(sf)
library(terra)
library(scales)
library(patchwork)

source(file.path("support_scripts", "figure_data_helpers.R"))

base_dir <- "HmscOutputs"
pattern <- "2026-03-13"
bandwidth_m <- 100000
min_species_per_group <- 5
atlas_from <- "1"
atlas_to <- "3"
period_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")

out_dir <- file.path("notebooks", "exploratory", "outputs", "trait-probability-change-gwr")
plot_dir <- file.path(out_dir, "plots")
model_dir <- file.path(out_dir, "models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

grid_shape_path <- path.expand("~/box/PhD/logistics/data/distributions/DK5km_ED50grid_approx_kvadrkod_DOF/DK5km_ED50grid_approx_kvadrkod_DOF.shp")
shape_sf <- if (file.exists(grid_shape_path)) st_as_sf(vect(grid_shape_path)) else NULL

predictor_labels <- c(
  "tmean_breeding" = "Temperature",
  "prec_breeding" = "Precipitation",
  "hh" = "Land-use heterogeneity",
  "perc_urban" = "Urban",
  "perc_cropland" = "Cropland",
  "perc_pasture" = "Pasture",
  "perc_forest" = "Forest",
  "perc_grass_shrub" = "Grass/Shrubland"
)

predictor_colours <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban" = "#4d4d4d",
  "Cropland" = "goldenrod1",
  "Pasture" = "darkorange",
  "Forest" = "springgreen4",
  "Grass/Shrubland" = "springgreen2",
  "No supported driver" = "#d9d9d9"
)

all_predictor_vars <- names(predictor_labels)
multivariate_predictor_vars <- c(
  "tmean_breeding",
  "prec_breeding",
  "perc_urban",
  "perc_cropland",
  "perc_forest"
)

scaled_name <- function(x) paste0("scaled_delta_", x)
delta_name <- function(x) paste0("delta_", x)

get_base_site <- function(x) sub("_[123]$", "", as.character(x))

message("Loading models, study designs, and site predictions...")
model_folders <- figure_model_folders(pattern = pattern, base_dir = base_dir)
atlas_ids <- as.character(atlas_numbers(model_folders))
mods <- load_hmsc_posteriors(model_folders, base_dir = base_dir)
designs <- load_hmsc_study_designs(mods)
preds_y <- load_or_compute_site_predictions(mods, model_folders, base_dir = base_dir)
names(mods) <- names(designs) <- names(preds_y) <- atlas_ids

focal_atlases <- c(atlas_from, atlas_to)
mods <- mods[focal_atlases]
designs <- designs[focal_atlases]
preds_y <- preds_y[focal_atlases]

load("Data/preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData")
traits <- Tr |>
  as.data.frame() |>
  rownames_to_column("species") |>
  select(species, Migration_a3_DOF, foraging_guild_consensus)

group_species_table <- function(traits, trait_col, trait_family, min_species = 1) {
  traits |>
    filter(!is.na(.data[[trait_col]]), .data[[trait_col]] != "") |>
    group_by(trait_family = trait_family, trait_value = .data[[trait_col]]) |>
    summarise(n_species = n_distinct(species), species = list(species), .groups = "drop") |>
    filter(n_species >= min_species)
}

trait_groups_all <- bind_rows(
  group_species_table(traits, "Migration_a3_DOF", "Migratory strategy", min_species = 1),
  group_species_table(traits, "foraging_guild_consensus", "Foraging guild", min_species = 1)
)
trait_groups_modelled <- trait_groups_all |> filter(n_species >= min_species_per_group)

write_csv(
  trait_groups_modelled |> select(trait_family, trait_value, n_species),
  file.path(out_dir, paste0(pattern, "-trait-probability-change-gwr-groups-n-ge-", min_species_per_group, ".csv"))
)

trait_probability_surface <- function(pred_y, design, trait_groups, traits) {
  model_species <- colnames(pred_y)
  denominator_species <- intersect(model_species, traits$species)
  expected_total_richness <- rowSums(pred_y[, denominator_species, drop = FALSE], na.rm = TRUE)

  pmap_dfr(trait_groups, function(trait_family, trait_value, n_species, species) {
    numerator_species <- intersect(model_species, species)
    expected_group_richness <- rowSums(pred_y[, numerator_species, drop = FALSE], na.rm = TRUE)

    tibble(
      survey = rownames(pred_y),
      trait_family = trait_family,
      trait_value = trait_value,
      n_species = n_species,
      trait_probability = ifelse(
        expected_total_richness > 0,
        expected_group_richness / expected_total_richness,
        NA_real_
      ),
      expected_group_richness = expected_group_richness,
      expected_total_richness = expected_total_richness
    )
  }) |>
    left_join(design, by = "survey")
}

message("Computing trait-probability deltas...")
trait_surfaces_all <- imap_dfr(preds_y, function(pred_y, atlas) {
  trait_probability_surface(pred_y, designs[[atlas]], trait_groups_all, traits) |>
    mutate(
      atlas = as.character(.env$atlas),
      period = period_lookup[[as.character(.env$atlas)]],
      base_site = get_base_site(site),
      .before = 1
    )
})

probability_sum_checks <- trait_surfaces_all |>
  group_by(atlas, period, trait_family, survey, site) |>
  summarise(probability_sum = sum(trait_probability, na.rm = TRUE), .groups = "drop") |>
  group_by(atlas, period, trait_family) |>
  summarise(
    n_sites = n(),
    min_sum = min(probability_sum, na.rm = TRUE),
    mean_sum = mean(probability_sum, na.rm = TRUE),
    max_sum = max(probability_sum, na.rm = TRUE),
    max_abs_deviation_from_one = max(abs(probability_sum - 1), na.rm = TRUE),
    .groups = "drop"
  )

trait_deltas <- trait_surfaces_all |>
  filter(atlas %in% focal_atlases) |>
  semi_join(trait_groups_modelled |> select(trait_family, trait_value), by = c("trait_family", "trait_value")) |>
  select(
    atlas, period, trait_family, trait_value, n_species, base_site, site, survey,
    X, Y, trait_probability, expected_group_richness, expected_total_richness
  ) |>
  pivot_wider(
    id_cols = c(trait_family, trait_value, n_species, base_site),
    names_from = atlas,
    values_from = c(site, survey, period, X, Y, trait_probability, expected_group_richness, expected_total_richness),
    names_sep = "_"
  ) |>
  filter(!is.na(.data[[paste0("trait_probability_", atlas_from)]]),
         !is.na(.data[[paste0("trait_probability_", atlas_to)]])) |>
  transmute(
    trait_family,
    trait_value,
    n_species,
    base_site,
    site_from = .data[[paste0("site_", atlas_from)]],
    site_to = .data[[paste0("site_", atlas_to)]],
    survey_from = .data[[paste0("survey_", atlas_from)]],
    survey_to = .data[[paste0("survey_", atlas_to)]],
    X = .data[[paste0("X_", atlas_from)]],
    Y = .data[[paste0("Y_", atlas_from)]],
    period_from = .data[[paste0("period_", atlas_from)]],
    period_to = .data[[paste0("period_", atlas_to)]],
    probability_from = .data[[paste0("trait_probability_", atlas_from)]],
    probability_to = .data[[paste0("trait_probability_", atlas_to)]],
    probability_delta = probability_to - probability_from,
    expected_group_richness_from = .data[[paste0("expected_group_richness_", atlas_from)]],
    expected_group_richness_to = .data[[paste0("expected_group_richness_", atlas_to)]],
    expected_group_richness_delta = expected_group_richness_to - expected_group_richness_from,
    expected_total_richness_from = .data[[paste0("expected_total_richness_", atlas_from)]],
    expected_total_richness_to = .data[[paste0("expected_total_richness_", atlas_to)]]
  )

message("Computing environmental deltas...")
environment_surfaces <- imap_dfr(mods, function(model, atlas) {
  model$XData |>
    as.data.frame() |>
    rownames_to_column("survey") |>
    left_join(designs[[atlas]] |> select(survey, site, X, Y), by = "survey") |>
    mutate(
      atlas = as.character(.env$atlas),
      period = period_lookup[[as.character(.env$atlas)]],
      base_site = get_base_site(site),
      .before = 1
    ) |>
    select(atlas, period, base_site, site, survey, X, Y, all_of(all_predictor_vars))
})

environment_wide <- environment_surfaces |>
  filter(atlas %in% focal_atlases) |>
  pivot_wider(
    id_cols = base_site,
    names_from = atlas,
    values_from = c(period, site, survey, X, Y, all_of(all_predictor_vars)),
    names_sep = "_"
  )

environment_deltas <- environment_wide |>
  transmute(
    base_site,
    site_from = .data[[paste0("site_", atlas_from)]],
    site_to = .data[[paste0("site_", atlas_to)]],
    survey_from = .data[[paste0("survey_", atlas_from)]],
    survey_to = .data[[paste0("survey_", atlas_to)]],
    X = .data[[paste0("X_", atlas_from)]],
    Y = .data[[paste0("Y_", atlas_from)]]
  )

for (var in all_predictor_vars) {
  environment_deltas[[delta_name(var)]] <-
    environment_wide[[paste0(var, "_", atlas_to)]] -
    environment_wide[[paste0(var, "_", atlas_from)]]
}

delta_vars <- delta_name(all_predictor_vars)
env_delta_summary <- map_dfr(all_predictor_vars, function(var) {
  delta_col <- delta_name(var)
  values <- environment_deltas[[delta_col]]
  tibble(
    variable = var,
    predictor_label = predictor_labels[[var]],
    n_sites = sum(!is.na(values)),
    min_delta = min(values, na.rm = TRUE),
    q10_delta = unname(quantile(values, 0.10, na.rm = TRUE)),
    median_delta = median(values, na.rm = TRUE),
    mean_delta = mean(values, na.rm = TRUE),
    q90_delta = unname(quantile(values, 0.90, na.rm = TRUE)),
    max_delta = max(values, na.rm = TRUE),
    sd_delta = sd(values, na.rm = TRUE),
    prop_positive = mean(values > 0, na.rm = TRUE),
    prop_zero = mean(values == 0, na.rm = TRUE),
    prop_negative = mean(values < 0, na.rm = TRUE)
  )
})

zero_variance <- env_delta_summary |> filter(!is.finite(sd_delta) | sd_delta <= 0)
if (nrow(zero_variance) > 0) {
  warning("Dropping zero-variance environmental deltas: ", paste(zero_variance$predictor_label, collapse = ", "))
}
valid_predictor_vars <- env_delta_summary |>
  filter(is.finite(sd_delta), sd_delta > 0) |>
  pull(variable)
valid_delta_vars <- delta_name(valid_predictor_vars)
valid_scaled_vars <- scaled_name(valid_predictor_vars)

environment_deltas <- environment_deltas |>
  mutate(across(all_of(valid_delta_vars), ~ as.numeric(scale(.x)), .names = "scaled_{.col}"))

analysis_input <- trait_deltas |>
  inner_join(
    environment_deltas |> select(base_site, all_of(valid_delta_vars), all_of(valid_scaled_vars)),
    by = "base_site"
  )

write_csv(probability_sum_checks, file.path(out_dir, paste0(pattern, "-trait-probability-change-gwr-probability-sum-checks.csv")))
write_csv(env_delta_summary, file.path(out_dir, paste0(pattern, "-environment-delta-summary-2010s-minus-1970s.csv")))
write_csv(analysis_input, file.path(out_dir, paste0(pattern, "-trait-probability-environment-delta-input-2010s-minus-1970s.csv")))

fit_gwr <- function(df, vars, model_set, predictor_label = NA_character_) {
  working <- df |>
    select(base_site, site_from, site_to, X, Y, probability_delta, all_of(vars)) |>
    drop_na()

  if (nrow(working) < 40) {
    return(NULL)
  }

  metadata <- working |> select(base_site, site_from, site_to)
  coordinates(working) <- c("X", "Y")
  formula <- reformulate(vars, response = "probability_delta")

  model <- gwr.basic(
    formula,
    data = working,
    bw = bandwidth_m,
    kernel = "bisquare",
    adaptive = FALSE
  )

  list(
    mod = model,
    bandwidth = bandwidth_m,
    n_rows = nrow(working),
    vars = vars,
    model_set = model_set,
    predictor_label = predictor_label,
    metadata = metadata
  )
}

safe_fit_gwr <- function(df, vars, model_set, predictor_label = NA_character_) {
  tryCatch(
    fit_gwr(df, vars, model_set, predictor_label),
    error = function(e) {
      warning(
        "GWR failed for ",
        paste(unique(df$trait_family), unique(df$trait_value), model_set, predictor_label, sep = " / "),
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
}

summarise_gwr <- function(fit, df, model_set, predictor_label = NA_character_) {
  if (is.null(fit)) {
    return(tibble(
      trait_family = unique(df$trait_family),
      trait_value = unique(df$trait_value),
      n_species = unique(df$n_species),
      model_set = model_set,
      predictor = predictor_label,
      n_rows = nrow(df),
      bandwidth_m = bandwidth_m,
      aicc = NA_real_,
      rss = NA_real_,
      local_r2_mean = NA_real_,
      local_r2_median = NA_real_,
      local_r2_max = NA_real_,
      max_abs_local_t = NA_real_,
      status = "failed"
    ))
  }

  sdf <- as.data.frame(fit$mod$SDF)
  tv_cols <- paste0(fit$vars, "_TV")
  tv_cols <- tv_cols[tv_cols %in% names(sdf)]
  diagnostics <- fit$mod$GW.diagnostic

  tibble(
    trait_family = unique(df$trait_family),
    trait_value = unique(df$trait_value),
    n_species = unique(df$n_species),
    model_set = model_set,
    predictor = predictor_label,
    n_rows = fit$n_rows,
    bandwidth_m = fit$bandwidth,
    aicc = diagnostics$AICc,
    rss = diagnostics$RSS.gw,
    local_r2_mean = mean(sdf$Local_R2, na.rm = TRUE),
    local_r2_median = median(sdf$Local_R2, na.rm = TRUE),
    local_r2_max = max(sdf$Local_R2, na.rm = TRUE),
    max_abs_local_t = if (length(tv_cols) == 0) NA_real_ else max(abs(as.matrix(sdf[, tv_cols, drop = FALSE])), na.rm = TRUE),
    status = "ok"
  )
}

extract_multivariate_local <- function(fit, df) {
  if (is.null(fit)) return(tibble())

  sdf <- as.data.frame(fit$mod$SDF)
  vars <- fit$vars
  tv_cols <- paste0(vars, "_TV")
  coef_cols <- vars
  tv_matrix <- as.matrix(sdf[, tv_cols, drop = FALSE])
  dominant_idx <- max.col(abs(tv_matrix), ties.method = "first")
  dominant_tv <- tv_matrix[cbind(seq_len(nrow(tv_matrix)), dominant_idx)]
  dominant_var <- vars[dominant_idx]

  local_long <- map_dfr(seq_along(vars), function(i) {
    var <- vars[[i]]
    raw_var <- str_remove(var, "^scaled_delta_")
    tibble(
      predictor = raw_var,
      predictor_label = predictor_labels[[raw_var]],
      coefficient = sdf[[coef_cols[[i]]]],
      local_t = sdf[[tv_cols[[i]]]],
      abs_local_t = abs(local_t),
      direction = ifelse(local_t >= 0, "Positive", "Negative")
    )
  }) |>
    mutate(local_row = rep(seq_len(nrow(sdf)), times = length(vars)))

  metadata <- fit$metadata |>
    mutate(local_row = row_number())

  local_long |>
    left_join(metadata, by = "local_row") |>
    mutate(
      X = coordinates(fit$mod$SDF)[local_row, 1],
      Y = coordinates(fit$mod$SDF)[local_row, 2],
      local_r2 = sdf$Local_R2[local_row],
      trait_family = unique(df$trait_family),
      trait_value = unique(df$trait_value),
      n_species = unique(df$n_species),
      model_set = "multivariate_reduced",
      dominant_predictor = str_remove(dominant_var[local_row], "^scaled_delta_"),
      dominant_predictor_label = predictor_labels[dominant_predictor],
      dominant_local_t = dominant_tv[local_row],
      dominant_abs_local_t = abs(dominant_local_t),
      supported = dominant_abs_local_t >= 1.96,
      dominant_driver = ifelse(supported, dominant_predictor_label, "No supported driver"),
      dominant_direction = ifelse(supported, ifelse(dominant_local_t >= 0, "Positive", "Negative"), "Not supported")
    ) |>
    select(-local_row)
}

extract_single_local <- function(fit, df) {
  if (is.null(fit)) return(tibble())

  var <- fit$vars[[1]]
  raw_var <- str_remove(var, "^scaled_delta_")
  sdf <- as.data.frame(fit$mod$SDF)
  metadata <- fit$metadata
  tv_col <- paste0(var, "_TV")

  tibble(
    base_site = metadata$base_site,
    site_from = metadata$site_from,
    site_to = metadata$site_to,
    X = coordinates(fit$mod$SDF)[, 1],
    Y = coordinates(fit$mod$SDF)[, 2],
    predictor = raw_var,
    predictor_label = predictor_labels[[raw_var]],
    coefficient = sdf[[var]],
    local_t = sdf[[tv_col]],
    abs_local_t = abs(local_t),
    direction = ifelse(local_t >= 0, "Positive", "Negative"),
    local_r2 = sdf$Local_R2,
    trait_family = unique(df$trait_family),
    trait_value = unique(df$trait_value),
    n_species = unique(df$n_species),
    model_set = "single_predictor"
  )
}

message("Fitting reduced multivariate GWRs...")
multivariate_vars <- scaled_name(multivariate_predictor_vars)
multivariate_vars <- multivariate_vars[multivariate_vars %in% names(analysis_input)]
trait_splits <- analysis_input |>
  group_by(trait_family, trait_value) |>
  group_split()
names(trait_splits) <- map_chr(trait_splits, ~ paste(unique(.x$trait_family), unique(.x$trait_value), sep = "__"))

multivariate_models <- map(trait_splits, ~ safe_fit_gwr(.x, multivariate_vars, "multivariate_reduced"))
multivariate_summary <- map2_dfr(multivariate_models, trait_splits, summarise_gwr, model_set = "multivariate_reduced")
multivariate_local <- map2_dfr(multivariate_models, trait_splits, extract_multivariate_local)
multivariate_dominant <- multivariate_local |>
  distinct(
    trait_family, trait_value, n_species, base_site, site_from, site_to, X, Y,
    local_r2, dominant_predictor, dominant_predictor_label, dominant_local_t,
    dominant_abs_local_t, supported, dominant_driver, dominant_direction
  )

message("Fitting single-predictor GWR sensitivity layer...")
single_specs <- crossing(
  split_name = names(trait_splits),
  predictor = valid_predictor_vars
) |>
  mutate(
    scaled_predictor = scaled_name(predictor),
    predictor_label = predictor_labels[predictor]
  ) |>
  filter(scaled_predictor %in% names(analysis_input))

single_models <- pmap(single_specs, function(split_name, predictor, scaled_predictor, predictor_label) {
  safe_fit_gwr(trait_splits[[split_name]], scaled_predictor, "single_predictor", predictor_label)
})

single_summary <- pmap_dfr(
  list(single_models, single_specs$split_name, single_specs$predictor_label),
  function(fit, split_name, predictor_label) {
    summarise_gwr(fit, trait_splits[[split_name]], model_set = "single_predictor", predictor_label = predictor_label)
  }
)

single_local <- pmap_dfr(
  list(single_models, single_specs$split_name),
  function(fit, split_name) extract_single_local(fit, trait_splits[[split_name]])
)

single_dominant <- single_local |>
  group_by(trait_family, trait_value, n_species, base_site) |>
  slice_max(abs_local_t, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    supported = abs_local_t >= 1.96,
    dominant_driver = ifelse(supported, predictor_label, "No supported driver"),
    dominant_direction = ifelse(supported, direction, "Not supported")
  )

write_csv(multivariate_summary, file.path(out_dir, paste0(pattern, "-multivariate-reduced-gwr-summary-100km.csv")))
write_csv(multivariate_local, file.path(out_dir, paste0(pattern, "-multivariate-reduced-gwr-local-results-100km.csv")))
write_csv(multivariate_dominant, file.path(out_dir, paste0(pattern, "-multivariate-reduced-gwr-dominant-local-drivers-100km.csv")))
write_csv(single_summary, file.path(out_dir, paste0(pattern, "-single-predictor-gwr-summary-100km.csv")))
write_csv(single_local, file.path(out_dir, paste0(pattern, "-single-predictor-gwr-local-results-100km.csv")))
write_csv(single_dominant, file.path(out_dir, paste0(pattern, "-single-predictor-gwr-dominant-local-drivers-100km.csv")))
saveRDS(multivariate_models, file.path(model_dir, paste0(pattern, "-multivariate-reduced-gwr-models-100km.rds")))
saveRDS(single_models, file.path(model_dir, paste0(pattern, "-single-predictor-gwr-models-100km.rds")))

dominant_bar_summary <- bind_rows(
  multivariate_dominant |> mutate(model_set = "Multivariate reduced"),
  single_dominant |>
    transmute(
      trait_family, trait_value, n_species, base_site, site_from, site_to, X, Y, local_r2,
      dominant_predictor = predictor,
      dominant_predictor_label = predictor_label,
      dominant_local_t = local_t,
      dominant_abs_local_t = abs_local_t,
      supported,
      dominant_driver,
      dominant_direction,
      model_set = "Single predictor"
    )
) |>
  filter(supported) |>
  count(model_set, trait_family, trait_value, n_species, dominant_driver, dominant_direction, name = "n_cells") |>
  group_by(model_set, trait_family, trait_value, n_species) |>
  mutate(
    prop_cells = n_cells / sum(n_cells),
    signed_prop = ifelse(dominant_direction == "Positive", prop_cells, -prop_cells)
  ) |>
  ungroup()

write_csv(dominant_bar_summary, file.path(out_dir, paste0(pattern, "-dominant-driver-signed-bar-summary-100km.csv")))

plot_dominant_bars <- function(df, model_set_name) {
  plot_df <- df |>
    filter(model_set == model_set_name) |>
    mutate(
      trait_label = paste0(trait_value, " (n=", n_species, ")"),
      trait_label = fct_reorder(trait_label, signed_prop, .fun = sum)
    )

  if (nrow(plot_df) == 0) return(NULL)

  ggplot(plot_df, aes(signed_prop, trait_label, fill = dominant_driver)) +
    geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.35) +
    geom_col(width = 0.72) +
    facet_grid(trait_family ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(labels = function(x) paste0(abs(round(x * 100)), "%")) +
    scale_fill_manual(values = predictor_colours, name = "Driver", drop = FALSE) +
    labs(
      title = paste0(model_set_name, ": dominant local environmental-change drivers"),
      subtitle = "Positive effects are to the right; negative effects are to the left. Only |local t| >= 1.96 cells are counted.",
      x = "Share of supported cells",
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey40", size = 9)
    )
}

p_multi_bars <- plot_dominant_bars(dominant_bar_summary, "Multivariate reduced")
p_single_bars <- plot_dominant_bars(dominant_bar_summary, "Single predictor")

if (!is.null(p_multi_bars)) {
  ggsave(
    file.path(plot_dir, paste0(pattern, "-multivariate-reduced-dominant-driver-bars-100km.png")),
    p_multi_bars,
    width = 10,
    height = 9,
    dpi = 300
  )
}

if (!is.null(p_single_bars)) {
  ggsave(
    file.path(plot_dir, paste0(pattern, "-single-predictor-dominant-driver-bars-100km.png")),
    p_single_bars,
    width = 10,
    height = 9,
    dpi = 300
  )
}

trait_change_summary <- analysis_input |>
  distinct(trait_family, trait_value, n_species, base_site, probability_delta) |>
  group_by(trait_family, trait_value, n_species) |>
  summarise(mean_probability_delta = mean(probability_delta, na.rm = TRUE), .groups = "drop")

selected_map_groups <- bind_rows(
  trait_change_summary |>
    filter(trait_family == "Foraging guild") |>
    slice_max(abs(mean_probability_delta), n = 8, with_ties = FALSE),
  trait_change_summary |>
    filter(trait_family == "Migratory strategy") |>
    slice_max(abs(mean_probability_delta), n = 5, with_ties = FALSE)
)

join_shape <- function(df) {
  if (is.null(shape_sf)) return(NULL)
  shape_sf |> left_join(df, by = c("kvadratkod" = "base_site"))
}

map_theme <- theme_minimal(base_size = 9) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
    plot.background = element_rect(fill = "white", colour = NA)
  )

plot_dominant_maps <- function(dominant_df, selected_groups, model_label, file_stub) {
  plot_df <- dominant_df |>
    semi_join(selected_groups |> select(trait_family, trait_value), by = c("trait_family", "trait_value")) |>
    mutate(
      dominant_driver = factor(dominant_driver, levels = names(predictor_colours)),
      dominant_direction = factor(dominant_direction, levels = c("Positive", "Negative", "Not supported")),
      trait_label = paste0(trait_value, " (", trait_family, ")")
    )

  if (nrow(plot_df) == 0) return(NULL)

  if (!is.null(shape_sf)) {
    map_df <- join_shape(plot_df)
    direction_points <- map_df |>
      filter(dominant_direction %in% c("Positive", "Negative")) |>
      st_point_on_surface()

    p <- ggplot(map_df) +
      geom_sf(aes(fill = dominant_driver), colour = "grey35", linewidth = 0.04) +
      geom_sf(data = direction_points, aes(shape = dominant_direction), size = 0.35, colour = "black", stroke = 0.2) +
      facet_wrap(~ trait_label, ncol = 4) +
      scale_fill_manual(values = predictor_colours, drop = FALSE, name = "Dominant driver") +
      scale_shape_manual(values = c("Positive" = 16, "Negative" = 4), name = "Direction") +
      labs(title = paste0(model_label, ": dominant local environmental-change drivers")) +
      map_theme
  } else {
    p <- ggplot(plot_df, aes(X, Y, colour = dominant_driver, shape = dominant_direction)) +
      geom_point(size = 0.7, stroke = 0.25) +
      coord_equal() +
      facet_wrap(~ trait_label, ncol = 4) +
      scale_colour_manual(values = predictor_colours, drop = FALSE, name = "Dominant driver") +
      scale_shape_manual(values = c("Positive" = 16, "Negative" = 4, "Not supported" = 3), name = "Direction") +
      labs(title = paste0(model_label, ": dominant local environmental-change drivers")) +
      map_theme
  }

  ggsave(
    file.path(plot_dir, paste0(pattern, "-", file_stub, "-selected-dominant-driver-maps-100km.png")),
    p,
    width = 12,
    height = 9,
    dpi = 300
  )
  p
}

plot_dominant_maps(multivariate_dominant, selected_map_groups, "Multivariate reduced", "multivariate-reduced")
plot_dominant_maps(single_dominant, selected_map_groups, "Single predictor", "single-predictor")

message("Finished trait-probability change GWR analysis.")
print(probability_sum_checks)
print(env_delta_summary)
print(multivariate_summary |> count(status))
print(single_summary |> count(status))
