rm(list = ls())

# Purpose:
# Build site-level community driver maps for the three atlas periods.
# For each fitted model grid cell, we combine species predicted probabilities,
# TjurR2-scaled variance partitioning, and certainty-weighted beta effects.
# The main mapped class is the environmental predictor with the strongest total
# absolute community contribution in that grid cell. The script also writes
# positive-only and negative-only versions to show which predictors dominate
# each side of the species response balance.

.libPaths(c("~/Rlibs", .libPaths()))

required_packages <- c(
  "tidyverse",
  "Hmsc",
  "readxl",
  "sf",
  "cowplot"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(Hmsc)
  library(readxl)
  library(sf)
  library(cowplot)
})

source(file.path("support_scripts", "figure_data_helpers.R"))

base_dir <- "HmscOutputs"
pattern <- "2026-03-13.*CoverageGoodAverage"
direction_threshold <- 0.25

out_dir <- file.path("notebooks", "exploratory", "outputs", "community-driver-maps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

period_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")

predictor_labels <- c(
  "tmean_breeding" = "Temperature",
  "prec_breeding" = "Precipitation",
  "hh" = "Land-use heterogeneity",
  "perc_urban" = "Urban (% coverage)",
  "perc_cropland" = "Cropland (% coverage)",
  "perc_pasture" = "Pasture (% coverage)",
  "perc_forest" = "Forest (% coverage)",
  "perc_grass_shrub" = "Grass/Shrubland (% coverage)"
)

driver_colours <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban (% coverage)" = "#4d4d4d",
  "Cropland (% coverage)" = "goldenrod1",
  "Pasture (% coverage)" = "darkorange",
  "Forest (% coverage)" = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2"
)

environment_order <- unname(predictor_labels)

grid_shape_path <- path.expand(
  "~/box/PhD/logistics/data/distributions/DK5km_ED50grid_approx_kvadrkod_DOF/DK5km_ED50grid_approx_kvadrkod_DOF.shp"
)

mainland_bbox <- st_bbox(
  c(xmin = 400000, xmax = 750000, ymin = 6000000, ymax = 6450000),
  crs = st_crs(25832)
)
bornholm_bbox <- st_bbox(
  c(xmin = 855000, xmax = 905000, ymin = 6100000, ymax = 6160000),
  crs = st_crs(25832)
)

mainland_width <- as.numeric(mainland_bbox["xmax"] - mainland_bbox["xmin"])
mainland_height <- as.numeric(mainland_bbox["ymax"] - mainland_bbox["ymin"])
bornholm_width <- as.numeric(bornholm_bbox["xmax"] - bornholm_bbox["xmin"])
bornholm_height <- as.numeric(bornholm_bbox["ymax"] - bornholm_bbox["ymin"])
inset_w <- bornholm_width / mainland_width
inset_h <- bornholm_height / mainland_height

read_scaled_vp_table <- function(model_folder, base_dir, variables) {
  path <- file.path(
    base_dir,
    model_folder,
    "Results",
    paste0(model_folder, "parameter_estimates_VP_.csv")
  )

  raw_vp <- read.csv(path, check.names = FALSE)
  names(raw_vp)[1] <- "Factor"
  raw_vp <- raw_vp |>
    tibble::column_to_rownames(var = "Factor")

  missing_rows <- setdiff(c(variables, "TjurR2"), rownames(raw_vp))
  if (length(missing_rows) > 0) {
    stop(
      "Missing VP row(s) in ",
      path,
      ": ",
      paste(missing_rows, collapse = ", "),
      call. = FALSE
    )
  }

  r2_values <- as.numeric(raw_vp["TjurR2", ])
  names(r2_values) <- colnames(raw_vp)

  raw_vp[variables, , drop = FALSE] |>
    mutate(across(everything(), as.numeric)) |>
    as.matrix() |>
    sweep(2, r2_values, `*`) |>
    as.data.frame(check.names = FALSE) |>
    tibble::rownames_to_column("variable")
}

read_beta_parameter_workbook <- function(model_folder, base_dir) {
  path <- file.path(
    base_dir,
    model_folder,
    "Results",
    paste0(model_folder, "parameter_estimates_Beta_.xlsx")
  )

  posterior_mean <- readxl::read_excel(path, sheet = "Posterior mean") |>
    pivot_longer(
      cols = -Species,
      names_to = "variable",
      values_to = "beta_mean"
    )

  posterior_positive <- readxl::read_excel(path, sheet = "Pr(x>0)") |>
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
      beta_mean = as.numeric(beta_mean),
      prob_positive = as.numeric(prob_positive),
      certainty_weight = 2 * abs(prob_positive - 0.5)
    )
}

predictor_sds_by_atlas <- function(models, variables) {
  imap_dfr(models, function(model, atlas) {
    missing_predictors <- setdiff(variables, colnames(model$XData))
    if (length(missing_predictors) > 0) {
      stop(
        "Missing predictor(s) in Atlas ",
        atlas,
        " XData: ",
        paste(missing_predictors, collapse = ", "),
        call. = FALSE
      )
    }

    tibble(
      atlas = as.character(atlas),
      variable = variables,
      predictor_sd = map_dbl(variables, ~ sd(model$XData[[.x]], na.rm = TRUE))
    )
  })
}

make_species_weights <- function(pred_y, design, atlas) {
  expected_richness <- rowSums(pred_y, na.rm = TRUE)
  if (any(!is.finite(expected_richness)) || any(expected_richness <= 0)) {
    stop(
      "Atlas ",
      atlas,
      " has site(s) with non-finite or zero expected richness; cannot normalize species weights.",
      call. = FALSE
    )
  }

  pred_y |>
    as.data.frame(check.names = FALSE) |>
    tibble::rownames_to_column("survey") |>
    pivot_longer(
      cols = -survey,
      names_to = "species",
      values_to = "predicted_probability"
    ) |>
    left_join(
      tibble(
        survey = rownames(pred_y),
        expected_richness = expected_richness
      ),
      by = "survey"
    ) |>
    mutate(
      species_weight = predicted_probability / expected_richness,
      atlas = as.character(atlas),
      period = period_lookup[atlas]
    ) |>
    left_join(
      design |> select(survey, site, X, Y),
      by = "survey"
    )
}

build_atlas_driver_scores <- function(atlas, model_folder, design, pred_y,
                                      predictor_sds, variables) {
  message("Scoring community drivers for Atlas ", atlas, " (", period_lookup[[atlas]], ")...")

  vp_scaled <- read_scaled_vp_table(model_folder, base_dir, variables)
  beta_long <- read_beta_parameter_workbook(model_folder, base_dir) |>
    filter(variable %in% variables) |>
    left_join(
      predictor_sds |> filter(atlas == .env$atlas),
      by = "variable"
    ) |>
    mutate(
      atlas = as.character(.env$atlas),
      period = period_lookup[atlas],
      standardized_beta = beta_mean * predictor_sd,
      weighted_standardized_beta = standardized_beta * certainty_weight
    )

  prediction_species <- colnames(pred_y)
  vp_species <- setdiff(colnames(vp_scaled), "variable")
  beta_species <- unique(beta_long$species)

  if (!setequal(prediction_species, vp_species) || !setequal(prediction_species, beta_species)) {
    stop(
      "Species names do not align for Atlas ",
      atlas,
      ". prediction species = ",
      length(prediction_species),
      "; VP species = ",
      length(vp_species),
      "; beta species = ",
      length(beta_species),
      call. = FALSE
    )
  }

  species_weights <- make_species_weights(pred_y, design, atlas)

  vp_long <- vp_scaled |>
    pivot_longer(
      cols = -variable,
      names_to = "species",
      values_to = "scaled_vp"
    ) |>
    mutate(
      atlas = as.character(.env$atlas),
      period = period_lookup[atlas],
      variable_label = factor(unname(predictor_labels[variable]), levels = environment_order)
    )

  contributions <- species_weights |>
    inner_join(
      vp_long,
      by = c("atlas", "period", "species"),
      relationship = "many-to-many"
    ) |>
    inner_join(
      beta_long |>
        select(
          atlas, period, species, variable, beta_mean, prob_positive,
          certainty_weight, predictor_sd, standardized_beta,
          weighted_standardized_beta
        ),
      by = c("atlas", "period", "species", "variable"),
      relationship = "many-to-one"
    ) |>
    mutate(
      contribution = species_weight * scaled_vp * weighted_standardized_beta
    )

  site_scores <- contributions |>
    group_by(atlas, period, survey, site, X, Y, variable, variable_label) |>
    summarise(
      n_species_contributing = n_distinct(species),
      expected_richness = first(expected_richness),
      positive_importance = sum(pmax(contribution, 0), na.rm = TRUE),
      negative_importance = sum(abs(pmin(contribution, 0)), na.rm = TRUE),
      total_importance = sum(abs(contribution), na.rm = TRUE),
      net_score = sum(contribution, na.rm = TRUE),
      direction_balance = if_else(
        total_importance > 0,
        net_score / total_importance,
        NA_real_
      ),
      .groups = "drop"
    ) |>
    group_by(atlas, period, survey, site) |>
    mutate(
      total_site_importance = sum(total_importance, na.rm = TRUE),
      relative_importance_at_site = if_else(
        total_site_importance > 0,
        total_importance / total_site_importance,
        NA_real_
      )
    ) |>
    ungroup()

  dominant_sites <- site_scores |>
    arrange(atlas, survey, desc(total_importance), variable) |>
    group_by(atlas, period, survey, site, X, Y) |>
    slice(1) |>
    ungroup() |>
    transmute(
      atlas,
      period,
      survey,
      site,
      X,
      Y,
      expected_richness,
      dominant_variable = variable,
      dominant_variable_label = as.character(variable_label),
      dominant_direction = case_when(
        direction_balance >= direction_threshold ~ "Mostly positive",
        direction_balance <= -direction_threshold ~ "Mostly negative",
        TRUE ~ "Mixed"
      ),
      dominant_total_importance = total_importance,
      dominant_net_score = net_score,
      dominant_direction_balance = direction_balance,
      total_site_importance,
      relative_dominance = relative_importance_at_site,
      n_species_contributing
    )

  list(site_scores = site_scores, dominant_sites = dominant_sites)
}

make_signed_dominant_sites <- function(site_scores, effect_subset) {
  score_column <- switch(
    effect_subset,
    positive = "positive_importance",
    negative = "negative_importance",
    stop("Unknown effect subset: ", effect_subset, call. = FALSE)
  )
  total_column <- paste0("total_site_", effect_subset, "_importance")
  relative_column <- paste0("relative_", effect_subset, "_dominance")
  dominant_score_column <- paste0("dominant_", effect_subset, "_importance")

  site_scores |>
    group_by(atlas, period, survey, site) |>
    mutate(
      "{total_column}" := sum(.data[[score_column]], na.rm = TRUE)
    ) |>
    ungroup() |>
    arrange(atlas, survey, desc(.data[[score_column]]), variable) |>
    group_by(atlas, period, survey, site, X, Y) |>
    slice(1) |>
    ungroup() |>
    transmute(
      effect_subset = effect_subset,
      atlas,
      period,
      survey,
      site,
      X,
      Y,
      expected_richness,
      dominant_variable = if_else(.data[[score_column]] > 0, variable, NA_character_),
      dominant_variable_label = if_else(
        .data[[score_column]] > 0,
        as.character(variable_label),
        NA_character_
      ),
      "{dominant_score_column}" := .data[[score_column]],
      "{total_column}" := .data[[total_column]],
      "{relative_column}" := if_else(
        .data[[total_column]] > 0,
        .data[[score_column]] / .data[[total_column]],
        NA_real_
      ),
      n_species_contributing
    )
}

make_signed_summary <- function(dominant_sites, effect_subset) {
  dominant_sites |>
    mutate(
      dominant_variable = coalesce(
        as.character(dominant_variable),
        paste0("no_", effect_subset, "_effect")
      ),
      dominant_variable_label = coalesce(
        as.character(dominant_variable_label),
        paste("No", effect_subset, "effect")
      )
    ) |>
    count(
      effect_subset,
      atlas,
      period,
      dominant_variable,
      dominant_variable_label,
      name = "n_sites"
    ) |>
    group_by(effect_subset, atlas, period) |>
    mutate(prop_sites = n_sites / sum(n_sites)) |>
    ungroup() |>
    arrange(effect_subset, atlas, desc(prop_sites), dominant_variable_label)
}

make_polygon_maps <- function(dominant_sites, shape_sf, legend_title = "Dominant driver") {
  shape_sf <- shape_sf |>
    mutate(kvadratkod = as.character(kvadratkod))

  plot_layer <- function(plot_data) {
    list(
      geom_sf(
        data = plot_data,
        aes(fill = dominant_variable_label),
        color = "grey75",
        linewidth = 0.08
      ),
      scale_fill_manual(
        values = driver_colours,
        breaks = environment_order,
        drop = FALSE,
        na.value = "#f2f2f2",
        name = legend_title
      )
    )
  }

  atlas_map <- function(atlas_id) {
    period_label <- period_lookup[[atlas_id]]
    plot_data <- shape_sf |>
      left_join(
        dominant_sites |> filter(atlas == atlas_id),
        by = c("kvadratkod" = "site")
      ) |>
      mutate(
        dominant_variable_label = factor(dominant_variable_label, levels = environment_order)
      )

    p_main <- ggplot() +
      plot_layer(plot_data) +
      coord_sf(
        xlim = c(mainland_bbox["xmin"], mainland_bbox["xmax"]),
        ylim = c(mainland_bbox["ymin"], mainland_bbox["ymax"]),
        expand = FALSE
      ) +
      labs(title = period_label, x = NULL, y = NULL) +
      theme_void(base_size = 11) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.background = element_rect(fill = "white", color = NA)
      )

    p_inset <- ggplot() +
      plot_layer(plot_data) +
      coord_sf(
        xlim = c(bornholm_bbox["xmin"], bornholm_bbox["xmax"]),
        ylim = c(bornholm_bbox["ymin"], bornholm_bbox["ymax"]),
        expand = FALSE
      ) +
      theme_void(base_size = 11) +
      theme(
        legend.position = "none",
        panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.35)
      )

    ggdraw(p_main) +
      draw_plot(
        p_inset,
        x = 1 - inset_w - 0.18,
        y = 1 - inset_h - 0.18,
        width = inset_w,
        height = inset_h
      )
  }

  map_plots <- map(names(period_lookup), atlas_map)

  legend <- cowplot::get_legend(
    tibble(
      x = seq_along(environment_order),
      y = 1,
      dominant_variable_label = factor(environment_order, levels = environment_order)
    ) |>
      ggplot(aes(x = x, y = y, fill = dominant_variable_label)) +
      geom_tile() +
      scale_fill_manual(
        values = driver_colours,
        breaks = environment_order,
        drop = FALSE,
        name = legend_title
      ) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
      theme_void(base_size = 10) +
      theme(
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 9)
      )
  )

  cowplot::plot_grid(
    cowplot::plot_grid(plotlist = map_plots, nrow = 1),
    legend,
    ncol = 1,
    rel_heights = c(1, 0.14)
  )
}

make_point_map <- function(dominant_sites, legend_title = "Dominant driver") {
  dominant_sites |>
    mutate(
      dominant_variable_label = factor(dominant_variable_label, levels = environment_order),
      period = factor(period, levels = unname(period_lookup))
    ) |>
    ggplot(aes(x = X, y = Y, colour = dominant_variable_label)) +
    geom_point(shape = 15, size = 1.1, alpha = 0.9) +
    coord_equal() +
    facet_wrap(~ period, nrow = 1) +
    scale_colour_manual(
      values = driver_colours,
      breaks = environment_order,
      drop = FALSE,
      name = legend_title
    ) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 13),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

message("Loading atlas model objects and fitted-site predictions...")
model_folders <- figure_model_folders(pattern = pattern, base_dir = base_dir)
atlas_ids <- as.character(atlas_numbers(model_folders))
names(model_folders) <- atlas_ids
model_folders <- model_folders[names(period_lookup)]

if (any(is.na(model_folders))) {
  stop("Expected Atlas 1, 2, and 3 model folders for pattern: ", pattern, call. = FALSE)
}

models <- load_hmsc_posteriors(model_folders, base_dir = base_dir)
designs <- load_hmsc_study_designs(models)
preds_y <- load_or_compute_site_predictions(models, model_folders, base_dir = base_dir)
names(models) <- names(designs) <- names(preds_y) <- names(period_lookup)

variables <- names(predictor_labels)
predictor_sds <- predictor_sds_by_atlas(models, variables)

atlas_outputs <- imap(model_folders, function(model_folder, atlas) {
  build_atlas_driver_scores(
    atlas = atlas,
    model_folder = model_folder,
    design = designs[[atlas]],
    pred_y = preds_y[[atlas]],
    predictor_sds = predictor_sds,
    variables = variables
  )
})

site_scores <- map_dfr(atlas_outputs, "site_scores") |>
  mutate(variable_label = as.character(variable_label))

dominant_sites <- map_dfr(atlas_outputs, "dominant_sites") |>
  mutate(
    period = factor(period, levels = unname(period_lookup)),
    dominant_variable_label = factor(dominant_variable_label, levels = environment_order),
    dominant_direction = factor(
      dominant_direction,
      levels = c("Mostly positive", "Mixed", "Mostly negative")
    )
  )

expected_dominant_rows <- sum(map_int(preds_y, nrow))
if (nrow(dominant_sites) != expected_dominant_rows) {
  stop(
    "Dominant-site output has ",
    nrow(dominant_sites),
    " rows, but expected ",
    expected_dominant_rows,
    " fitted sites across atlases.",
    call. = FALSE
  )
}

if (any(!is.finite(dominant_sites$relative_dominance))) {
  stop("Some sites have non-finite relative dominance values.", call. = FALSE)
}

atlas_summary <- dominant_sites |>
  count(
    atlas,
    period,
    dominant_variable,
    dominant_variable_label,
    dominant_direction,
    name = "n_sites"
  ) |>
  group_by(atlas, period) |>
  mutate(prop_sites = n_sites / sum(n_sites)) |>
  ungroup() |>
  arrange(atlas, desc(prop_sites), dominant_variable_label, dominant_direction)

prop_checks <- atlas_summary |>
  group_by(atlas) |>
  summarise(prop_sum = sum(prop_sites), .groups = "drop")

if (any(abs(prop_checks$prop_sum - 1) > 1e-8)) {
  stop("Atlas summary proportions do not sum to one within atlas.", call. = FALSE)
}

positive_dominant_sites <- make_signed_dominant_sites(site_scores, "positive") |>
  mutate(
    period = factor(period, levels = unname(period_lookup)),
    dominant_variable_label = factor(dominant_variable_label, levels = environment_order)
  )

negative_dominant_sites <- make_signed_dominant_sites(site_scores, "negative") |>
  mutate(
    period = factor(period, levels = unname(period_lookup)),
    dominant_variable_label = factor(dominant_variable_label, levels = environment_order)
  )

positive_summary <- make_signed_summary(positive_dominant_sites, "positive")
negative_summary <- make_signed_summary(negative_dominant_sites, "negative")

for (split_name in c("positive", "negative")) {
  split_sites <- if (split_name == "positive") positive_dominant_sites else negative_dominant_sites
  split_summary <- if (split_name == "positive") positive_summary else negative_summary

  if (nrow(split_sites) != expected_dominant_rows) {
    stop(
      "The ",
      split_name,
      "-only dominant-site output has ",
      nrow(split_sites),
      " rows, but expected ",
      expected_dominant_rows,
      ".",
      call. = FALSE
    )
  }

  split_prop_checks <- split_summary |>
    group_by(atlas) |>
    summarise(prop_sum = sum(prop_sites), .groups = "drop")

  if (any(abs(split_prop_checks$prop_sum - 1) > 1e-8)) {
    stop(
      "The ",
      split_name,
      "-only summary proportions do not sum to one within atlas.",
      call. = FALSE
    )
  }
}

write_csv(
  site_scores,
  file.path(out_dir, "community-driver-site-scores.csv"),
  na = ""
)
write_csv(
  dominant_sites,
  file.path(out_dir, "community-driver-dominant-sites.csv"),
  na = ""
)
write_csv(
  atlas_summary,
  file.path(out_dir, "community-driver-atlas-summary.csv"),
  na = ""
)
write_csv(
  positive_dominant_sites,
  file.path(out_dir, "community-driver-positive-dominant-sites.csv"),
  na = ""
)
write_csv(
  negative_dominant_sites,
  file.path(out_dir, "community-driver-negative-dominant-sites.csv"),
  na = ""
)
write_csv(
  positive_summary,
  file.path(out_dir, "community-driver-positive-atlas-summary.csv"),
  na = ""
)
write_csv(
  negative_summary,
  file.path(out_dir, "community-driver-negative-atlas-summary.csv"),
  na = ""
)

message("Building the dominant-driver map...")
if (file.exists(grid_shape_path)) {
  shape_sf <- sf::st_read(grid_shape_path, quiet = TRUE)
  if (is.na(sf::st_crs(shape_sf))) {
    sf::st_crs(shape_sf) <- 25832
  }
  driver_map <- make_polygon_maps(dominant_sites, shape_sf, legend_title = "Dominant driver")
  positive_driver_map <- make_polygon_maps(
    positive_dominant_sites,
    shape_sf,
    legend_title = "Positive driver"
  )
  negative_driver_map <- make_polygon_maps(
    negative_dominant_sites,
    shape_sf,
    legend_title = "Negative driver"
  )
} else {
  warning(
    "Grid shapefile not found at ",
    grid_shape_path,
    ". Falling back to a coordinate point map."
  )
  driver_map <- make_point_map(dominant_sites, legend_title = "Dominant driver")
  positive_driver_map <- make_point_map(
    positive_dominant_sites,
    legend_title = "Positive driver"
  )
  negative_driver_map <- make_point_map(
    negative_dominant_sites,
    legend_title = "Negative driver"
  )
}

plot_path <- file.path(out_dir, "community-driver-dominant-maps.png")
positive_plot_path <- file.path(out_dir, "community-driver-positive-dominant-maps.png")
negative_plot_path <- file.path(out_dir, "community-driver-negative-dominant-maps.png")
ggsave(
  filename = plot_path,
  plot = driver_map,
  width = 12,
  height = 5.6,
  dpi = 300,
  bg = "white"
)
ggsave(
  filename = positive_plot_path,
  plot = positive_driver_map,
  width = 12,
  height = 5.6,
  dpi = 300,
  bg = "white"
)
ggsave(
  filename = negative_plot_path,
  plot = negative_driver_map,
  width = 12,
  height = 5.6,
  dpi = 300,
  bg = "white"
)

non_png_plots <- list.files(
  out_dir,
  pattern = "\\.(pdf|svg|jpg|jpeg|tif|tiff)$",
  ignore.case = TRUE,
  full.names = TRUE
)
if (length(non_png_plots) > 0) {
  stop(
    "Found non-PNG plot output(s) in ",
    out_dir,
    ": ",
    paste(basename(non_png_plots), collapse = ", "),
    call. = FALSE
  )
}

message("Finished community driver maps.")
message("Outputs written to: ", out_dir)
message("Map written to: ", plot_path)
message("Positive-only map written to: ", positive_plot_path)
message("Negative-only map written to: ", negative_plot_path)
