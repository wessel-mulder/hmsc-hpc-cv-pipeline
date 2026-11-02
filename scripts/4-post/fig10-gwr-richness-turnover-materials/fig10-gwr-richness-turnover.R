# Purpose:
# Build the publication-facing geographically weighted regression figure for
# modelled richness and all-site beta turnover. The script owns the GWR fitting
# workflow, stores reusable model/result objects locally, and reuses those cached
# outputs on later runs so the expensive regressions do not run every time.

rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
force_refit <- "--force-refit" %in% args

base_dir <- "./HmscOutputs"
pattern <- "2026-03-13"
bandwidth_label <- "fixed_100km"
bandwidth_m <- 100000
turnover_definition <- "mean_bray_curtis_to_all_sites"
significance_t <- 1.96

gwr_out_dir <- file.path("intermediate-outputs", "gwr-output", "fig10-gwr-richness-turnover", "gwr")
gwr_model_dir <- file.path(gwr_out_dir, "models")
figure_out_dir <- file.path("manuscript-outputs", "figures", "main")

gwr_fit_path <- file.path(gwr_model_dir, paste0(pattern, "-gwr-fits-richness-beta-100km.rds"))
gwr_inputs_path <- file.path(gwr_model_dir, paste0(pattern, "-gwr-analysis-inputs-100km.rds"))
gwr_summary_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-100km-summary.csv"))
gwr_local_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-local-coefficients-100km.csv"))
gwr_dominant_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-dominant-drivers-100km.csv"))
png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10-gwr-richness-turnover.png"))

dir.create(gwr_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gwr_model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_out_dir, recursive = TRUE, showWarnings = FALSE)

cache_uses_current_turnover_definition <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }

  inputs <- tryCatch(readRDS(path), error = function(e) NULL)
  identical(inputs$turnover_definition, turnover_definition)
}

needs_refit <- force_refit ||
  !file.exists(gwr_fit_path) ||
  !file.exists(gwr_inputs_path) ||
  !cache_uses_current_turnover_definition(gwr_inputs_path)

load_required_packages <- function(needs_refit) {
  core_packages <- c("tidyverse", "cowplot", "scales", "sp", "sf", "terra")
  refit_packages <- c("GWmodel", "Hmsc", "vegan")
  required_packages <- unique(c(core_packages, if (needs_refit) refit_packages else character()))

  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    stop(
      "Install these R packages before running this script: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  suppressPackageStartupMessages({
    library(tidyverse)
    library(cowplot)
    library(scales)
    library(sp)
    library(sf)
    library(terra)
    if (needs_refit) {
      library(GWmodel)
      library(Hmsc)
      library(vegan)
    }
  })
}

load_required_packages(needs_refit)

period_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")
period_levels <- unname(period_lookup)

response_labels <- c(
  "predicted_richness" = "Predicted richness",
  "beta_turnover" = "All-site beta turnover"
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

fig1_variable_colours <- c(
  "Site-level random effect" = "ivory2",
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban (% coverage)" = "snow3",
  "Cropland (% coverage)" = "goldenrod1",
  "Pasture (% coverage)" = "darkorange",
  "Forest (% coverage)" = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2"
)

figure1_driver_colours <- fig1_variable_colours[
  names(fig1_variable_colours) != "Site-level random effect"
]

driver_colours <- c(
  figure1_driver_colours,
  "Multiple associations" = "#7f3b08",
  "No significant driver" = "white"
)

no_significant_driver <- "No significant driver"
multiple_associations_driver <- "Multiple associations"
no_significant_outline <- "black"
dominance_ratio_threshold <- 2
driver_levels <- names(driver_colours)
direction_levels <- c("Positive", "Negative", "Not significant")
direction_shapes <- c("Positive" = 19, "Negative" = 4, "Not significant" = 19)
bandwidths_m <- c("fixed_100km" = bandwidth_m)

predictor_names <- function(df) {
  predictors <- c(
    "tmean_breeding",
    "prec_breeding",
    "hh",
    grep("^perc_", colnames(df), value = TRUE)
  )

  missing_predictors <- setdiff(predictors, colnames(df))
  if (length(missing_predictors) > 0) {
    stop("Missing GWR predictors: ", paste(missing_predictors, collapse = ", "), call. = FALSE)
  }

  predictors
}

prepare_gwr_data <- function(df, response) {
  vars <- predictor_names(df)
  required_cols <- c("survey", "site", "X", "Y", response, vars)
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop("Missing GWR input columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  out <- df |>
    select(all_of(required_cols)) |>
    drop_na()

  # GWR coefficients are easier to compare across predictors when each
  # environmental predictor is standardized before fitting.
  out[vars] <- as.data.frame(scale(out[vars]))
  out <- as.data.frame(out)
  rownames(out) <- out$survey
  coordinates(out) <- c("X", "Y")
  out
}

site_beta_turnover <- function(pred_y, design) {
  design <- design[match(rownames(pred_y), design$survey), ]

  community_dissimilarity <- as.matrix(vegan::vegdist(pred_y, method = "bray"))
  diag(community_dissimilarity) <- NA_real_

  tibble(
    survey = rownames(pred_y),
    beta_turnover = rowMeans(community_dissimilarity, na.rm = TRUE)
  ) |>
    left_join(design, by = "survey")
}

fit_gwr_model <- function(df, response, bandwidth_m = 100000) {
  vars <- predictor_names(df)
  gwr_data <- prepare_gwr_data(df, response)
  formula <- reformulate(vars, response = response)

  gwr_model <- GWmodel::gwr.basic(
    formula,
    data = gwr_data,
    bw = bandwidth_m,
    kernel = "bisquare",
    adaptive = FALSE
  )

  collin_diag <- GWmodel::gwr.collin.diagno(
    formula,
    data = gwr_data,
    bw = bandwidth_m,
    kernel = "bisquare",
    adaptive = FALSE
  )

  list(
    mod = gwr_model,
    collin = collin_diag,
    bandwidth = bandwidth_m,
    adaptive = FALSE
  )
}

fit_response_bandwidths <- function(dfs, response) {
  imap(dfs, function(df, atlas) {
    setNames(
      list(fit_gwr_model(df, response = response, bandwidth_m = bandwidth_m)),
      bandwidth_label
    )
  })
}

summarise_gwr_model <- function(fit, response, atlas, label) {
  sdf <- as.data.frame(fit$mod$SDF)
  diagnostics <- fit$mod$GW.diagnostic
  local_cn <- fit$collin$SDF$local_CN

  tibble(
    response = response,
    atlas = as.integer(atlas),
    period = period_lookup[[as.character(atlas)]],
    bandwidth = label,
    adaptive = fit$adaptive,
    bandwidth_value = fit$bandwidth,
    aicc = diagnostics$AICc,
    rss = diagnostics$RSS.gw,
    local_r2_mean = mean(sdf$Local_R2, na.rm = TRUE),
    local_r2_median = median(sdf$Local_R2, na.rm = TRUE),
    local_r2_min = min(sdf$Local_R2, na.rm = TRUE),
    local_r2_max = max(sdf$Local_R2, na.rm = TRUE),
    local_cn_mean = mean(local_cn, na.rm = TRUE),
    local_cn_max = max(local_cn, na.rm = TRUE),
    local_cn_over_30 = mean(local_cn > 30, na.rm = TRUE)
  )
}

extract_gwr_local_coefficients <- function(fit, response, atlas, label) {
  sdf <- as.data.frame(fit$mod$SDF)
  surveys <- rownames(sdf)
  vars <- names(var_labels)
  vars <- vars[vars %in% colnames(sdf)]

  map_dfr(vars, function(var) {
    t_col <- paste0(var, "_TV")
    se_col <- paste0(var, "_SE")
    local_t <- if (t_col %in% colnames(sdf)) sdf[[t_col]] else rep(NA_real_, nrow(sdf))

    tibble(
      response = response,
      response_label = response_labels[[response]],
      atlas = as.integer(atlas),
      period = period_lookup[[as.character(atlas)]],
      bandwidth = label,
      bandwidth_value = fit$bandwidth,
      survey = surveys,
      site = sub("_[123]$", "", surveys),
      X = sdf$X,
      Y = sdf$Y,
      predictor = var,
      driver_label = unname(var_labels[[var]]),
      coefficient = sdf[[var]],
      local_se = if (se_col %in% colnames(sdf)) sdf[[se_col]] else NA_real_,
      local_t = local_t,
      abs_local_t = abs(local_t),
      local_r2 = sdf$Local_R2,
      direction = if_else(local_t >= 0, "Positive", "Negative")
    )
  })
}

extract_all_gwr_local_coefficients <- function(gwr_fits) {
  imap_dfr(gwr_fits, function(response_fits, response) {
    imap_dfr(response_fits, function(atlas_fits, atlas) {
      extract_gwr_local_coefficients(
        fit = atlas_fits[[bandwidth_label]],
        response = response,
        atlas = atlas,
        label = bandwidth_label
      )
    })
  })
}

summarise_all_gwr_models <- function(gwr_fits) {
  imap_dfr(gwr_fits, function(response_fits, response) {
    imap_dfr(response_fits, function(atlas_fits, atlas) {
      summarise_gwr_model(
        fit = atlas_fits[[bandwidth_label]],
        response = response,
        atlas = atlas,
        label = bandwidth_label
      )
    })
  })
}

dominant_driver_table <- function(local_coefficients) {
  local_coefficients |>
    mutate(abs_coefficient = abs(.data$coefficient)) |>
    group_by(.data$response, .data$response_label, .data$atlas, .data$period, .data$survey) |>
    arrange(desc(.data$abs_coefficient), desc(.data$abs_local_t), .by_group = TRUE) |>
    mutate(
      runner_up_driver_label = dplyr::nth(.data$driver_label, 2, default = NA_character_),
      runner_up_abs_coefficient = dplyr::nth(.data$abs_coefficient, 2, default = NA_real_),
      top_to_runner_up_ratio = case_when(
        is.na(.data$runner_up_abs_coefficient) ~ Inf,
        .data$runner_up_abs_coefficient == 0 & .data$abs_coefficient > 0 ~ Inf,
        .data$runner_up_abs_coefficient == 0 ~ NA_real_,
        TRUE ~ .data$abs_coefficient / .data$runner_up_abs_coefficient
      ),
      clear_dominant = .data$top_to_runner_up_ratio >= dominance_ratio_threshold
    ) |>
    slice_head(n = 1) |>
    ungroup() |>
    mutate(
      supported = .data$abs_local_t >= significance_t,
      dominant_driver = case_when(
        !.data$supported ~ no_significant_driver,
        .data$clear_dominant ~ .data$driver_label,
        TRUE ~ multiple_associations_driver
      ),
      direction_plot = if_else(
        .data$supported & .data$clear_dominant,
        .data$direction,
        "Not significant"
      ),
      dominant_driver = factor(.data$dominant_driver, levels = driver_levels),
      direction_plot = factor(.data$direction_plot, levels = direction_levels),
      period = factor(.data$period, levels = period_levels),
      response_label = factor(.data$response_label, levels = unname(response_labels))
    )
}

if (needs_refit) {
  message("Refitting GWR models and refreshing local caches.")
  source(file.path("scripts", "4-data-visualisation", "fig10-gwr-richness-turnover-materials", "figure_data_helpers.R"))

  matching_folders <- figure_model_folders(pattern = pattern, base_dir = base_dir)
  model_nums <- atlas_numbers(matching_folders)

  mods <- load_hmsc_posteriors(matching_folders, base_dir = base_dir)
  designs <- load_hmsc_study_designs(mods)
  preds_y <- load_or_compute_site_predictions(mods, matching_folders, base_dir = base_dir)

  names(mods) <- names(designs) <- names(preds_y) <- model_nums

  add_environment <- function(response_df, model) {
    response_df |>
      left_join(model$XData |> rownames_to_column(var = "survey"), by = "survey")
  }

  predicted_richness <- predicted_richness_frames(preds_y, designs)
  richness_dfs <- map2(predicted_richness, mods, add_environment)
  beta_turnover <- map2(preds_y, designs, site_beta_turnover)
  beta_turnover_dfs <- map2(beta_turnover, mods, add_environment)

  gwr_fits <- list(
    predicted_richness = fit_response_bandwidths(richness_dfs, "richness"),
    beta_turnover = fit_response_bandwidths(beta_turnover_dfs, "beta_turnover")
  )

  saveRDS(gwr_fits, gwr_fit_path)
  saveRDS(
    list(
      richness_dfs = richness_dfs,
      beta_turnover_dfs = beta_turnover_dfs,
      var_labels = var_labels,
      driver_colours = driver_colours,
      bandwidths_m = bandwidths_m,
      turnover_definition = turnover_definition,
      significance_t = significance_t
    ),
    gwr_inputs_path
  )
} else {
  message("Using cached GWR fits: ", gwr_fit_path)
  gwr_fits <- readRDS(gwr_fit_path)
}

gwr_summary <- summarise_all_gwr_models(gwr_fits)
gwr_local <- extract_all_gwr_local_coefficients(gwr_fits)
gwr_dominant <- dominant_driver_table(gwr_local)

readr::write_csv(gwr_summary, gwr_summary_path)
readr::write_csv(gwr_local, gwr_local_path)
readr::write_csv(gwr_dominant, gwr_dominant_path)

shape <- vect(file.path("data-new", "derived", "atlas-grid", "dof-5km-grid", "DK5km_ED50grid_approx_kvadrkod_DOF.shp"))
shape_sf <- st_as_sf(shape)

mainland_bbox <- st_bbox(c(xmin = 400000,
                           xmax = 750000,
                           ymin = 6000000,
                           ymax = 6450000),
                         crs = st_crs(25832))
bornholm_bbox <- st_bbox(c(xmin = 855000, xmax = 905000,
                           ymin = 6100000, ymax = 6160000),
                         crs = st_crs(25832))

mainland_width <- as.numeric(mainland_bbox["xmax"] - mainland_bbox["xmin"])
mainland_height <- as.numeric(mainland_bbox["ymax"] - mainland_bbox["ymin"])
bornholm_width <- as.numeric(bornholm_bbox["xmax"] - bornholm_bbox["xmin"])
bornholm_height <- as.numeric(bornholm_bbox["ymax"] - bornholm_bbox["ymin"])

inset_w <- bornholm_width / mainland_width
inset_h <- bornholm_height / mainland_height

map_theme <- function(legend_position = "none") {
  theme_minimal(base_size = 10) +
    theme(
      legend.position = legend_position,
      legend.box = "vertical",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 8),
      legend.spacing.y = unit(0.08, "cm"),
      legend.key.height = unit(0.32, "cm"),
      legend.key.width = unit(0.42, "cm"),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(0, 0, 0, 0)
    )
}

driver_scales <- function() {
  list(
    scale_colour_manual(
      values = driver_colours,
      limits = driver_levels,
      drop = FALSE,
      guide = "none"
    ),
    scale_shape_manual(
      values = direction_shapes,
      limits = direction_levels,
      breaks = c("Positive", "Negative"),
      drop = FALSE,
      guide = "none"
    )
  )
}

legend_theme <- theme_minimal(base_size = 10) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.32, "cm"),
    legend.key.width = unit(0.42, "cm"),
    legend.spacing.y = unit(0.08, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA)
  )

figure1_style_legend_plot <- function() {
  legend_data <- tibble(
    dominant_driver = factor(driver_levels, levels = driver_levels),
    legend_outline = if_else(
      driver_levels == no_significant_driver,
      no_significant_driver,
      "Other driver"
    )
  )

  ggplot(legend_data, aes(x = .data$dominant_driver, y = 1)) +
    geom_col(
      aes(fill = .data$dominant_driver, colour = .data$legend_outline),
      width = 0.8
    ) +
    geom_point(
      data = tibble(direction_plot = factor(c("Positive", "Negative"), levels = direction_levels), x = 1:2, y = 1),
      aes(x = .data$x, y = .data$y, shape = .data$direction_plot),
      inherit.aes = FALSE,
      colour = "black",
      size = 2.7
    ) +
    scale_fill_manual(
      values = driver_colours,
      limits = driver_levels,
      drop = FALSE,
      name = "Variable"
    ) +
    scale_colour_manual(
      values = c("Other driver" = NA, "No significant driver" = no_significant_outline),
      guide = "none"
    ) +
    scale_shape_manual(
      values = direction_shapes,
      limits = direction_levels,
      breaks = c("Positive", "Negative"),
      drop = FALSE,
      name = "Direction"
    ) +
    guides(
      fill = guide_legend(order = 1, ncol = 1, byrow = TRUE),
      shape = guide_legend(order = 2, ncol = 1, override.aes = list(size = 3))
    ) +
    legend_theme
}

driver_point_data <- function(df) {
  plot_data <- shape_sf |>
    left_join(
      df |>
        select(all_of(c("site", "dominant_driver", "direction_plot"))),
      by = c("kvadratkod" = "site")
    ) |>
    filter(!is.na(.data$dominant_driver))

  suppressWarnings(st_point_on_surface(plot_data))
}

plot_dominant_driver_base <- function(df, bbox = mainland_bbox,
                                      show_legend = FALSE, border = FALSE,
                                      point_size = 1.34) {
  plot_data <- driver_point_data(df)
  no_driver_points <- plot_data |>
    filter(as.character(.data$dominant_driver) == no_significant_driver)
  driver_points <- plot_data |>
    filter(as.character(.data$dominant_driver) != no_significant_driver)

  ggplot() +
    geom_sf(
      data = no_driver_points,
      shape = 21,
      fill = driver_colours[[no_significant_driver]],
      colour = no_significant_outline,
      size = point_size,
      alpha = 0.95,
      stroke = 0.42
    ) +
    geom_sf(
      data = driver_points,
      aes(colour = .data$dominant_driver, shape = .data$direction_plot),
      size = point_size,
      alpha = 0.95,
      stroke = 0.55
    ) +
    driver_scales() +
    coord_sf(
      xlim = c(bbox["xmin"], bbox["xmax"]),
      ylim = c(bbox["ymin"], bbox["ymax"]),
      expand = FALSE
    ) +
    labs(x = NULL, y = NULL) +
    map_theme(if (show_legend) "right" else "none") +
    theme(
      panel.border = if (border) {
        element_rect(colour = "grey30", fill = NA, linewidth = 0.5)
      } else {
        element_blank()
      }
    )
}

facet_strip <- function(label) {
  ggdraw() +
    draw_label(
      label,
      x = 0.5,
      y = 0.52,
      hjust = 0.5,
      vjust = 0.5,
      fontface = "bold",
      size = 9,
      colour = "black"
    ) +
    theme(
      plot.background = element_rect(fill = "grey92", colour = "white"),
      plot.margin = margin(0, 0, 1.2, 0)
    )
}

plot_map_with_inset <- function(df) {
  p_main <- plot_dominant_driver_base(df, bbox = mainland_bbox)
  p_inset <- plot_dominant_driver_base(
    df,
    bbox = bornholm_bbox,
    border = TRUE,
    point_size = 1.45
  ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.7),
      plot.margin = margin(0, 0, 0, 0)
    )

  ggdraw(p_main) +
    draw_plot(
      p_inset,
      x = 1 - inset_w - 0.13,
      y = 1 - inset_h - 0.13,
      width = inset_w,
      height = inset_h
    ) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

period_panel <- function(df, period_name) {
  row_df <- df |>
    filter(.data$period == period_name)

  plot_grid(
    facet_strip(period_name),
    plot_map_with_inset(row_df),
    ncol = 1,
    rel_heights = c(0.072, 1)
  )
}

panel_title <- function(label) {
  ggdraw() +
    draw_label(
      label,
      x = 0,
      hjust = 0,
      fontface = "bold",
      size = 12,
      colour = "black"
    ) +
    theme(plot.margin = margin(0, 0, 2, 0))
}

response_map_row <- function(response_name) {
  response_df <- gwr_dominant |>
    filter(.data$response == response_name)

  response_maps <- map(period_levels, ~ period_panel(response_df, .x))

  plot_grid(
    plotlist = response_maps,
    nrow = 1,
    align = "hv",
    axis = "tblr",
    rel_widths = c(1, 1, 1)
  )
}

figure_legend <- get_legend(figure1_style_legend_plot())

panel_a <- plot_grid(
  panel_title("A  Species richness"),
  response_map_row("predicted_richness"),
  ncol = 1,
  rel_heights = c(0.07, 1)
)

panel_b <- plot_grid(
  panel_title("B  Beta turnover"),
  response_map_row("beta_turnover"),
  ncol = 1,
  rel_heights = c(0.07, 1)
)

map_panels <- plot_grid(
  panel_a,
  panel_b,
  ncol = 1,
  rel_heights = c(1, 1),
  align = "v",
  axis = "lr"
)

fig10_gwr_richness_turnover <- plot_grid(
  map_panels,
  figure_legend,
  ncol = 2,
  rel_widths = c(1, 0.18),
  align = "h",
  axis = "tb"
) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(png_path, fig10_gwr_richness_turnover, width = 12.8, height = 9.8, units = "in", dpi = 300, bg = "white")

message("Saved: ", gwr_summary_path)
message("Saved: ", gwr_local_path)
message("Saved: ", gwr_dominant_path)
message("Saved: ", png_path)
