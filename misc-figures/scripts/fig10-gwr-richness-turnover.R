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

gwr_out_dir <- file.path("notebooks", "exploratory", "outputs", "gwr")
gwr_model_dir <- file.path(gwr_out_dir, "models")
figure_out_dir <- file.path("misc-figures", "outputs", "main")

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
  core_packages <- c("tidyverse", "cowplot", "scales", "sp")
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

driver_colours <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban (% coverage)" = "#4d4d4d",
  "Cropland (% coverage)" = "goldenrod1",
  "Pasture (% coverage)" = "darkorange",
  "Forest (% coverage)" = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2",
  "No significant driver" = "#d9d9d9"
)

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
    group_by(.data$response, .data$response_label, .data$atlas, .data$period, .data$survey) |>
    slice_max(.data$abs_local_t, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      supported = .data$abs_local_t >= significance_t,
      dominant_driver = if_else(.data$supported, .data$driver_label, "No significant driver"),
      direction_plot = if_else(.data$supported, .data$direction, "Not significant"),
      dominant_driver = factor(.data$dominant_driver, levels = driver_levels),
      direction_plot = factor(.data$direction_plot, levels = direction_levels),
      period = factor(.data$period, levels = period_levels),
      response_label = factor(.data$response_label, levels = unname(response_labels))
    )
}

if (needs_refit) {
  message("Refitting GWR models and refreshing local caches.")
  source(file.path("support_scripts", "figure_data_helpers.R"))

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

mainland_bbox <- c(xmin = 400000, xmax = 750000, ymin = 6000000, ymax = 6450000)
bornholm_bbox <- c(xmin = 855000, xmax = 905000, ymin = 6100000, ymax = 6160000)
mainland_width <- mainland_bbox[["xmax"]] - mainland_bbox[["xmin"]]
mainland_height <- mainland_bbox[["ymax"]] - mainland_bbox[["ymin"]]
bornholm_width <- bornholm_bbox[["xmax"]] - bornholm_bbox[["xmin"]]
bornholm_height <- bornholm_bbox[["ymax"]] - bornholm_bbox[["ymin"]]
bornholm_inset_width <- bornholm_width / mainland_width
bornholm_inset_height <- bornholm_height / mainland_height

map_theme <- function(legend_position = "none") {
  theme_minimal(base_size = 10) +
    theme(
      legend.position = legend_position,
      legend.box = "vertical",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 8),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
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
      name = "Dominant driver",
      guide = guide_legend(
        nrow = 3,
        byrow = TRUE,
        override.aes = list(size = 3.2, shape = 19)
      )
    ),
    scale_shape_manual(
      values = direction_shapes,
      limits = direction_levels,
      breaks = c("Positive", "Negative"),
      drop = FALSE,
      name = "Direction",
      guide = guide_legend(nrow = 1, override.aes = list(size = 3))
    )
  )
}

plot_dominant_driver_base <- function(df, title = NULL, bbox = mainland_bbox,
                                      show_legend = FALSE, border = FALSE,
                                      point_size = 1.25) {
  ggplot(df) +
    geom_point(
      aes(x = .data$X, y = .data$Y, colour = .data$dominant_driver, shape = .data$direction_plot),
      size = point_size,
      alpha = 0.95,
      stroke = 0.55
    ) +
    driver_scales() +
    coord_fixed(
      xlim = c(bbox[["xmin"]], bbox[["xmax"]]),
      ylim = c(bbox[["ymin"]], bbox[["ymax"]]),
      expand = FALSE
    ) +
    labs(title = title) +
    map_theme(if (show_legend) "bottom" else "none") +
    theme(
      panel.border = if (border) {
        element_rect(colour = "grey35", fill = NA, linewidth = 0.45)
      } else {
        element_blank()
      }
    )
}

plot_map_with_inset <- function(df, title) {
  p_main <- plot_dominant_driver_base(df, title = title, bbox = mainland_bbox)
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

panel_title <- function(label) {
  ggdraw() +
    draw_label(
      label,
      x = 0,
      hjust = 0,
      fontface = "bold",
      size = 15
    ) +
    theme(plot.margin = margin(0, 0, 2, 0))
}

response_map_row <- function(response_name) {
  response_maps <- map(period_levels, function(period_name) {
    row_df <- gwr_dominant |>
      filter(.data$response == response_name, .data$period == period_name)

    plot_map_with_inset(row_df, title = period_name)
  })

  plot_grid(plotlist = response_maps, nrow = 1, align = "h", rel_widths = c(1, 1, 1))
}

figure_legend <- get_legend(
  plot_dominant_driver_base(gwr_dominant, show_legend = TRUE) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 7.8),
      legend.title = element_text(size = 9, face = "bold"),
      legend.margin = margin(0, 0, 0, 0)
    )
)

panel_a <- plot_grid(
  panel_title("A  Predicted richness"),
  response_map_row("predicted_richness"),
  ncol = 1,
  rel_heights = c(0.08, 1)
)

panel_b <- plot_grid(
  panel_title("B  All-site beta turnover"),
  response_map_row("beta_turnover"),
  ncol = 1,
  rel_heights = c(0.08, 1)
)

fig10_gwr_richness_turnover <- plot_grid(
  panel_a,
  panel_b,
  figure_legend,
  ncol = 1,
  rel_heights = c(1, 1, 0.28)
) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(png_path, fig10_gwr_richness_turnover, width = 13.5, height = 11.8, units = "in", dpi = 300, bg = "white")

message("Saved: ", gwr_summary_path)
message("Saved: ", gwr_local_path)
message("Saved: ", gwr_dominant_path)
message("Saved: ", png_path)
