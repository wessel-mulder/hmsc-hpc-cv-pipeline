# Purpose:
# Fit geographically weighted regressions for modelled species richness, one
# predictor at a time, across the three atlas periods. The script stores the
# fitted objects and tabular summaries locally so the exploratory report and
# plotting workflow can be refreshed without re-running the expensive GWR fits.

rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
force_refit <- "--force-refit" %in% args

base_dir <- "./HmscOutputs"
pattern <- "2026-03-13"
bandwidth_label <- "fixed_100km"
bandwidth_m <- 100000
significance_t <- 1.96
local_t_limit <- 4
site_filter_label <- "common-sites"
site_filter_description <- "Grid cells present in all three atlas periods"

gwr_out_dir <- file.path("notebooks", "exploratory", "outputs", "gwr")
gwr_model_dir <- file.path(gwr_out_dir, "models")
figure_out_dir <- file.path("misc-figures", "outputs", "main")

shared_inputs_path <- file.path(gwr_model_dir, paste0(pattern, "-gwr-analysis-inputs-100km.rds"))
gwr_fit_path <- file.path(gwr_model_dir, paste0(pattern, "-richness-univariate-gwr-", site_filter_label, "-fits-100km.rds"))
gwr_inputs_path <- file.path(gwr_model_dir, paste0(pattern, "-richness-univariate-gwr-", site_filter_label, "-analysis-inputs-100km.rds"))
gwr_summary_path <- file.path(gwr_out_dir, paste0(pattern, "-richness-univariate-gwr-", site_filter_label, "-summary-100km.csv"))
gwr_local_path <- file.path(gwr_out_dir, paste0(pattern, "-richness-univariate-gwr-", site_filter_label, "-local-results-100km.csv"))
gwr_comparison_path <- file.path(gwr_out_dir, paste0(pattern, "-richness-univariate-gwr-", site_filter_label, "-model-comparison-100km.csv"))
gwr_strongest_path <- file.path(gwr_out_dir, paste0(pattern, "-richness-univariate-gwr-", site_filter_label, "-strongest-associations-100km.csv"))

model_comparison_png_path <- file.path(
  figure_out_dir,
  paste0(pattern, "-fig11-gwr-richness-univariate-", site_filter_label, "-model-comparison.png")
)
strongest_map_png_path <- file.path(
  figure_out_dir,
  paste0(pattern, "-fig11-gwr-richness-univariate-", site_filter_label, "-strongest-associations.png")
)
predictor_map_prefix <- file.path(
  figure_out_dir,
  paste0(pattern, "-fig11-gwr-richness-univariate-", site_filter_label, "-map-")
)

dir.create(gwr_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gwr_model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_out_dir, recursive = TRUE, showWarnings = FALSE)

needs_model_fit <- force_refit || !file.exists(gwr_fit_path)
needs_input_refresh <- force_refit || !file.exists(gwr_inputs_path)
needs_hmsc_inputs <- needs_input_refresh && !file.exists(shared_inputs_path)

load_required_packages <- function(needs_model_fit, needs_hmsc_inputs) {
  core_packages <- c("tidyverse", "cowplot", "scales", "sp")
  refit_packages <- if (needs_model_fit) "GWmodel" else character()
  hmsc_packages <- if (needs_hmsc_inputs) "Hmsc" else character()
  required_packages <- unique(c(core_packages, refit_packages, hmsc_packages))

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
    if (needs_model_fit) {
      library(GWmodel)
    }
    if (needs_hmsc_inputs) {
      library(Hmsc)
    }
  })
}

load_required_packages(needs_model_fit, needs_hmsc_inputs)

period_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")
period_levels <- unname(period_lookup)

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
  "No supported association" = "#d9d9d9"
)

direction_levels <- c("Positive", "Negative", "Not supported")
direction_shapes <- c("Positive" = 19, "Negative" = 4, "Not supported" = 19)

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

model_specs <- function(df) {
  predictors <- predictor_names(df)

  bind_rows(
    tibble(
      model_id = "intercept_only",
      model_type = "intercept_only",
      predictor = NA_character_,
      model_label = "Intercept only"
    ),
    tibble(
      model_id = predictors,
      model_type = "single_predictor",
      predictor = predictors,
      model_label = unname(var_labels[predictors])
    )
  )
}

filter_to_common_sites <- function(richness_dfs) {
  common_sites <- richness_dfs |>
    map(~ unique(.x$site)) |>
    reduce(intersect) |>
    sort()

  if (length(common_sites) == 0) {
    stop("No grid-cell sites are shared by all atlas periods.", call. = FALSE)
  }

  site_summary <- imap_dfr(richness_dfs, function(df, atlas) {
    tibble(
      atlas = as.integer(atlas),
      period = period_lookup[[as.character(atlas)]],
      original_n_sites = n_distinct(df$site),
      common_n_sites = length(common_sites),
      dropped_n_sites = n_distinct(df$site) - length(common_sites)
    )
  })

  filtered_dfs <- map(richness_dfs, function(df) {
    df |>
      filter(.data$site %in% common_sites) |>
      arrange(.data$site)
  })

  list(
    richness_dfs = filtered_dfs,
    common_sites = common_sites,
    site_summary = site_summary
  )
}

make_common_site_inputs <- function(richness_dfs, source) {
  common_site_data <- filter_to_common_sites(richness_dfs)

  message(
    "Keeping ",
    length(common_site_data$common_sites),
    " grid cells present in all three atlas periods."
  )

  list(
    richness_dfs = common_site_data$richness_dfs,
    common_sites = common_site_data$common_sites,
    common_site_summary = common_site_data$site_summary,
    site_filter = site_filter_label,
    site_filter_description = site_filter_description,
    var_labels = var_labels,
    bandwidth_m = bandwidth_m,
    bandwidth_label = bandwidth_label,
    significance_t = significance_t,
    source = source
  )
}

load_or_create_richness_inputs <- function(force_refresh = FALSE) {
  if (!force_refresh && file.exists(gwr_inputs_path)) {
    message("Using cached richness-univariate inputs: ", gwr_inputs_path)
    return(readRDS(gwr_inputs_path))
  }

  if (file.exists(shared_inputs_path)) {
    message("Using existing Fig. 10 richness inputs: ", shared_inputs_path)
    shared_inputs <- readRDS(shared_inputs_path)
    if (!"richness_dfs" %in% names(shared_inputs)) {
      stop("Shared GWR input cache does not contain `richness_dfs`: ", shared_inputs_path, call. = FALSE)
    }

    richness_inputs <- make_common_site_inputs(
      richness_dfs = shared_inputs$richness_dfs,
      source = shared_inputs_path
    )
    saveRDS(richness_inputs, gwr_inputs_path)
    return(richness_inputs)
  }

  message("Shared richness inputs are missing, rebuilding them from HMSC predictions.")
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

  richness_inputs <- make_common_site_inputs(
    richness_dfs = map2(predicted_richness_frames(preds_y, designs), mods, add_environment),
    source = "HMSC predictions"
  )

  saveRDS(richness_inputs, gwr_inputs_path)
  richness_inputs
}

prepare_gwr_data <- function(df, response, predictor = NULL) {
  required_cols <- c("survey", "site", "X", "Y", response, predictor)
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop("Missing GWR input columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  out <- df |>
    select(all_of(required_cols)) |>
    drop_na()

  # Each univariate predictor is standardized within atlas period. That makes
  # coefficients comparable as richness change per one local-period SD.
  if (!is.null(predictor)) {
    out[[predictor]] <- as.numeric(scale(out[[predictor]]))
  }

  out <- as.data.frame(out)
  rownames(out) <- out$survey
  coordinates(out) <- c("X", "Y")
  out
}

fit_single_gwr_model <- function(df, response, spec) {
  predictor <- spec$predictor[[1]]
  formula <- if (is.na(predictor)) {
    as.formula(paste(response, "~ 1"))
  } else {
    reformulate(predictor, response = response)
  }

  gwr_data <- prepare_gwr_data(
    df = df,
    response = response,
    predictor = if (is.na(predictor)) NULL else predictor
  )

  gwr_model <- GWmodel::gwr.basic(
    formula,
    data = gwr_data,
    bw = bandwidth_m,
    kernel = "bisquare",
    adaptive = FALSE
  )

  list(
    mod = gwr_model,
    model_id = spec$model_id[[1]],
    model_type = spec$model_type[[1]],
    predictor = predictor,
    model_label = spec$model_label[[1]],
    formula = formula,
    bandwidth = bandwidth_m,
    adaptive = FALSE
  )
}

fit_all_richness_models <- function(richness_dfs) {
  imap(richness_dfs, function(df, atlas) {
    specs <- model_specs(df)

    specs |>
      split(specs$model_id) |>
      imap(function(spec, model_id) {
        message("Fitting richness GWR: atlas ", atlas, ", model ", model_id)
        fit_single_gwr_model(df = df, response = "richness", spec = spec)
      })
  })
}

safe_stat <- function(x, fun) {
  if (all(is.na(x))) return(NA_real_)
  fun(x, na.rm = TRUE)
}

summarise_gwr_model <- function(fit, atlas) {
  sdf <- as.data.frame(fit$mod$SDF)
  diagnostics <- fit$mod$GW.diagnostic
  local_r2 <- if ("Local_R2" %in% colnames(sdf)) sdf$Local_R2 else rep(NA_real_, nrow(sdf))

  tibble(
    response = "richness",
    site_filter = site_filter_label,
    atlas = as.integer(atlas),
    period = period_lookup[[as.character(atlas)]],
    model_id = fit$model_id,
    model_type = fit$model_type,
    predictor = fit$predictor,
    model_label = fit$model_label,
    formula = paste(deparse(fit$formula), collapse = " "),
    bandwidth = bandwidth_label,
    adaptive = fit$adaptive,
    bandwidth_value = fit$bandwidth,
    n_sites = nrow(sdf),
    aicc = diagnostics$AICc,
    aic = diagnostics$AIC,
    bic = diagnostics$BIC,
    rss = diagnostics$RSS.gw,
    enp = diagnostics$enp,
    edf = diagnostics$edf,
    gw_r2 = diagnostics$gw.R2,
    gw_r2_adj = diagnostics$gwR2.adj,
    local_r2_mean = safe_stat(local_r2, mean),
    local_r2_median = safe_stat(local_r2, median),
    local_r2_min = safe_stat(local_r2, min),
    local_r2_max = safe_stat(local_r2, max)
  )
}

summarise_all_gwr_models <- function(gwr_fits) {
  imap_dfr(gwr_fits, function(atlas_fits, atlas) {
    imap_dfr(atlas_fits, function(fit, model_id) {
      summarise_gwr_model(fit = fit, atlas = atlas)
    })
  })
}

extract_local_results <- function(fit, atlas) {
  sdf <- as.data.frame(fit$mod$SDF)
  surveys <- rownames(sdf)
  coefficient_col <- if (fit$model_type == "intercept_only") "Intercept" else fit$predictor
  se_col <- paste0(coefficient_col, "_SE")
  t_col <- paste0(coefficient_col, "_TV")

  local_t <- if (t_col %in% colnames(sdf)) sdf[[t_col]] else rep(NA_real_, nrow(sdf))
  local_r2 <- if ("Local_R2" %in% colnames(sdf)) sdf$Local_R2 else rep(NA_real_, nrow(sdf))
  supported <- fit$model_type == "single_predictor" & abs(local_t) >= significance_t

  tibble(
    response = "richness",
    site_filter = site_filter_label,
    atlas = as.integer(atlas),
    period = period_lookup[[as.character(atlas)]],
    bandwidth = bandwidth_label,
    bandwidth_value = fit$bandwidth,
    model_id = fit$model_id,
    model_type = fit$model_type,
    predictor = fit$predictor,
    model_label = fit$model_label,
    survey = surveys,
    site = sub("_[123]$", "", surveys),
    X = sdf$X,
    Y = sdf$Y,
    observed = sdf$y,
    fitted = sdf$yhat,
    residual = sdf$residual,
    coefficient = sdf[[coefficient_col]],
    local_se = if (se_col %in% colnames(sdf)) sdf[[se_col]] else NA_real_,
    local_t = local_t,
    abs_local_t = abs(local_t),
    local_r2 = local_r2,
    supported = supported,
    direction = case_when(
      !supported ~ "Not supported",
      local_t >= 0 ~ "Positive",
      TRUE ~ "Negative"
    )
  )
}

extract_all_local_results <- function(gwr_fits) {
  imap_dfr(gwr_fits, function(atlas_fits, atlas) {
    imap_dfr(atlas_fits, function(fit, model_id) {
      extract_local_results(fit = fit, atlas = atlas)
    })
  })
}

model_comparison_table <- function(gwr_summary) {
  gwr_summary |>
    group_by(.data$atlas, .data$period) |>
    mutate(
      delta_aicc = .data$aicc - min(.data$aicc, na.rm = TRUE),
      aicc_weight = exp(-0.5 * .data$delta_aicc) / sum(exp(-0.5 * .data$delta_aicc), na.rm = TRUE),
      period_rank = min_rank(.data$aicc),
      beats_intercept = .data$aicc < .data$aicc[.data$model_id == "intercept_only"][[1]]
    ) |>
    ungroup() |>
    arrange(.data$atlas, .data$period_rank, .data$model_id)
}

strongest_univariate_association <- function(local_results) {
  local_results |>
    filter(.data$model_type == "single_predictor") |>
    group_by(.data$atlas, .data$period, .data$survey) |>
    slice_max(.data$abs_local_t, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      supported = .data$abs_local_t >= significance_t,
      dominant_driver = if_else(.data$supported, .data$model_label, "No supported association"),
      direction_plot = if_else(.data$supported, .data$direction, "Not supported"),
      dominant_driver = factor(.data$dominant_driver, levels = names(driver_colours)),
      direction_plot = factor(.data$direction_plot, levels = direction_levels),
      period = factor(.data$period, levels = period_levels)
    )
}

if (needs_input_refresh || needs_model_fit) {
  richness_inputs <- load_or_create_richness_inputs(force_refresh = needs_input_refresh)
} else {
  richness_inputs <- readRDS(gwr_inputs_path)
}

common_n_sites <- length(richness_inputs$common_sites)

if (needs_model_fit) {
  message("Refitting richness-only univariate GWR models.")
  gwr_fits <- fit_all_richness_models(richness_inputs$richness_dfs)
  saveRDS(gwr_fits, gwr_fit_path)
} else {
  message("Using cached richness-only univariate GWR fits: ", gwr_fit_path)
  gwr_fits <- readRDS(gwr_fit_path)
}

gwr_summary <- summarise_all_gwr_models(gwr_fits)
gwr_local <- extract_all_local_results(gwr_fits)
gwr_comparison <- model_comparison_table(gwr_summary)
gwr_strongest <- strongest_univariate_association(gwr_local)

readr::write_csv(gwr_summary, gwr_summary_path)
readr::write_csv(gwr_local, gwr_local_path)
readr::write_csv(gwr_comparison, gwr_comparison_path)
readr::write_csv(gwr_strongest, gwr_strongest_path)

mainland_bbox <- c(xmin = 400000, xmax = 750000, ymin = 6000000, ymax = 6450000)
bornholm_bbox <- c(xmin = 855000, xmax = 905000, ymin = 6100000, ymax = 6160000)
mainland_width <- mainland_bbox[["xmax"]] - mainland_bbox[["xmin"]]
mainland_height <- mainland_bbox[["ymax"]] - mainland_bbox[["ymin"]]
bornholm_width <- bornholm_bbox[["xmax"]] - bornholm_bbox[["xmin"]]
bornholm_height <- bornholm_bbox[["ymax"]] - bornholm_bbox[["ymin"]]
bornholm_inset_width <- bornholm_width / mainland_width
bornholm_inset_height <- bornholm_height / mainland_height

slugify <- function(x) {
  x |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]+", "-") |>
    str_replace_all("^-|-$", "")
}

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

plot_local_t_base <- function(df, title = NULL, bbox = mainland_bbox,
                              show_legend = FALSE, border = FALSE,
                              point_size = 1.25) {
  df <- df |>
    mutate(local_t_plot = pmax(pmin(.data$local_t, local_t_limit), -local_t_limit))

  ggplot() +
    geom_point(
      data = filter(df, !.data$supported),
      aes(x = .data$X, y = .data$Y),
      colour = "#d9d9d9",
      size = point_size,
      alpha = 0.9
    ) +
    geom_point(
      data = filter(df, .data$supported),
      aes(x = .data$X, y = .data$Y, colour = .data$local_t_plot),
      size = point_size,
      alpha = 0.95
    ) +
    scale_colour_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-local_t_limit, local_t_limit),
      oob = scales::squish,
      name = "Local t-value\n(capped)"
    ) +
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

plot_local_t_map_with_inset <- function(df, title) {
  p_main <- plot_local_t_base(df, title = title, bbox = mainland_bbox)
  p_inset <- plot_local_t_base(
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

plot_predictor_map_panel <- function(local_results, predictor_name) {
  predictor_label <- unname(var_labels[[predictor_name]])
  predictor_df <- local_results |>
    filter(.data$predictor == predictor_name) |>
    mutate(period = factor(.data$period, levels = period_levels))

  period_maps <- map(period_levels, function(period_name) {
    plot_local_t_map_with_inset(
      filter(predictor_df, .data$period == period_name),
      title = period_name
    )
  })

  legend_plot <- plot_local_t_base(predictor_df, show_legend = TRUE) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold")
    )

  plot_grid(
    panel_title(paste0("Richness association: ", predictor_label, " (common grid cells)")),
    plot_grid(plotlist = period_maps, nrow = 1, align = "h", rel_widths = c(1, 1, 1)),
    get_legend(legend_plot),
    ncol = 1,
    rel_heights = c(0.11, 1, 0.16)
  ) +
    theme(plot.background = element_rect(fill = "white", colour = NA))
}

strongest_driver_scales <- function() {
  list(
    scale_colour_manual(
      values = driver_colours,
      limits = names(driver_colours),
      drop = FALSE,
      name = "Strongest univariate\nassociation",
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

plot_strongest_base <- function(df, title = NULL, bbox = mainland_bbox,
                                show_legend = FALSE, border = FALSE,
                                point_size = 1.25) {
  ggplot(df) +
    geom_point(
      aes(x = .data$X, y = .data$Y, colour = .data$dominant_driver, shape = .data$direction_plot),
      size = point_size,
      alpha = 0.95,
      stroke = 0.55
    ) +
    strongest_driver_scales() +
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

plot_strongest_map_with_inset <- function(df, title) {
  p_main <- plot_strongest_base(df, title = title, bbox = mainland_bbox)
  p_inset <- plot_strongest_base(
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

plot_strongest_panel <- function(strongest_df) {
  period_maps <- map(period_levels, function(period_name) {
    plot_strongest_map_with_inset(
      filter(strongest_df, .data$period == period_name),
      title = period_name
    )
  })

  legend_plot <- plot_strongest_base(strongest_df, show_legend = TRUE) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 7.8),
      legend.title = element_text(size = 9, face = "bold")
    )

  plot_grid(
    panel_title(paste0("Strongest richness association from univariate GWRs (", common_n_sites, " common grid cells)")),
    plot_grid(plotlist = period_maps, nrow = 1, align = "h", rel_widths = c(1, 1, 1)),
    get_legend(legend_plot),
    ncol = 1,
    rel_heights = c(0.11, 1, 0.24)
  ) +
    theme(plot.background = element_rect(fill = "white", colour = NA))
}

plot_model_comparison <- function(comparison_df) {
  plot_df <- comparison_df |>
    mutate(
      period = factor(.data$period, levels = period_levels),
      model_label_wrapped = str_wrap(.data$model_label, width = 24),
      model_label_wrapped = factor(
        .data$model_label_wrapped,
        levels = rev(unique(.data$model_label_wrapped[order(.data$model_type, .data$model_label)]))
      ),
      delta_aicc_capped = pmin(.data$delta_aicc, 50),
      tile_label = if_else(
        .data$delta_aicc == 0,
        paste0("best\nR2=", sprintf("%.2f", .data$local_r2_mean)),
        paste0("dAICc=", sprintf("%.1f", .data$delta_aicc), "\nR2=", sprintf("%.2f", .data$local_r2_mean))
      ),
      label_colour = if_else(.data$delta_aicc_capped <= 12, "white", "grey15")
    )

  ggplot(plot_df, aes(x = .data$period, y = .data$model_label_wrapped)) +
    geom_tile(aes(fill = .data$delta_aicc_capped), colour = "white", linewidth = 0.8) +
    geom_text(aes(label = .data$tile_label, colour = .data$label_colour), size = 2.8, lineheight = 0.9) +
    scale_fill_gradient(
      low = "#2166ac",
      high = "#f7f7f7",
      limits = c(0, 50),
      name = "Delta AICc\n(capped at 50)"
    ) +
    scale_colour_identity() +
    labs(
      title = "Richness-only single-predictor GWR support",
      subtitle = paste0(
        "Common grid cells only (n = ", common_n_sites,
        " per period); one fixed 100 km bisquare GWR per atlas period and predictor."
      ),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

model_comparison_plot <- plot_model_comparison(gwr_comparison)
strongest_map_plot <- plot_strongest_panel(gwr_strongest)

ggsave(model_comparison_png_path, model_comparison_plot, width = 8.5, height = 5.6, units = "in", dpi = 300, bg = "white")
ggsave(strongest_map_png_path, strongest_map_plot, width = 13.5, height = 7.2, units = "in", dpi = 300, bg = "white")

predictor_map_paths <- map_chr(names(var_labels), function(predictor_name) {
  map_path <- paste0(predictor_map_prefix, slugify(predictor_name), ".png")
  predictor_map <- plot_predictor_map_panel(gwr_local, predictor_name)
  ggsave(map_path, predictor_map, width = 13.5, height = 7.2, units = "in", dpi = 300, bg = "white")
  map_path
})

message("Saved: ", gwr_fit_path)
message("Saved: ", gwr_inputs_path)
message("Saved: ", gwr_summary_path)
message("Saved: ", gwr_local_path)
message("Saved: ", gwr_comparison_path)
message("Saved: ", gwr_strongest_path)
message("Saved: ", model_comparison_png_path)
message("Saved: ", strongest_map_png_path)
walk(predictor_map_paths, ~ message("Saved: ", .x))
