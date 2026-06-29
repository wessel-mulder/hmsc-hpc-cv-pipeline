# Visualise HMSC Eta spatial latent factors for the March 2026 atlas models.
#
# The fitted HMSC objects store Eta as site-by-latent-factor posterior samples.
# Eta row names are blank after fitting, so this script intentionally joins Eta
# values to the random-level coordinates by row order and uses those coordinate
# row names as the grid-cell/site IDs.

remove(list = ls())

# Keep the project-specific library path first, matching the main HMSC scripts.
.libPaths(c("~/Rlibs", .libPaths()))

suppressPackageStartupMessages({
  library(Hmsc)
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(scales)
  library(sf)
  library(tibble)
  library(tidyr)
})

model_prefix <- "2026-03-13_06-58-56"
output_prefix <- "2026-03-13"

model_specs <- tibble(
  atlas = c("1", "2", "3"),
  period = c("1970s", "1990s", "2010s"),
  model_folder = paste0(
    model_prefix,
    "_Atlas",
    atlas,
    "_MinOccs5_CoverageGoodAverage"
  ),
  expected_sites = c(1776L, 1850L, 1904L)
)

fitted_file <- function(model_folder) {
  file.path(
    "HmscOutputs",
    model_folder,
    "Models",
    "Fitted",
    "HPC_samples_0250_thin_100_chains_4.Rdata"
  )
}

grid_path <- file.path(
  "Data",
  "data",
  "1_preprocessing",
  "atlas-grids",
  "grids-ocean-thresholds",
  "grids_ocean_thresholds.shp"
)

plot_dir <- file.path(
  "notebooks",
  "exploratory",
  "outputs",
  "eta-spatial-loadings",
  "plots"
)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

mainland_bbox <- c(xmin = 400000, xmax = 750000, ymin = 6000000, ymax = 6450000)
bornholm_bbox <- c(xmin = 855000, xmax = 905000, ymin = 6100000, ymax = 6160000)

eta_levels <- paste("Eta", seq_len(5))

eta_fill_scale <- function(eta_limit) {
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-eta_limit, eta_limit),
    oob = squish,
    name = "Posterior\nmean Eta"
  )
}

eta_alpha_scale <- function() {
  scale_alpha_continuous(
    range = c(0.25, 1),
    limits = c(0.5, 1),
    oob = squish,
    breaks = c(0.5, 0.75, 0.95, 1),
    labels = percent_format(accuracy = 1),
    name = "Sign\nsupport"
  )
}

eta_map_theme <- function(legend_position = "none") {
  theme_void(base_size = 9) +
    theme(
      legend.position = legend_position,
      legend.box = "horizontal",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(0, 1, 0, 1)
    )
}

load_fitted_hmsc_model <- function(model_file) {
  if (!file.exists(model_file)) {
    stop("Missing fitted model file: ", model_file, call. = FALSE)
  }

  # Load into a private environment so this script never depends on global names
  # left behind by a previous model.
  model_env <- new.env(parent = emptyenv())
  load(model_file, envir = model_env)

  if (!exists("fitted_model", envir = model_env, inherits = FALSE)) {
    stop("Expected object `fitted_model` in: ", model_file, call. = FALSE)
  }

  fitted_model <- get("fitted_model", envir = model_env, inherits = FALSE)
  if (!"posteriors" %in% names(fitted_model)) {
    stop("Expected `fitted_model$posteriors` in: ", model_file, call. = FALSE)
  }

  fitted_model$posteriors
}

matrix_to_eta_long <- function(mat, site_ids, value_name) {
  if (ncol(mat) != length(eta_levels)) {
    stop(
      "Expected ",
      length(eta_levels),
      " Eta columns, found ",
      ncol(mat),
      ".",
      call. = FALSE
    )
  }

  colnames(mat) <- eta_levels

  as_tibble(mat) |>
    mutate(site = site_ids, .before = 1) |>
    pivot_longer(
      cols = all_of(eta_levels),
      names_to = "eta_factor",
      values_to = value_name
    )
}

extract_eta_table <- function(model, atlas, period, expected_sites) {
  eta_est <- getPostEstimate(model, parName = "Eta")
  coords <- as.data.frame(model$ranLevels$site$s)
  coords <- coords |>
    rownames_to_column("site") |>
    as_tibble()

  if (nrow(coords) != nrow(eta_est$mean)) {
    stop(
      "Eta rows and random-level coordinate rows differ for Atlas ",
      atlas,
      ": ",
      nrow(eta_est$mean),
      " Eta rows vs ",
      nrow(coords),
      " coordinate rows.",
      call. = FALSE
    )
  }

  if (nrow(coords) != expected_sites) {
    stop(
      "Atlas ",
      atlas,
      " expected ",
      expected_sites,
      " modelled sites, found ",
      nrow(coords),
      ".",
      call. = FALSE
    )
  }

  eta_mean <- matrix_to_eta_long(eta_est$mean, coords$site, "eta_mean")
  eta_support <- matrix_to_eta_long(eta_est$support, coords$site, "support_positive")
  eta_support_neg <- matrix_to_eta_long(eta_est$supportNeg, coords$site, "support_negative")

  eta_mean |>
    left_join(eta_support, by = c("site", "eta_factor")) |>
    left_join(eta_support_neg, by = c("site", "eta_factor")) |>
    left_join(coords, by = "site") |>
    mutate(
      atlas = atlas,
      period = period,
      eta_factor = factor(.data$eta_factor, levels = eta_levels),
      # This is the posterior support for the sign implied by the posterior mean.
      sign_support = if_else(
        .data$eta_mean >= 0,
        .data$support_positive,
        .data$support_negative
      )
    )
}

plot_eta_base <- function(df, eta_limit, bbox, title = NULL,
                          show_legend = FALSE, show_border = FALSE) {
  ggplot(df) +
    geom_sf(
      aes(fill = .data$eta_mean, alpha = .data$sign_support),
      colour = NA
    ) +
    eta_fill_scale(eta_limit) +
    eta_alpha_scale() +
    coord_sf(
      xlim = c(bbox[["xmin"]], bbox[["xmax"]]),
      ylim = c(bbox[["ymin"]], bbox[["ymax"]]),
      expand = FALSE
    ) +
    guides(
      fill = guide_colourbar(order = 1, barwidth = grid::unit(5, "cm")),
      alpha = guide_legend(order = 2, override.aes = list(fill = "grey35"))
    ) +
    labs(title = title) +
    eta_map_theme(if (show_legend) "bottom" else "none") +
    theme(
      panel.border = if (show_border) {
        element_rect(colour = "grey35", fill = NA, linewidth = 0.45)
      } else {
        element_blank()
      }
    )
}

plot_eta_inset_map <- function(df, eta_factor, eta_limit) {
  eta_df <- df |> filter(.data$eta_factor == .env$eta_factor)

  p_main <- plot_eta_base(
    eta_df,
    eta_limit = eta_limit,
    bbox = mainland_bbox,
    title = as.character(eta_factor)
  )

  p_inset <- plot_eta_base(
    eta_df,
    eta_limit = eta_limit,
    bbox = bornholm_bbox,
    show_border = TRUE
  ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.55)
    )

  ggdraw(p_main) +
    draw_plot(
      p_inset,
      x = 0.64,
      y = 0.58,
      width = 0.28,
      height = 0.28
    )
}

make_atlas_panel <- function(df, atlas, period) {
  eta_limit <- max(abs(df$eta_mean), na.rm = TRUE)
  eta_limit <- ceiling(eta_limit * 10) / 10

  factor_maps <- map(
    eta_levels,
    \(eta_factor) plot_eta_inset_map(df, eta_factor, eta_limit)
  )

  legend_plot <- plot_eta_base(
    df |> filter(.data$eta_factor == eta_levels[[1]]),
    eta_limit = eta_limit,
    bbox = mainland_bbox,
    show_legend = TRUE
  )
  shared_legend <- get_legend(legend_plot)

  title <- ggdraw() +
    draw_label(
      paste0("Atlas ", atlas, " (", period, ") Eta spatial loadings"),
      fontface = "bold",
      size = 15,
      x = 0.5,
      hjust = 0.5
    )

  subtitle <- ggdraw() +
    draw_label(
      "Within-model latent spatial gradients; transparency shows posterior support for each cell's sign.",
      size = 10,
      x = 0.5,
      hjust = 0.5
    )

  plot_grid(
    title,
    subtitle,
    plot_grid(plotlist = factor_maps, nrow = 1, align = "h"),
    shared_legend,
    ncol = 1,
    rel_heights = c(0.10, 0.07, 1, 0.17)
  )
}

grid_sf <- st_read(grid_path, quiet = TRUE) |>
  select(site = kvdrtkd, geometry)

message("Loaded grid cells: ", nrow(grid_sf))

walk(seq_len(nrow(model_specs)), function(i) {
  spec <- model_specs[i, ]
  model_file <- fitted_file(spec$model_folder)

  message("Loading Atlas ", spec$atlas, " model: ", model_file)
  model <- load_fitted_hmsc_model(model_file)

  eta_table <- extract_eta_table(
    model = model,
    atlas = spec$atlas,
    period = spec$period,
    expected_sites = spec$expected_sites
  )

  eta_sf <- grid_sf |>
    inner_join(eta_table, by = "site")

  joined_sites <- n_distinct(eta_sf$site)
  if (joined_sites != spec$expected_sites) {
    stop(
      "Atlas ",
      spec$atlas,
      " expected ",
      spec$expected_sites,
      " joined grid cells, found ",
      joined_sites,
      ".",
      call. = FALSE
    )
  }

  message(
    "Atlas ",
    spec$atlas,
    ": joined ",
    joined_sites,
    " sites and ",
    nlevels(eta_sf$eta_factor),
    " Eta factors."
  )

  atlas_panel <- make_atlas_panel(
    eta_sf,
    atlas = spec$atlas,
    period = spec$period
  )

  out_file <- file.path(
    plot_dir,
    paste0(output_prefix, "-eta-spatial-loadings-atlas-", spec$atlas, ".png")
  )

  ggsave(
    filename = out_file,
    plot = atlas_panel,
    width = 16,
    height = 5.2,
    dpi = 300,
    bg = "white"
  )

  message("Saved: ", out_file)
})

message("Done. Eta spatial loading PNGs are in: ", plot_dir)
