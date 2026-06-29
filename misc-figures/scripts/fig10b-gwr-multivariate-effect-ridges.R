# Purpose:
# Make a ridge-style summary of the multivariate GWR environmental effects.
# The script reads the cached local coefficient/t-value table from the Fig. 10
# multivariate GWR workflow, so it never refits the GWR models. Ridges show the
# distribution of local t-values across grid cells for each environmental layer,
# response, and atlas period. Side labels count supported negative and positive
# cells using the same |t| >= 1.96 rule as the maps.

rm(list = ls())

pattern <- "2026-03-13"
bandwidth_label <- "fixed_100km"
significance_t <- 1.96
local_t_xlim <- c(-8, 8)
ridge_height <- 0.78

gwr_out_dir <- file.path("notebooks", "exploratory", "outputs", "gwr")
figure_out_dir <- file.path("misc-figures", "outputs", "main")

gwr_local_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-local-coefficients-100km.csv"))
ridge_summary_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-multivariate-effect-ridge-summary-100km.csv"))
ridge_density_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-multivariate-effect-ridge-density-100km.csv"))
richness_line_summary_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-richness-multivariate-effect-ridge-line-summary-100km.csv"))
richness_line_density_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-richness-multivariate-effect-ridge-line-density-100km.csv"))
ridge_png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10b-gwr-multivariate-effect-ridges.png"))
richness_line_png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10c-gwr-richness-multivariate-effect-ridge-lines.png"))

dir.create(gwr_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_out_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("tidyverse", "scales")
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
  library(scales)
})

if (!file.exists(gwr_local_path)) {
  stop(
    "The multivariate GWR local coefficient cache is missing. Run ",
    "`Rscript misc-figures/scripts/fig10-gwr-richness-turnover.R` first.",
    call. = FALSE
  )
}

period_levels <- c("1970s", "1990s", "2010s")
response_levels <- c("Predicted richness", "All-site beta turnover")
driver_levels <- c(
  "Temperature",
  "Precipitation",
  "Land-use heterogeneity",
  "Urban (% coverage)",
  "Cropland (% coverage)",
  "Pasture (% coverage)",
  "Forest (% coverage)",
  "Grass/Shrubland (% coverage)"
)

effect_colours <- c(
  "Negative local effect" = "#2166ac",
  "Positive local effect" = "#b2182b"
)

# This is the Figure 2 variable-importance Sankey/bump palette, reused here so
# each environmental layer keeps the same colour identity across figures.
variable_colours <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban (% coverage)" = "snow3",
  "Cropland (% coverage)" = "goldenrod1",
  "Pasture (% coverage)" = "darkorange",
  "Forest (% coverage)" = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2"
)

period_alphas <- c("1970s" = 0.34, "1990s" = 0.58, "2010s" = 0.92)
period_linewidths <- c("1970s" = 0.55, "1990s" = 0.72, "2010s" = 0.9)

gwr_local <- readr::read_csv(
  gwr_local_path,
  col_types = cols(period = col_character(), .default = col_guess()),
  show_col_types = FALSE
) |>
  filter(.data$bandwidth == bandwidth_label) |>
  mutate(
    response_label = factor(.data$response_label, levels = response_levels),
    period = factor(.data$period, levels = period_levels),
    driver_label = factor(.data$driver_label, levels = driver_levels)
  )

required_columns <- c(
  "response_label", "period", "driver_label", "coefficient", "local_t"
)
missing_columns <- setdiff(required_columns, colnames(gwr_local))

if (length(missing_columns) > 0) {
  stop(
    "The local coefficient cache is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

driver_y <- tibble(
  driver_label = factor(driver_levels, levels = driver_levels),
  base_y = length(driver_levels):1
)

safe_quantile <- function(x, prob) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  as.numeric(quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}

ridge_summary <- gwr_local |>
  group_by(.data$response_label, .data$period, .data$driver_label) |>
  summarise(
    n_cells = n(),
    n_supported_negative = sum(.data$local_t <= -significance_t, na.rm = TRUE),
    n_supported_positive = sum(.data$local_t >= significance_t, na.rm = TRUE),
    n_not_supported = sum(abs(.data$local_t) < significance_t, na.rm = TRUE),
    percent_supported_negative = 100 * .data$n_supported_negative / .data$n_cells,
    percent_supported_positive = 100 * .data$n_supported_positive / .data$n_cells,
    percent_not_supported = 100 * .data$n_not_supported / .data$n_cells,
    coefficient_q05 = safe_quantile(.data$coefficient, 0.05),
    coefficient_median = median(.data$coefficient, na.rm = TRUE),
    coefficient_q95 = safe_quantile(.data$coefficient, 0.95),
    local_t_q05 = safe_quantile(.data$local_t, 0.05),
    local_t_median = median(.data$local_t, na.rm = TRUE),
    local_t_q95 = safe_quantile(.data$local_t, 0.95),
    .groups = "drop"
  ) |>
  left_join(driver_y, by = "driver_label") |>
  mutate(
    negative_count_label = if_else(
      .data$n_supported_negative > 0,
      paste0("-", comma(.data$n_supported_negative)),
      ""
    ),
    positive_count_label = if_else(
      .data$n_supported_positive > 0,
      paste0("+", comma(.data$n_supported_positive)),
      ""
    )
  )

make_density <- function(x, n = 320) {
  x <- x[is.finite(x)]

  if (length(x) < 2 || length(unique(x)) < 2) {
    return(tibble(local_t_plot = seq(local_t_xlim[1], local_t_xlim[2], length.out = n), density = 0))
  }

  density_estimate <- density(
    x,
    from = local_t_xlim[1],
    to = local_t_xlim[2],
    n = n,
    adjust = 1.05,
    na.rm = TRUE
  )

  tibble(local_t_plot = density_estimate$x, density = density_estimate$y)
}

ridge_density <- gwr_local |>
  group_by(.data$response_label, .data$period, .data$driver_label) |>
  group_modify(~ make_density(.x$local_t)) |>
  ungroup() |>
  left_join(driver_y, by = "driver_label") |>
  group_by(.data$response_label, .data$period, .data$driver_label) |>
  mutate(
    max_density = max(.data$density, na.rm = TRUE),
    density_scaled = if_else(.data$max_density > 0, .data$density / .data$max_density, 0),
    ridge_y = .data$base_y + ridge_height * .data$density_scaled,
    effect_sign = if_else(.data$local_t_plot < 0, "Negative local effect", "Positive local effect"),
    effect_sign = factor(.data$effect_sign, levels = names(effect_colours))
  ) |>
  ungroup()

readr::write_csv(ridge_summary, ridge_summary_path)
readr::write_csv(ridge_density, ridge_density_path)

richness_line_summary <- ridge_summary |>
  filter(.data$response_label == "Predicted richness")

richness_line_density <- ridge_density |>
  filter(.data$response_label == "Predicted richness") |>
  mutate(
    period = factor(.data$period, levels = period_levels),
    driver_label = factor(.data$driver_label, levels = driver_levels)
  )

readr::write_csv(richness_line_summary, richness_line_summary_path)
readr::write_csv(richness_line_density, richness_line_density_path)

plot_ridge <- function() {
  ggplot() +
    geom_segment(
      data = ridge_summary,
      aes(
        x = local_t_xlim[1],
        xend = local_t_xlim[2],
        y = .data$base_y,
        yend = .data$base_y
      ),
      colour = "grey86",
      linewidth = 0.25
    ) +
    geom_ribbon(
      data = filter(ridge_density, .data$local_t_plot <= 0),
      aes(
        x = .data$local_t_plot,
        ymin = .data$base_y,
        ymax = .data$ridge_y,
        group = interaction(.data$response_label, .data$period, .data$driver_label, .data$effect_sign),
        fill = .data$effect_sign
      ),
      alpha = 0.62,
      colour = NA
    ) +
    geom_ribbon(
      data = filter(ridge_density, .data$local_t_plot >= 0),
      aes(
        x = .data$local_t_plot,
        ymin = .data$base_y,
        ymax = .data$ridge_y,
        group = interaction(.data$response_label, .data$period, .data$driver_label, .data$effect_sign),
        fill = .data$effect_sign
      ),
      alpha = 0.62,
      colour = NA
    ) +
    geom_line(
      data = ridge_density,
      aes(
        x = .data$local_t_plot,
        y = .data$ridge_y,
        group = interaction(.data$response_label, .data$period, .data$driver_label)
      ),
      colour = "grey18",
      linewidth = 0.22
    ) +
    geom_vline(xintercept = 0, colour = "grey20", linewidth = 0.35) +
    geom_vline(xintercept = c(-significance_t, significance_t), colour = "grey45", linetype = "22", linewidth = 0.35) +
    geom_text(
      data = ridge_summary,
      aes(x = local_t_xlim[1] + 0.22, y = .data$base_y + 0.52, label = .data$negative_count_label),
      colour = effect_colours[["Negative local effect"]],
      hjust = 0,
      size = 2.25,
      fontface = "bold"
    ) +
    geom_text(
      data = ridge_summary,
      aes(x = local_t_xlim[2] - 0.22, y = .data$base_y + 0.52, label = .data$positive_count_label),
      colour = effect_colours[["Positive local effect"]],
      hjust = 1,
      size = 2.25,
      fontface = "bold"
    ) +
    facet_grid(.data$response_label ~ .data$period) +
    scale_fill_manual(values = effect_colours, name = NULL) +
    scale_x_continuous(
      limits = local_t_xlim,
      breaks = c(-8, -4, -significance_t, 0, significance_t, 4, 8),
      labels = c("-8", "-4", "-1.96", "0", "1.96", "4", "8"),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = driver_y$base_y,
      labels = as.character(driver_y$driver_label),
      limits = c(0.65, length(driver_levels) + 1.05),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = "Multivariate GWR environmental effects across grid cells",
      subtitle = "Ridges are local t-value distributions; labels count supported cells at t <= -1.96 and t >= 1.96.",
      x = "Local t-value from multivariate GWR",
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 8.5),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9.5, colour = "grey30"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

plot_richness_line_ridge <- function() {
  ggplot() +
    geom_segment(
      data = richness_line_summary,
      aes(
        x = local_t_xlim[1],
        xend = local_t_xlim[2],
        y = .data$base_y,
        yend = .data$base_y
      ),
      colour = "grey88",
      linewidth = 0.25
    ) +
    geom_vline(xintercept = 0, colour = "grey22", linewidth = 0.35) +
    geom_vline(xintercept = c(-significance_t, significance_t), colour = "grey50", linetype = "22", linewidth = 0.35) +
    geom_line(
      data = richness_line_density,
      aes(
        x = .data$local_t_plot,
        y = .data$ridge_y,
        colour = .data$driver_label,
        alpha = .data$period,
        linewidth = .data$period,
        group = interaction(.data$period, .data$driver_label)
      ),
      lineend = "round"
    ) +
    scale_colour_manual(values = variable_colours, name = "Environmental layer", drop = FALSE) +
    scale_alpha_manual(values = period_alphas, name = "Atlas period", drop = FALSE) +
    scale_linewidth_manual(values = period_linewidths, name = "Atlas period", drop = FALSE) +
    scale_x_continuous(
      limits = local_t_xlim,
      breaks = c(-8, -4, -significance_t, 0, significance_t, 4, 8),
      labels = c("-8", "-4", "-1.96", "0", "1.96", "4", "8"),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = driver_y$base_y,
      labels = as.character(driver_y$driver_label),
      limits = c(0.65, length(driver_levels) + 1.08),
      expand = expansion(mult = c(0, 0))
    ) +
    guides(
      colour = guide_legend(override.aes = list(alpha = 0.95, linewidth = 1.1), nrow = 2, byrow = TRUE),
      alpha = guide_legend(override.aes = list(colour = "grey20")),
      linewidth = "none"
    ) +
    labs(
      title = "Richness GWR environmental effects through time",
      subtitle = "One ridge line per atlas period; variables use the Figure 2 palette and period is shown with alpha.",
      x = "Local t-value from multivariate richness GWR",
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 8),
      legend.title = element_text(face = "bold", size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey91", linewidth = 0.25),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9.5, colour = "grey30"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

ridge_plot <- plot_ridge()
richness_line_plot <- plot_richness_line_ridge()

ggsave(
  ridge_png_path,
  ridge_plot,
  width = 13.5,
  height = 8.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  richness_line_png_path,
  richness_line_plot,
  width = 11.2,
  height = 6.4,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Saved: ", ridge_summary_path)
message("Saved: ", ridge_density_path)
message("Saved: ", richness_line_summary_path)
message("Saved: ", richness_line_density_path)
message("Saved: ", ridge_png_path)
message("Saved: ", richness_line_png_path)
