# Purpose:
# Make a forest-plot style summary of the multivariate GWR environmental-layer
# effects for both predicted richness and all-site beta turnover. This is the
# polished version of the "Draft 2" exploratory plot: thin intervals show the
# 5-95% spread of local t-values across grid cells, thick intervals show the
# 25-75% spread, and points show the median. The script reads the cached Fig. 10
# local coefficient table and does not refit any GWR models.

rm(list = ls())

pattern <- "2026-03-13"
bandwidth_label <- "fixed_100km"
significance_t <- 1.96
local_t_xlim <- c(-8, 8)

gwr_out_dir <- file.path("notebooks", "exploratory", "outputs", "gwr")
figure_out_dir <- file.path("misc-figures", "outputs", "main")

gwr_local_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-local-coefficients-100km.csv"))
forest_summary_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-multivariate-effect-forest-interval-summary-100km.csv"))

richness_png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10e1-gwr-richness-effect-forest-intervals.png"))
turnover_png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10e2-gwr-turnover-effect-forest-intervals.png"))
combined_png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10e-gwr-richness-turnover-effect-forest-intervals.png"))

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

# Figure 2 variable-importance Sankey/bump palette.
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

period_alphas <- c("1970s" = 0.35, "1990s" = 0.62, "2010s" = 0.95)
period_shapes <- c("1970s" = 21, "1990s" = 22, "2010s" = 24)

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

required_columns <- c("response_label", "period", "driver_label", "coefficient", "local_t")
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

forest_summary <- gwr_local |>
  group_by(.data$response_label, .data$period, .data$driver_label) |>
  summarise(
    n_cells = n(),
    n_supported_negative = sum(.data$local_t <= -significance_t, na.rm = TRUE),
    n_supported_positive = sum(.data$local_t >= significance_t, na.rm = TRUE),
    percent_supported_negative = 100 * .data$n_supported_negative / .data$n_cells,
    percent_supported_positive = 100 * .data$n_supported_positive / .data$n_cells,
    local_t_q05 = safe_quantile(.data$local_t, 0.05),
    local_t_q25 = safe_quantile(.data$local_t, 0.25),
    local_t_median = median(.data$local_t, na.rm = TRUE),
    local_t_q75 = safe_quantile(.data$local_t, 0.75),
    local_t_q95 = safe_quantile(.data$local_t, 0.95),
    coefficient_q05 = safe_quantile(.data$coefficient, 0.05),
    coefficient_median = median(.data$coefficient, na.rm = TRUE),
    coefficient_q95 = safe_quantile(.data$coefficient, 0.95),
    .groups = "drop"
  ) |>
  left_join(driver_y, by = "driver_label") |>
  mutate(
    y_period = .data$base_y + recode(
      as.character(.data$period),
      "1970s" = -0.22,
      "1990s" = 0,
      "2010s" = 0.22
    )
  )

readr::write_csv(forest_summary, forest_summary_path)

effect_theme <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9.5, colour = "grey30"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

plot_forest_interval <- function(df, title, subtitle = NULL, facet_response = FALSE) {
  p <- ggplot(df) +
    geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.35) +
    geom_vline(xintercept = c(-significance_t, significance_t), colour = "grey55", linetype = "22", linewidth = 0.35) +
    geom_segment(
      aes(
        x = .data$local_t_q05,
        xend = .data$local_t_q95,
        y = .data$y_period,
        yend = .data$y_period,
        colour = .data$driver_label,
        alpha = .data$period
      ),
      linewidth = 0.8,
      lineend = "round"
    ) +
    geom_segment(
      aes(
        x = .data$local_t_q25,
        xend = .data$local_t_q75,
        y = .data$y_period,
        yend = .data$y_period,
        colour = .data$driver_label,
        alpha = .data$period
      ),
      linewidth = 2.7,
      lineend = "round"
    ) +
    geom_point(
      aes(
        x = .data$local_t_median,
        y = .data$y_period,
        fill = .data$driver_label,
        shape = .data$period,
        alpha = .data$period
      ),
      colour = "grey15",
      size = 2.4,
      stroke = 0.25
    ) +
    scale_colour_manual(values = variable_colours, guide = "none") +
    scale_fill_manual(values = variable_colours, name = "Environmental layer", drop = FALSE) +
    scale_alpha_manual(values = period_alphas, name = "Atlas period", drop = FALSE) +
    scale_shape_manual(values = period_shapes, name = "Atlas period", drop = FALSE) +
    scale_x_continuous() +
    scale_y_continuous(
      breaks = driver_y$base_y,
      labels = as.character(driver_y$driver_label),
      expand = expansion(mult = c(0.06, 0.1))
    ) +
    coord_cartesian(xlim = local_t_xlim) +
    guides(
      fill = guide_legend(override.aes = list(alpha = 0.95, shape = 21), nrow = 2, byrow = TRUE),
      alpha = "none",
      shape = guide_legend(override.aes = list(fill = "grey45"))
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Local t-value from multivariate GWR",
      y = NULL
    ) +
    effect_theme()

  if (facet_response) {
    p <- p +
      facet_wrap(~ response_label, ncol = 1) +
      theme(strip.text = element_text(face = "bold", size = 10))
  }

  p
}

common_subtitle <- "Thin intervals are 5-95%; thick intervals are 25-75%; points are medians."

richness_plot <- forest_summary |>
  filter(.data$response_label == "Predicted richness") |>
  plot_forest_interval(
    title = "Richness GWR environmental effects",
    subtitle = common_subtitle
  )

turnover_plot <- forest_summary |>
  filter(.data$response_label == "All-site beta turnover") |>
  plot_forest_interval(
    title = "Turnover GWR environmental effects",
    subtitle = common_subtitle
  )

combined_plot <- forest_summary |>
  plot_forest_interval(
    title = "Richness and turnover GWR environmental effects",
    subtitle = common_subtitle,
    facet_response = TRUE
  )

ggsave(richness_png_path, richness_plot, width = 10.5, height = 6.5, units = "in", dpi = 300, bg = "white")
ggsave(turnover_png_path, turnover_plot, width = 10.5, height = 6.5, units = "in", dpi = 300, bg = "white")
ggsave(combined_png_path, combined_plot, width = 10.5, height = 10.2, units = "in", dpi = 300, bg = "white")

message("Saved: ", forest_summary_path)
message("Saved: ", richness_png_path)
message("Saved: ", turnover_png_path)
message("Saved: ", combined_png_path)
