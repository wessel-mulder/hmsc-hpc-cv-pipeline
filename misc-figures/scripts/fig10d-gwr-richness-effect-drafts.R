# Purpose:
# Make quick draft alternatives for visualising the environmental-layer effects
# from the multivariate richness GWR. These are intentionally exploratory PNGs:
# they all read the cached Fig. 10 local coefficient table and do not refit GWR
# models. The local t-value is used as the common effect scale because it is
# directly comparable across scaled environmental predictors and matches the
# support rule used in the maps.

rm(list = ls())

pattern <- "2026-03-13"
bandwidth_label <- "fixed_100km"
significance_t <- 1.96
local_t_xlim <- c(-8, 8)
ridge_height <- 0.72

gwr_out_dir <- file.path("notebooks", "exploratory", "outputs", "gwr")
figure_out_dir <- file.path("misc-figures", "outputs", "main")

gwr_local_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-local-coefficients-100km.csv"))
draft_summary_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-richness-multivariate-effect-draft-summary-100km.csv"))
draft_density_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-richness-multivariate-effect-draft-density-100km.csv"))

signed_support_png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10d1-gwr-richness-signed-support-draft.png"))
forest_interval_png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10d2-gwr-richness-forest-interval-draft.png"))
support_heatmap_png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10d3-gwr-richness-support-heatmap-draft.png"))
filled_ridges_png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10d4-gwr-richness-filled-ridges-draft.png"))
dumbbell_png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10d5-gwr-richness-support-dumbbell-draft.png"))

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

richness_local <- readr::read_csv(
  gwr_local_path,
  col_types = cols(period = col_character(), .default = col_guess()),
  show_col_types = FALSE
) |>
  filter(
    .data$bandwidth == bandwidth_label,
    .data$response_label == "Predicted richness"
  ) |>
  mutate(
    period = factor(.data$period, levels = period_levels),
    driver_label = factor(.data$driver_label, levels = driver_levels)
  )

required_columns <- c("period", "driver_label", "coefficient", "local_t")
missing_columns <- setdiff(required_columns, colnames(richness_local))

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

summary_df <- richness_local |>
  group_by(.data$period, .data$driver_label) |>
  summarise(
    n_cells = n(),
    n_supported_negative = sum(.data$local_t <= -significance_t, na.rm = TRUE),
    n_supported_positive = sum(.data$local_t >= significance_t, na.rm = TRUE),
    n_not_supported = sum(abs(.data$local_t) < significance_t, na.rm = TRUE),
    percent_supported_negative = 100 * .data$n_supported_negative / .data$n_cells,
    percent_supported_positive = 100 * .data$n_supported_positive / .data$n_cells,
    percent_not_supported = 100 * .data$n_not_supported / .data$n_cells,
    net_supported_percent = .data$percent_supported_positive - .data$percent_supported_negative,
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
  left_join(driver_y, by = "driver_label")

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

density_df <- richness_local |>
  group_by(.data$period, .data$driver_label) |>
  group_modify(~ make_density(.x$local_t)) |>
  ungroup() |>
  left_join(driver_y, by = "driver_label") |>
  group_by(.data$period, .data$driver_label) |>
  mutate(
    max_density = max(.data$density, na.rm = TRUE),
    density_scaled = if_else(.data$max_density > 0, .data$density / .data$max_density, 0),
    ridge_y = .data$base_y + ridge_height * .data$density_scaled,
    effect_sign = if_else(.data$local_t_plot < 0, "Negative", "Positive")
  ) |>
  ungroup()

readr::write_csv(summary_df, draft_summary_path)
readr::write_csv(density_df, draft_density_path)

base_effect_theme <- function(base_size = 10) {
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

plot_signed_support <- function() {
  support_long <- summary_df |>
    select(
      period,
      driver_label,
      base_y,
      percent_supported_negative,
      percent_supported_positive
    ) |>
    pivot_longer(
      cols = c(percent_supported_negative, percent_supported_positive),
      names_to = "direction",
      values_to = "percent_supported"
    ) |>
    mutate(
      signed_percent = if_else(.data$direction == "percent_supported_negative", -.data$percent_supported, .data$percent_supported),
      direction = if_else(.data$direction == "percent_supported_negative", "Negative", "Positive"),
      y_period = .data$base_y + recode(as.character(.data$period), "1970s" = -0.22, "1990s" = 0, "2010s" = 0.22)
    )

  ggplot(support_long) +
    geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.35) +
    geom_segment(
      aes(
        x = 0,
        xend = .data$signed_percent,
        y = .data$y_period,
        yend = .data$y_period,
        colour = .data$driver_label,
        alpha = .data$period
      ),
      linewidth = 2.4,
      lineend = "round"
    ) +
    geom_point(
      aes(
        x = .data$signed_percent,
        y = .data$y_period,
        fill = .data$driver_label,
        shape = .data$period,
        alpha = .data$period
      ),
      colour = "grey20",
      size = 2.7,
      stroke = 0.25
    ) +
    scale_colour_manual(values = variable_colours, guide = "none") +
    scale_fill_manual(values = variable_colours, name = "Environmental layer", drop = FALSE) +
    scale_alpha_manual(values = period_alphas, name = "Atlas period", drop = FALSE) +
    scale_shape_manual(values = period_shapes, name = "Atlas period", drop = FALSE) +
    scale_x_continuous(labels = label_percent(scale = 1)) +
    scale_y_continuous(
      breaks = driver_y$base_y,
      labels = as.character(driver_y$driver_label),
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    guides(
      fill = guide_legend(override.aes = list(alpha = 0.95, shape = 21), nrow = 2, byrow = TRUE),
      alpha = "none",
      shape = guide_legend(override.aes = list(fill = "grey45"))
    ) +
    labs(
      title = "Draft 1: signed support",
      subtitle = "Percent of grid cells with supported negative effects to the left and positive effects to the right.",
      x = "Supported grid cells",
      y = NULL
    ) +
    base_effect_theme()
}

plot_forest_interval <- function() {
  interval_df <- summary_df |>
    mutate(
      y_period = .data$base_y + recode(as.character(.data$period), "1970s" = -0.22, "1990s" = 0, "2010s" = 0.22)
    )

  ggplot(interval_df) +
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
    coord_cartesian(xlim = local_t_xlim) +
    scale_y_continuous(
      breaks = driver_y$base_y,
      labels = as.character(driver_y$driver_label),
      expand = expansion(mult = c(0.05, 0.08))
    ) +
    guides(
      fill = guide_legend(override.aes = list(alpha = 0.95, shape = 21), nrow = 2, byrow = TRUE),
      alpha = "none",
      shape = guide_legend(override.aes = list(fill = "grey45"))
    ) +
    labs(
      title = "Draft 2: forest-style distribution intervals",
      subtitle = "Thin intervals are 5-95%; thick intervals are 25-75%; points are medians.",
      x = "Local t-value from multivariate richness GWR",
      y = NULL
    ) +
    base_effect_theme()
}

plot_support_heatmap <- function() {
  heatmap_df <- summary_df |>
    mutate(
      label = paste0(
        "+", comma(.data$n_supported_positive),
        "\n",
        if_else(.data$n_supported_negative > 0, paste0("-", comma(.data$n_supported_negative)), "0")
      )
    )

  ggplot(heatmap_df, aes(x = .data$period, y = .data$driver_label)) +
    geom_tile(aes(fill = .data$net_supported_percent), colour = "white", linewidth = 0.7) +
    geom_text(aes(label = .data$label), size = 3, lineheight = 0.85, fontface = "bold") +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      name = "Positive minus\nnegative support (%)"
    ) +
    labs(
      title = "Draft 3: net-support heatmap",
      subtitle = "Cell fill is positive minus negative supported grid-cell percentage; text gives +cells / -cells.",
      x = NULL,
      y = NULL
    ) +
    base_effect_theme() +
    scale_y_discrete(limits = rev(driver_levels)) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(size = 9),
      legend.position = "right"
    )
}

plot_filled_ridges <- function() {
  ggplot() +
    geom_segment(
      data = summary_df,
      aes(x = local_t_xlim[1], xend = local_t_xlim[2], y = .data$base_y, yend = .data$base_y),
      colour = "grey86",
      linewidth = 0.25
    ) +
    geom_ribbon(
      data = filter(density_df, .data$local_t_plot <= 0),
      aes(
        x = .data$local_t_plot,
        ymin = .data$base_y,
        ymax = .data$ridge_y,
        group = interaction(.data$period, .data$driver_label)
      ),
      fill = "#2166ac",
      alpha = 0.58
    ) +
    geom_ribbon(
      data = filter(density_df, .data$local_t_plot >= 0),
      aes(
        x = .data$local_t_plot,
        ymin = .data$base_y,
        ymax = .data$ridge_y,
        group = interaction(.data$period, .data$driver_label)
      ),
      fill = "#b2182b",
      alpha = 0.58
    ) +
    geom_line(
      data = density_df,
      aes(
        x = .data$local_t_plot,
        y = .data$ridge_y,
        group = interaction(.data$period, .data$driver_label)
      ),
      colour = "grey20",
      linewidth = 0.25
    ) +
    geom_vline(xintercept = 0, colour = "grey22", linewidth = 0.35) +
    geom_vline(xintercept = c(-significance_t, significance_t), colour = "grey50", linetype = "22", linewidth = 0.35) +
    facet_wrap(~ period, nrow = 1) +
    scale_x_continuous(limits = local_t_xlim) +
    scale_y_continuous(
      breaks = driver_y$base_y,
      labels = as.character(driver_y$driver_label),
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    labs(
      title = "Draft 4: small-multiple filled ridges",
      subtitle = "Richness-only version of the filled ridge plot; periods are separated rather than overlaid.",
      x = "Local t-value from multivariate richness GWR",
      y = NULL
    ) +
    base_effect_theme() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )
}

plot_dumbbell <- function() {
  dumbbell_df <- summary_df |>
    select(
      period,
      driver_label,
      percent_supported_negative,
      percent_supported_positive
    ) |>
    pivot_longer(
      cols = c(percent_supported_negative, percent_supported_positive),
      names_to = "direction",
      values_to = "percent_supported"
    ) |>
    mutate(
      direction = recode(
        .data$direction,
        percent_supported_negative = "Negative support",
        percent_supported_positive = "Positive support"
      ),
      direction = factor(.data$direction, levels = c("Negative support", "Positive support"))
    )

  ggplot(dumbbell_df, aes(x = .data$period, y = .data$percent_supported, group = .data$driver_label)) +
    geom_line(aes(colour = .data$driver_label), linewidth = 0.75, alpha = 0.55) +
    geom_point(aes(fill = .data$driver_label), shape = 21, colour = "grey15", size = 2.5, stroke = 0.25) +
    facet_wrap(~ direction, nrow = 1) +
    scale_colour_manual(values = variable_colours, guide = "none") +
    scale_fill_manual(values = variable_colours, name = "Environmental layer", drop = FALSE) +
    scale_y_continuous(labels = label_percent(scale = 1)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(
      title = "Draft 5: support through time",
      subtitle = "Each line follows one environmental layer across atlas periods.",
      x = NULL,
      y = "Supported grid cells"
    ) +
    base_effect_theme() +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold")
    )
}

ggsave(signed_support_png_path, plot_signed_support(), width = 10.5, height = 6.5, units = "in", dpi = 300, bg = "white")
ggsave(forest_interval_png_path, plot_forest_interval(), width = 10.5, height = 6.5, units = "in", dpi = 300, bg = "white")
ggsave(support_heatmap_png_path, plot_support_heatmap(), width = 7.4, height = 5.8, units = "in", dpi = 300, bg = "white")
ggsave(filled_ridges_png_path, plot_filled_ridges(), width = 12.5, height = 6.2, units = "in", dpi = 300, bg = "white")
ggsave(dumbbell_png_path, plot_dumbbell(), width = 10.8, height = 5.8, units = "in", dpi = 300, bg = "white")

message("Saved: ", draft_summary_path)
message("Saved: ", draft_density_path)
message("Saved: ", signed_support_png_path)
message("Saved: ", forest_interval_png_path)
message("Saved: ", support_heatmap_png_path)
message("Saved: ", filled_ridges_png_path)
message("Saved: ", dumbbell_png_path)
