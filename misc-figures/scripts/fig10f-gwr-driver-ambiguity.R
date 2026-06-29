# Purpose:
# Diagnose whether the mapped dominant GWR driver is clearly separated from
# other locally supported predictors. This script reads the cached multivariate
# GWR local coefficient table, derives three ambiguity metrics, and maps them
# using the same mainland plus Bornholm inset layout as Fig. 10:
#   1. dominance margin: top absolute local t-value minus runner-up
#   2. number of predictors with |local t| >= 1.96
#   3. categorical single-driver versus multi-driver support
#
# The script does not refit any GWR models.

rm(list = ls())

pattern <- "2026-03-13"
bandwidth_label <- "fixed_100km"
significance_t <- 1.96

gwr_out_dir <- file.path("notebooks", "exploratory", "outputs", "gwr")
figure_out_dir <- file.path("misc-figures", "outputs", "main")

gwr_local_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-local-coefficients-100km.csv"))
ambiguity_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-driver-ambiguity-100km.csv"))
ambiguity_summary_path <- file.path(gwr_out_dir, paste0(pattern, "-gwr-driver-ambiguity-summary-100km.csv"))
margin_png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10f1-gwr-dominance-margin.png"))
count_png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10f2-gwr-supported-predictor-count.png"))
class_png_path <- file.path(figure_out_dir, paste0(pattern, "-fig10f3-gwr-driver-ambiguity-class.png"))

dir.create(gwr_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_out_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("tidyverse", "cowplot", "scales")
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
count_levels <- as.character(0:8)
class_levels <- c(
  "No supported predictors",
  "One supported predictor",
  "Two supported predictors",
  "Three supported predictors",
  "Four or more supported predictors"
)

count_colours <- c(
  "0" = "#d9d9d9",
  "1" = "#f7fbff",
  "2" = "#deebf7",
  "3" = "#c6dbef",
  "4" = "#9ecae1",
  "5" = "#6baed6",
  "6" = "#3182bd",
  "7" = "#08519c",
  "8" = "#08306b"
)

class_colours <- c(
  "No supported predictors" = "#d9d9d9",
  "One supported predictor" = "#2ca25f",
  "Two supported predictors" = "#fec44f",
  "Three supported predictors" = "#f03b20",
  "Four or more supported predictors" = "#7a0177"
)

gwr_local <- readr::read_csv(
  gwr_local_path,
  col_types = cols(period = col_character(), .default = col_guess()),
  show_col_types = FALSE
) |>
  filter(.data$bandwidth == bandwidth_label) |>
  mutate(
    response_label = factor(.data$response_label, levels = response_levels),
    period = factor(.data$period, levels = period_levels)
  )

required_columns <- c(
  "response", "response_label", "atlas", "period", "survey", "site", "X", "Y",
  "predictor", "driver_label", "local_t", "abs_local_t"
)
missing_columns <- setdiff(required_columns, colnames(gwr_local))

if (length(missing_columns) > 0) {
  stop(
    "The local coefficient cache is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

gwr_ambiguity <- gwr_local |>
  group_by(
    .data$response,
    .data$response_label,
    .data$atlas,
    .data$period,
    .data$survey,
    .data$site,
    .data$X,
    .data$Y
  ) |>
  arrange(desc(.data$abs_local_t), .by_group = TRUE) |>
  summarise(
    top_predictor = first(.data$predictor),
    top_driver_label = first(.data$driver_label),
    top_local_t = first(.data$local_t),
    top_abs_local_t = first(.data$abs_local_t),
    second_predictor = nth(.data$predictor, 2),
    second_driver_label = nth(.data$driver_label, 2),
    second_local_t = nth(.data$local_t, 2),
    second_abs_local_t = nth(.data$abs_local_t, 2),
    t_margin = .data$top_abs_local_t - .data$second_abs_local_t,
    n_supported_predictors = sum(.data$abs_local_t >= significance_t, na.rm = TRUE),
    supported_predictors = paste(.data$driver_label[.data$abs_local_t >= significance_t], collapse = "; "),
    .groups = "drop"
  ) |>
  mutate(
    top_supported = .data$top_abs_local_t >= significance_t,
    margin_for_plot = if_else(.data$top_supported, .data$t_margin, NA_real_),
    supported_count_label = factor(
      as.character(.data$n_supported_predictors),
      levels = count_levels
    ),
    ambiguity_class = case_when(
      .data$n_supported_predictors == 0 ~ "No supported predictors",
      .data$n_supported_predictors == 1 ~ "One supported predictor",
      .data$n_supported_predictors == 2 ~ "Two supported predictors",
      .data$n_supported_predictors == 3 ~ "Three supported predictors",
      .data$n_supported_predictors >= 4 ~ "Four or more supported predictors"
    ),
    ambiguity_class = factor(.data$ambiguity_class, levels = class_levels),
    response_label = factor(.data$response_label, levels = response_levels),
    period = factor(.data$period, levels = period_levels)
  )

ambiguity_summary <- gwr_ambiguity |>
  group_by(.data$response, .data$response_label, .data$period) |>
  summarise(
    n_cells = n(),
    n_no_supported = sum(.data$n_supported_predictors == 0, na.rm = TRUE),
    n_one_supported = sum(.data$n_supported_predictors == 1, na.rm = TRUE),
    n_multi_supported = sum(.data$n_supported_predictors >= 2, na.rm = TRUE),
    percent_no_supported = 100 * .data$n_no_supported / .data$n_cells,
    percent_one_supported = 100 * .data$n_one_supported / .data$n_cells,
    percent_multi_supported = 100 * .data$n_multi_supported / .data$n_cells,
    supported_count_mean = mean(.data$n_supported_predictors, na.rm = TRUE),
    supported_count_median = median(.data$n_supported_predictors, na.rm = TRUE),
    top_abs_t_median = median(.data$top_abs_local_t, na.rm = TRUE),
    margin_median_supported = median(.data$margin_for_plot, na.rm = TRUE),
    margin_q25_supported = quantile(.data$margin_for_plot, 0.25, na.rm = TRUE, names = FALSE),
    margin_q75_supported = quantile(.data$margin_for_plot, 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )

readr::write_csv(gwr_ambiguity, ambiguity_path)
readr::write_csv(ambiguity_summary, ambiguity_summary_path)

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

metric_coord <- function(bbox) {
  coord_fixed(
    xlim = c(bbox[["xmin"]], bbox[["xmax"]]),
    ylim = c(bbox[["ymin"]], bbox[["ymax"]]),
    expand = FALSE
  )
}

margin_cap <- quantile(gwr_ambiguity$margin_for_plot, 0.95, na.rm = TRUE, names = FALSE)
margin_cap <- max(1, margin_cap)

margin_scale <- function() {
  scale_colour_gradientn(
    colours = c("#ffffcc", "#a1dab4", "#41b6c4", "#225ea8"),
    limits = c(0, margin_cap),
    oob = scales::squish,
    na.value = "#d9d9d9",
    name = "Top - second |t|"
  )
}

count_scale <- function() {
  scale_colour_manual(
    values = count_colours,
    limits = count_levels,
    drop = FALSE,
    name = "Supported predictors"
  )
}

class_scale <- function() {
  scale_colour_manual(
    values = class_colours,
    limits = class_levels,
    drop = FALSE,
    name = "Support class",
    guide = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3.2))
  )
}

plot_margin_base <- function(df, title = NULL, bbox = mainland_bbox,
                             show_legend = FALSE, border = FALSE,
                             point_size = 1.25) {
  ggplot(df) +
    geom_point(
      aes(x = .data$X, y = .data$Y, colour = .data$margin_for_plot),
      size = point_size,
      alpha = 0.95
    ) +
    margin_scale() +
    metric_coord(bbox) +
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

plot_count_base <- function(df, title = NULL, bbox = mainland_bbox,
                            show_legend = FALSE, border = FALSE,
                            point_size = 1.25) {
  ggplot(df) +
    geom_point(
      aes(x = .data$X, y = .data$Y, colour = .data$supported_count_label),
      size = point_size,
      alpha = 0.95
    ) +
    count_scale() +
    metric_coord(bbox) +
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

plot_class_base <- function(df, title = NULL, bbox = mainland_bbox,
                            show_legend = FALSE, border = FALSE,
                            point_size = 1.25) {
  ggplot(df) +
    geom_point(
      aes(x = .data$X, y = .data$Y, colour = .data$ambiguity_class),
      size = point_size,
      alpha = 0.95
    ) +
    class_scale() +
    metric_coord(bbox) +
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

plot_map_with_inset <- function(df, title, base_plotter) {
  p_main <- base_plotter(df, title = title, bbox = mainland_bbox)
  p_inset <- base_plotter(df, bbox = bornholm_bbox, border = TRUE, point_size = 1.45) +
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

response_map_row <- function(response_name, base_plotter) {
  response_maps <- map(period_levels, function(period_name) {
    row_df <- gwr_ambiguity |>
      filter(.data$response == response_name, .data$period == period_name)

    plot_map_with_inset(row_df, title = period_name, base_plotter = base_plotter)
  })

  plot_grid(plotlist = response_maps, nrow = 1, align = "h", rel_widths = c(1, 1, 1))
}

build_ambiguity_figure <- function(title_a, title_b, base_plotter, legend_plot) {
  figure_legend <- get_legend(
    legend_plot(gwr_ambiguity, show_legend = TRUE) +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        legend.margin = margin(0, 0, 0, 0)
      )
  )

  panel_a <- plot_grid(
    panel_title(title_a),
    response_map_row("predicted_richness", base_plotter),
    ncol = 1,
    rel_heights = c(0.08, 1)
  )

  panel_b <- plot_grid(
    panel_title(title_b),
    response_map_row("beta_turnover", base_plotter),
    ncol = 1,
    rel_heights = c(0.08, 1)
  )

  plot_grid(
    panel_a,
    panel_b,
    figure_legend,
    ncol = 1,
    rel_heights = c(1, 1, 0.22)
  ) +
    theme(plot.background = element_rect(fill = "white", colour = NA))
}

margin_figure <- build_ambiguity_figure(
  title_a = "A  Richness dominance margin",
  title_b = "B  Turnover dominance margin",
  base_plotter = plot_margin_base,
  legend_plot = plot_margin_base
)

count_figure <- build_ambiguity_figure(
  title_a = "A  Richness supported-predictor count",
  title_b = "B  Turnover supported-predictor count",
  base_plotter = plot_count_base,
  legend_plot = plot_count_base
)

class_figure <- build_ambiguity_figure(
  title_a = "A  Richness driver ambiguity",
  title_b = "B  Turnover driver ambiguity",
  base_plotter = plot_class_base,
  legend_plot = plot_class_base
)

ggsave(margin_png_path, margin_figure, width = 13.5, height = 11.6, units = "in", dpi = 300, bg = "white")
ggsave(count_png_path, count_figure, width = 13.5, height = 11.6, units = "in", dpi = 300, bg = "white")
ggsave(class_png_path, class_figure, width = 13.5, height = 11.6, units = "in", dpi = 300, bg = "white")

message("Saved: ", ambiguity_path)
message("Saved: ", ambiguity_summary_path)
message("Saved: ", margin_png_path)
message("Saved: ", count_png_path)
message("Saved: ", class_png_path)
