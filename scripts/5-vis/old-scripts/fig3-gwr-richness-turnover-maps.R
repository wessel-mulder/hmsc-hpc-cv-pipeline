#!/usr/bin/env Rscript

# Main Figure 3 GWR figure:
# Dominant environmental associations for predicted richness and all-site beta
# turnover across the three atlas periods, using SFig. 1-style map panels and
# Figure 1-style legend keys.

rm(list = ls())

required_packages <- c(
  "cowplot",
  "dplyr",
  "ggplot2",
  "readr",
  "sf",
  "tibble"
)
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
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(sf)
  library(tibble)
})

pattern <- "2026-03-13"

input_path <- file.path(
  "intermediate-outputs",
  "gwr-output",
  "fig10-gwr-richness-turnover",
  "gwr",
  paste0(pattern, "-gwr-local-coefficients-100km.csv")
)
grid_path <- file.path(
  "data-new",
  "derived",
  "atlas-grid",
  "dof-5km-grid",
  "DK5km_ED50grid_approx_kvadrkod_DOF.shp"
)

figure_dir <- file.path("manuscript-outputs", "figures", "main")
table_dir <- file.path(
  "intermediate-outputs",
  "gwr-output",
  "fig3-gwr-richness-turnover-maps"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

figure_path <- file.path(
  figure_dir,
  paste0(pattern, "-fig3-gwr-richness-turnover-maps.png")
)
map_values_path <- file.path(
  table_dir,
  paste0(pattern, "-fig3-gwr-map-values.csv")
)
map_summary_path <- file.path(
  table_dir,
  paste0(pattern, "-fig3-gwr-dominant-driver-counts.csv")
)

if (!file.exists(input_path)) {
  stop("Missing cached GWR local-coefficient table: ", input_path, call. = FALSE)
}
if (!file.exists(grid_path)) {
  stop("Missing atlas grid shapefile: ", grid_path, call. = FALSE)
}

period_levels <- c("1970s", "1990s", "2010s")
response_labels <- c(
  "Predicted richness" = "Species richness",
  "All-site beta turnover" = "Beta turnover"
)
response_levels <- unname(response_labels)

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

driver_colours <- c(
  fig1_variable_colours[names(fig1_variable_colours) != "Site-level random effect"],
  "No significant driver" = "white"
)
driver_levels <- names(driver_colours)
no_significant_driver <- "No significant driver"
no_significant_outline <- "black"
significance_t <- 1.96
direction_levels <- c("Positive", "Negative", "Not significant")
direction_shapes <- c("Positive" = 19, "Negative" = 4, "Not significant" = 19)
environmental_driver_levels <- setdiff(driver_levels, no_significant_driver)

local_coefficients <- read_csv(
  input_path,
  col_types = cols(
    response_label = col_character(),
    period = col_character(),
    .default = col_guess()
  )
) |>
  filter(.data$response_label %in% names(response_labels))

dominant_drivers <- local_coefficients |>
  group_by(.data$response, .data$response_label, .data$atlas, .data$period, .data$survey) |>
  slice_max(.data$abs_local_t, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    supported = .data$abs_local_t >= significance_t,
    dominant_driver = if_else(.data$supported, .data$driver_label, no_significant_driver),
    direction_plot = if_else(.data$supported, .data$direction, "Not significant")
) |>
  mutate(
    response_panel = unname(response_labels[.data$response_label]),
    response_panel = factor(.data$response_panel, levels = response_levels),
    period = factor(.data$period, levels = period_levels),
    dominant_driver = factor(.data$dominant_driver, levels = driver_levels)
  )

if (any(is.na(dominant_drivers$period))) {
  stop("Unexpected period labels in: ", input_path, call. = FALSE)
}
if (any(is.na(dominant_drivers$dominant_driver))) {
  unknown_drivers <- unique(dominant_drivers$dominant_driver[is.na(dominant_drivers$dominant_driver)])
  stop(
    "Unexpected dominant-driver labels in ",
    input_path,
    ": ",
    paste(unknown_drivers, collapse = ", "),
    call. = FALSE
  )
}

map_values <- dominant_drivers |>
  select(
    response,
    response_label,
    response_panel,
    atlas,
    period,
    survey,
    site,
    X,
    Y,
    dominant_driver,
    direction_plot
  ) |>
  arrange(.data$period, .data$response_panel, .data$site)

map_summary <- dominant_drivers |>
  count(.data$response_panel, .data$period, .data$dominant_driver, name = "n_grid_cells") |>
  group_by(.data$response_panel, .data$period) |>
  mutate(
    total_grid_cells = sum(.data$n_grid_cells),
    percent_grid_cells = .data$n_grid_cells / .data$total_grid_cells * 100
  ) |>
  ungroup() |>
  arrange(.data$response_panel, .data$period, .data$dominant_driver)

write_csv(map_values, map_values_path)
write_csv(map_summary, map_summary_path)

shape_sf <- st_read(grid_path, quiet = TRUE) |>
  st_transform(25832)

mainland_bbox <- st_bbox(
  c(xmin = 430000, xmax = 745000, ymin = 6040000, ymax = 6410000),
  crs = st_crs(25832)
)
bornholm_bbox <- st_bbox(
  c(xmin = 858000, xmax = 897000, ymin = 6107000, ymax = 6145000),
  crs = st_crs(25832)
)

mainland_width <- as.numeric(mainland_bbox["xmax"] - mainland_bbox["xmin"])
mainland_height <- as.numeric(mainland_bbox["ymax"] - mainland_bbox["ymin"])
bornholm_width <- as.numeric(bornholm_bbox["xmax"] - bornholm_bbox["xmin"])
bornholm_height <- as.numeric(bornholm_bbox["ymax"] - bornholm_bbox["ymin"])

inset_w <- bornholm_width / mainland_width
inset_h <- bornholm_height / mainland_height
inset_margin_x <- 0.2
inset_margin_y <- 0.2

map_theme <- function() {
  theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
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
      size = 12,
      colour = "black"
    ) +
    theme(
      plot.background = element_rect(fill = "grey92", colour = "white"),
      plot.margin = margin(0, 0, 1.2, 0)
    )
}

period_strip <- function(label) {
  ggdraw() +
    draw_label(
      label,
      x = 0.5,
      y = 0.5,
      hjust = 0.5,
      vjust = 0.5,
      angle = 90,
      fontface = "bold",
      size = 12,
      colour = "black"
    ) +
    theme(
      plot.background = element_rect(fill = "grey92", colour = "white"),
      plot.margin = margin(0, 1.2, 0, 0)
    )
}

map_cell_data <- function(df) {
  plot_data <- shape_sf |>
    left_join(
      df |>
        select(site, dominant_driver, direction_plot),
      by = c("kvadratkod" = "site")
    ) |>
    filter(!is.na(.data$dominant_driver))

  suppressWarnings(st_point_on_surface(plot_data))
}

plot_dominant_map_base <- function(df, bbox, border = FALSE, point_size = 1.75) {
  plot_data <- map_cell_data(df)
  no_driver_cells <- plot_data |>
    filter(as.character(.data$dominant_driver) == no_significant_driver)
  driver_cells <- plot_data |>
    filter(as.character(.data$dominant_driver) != no_significant_driver)

  ggplot() +
    geom_sf(
      data = no_driver_cells,
      shape = 21,
      fill = driver_colours[[no_significant_driver]],
      colour = no_significant_outline,
      size = point_size,
      alpha = 0.95,
      stroke = 0.52
    ) +
    geom_sf(
      data = driver_cells,
      aes(colour = .data$dominant_driver, shape = .data$direction_plot),
      size = point_size,
      alpha = 0.95,
      stroke = 0.72
    ) +
    scale_colour_manual(
      values = driver_colours,
      limits = driver_levels,
      drop = FALSE,
      guide = "none"
    ) +
    scale_shape_manual(
      values = direction_shapes,
      limits = direction_levels,
      breaks = c("Positive", "Negative"),
      drop = FALSE,
      guide = "none"
    ) +
    coord_sf(
      xlim = c(bbox["xmin"], bbox["xmax"]),
      ylim = c(bbox["ymin"], bbox["ymax"]),
      expand = FALSE
    ) +
    labs(x = NULL, y = NULL) +
    map_theme() +
    theme(
      panel.border = if (border) {
        element_rect(colour = "grey30", fill = NA, linewidth = 0.5)
      } else {
        element_blank()
      }
    )
}

plot_map_with_inset <- function(df) {
  p_main <- plot_dominant_map_base(df, mainland_bbox)
  p_inset <- plot_dominant_map_base(
    df,
    bornholm_bbox,
    border = TRUE,
    point_size = 2.05
  ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5),
      plot.margin = margin(0, 0, 0, 0)
    )

  ggdraw(p_main) +
    draw_plot(
      p_inset,
      x = 1 - inset_w - inset_margin_x,
      y = 1 - inset_h - inset_margin_y,
      width = inset_w,
      height = inset_h
    ) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

map_panel <- function(response_panel_name, period_name) {
  panel_df <- map_values |>
    filter(
      .data$response_panel == response_panel_name,
      .data$period == period_name
    )

  plot_map_with_inset(panel_df)
}

figure1_style_combined_legend_plot <- function() {
  variable_legend_data <- tibble(
    dominant_driver = factor(environmental_driver_levels, levels = environmental_driver_levels),
    x = seq_along(environmental_driver_levels),
    y = 1
  )
  direction_legend_data <- tibble(
    association = factor(
      c("Positive", "Negative", no_significant_driver),
      levels = c("Positive", "Negative", no_significant_driver)
    ),
    x = seq_len(3),
    y = 1
  )

  ggplot() +
    geom_col(
      data = variable_legend_data,
      aes(x = .data$x, y = .data$y, fill = .data$dominant_driver),
      colour = "white",
      linewidth = 0.1,
      width = 0.8
    ) +
    geom_point(
      data = direction_legend_data,
      aes(x = .data$x, y = .data$y, shape = .data$association),
      colour = "black",
      fill = "white",
      size = 2.7,
      stroke = 0.5
    ) +
    scale_fill_manual(
      values = driver_colours,
      limits = environmental_driver_levels,
      drop = FALSE,
      name = "Variable"
    ) +
    scale_shape_manual(
      values = c("Positive" = 19, "Negative" = 4, "No significant driver" = 21),
      limits = c("Positive", "Negative", no_significant_driver),
      breaks = c("Positive", "Negative", no_significant_driver),
      labels = c("Positive", "Negative", "No significant association"),
      drop = FALSE,
      name = "Direction"
    ) +
    guides(
      shape = guide_legend(
        order = 1,
        ncol = 1,
        override.aes = list(size = 3, colour = "black", fill = "white")
      ),
      fill = guide_legend(order = 2, ncol = 1, byrow = TRUE)
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.spacing.y = grid::unit(0, "pt"),
      legend.box.spacing = grid::unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.25),
      axis.text = element_blank()
    )
}


figure_legend <- get_legend(
  figure1_style_combined_legend_plot() +
    theme(
      legend.position = "right",
      legend.direction = "vertical"
    )
)

facet_title_row <- plot_grid(
  NULL,
  facet_strip("1970s"),
  facet_strip("1990s"),
  facet_strip("2010s"),
  nrow = 1,
  align = "h",
  axis = "tb",
  rel_widths = c(0.075, 1, 1, 1)
)

response_map_row <- function(response_name) {
  plot_grid(
    period_strip(response_name),
    map_panel(response_name, "1970s"),
    map_panel(response_name, "1990s"),
    map_panel(response_name, "2010s"),
    nrow = 1,
    align = "hv",
    axis = "tblr",
    rel_widths = c(0.075, 1, 1, 1)
  )
}

map_rows <- plot_grid(
  response_map_row("Species richness"),
  response_map_row("Beta turnover"),
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

main_panels <- plot_grid(
  facet_title_row,
  map_rows,
  ncol = 1,
  rel_heights = c(0.045, 1),
  align = "v",
  axis = "lr"
)

fig3_gwr_richness_turnover_maps <- plot_grid(
  main_panels,
  figure_legend,
  nrow = 1,
  rel_widths = c(1, 0.25),
  align = "h",
  axis = "tb"
) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(
  figure_path,
  fig3_gwr_richness_turnover_maps,
  width = 16,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Saved figure: ", figure_path)
message("Saved map values: ", map_values_path)
message("Saved map summary: ", map_summary_path)
