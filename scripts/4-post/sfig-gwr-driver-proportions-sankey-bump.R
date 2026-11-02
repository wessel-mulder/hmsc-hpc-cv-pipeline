#!/usr/bin/env Rscript

# Supplementary diagnostic: proportion of GWR grid cells assigned to each
# dominant driver through time for species richness and beta turnover.

rm(list = ls())

required_packages <- c("dplyr", "ggplot2", "ggsankey", "readr", "scales", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install required package(s) before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggsankey)
  library(readr)
  library(scales)
  library(tidyr)
})

pattern <- "2026-03-13"
gwr_driver_path <- file.path(
  "intermediate-outputs",
  "gwr-output",
  "fig10-gwr-richness-turnover",
  "gwr",
  paste0(pattern, "-gwr-dominant-drivers-100km.csv")
)

figure_out_dir <- file.path("manuscript-outputs", "figures", "supplement")
table_out_dir <- file.path("manuscript-outputs", "tables", "supplement")

figure_path <- file.path(
  figure_out_dir,
  paste0(pattern, "-sfig-gwr-driver-proportions-sankey-bump.png")
)
table_path <- file.path(
  table_out_dir,
  paste0(pattern, "-sfig-gwr-driver-proportions-sankey-bump-data.csv")
)

if (!file.exists(gwr_driver_path)) {
  stop("Missing GWR dominant-driver table: ", gwr_driver_path, call. = FALSE)
}

dir.create(figure_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_out_dir, recursive = TRUE, showWarnings = FALSE)

period_levels <- c("1970s", "1990s", "2010s")
response_labels <- c(
  "Predicted richness" = "Species richness",
  "All-site beta turnover" = "Beta turnover"
)

driver_colours <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban (% coverage)" = "snow3",
  "Cropland (% coverage)" = "goldenrod1",
  "Pasture (% coverage)" = "darkorange",
  "Forest (% coverage)" = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2",
  "No significant driver" = "#d9d9d9"
)

driver_levels <- names(driver_colours)

dominant_drivers <- read_csv(
  gwr_driver_path,
  col_types = cols(
    period = col_character(),
    .default = col_guess()
  )
) |>
  filter(.data$response_label %in% names(response_labels)) |>
  mutate(
    response_panel = unname(response_labels[.data$response_label]),
    response_panel = factor(
      .data$response_panel,
      levels = unname(response_labels)
    ),
    period = factor(.data$period, levels = period_levels),
    dominant_driver = factor(.data$dominant_driver, levels = driver_levels)
  )

if (any(is.na(dominant_drivers$period))) {
  stop("Unexpected period labels in: ", gwr_driver_path, call. = FALSE)
}
if (any(is.na(dominant_drivers$dominant_driver))) {
  unknown_drivers <- dominant_drivers |>
    filter(is.na(.data$dominant_driver)) |>
    distinct(.data$dominant_driver) |>
    pull(.data$dominant_driver)
  stop(
    "Unexpected dominant-driver labels in ",
    gwr_driver_path,
    ": ",
    paste(unknown_drivers, collapse = ", "),
    call. = FALSE
  )
}

panel_totals <- dominant_drivers |>
  count(.data$response_panel, .data$period, name = "total_grid_cells")

driver_proportions <- dominant_drivers |>
  count(.data$response_panel, .data$period, .data$dominant_driver, name = "n_grid_cells") |>
  complete(
    response_panel,
    period,
    dominant_driver = factor(driver_levels, levels = driver_levels),
    fill = list(n_grid_cells = 0L)
  ) |>
  left_join(panel_totals, by = c("response_panel", "period")) |>
  mutate(
    proportion = .data$n_grid_cells / .data$total_grid_cells,
    percent_grid_cells = .data$proportion * 100
  ) |>
  arrange(.data$response_panel, .data$period, .data$dominant_driver)

write_csv(driver_proportions, table_path)

percent_decimal_labels <- function(x) {
  vals <- format(round(x * 100, 1), trim = TRUE, scientific = FALSE)
  paste0(sub("\\.0$", "", vals), "%")
}

gwr_driver_sankey <- ggplot(
  driver_proportions,
  aes(
    x = period,
    node = dominant_driver,
    fill = dominant_driver,
    value = proportion,
    group = dominant_driver
  )
) +
  geom_sankey_bump(space = 0, type = "alluvial", alpha = 0.86) +
  facet_wrap(vars(response_panel), ncol = 1) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    minor_breaks = seq(0, 1, by = 0.125),
    labels = percent_decimal_labels,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_fill_manual(
    values = driver_colours,
    limits = driver_levels,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Grid cells (%)",
    fill = "Dominant driver"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.32, "cm"),
    legend.key.width = unit(0.42, "cm"),
    legend.spacing.y = unit(0.08, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.25),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(face = "bold", size = 11),
    plot.background = element_rect(fill = "grey99", color = NA)
  )


ggsave(
  figure_path,
  gwr_driver_sankey,
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Saved figure: ", figure_path)
message("Saved plotted data: ", table_path)
