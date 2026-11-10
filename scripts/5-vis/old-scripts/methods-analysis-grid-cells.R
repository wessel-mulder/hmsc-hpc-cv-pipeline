rm(list = ls())

# Three-panel methods figure showing the atlas grid cells used in the analysis.
# The map framing and Bornholm inset follow the richness-map figure style.

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(sf)
  library(stringr)
})

base_dir <- "HmscOutputs"
model_pattern <- "2026-07-27"
grid_path <- file.path(
  "data-new",
  "derived",
  "atlas-grid",
  "dof-5km-grid",
  "DK5km_ED50grid_approx_kvadrkod_DOF.shp"
)
output_path <- file.path(
  "manuscript-outputs",
  "figures",
  "methods",
  paste0(model_pattern, "-methods-analysis-grid-cells.png")
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

period_lookup <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

model_folders <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE) |>
  keep(~ str_detect(.x, fixed(model_pattern))) |>
  sort()

if (length(model_folders) == 0) {
  stop("No model folders found for pattern: ", model_pattern, call. = FALSE)
}
if (!file.exists(grid_path)) {
  stop("Missing grid shapefile: ", grid_path, call. = FALSE)
}

read_analysis_sites <- function(model_folder) {
  fit_path <- file.path(
    base_dir,
    model_folder,
    "Models",
    "Fitted",
    "HPC_samples_0250_thin_100_chains_4.Rdata"
  )
  if (!file.exists(fit_path)) {
    stop("Missing fitted model: ", fit_path, call. = FALSE)
  }

  fit_env <- new.env(parent = emptyenv())
  load(fit_path, envir = fit_env)

  atlas <- str_extract(model_folder, "(?<=Atlas)\\d+")
  sites <- fit_env$fitted_model$posteriors$studyDesign |>
    pull(.data$site) |>
    as.character() |>
    unique() |>
    sort()

  tibble(
    atlas = atlas,
    period = period_lookup[[atlas]],
    kvadratkod = sites
  )
}

analysis_sites <- map_dfr(model_folders, read_analysis_sites) |>
  mutate(period = factor(.data$period, levels = unname(period_lookup)))

grid_sf <- st_read(grid_path, quiet = TRUE) |>
  st_transform(25832)

mainland_bbox <- st_bbox(
  c(xmin = 400000, xmax = 750000, ymin = 6000000, ymax = 6450000),
  crs = st_crs(25832)
)
bornholm_bbox <- st_bbox(
  c(xmin = 855000, xmax = 905000, ymin = 6100000, ymax = 6160000),
  crs = st_crs(25832)
)

mainland_width <- as.numeric(mainland_bbox["xmax"] - mainland_bbox["xmin"])
mainland_height <- as.numeric(mainland_bbox["ymax"] - mainland_bbox["ymin"])
bornholm_width <- as.numeric(bornholm_bbox["xmax"] - bornholm_bbox["xmin"])
bornholm_height <- as.numeric(bornholm_bbox["ymax"] - bornholm_bbox["ymin"])

inset_w <- bornholm_width / mainland_width
inset_h <- bornholm_height / mainland_height

plot_grid_panel <- function(period_label) {
  used_sites <- analysis_sites |>
    filter(.data$period == period_label) |>
    pull(.data$kvadratkod)

  plot_data <- grid_sf |>
    mutate(used = .data$kvadratkod %in% used_sites)

  n_cells <- sum(plot_data$used, na.rm = TRUE)
  panel_title <- paste0(period_label, " (n = ", format(n_cells, big.mark = ","), ")")

  draw_map <- function(bbox, title = NULL) {
    ggplot() +
      geom_sf(
        data = plot_data,
        aes(fill = .data$used),
        color = "grey35",
        linewidth = 0.08
      ) +
      scale_fill_manual(
        values = c("FALSE" = "grey96", "TRUE" = "#4C78A8"),
        guide = "none"
      ) +
      coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]),
        expand = FALSE
      ) +
      labs(x = NULL, y = NULL, title = title) +
      theme_minimal() +
      theme(
        axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(6, 0, 0, 0),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11)
      )
  }

  p_main <- draw_map(mainland_bbox, panel_title)
  p_inset <- draw_map(bornholm_bbox) +
    theme_void() +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.5)
    )

  ggdraw(p_main) +
    draw_plot(
      p_inset,
      x = 1 - inset_w - 0.2,
      y = 1 - inset_h - 0.2,
      width = inset_w,
      height = inset_h
    )
}

period_panels <- map(unname(period_lookup), plot_grid_panel)

methods_grid_cells <- plot_grid(
  plotlist = period_panels,
  nrow = 1,
  labels = c("A", "B", "C"),
  label_fontface = "bold",
  label_size = 12
)

ggsave(
  output_path,
  methods_grid_cells,
  width = 10.5,
  height = 4.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Wrote: ", output_path)
