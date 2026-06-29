# Atlas species rarity hotspot maps.
#
# Purpose:
#   Map where each Danish atlas period holds the strongest concentrations of
#   rare species.
#
# Rarity definition:
#   This script uses a weighted-endemism style score separately within each
#   atlas period. For a given atlas, each species receives a rarity weight of
#   1 / occupancy, where occupancy is the number of grid cells where that
#   species is present in that atlas. Each grid cell's rarity score is the sum
#   of those weights for all species present there.
#
# Data boundary:
#   The occurrence matrix comes from the filtered HMSC model input. That means
#   it uses the same Good/Average coverage, land threshold, environmental
#   completeness, and minimum-occurrence filtering as the current model family.
#
# Output:
#   Plot output is PNG only. A CSV of the mapped scores is also written for each
#   atlas so the highest-scoring cells can be inspected without rerunning the
#   script.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(cowplot)
})

#### PATHS ####

model_input_path <- file.path(
  "Data",
  "preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData"
)

grid_path <- file.path(
  "Data",
  "data",
  "1_preprocessing",
  "atlas-grids",
  "DOF_Shapefiles_",
  "DK5km_ED50grid_approx_kvadrkod_DOF.shp"
)

#### CONSTANTS ####

atlas_specs <- data.frame(
  atlas = c("1", "2", "3"),
  period = c("1970s", "1990s", "2010s"),
  expected_rows = c(1776, 1850, 1904),
  stringsAsFactors = FALSE
)

hotspot_probability <- 0.90

# These bounding boxes match the publication-facing richness-map workflow. The
# Bornholm inset keeps the island readable without stretching the mainland map.
mainland_bbox <- sf::st_bbox(
  c(xmin = 400000, xmax = 750000, ymin = 6000000, ymax = 6450000),
  crs = sf::st_crs(25832)
)

bornholm_bbox <- sf::st_bbox(
  c(xmin = 855000, xmax = 905000, ymin = 6100000, ymax = 6160000),
  crs = sf::st_crs(25832)
)

rarity_palette <- c(
  "#f7fcf5",
  "#d9f0d3",
  "#a6dba0",
  "#5aae61",
  "#1b7837",
  "#00441b"
)

#### HELPERS ####

required_objects_present <- function(env, object_names) {
  all(vapply(object_names, exists, logical(1), envir = env, inherits = FALSE))
}

as_occurrence_matrix <- function(y, atlas_id) {
  # Keep matrix handling explicit because the rarity calculation is a matrix
  # multiplication: site-by-species presence values times species rarity weights.
  y_matrix <- as.matrix(y)

  if (!is.numeric(y_matrix)) {
    stop("Expected numeric 0/1 occurrence values in Y.", call. = FALSE)
  }

  if (anyNA(y_matrix)) {
    stop("Atlas ", atlas_id, " occurrence matrix contains NA values.", call. = FALSE)
  }

  if (!all(y_matrix %in% c(0, 1))) {
    stop(
      "Expected Atlas ",
      atlas_id,
      " occurrence matrix to contain only 0/1 values.",
      call. = FALSE
    )
  }

  y_matrix
}

percentile_rank <- function(x) {
  # A compact base-R percentile rank for visualizing the top decile. The lowest
  # score is 0 and the highest score is 1.
  if (length(x) == 1) {
    return(1)
  }

  (rank(x, ties.method = "average") - 1) / (length(x) - 1)
}

atlas_output_paths <- function(atlas_id) {
  out_dir <- file.path(
    "notebooks",
    "exploratory",
    "outputs",
    paste0("atlas", atlas_id, "-species-rarity-map")
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  list(
    score_csv = file.path(
      out_dir,
      paste0("atlas", atlas_id, "-inverse-occupancy-rarity-scores.csv")
    ),
    png = file.path(
      out_dir,
      paste0("atlas", atlas_id, "-inverse-occupancy-rarity-map.png")
    )
  )
}

calculate_rarity_scores <- function(Y, Design, atlas_id, period_label, expected_rows) {
  atlas_surveys <- rownames(Design)[as.character(Design$atlas) == atlas_id]
  atlas_y <- Y[atlas_surveys, , drop = FALSE]
  atlas_design <- Design[atlas_surveys, , drop = FALSE]

  if (!identical(rownames(atlas_y), rownames(atlas_design))) {
    stop("Y and Design row names do not align after Atlas ", atlas_id, " subsetting.", call. = FALSE)
  }

  if (nrow(atlas_y) != expected_rows) {
    stop(
      "Expected ",
      expected_rows,
      " Atlas ",
      atlas_id,
      " rows, but found ",
      nrow(atlas_y),
      ".",
      call. = FALSE
    )
  }

  if (anyDuplicated(as.character(atlas_design$site)) > 0) {
    stop("Atlas ", atlas_id, " contains duplicate site codes after subsetting.", call. = FALSE)
  }

  atlas_y_matrix <- as_occurrence_matrix(atlas_y, atlas_id)
  species_occupancy <- colSums(atlas_y_matrix)

  if (any(species_occupancy == 0)) {
    zero_occupancy_species <- names(species_occupancy)[species_occupancy == 0]
    stop(
      "Cannot calculate inverse occupancy for zero-occupancy Atlas ",
      atlas_id,
      " species: ",
      paste(zero_occupancy_species, collapse = ", "),
      call. = FALSE
    )
  }

  rarity_weights <- 1 / species_occupancy
  rarity_score <- as.numeric(atlas_y_matrix %*% rarity_weights)
  species_richness <- rowSums(atlas_y_matrix)

  rarity_scores <- data.frame(
    survey = rownames(atlas_y),
    site = as.character(atlas_design$site),
    atlas = atlas_id,
    period = period_label,
    X = atlas_design$lon,
    Y = atlas_design$lat,
    species_richness = species_richness,
    rarity_score = rarity_score,
    rarity_percentile = percentile_rank(rarity_score),
    stringsAsFactors = FALSE
  )

  hotspot_cutoff <- as.numeric(
    stats::quantile(rarity_scores$rarity_score, hotspot_probability, na.rm = TRUE)
  )

  rarity_scores <- rarity_scores |>
    mutate(
      rarity_rank = rank(-rarity_score, ties.method = "min"),
      is_top_decile = rarity_score >= hotspot_cutoff
    ) |>
    arrange(rarity_rank)

  list(
    scores = rarity_scores,
    species_occupancy = species_occupancy,
    hotspot_cutoff = hotspot_cutoff
  )
}

join_scores_to_grid <- function(grid_sf, rarity_scores, atlas_id) {
  missing_grid_sites <- setdiff(rarity_scores$site, grid_sf$kvadratkod)

  if (length(missing_grid_sites) > 0) {
    stop(
      "Some Atlas ",
      atlas_id,
      " sites were not found in the grid shapefile: ",
      paste(head(missing_grid_sites, 10), collapse = ", "),
      call. = FALSE
    )
  }

  plot_sf <- grid_sf |>
    left_join(rarity_scores, by = c("kvadratkod" = "site"))

  mapped_sites <- plot_sf |>
    sf::st_drop_geometry() |>
    filter(!is.na(rarity_score)) |>
    pull(kvadratkod)

  unmapped_sites <- setdiff(rarity_scores$site, mapped_sites)

  if (length(unmapped_sites) > 0) {
    stop(
      "Some Atlas ",
      atlas_id,
      " rarity scores did not join to grid polygons: ",
      paste(head(unmapped_sites, 10), collapse = ", "),
      call. = FALSE
    )
  }

  if (length(mapped_sites) != nrow(rarity_scores)) {
    stop("Grid join produced an unexpected number of mapped Atlas ", atlas_id, " cells.", call. = FALSE)
  }

  plot_sf
}

make_rarity_map <- function(grid_sf, mapped_plot_sf, hotspot_sf, xlim, ylim,
                            show_legend = TRUE, title = NULL,
                            subtitle = NULL) {
  ggplot() +
    geom_sf(
      data = grid_sf,
      fill = "grey94",
      color = "white",
      linewidth = 0.04
    ) +
    geom_sf(
      data = mapped_plot_sf,
      aes(fill = rarity_score),
      color = "grey65",
      linewidth = 0.04
    ) +
    geom_sf(
      data = hotspot_sf,
      fill = NA,
      color = "grey10",
      linewidth = 0.18
    ) +
    scale_fill_gradientn(
      colors = rarity_palette,
      name = "Rarity\nscore",
      labels = function(x) sprintf("%.2f", x)
    ) +
    coord_sf(
      xlim = c(xlim["xmin"], xlim["xmax"]),
      ylim = c(ylim["ymin"], ylim["ymax"]),
      expand = FALSE
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_blank(),
      legend.position = ifelse(show_legend, "right", "none"),
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey25")
    )
}

write_rarity_outputs <- function(grid_sf, plot_sf, rarity_result, atlas_id, period_label) {
  paths <- atlas_output_paths(atlas_id)
  rarity_scores <- rarity_result$scores

  utils::write.csv(rarity_scores, paths$score_csv, row.names = FALSE, na = "")

  mainland_width <- as.numeric(mainland_bbox["xmax"] - mainland_bbox["xmin"])
  mainland_height <- as.numeric(mainland_bbox["ymax"] - mainland_bbox["ymin"])
  bornholm_width <- as.numeric(bornholm_bbox["xmax"] - bornholm_bbox["xmin"])
  bornholm_height <- as.numeric(bornholm_bbox["ymax"] - bornholm_bbox["ymin"])

  inset_width <- bornholm_width / mainland_width
  inset_height <- bornholm_height / mainland_height

  mapped_plot_sf <- plot_sf |>
    filter(!is.na(rarity_score))

  hotspot_sf <- mapped_plot_sf |>
    filter(is_top_decile)

  subtitle <- paste0(
    "Atlas ",
    atlas_id,
    " (",
    period_label,
    ") raw presences; score = sum of inverse atlas-specific occupancy"
  )

  main_map <- make_rarity_map(
    grid_sf = grid_sf,
    mapped_plot_sf = mapped_plot_sf,
    hotspot_sf = hotspot_sf,
    xlim = mainland_bbox,
    ylim = mainland_bbox,
    show_legend = TRUE,
    title = paste0("Atlas ", atlas_id, " species rarity hotspots"),
    subtitle = subtitle
  )

  bornholm_map <- make_rarity_map(
    grid_sf = grid_sf,
    mapped_plot_sf = mapped_plot_sf,
    hotspot_sf = hotspot_sf,
    xlim = bornholm_bbox,
    ylim = bornholm_bbox,
    show_legend = FALSE
  ) +
    theme_void() +
    theme(
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.5),
      plot.background = element_rect(fill = "white", color = NA)
    )

  rarity_map <- cowplot::ggdraw(main_map) +
    cowplot::draw_plot(
      bornholm_map,
      x = 1 - inset_width - 0.18,
      y = 1 - inset_height - 0.18,
      width = inset_width,
      height = inset_height
    )

  ggplot2::ggsave(
    filename = paths$png,
    plot = rarity_map,
    width = 7.5,
    height = 8.2,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  paths
}

print_atlas_summary <- function(rarity_result, atlas_id, paths) {
  rarity_scores <- rarity_result$scores

  message("Atlas ", atlas_id, " outputs")
  message("  Wrote score table: ", paths$score_csv)
  message("  Wrote PNG map: ", paths$png)
  message("  Rows: ", nrow(rarity_scores))
  message("  Species occupancy minimum: ", min(rarity_result$species_occupancy))
  message("  Maximum rarity score: ", round(max(rarity_scores$rarity_score), 6))
  message("  Top-decile cutoff: ", round(rarity_result$hotspot_cutoff, 6))
  message("  Top 10 rarity cells:")
  print(
    rarity_scores |>
      select(site, survey, rarity_score, species_richness, rarity_rank) |>
      head(10)
  )
}

#### LOAD MODEL-INPUT DATA ####

model_env <- new.env(parent = emptyenv())
load(model_input_path, envir = model_env)

if (!required_objects_present(model_env, c("Y", "Design"))) {
  stop(
    "Expected objects Y and Design in ",
    model_input_path,
    ".",
    call. = FALSE
  )
}

Y <- model_env$Y
Design <- model_env$Design

grid_sf <- sf::st_read(grid_path, quiet = TRUE)

if (!"kvadratkod" %in% names(grid_sf)) {
  stop("Expected a kvadratkod column in the atlas grid shapefile.", call. = FALSE)
}

#### CALCULATE, MAP, AND WRITE EACH ATLAS ####

for (i in seq_len(nrow(atlas_specs))) {
  atlas_id <- atlas_specs$atlas[[i]]
  period_label <- atlas_specs$period[[i]]
  expected_rows <- atlas_specs$expected_rows[[i]]

  rarity_result <- calculate_rarity_scores(
    Y = Y,
    Design = Design,
    atlas_id = atlas_id,
    period_label = period_label,
    expected_rows = expected_rows
  )

  plot_sf <- join_scores_to_grid(
    grid_sf = grid_sf,
    rarity_scores = rarity_result$scores,
    atlas_id = atlas_id
  )

  output_paths <- write_rarity_outputs(
    grid_sf = grid_sf,
    plot_sf = plot_sf,
    rarity_result = rarity_result,
    atlas_id = atlas_id,
    period_label = period_label
  )

  print_atlas_summary(
    rarity_result = rarity_result,
    atlas_id = atlas_id,
    paths = output_paths
  )
}
