# Atlas 3 rarity-hotspot environment characterization.
#
# Purpose:
#   Ask what distinguishes the Atlas 3 species-rarity hotspots from other
#   Atlas 3 grid cells. This script uses the top-decile hotspot flag from the
#   existing Atlas 3 inverse-occupancy rarity map and compares those cells
#   against the remaining Atlas 3 cells using local environmental covariates.
#
# Important boundary:
#   This script does not estimate protected-area coverage because this repo does
#   not currently include a conservation-area layer such as WDPA, Natura 2000,
#   national protected areas, or IBA polygons. It writes a short note documenting
#   that limitation so the analysis is not mistaken for a conservation overlay.
#
# Output:
#   Plot output is PNG only.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(tibble)
  library(tidyr)
})

source(file.path("support_scripts", "data_helpers.R"))

#### PATHS ####

scores_path <- file.path(
  "notebooks",
  "exploratory",
  "outputs",
  "atlas3-species-rarity-map",
  "atlas3-inverse-occupancy-rarity-scores.csv"
)

environment_path <- file.path(
  "Data",
  "data",
  "1_preprocessing",
  "X_environmental",
  "X_Environmental.csv"
)

grid_path <- file.path(
  "Data",
  "data",
  "1_preprocessing",
  "atlas-grids",
  "DOF_Shapefiles_",
  "DK5km_ED50grid_approx_kvadrkod_DOF.shp"
)

out_dir <- file.path(
  "notebooks",
  "exploratory",
  "outputs",
  "atlas3-species-rarity-map",
  "environment-characterization"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_data_path <- file.path(out_dir, "atlas3-rarity-hotspot-environment-data.csv")
variable_summary_path <- file.path(out_dir, "atlas3-rarity-hotspot-environment-summary.csv")
variable_difference_path <- file.path(out_dir, "atlas3-rarity-hotspot-variable-differences.csv")
dominant_landcover_path <- file.path(out_dir, "atlas3-rarity-hotspot-dominant-landcover-summary.csv")
forest_threshold_path <- file.path(out_dir, "atlas3-rarity-hotspot-forest-threshold-summary.csv")
conservation_note_path <- file.path(out_dir, "atlas3-rarity-hotspot-conservation-overlay-note.txt")

landcover_png_path <- file.path(out_dir, "atlas3-rarity-hotspot-landcover-comparison.png")
difference_png_path <- file.path(out_dir, "atlas3-rarity-hotspot-variable-differences.png")
dominant_map_png_path <- file.path(out_dir, "atlas3-rarity-hotspot-dominant-landcover-map.png")

#### CONSTANTS ####

hotspot_label <- c(
  "FALSE" = "Other Atlas 3 cells",
  "TRUE" = "Top-decile rarity hotspots"
)

landcover_labels <- c(
  perc_urban = "Urban",
  perc_cropland = "Cropland",
  perc_pasture = "Pasture",
  perc_forest = "Forest",
  perc_grass_shrub = "Grass/shrub",
  perc_fresh_saltwater = "Fresh/saltwater"
)

variable_labels <- c(
  tmean_year = "Annual temperature",
  tmean_winter = "Winter temperature",
  tmean_breeding = "Breeding temperature",
  prec_year = "Annual precipitation",
  prec_winter = "Winter precipitation",
  prec_breeding = "Breeding precipitation",
  hh = "Habitat heterogeneity",
  unique = "Land-cover class richness",
  perc_urban = "Urban",
  perc_cropland = "Cropland",
  perc_pasture = "Pasture",
  perc_forest = "Forest",
  perc_grass_shrub = "Grass/shrub",
  perc_fresh_saltwater = "Fresh/saltwater"
)

landcover_palette <- c(
  "Urban" = "#8c510a",
  "Cropland" = "#dfc27d",
  "Pasture" = "#bf812d",
  "Forest" = "#1b7837",
  "Grass/shrub" = "#80cdc1",
  "Fresh/saltwater" = "#2166ac"
)

#### HELPERS ####

safe_pooled_sd <- function(x, group) {
  hotspot_values <- x[group]
  other_values <- x[!group]
  hotspot_sd <- stats::sd(hotspot_values, na.rm = TRUE)
  other_sd <- stats::sd(other_values, na.rm = TRUE)
  hotspot_n <- sum(!is.na(hotspot_values))
  other_n <- sum(!is.na(other_values))

  if (hotspot_n < 2 || other_n < 2) {
    return(NA_real_)
  }

  sqrt(((hotspot_n - 1) * hotspot_sd^2 + (other_n - 1) * other_sd^2) / (hotspot_n + other_n - 2))
}

summarise_one_variable <- function(df, variable) {
  grouped <- df |>
    group_by(.data$hotspot_group) |>
    summarise(
      variable = variable,
      n_cells = n(),
      mean = mean(.data[[variable]], na.rm = TRUE),
      median = median(.data[[variable]], na.rm = TRUE),
      sd = sd(.data[[variable]], na.rm = TRUE),
      q25 = unname(stats::quantile(.data[[variable]], 0.25, na.rm = TRUE)),
      q75 = unname(stats::quantile(.data[[variable]], 0.75, na.rm = TRUE)),
      min = min(.data[[variable]], na.rm = TRUE),
      max = max(.data[[variable]], na.rm = TRUE),
      .groups = "drop"
    )

  grouped
}

variable_difference <- function(df, variable) {
  hotspot_values <- df[[variable]][df$is_top_decile]
  other_values <- df[[variable]][!df$is_top_decile]
  pooled_sd <- safe_pooled_sd(df[[variable]], df$is_top_decile)

  data.frame(
    variable = variable,
    variable_label = unname(variable_labels[[variable]]),
    hotspot_mean = mean(hotspot_values, na.rm = TRUE),
    other_mean = mean(other_values, na.rm = TRUE),
    hotspot_median = median(hotspot_values, na.rm = TRUE),
    other_median = median(other_values, na.rm = TRUE),
    mean_difference = mean(hotspot_values, na.rm = TRUE) - mean(other_values, na.rm = TRUE),
    median_difference = median(hotspot_values, na.rm = TRUE) - median(other_values, na.rm = TRUE),
    standardized_mean_difference = ifelse(
      is.na(pooled_sd) || pooled_sd == 0,
      NA_real_,
      (mean(hotspot_values, na.rm = TRUE) - mean(other_values, na.rm = TRUE)) / pooled_sd
    ),
    stringsAsFactors = FALSE
  )
}

dominant_landcover <- function(df) {
  landcover_cols <- names(landcover_labels)
  landcover_matrix <- as.matrix(df[, landcover_cols, drop = FALSE])
  dominant_col <- landcover_cols[max.col(landcover_matrix, ties.method = "first")]
  unname(landcover_labels[dominant_col])
}

#### LOAD AND PREPARE DATA ####

scores <- read.csv(scores_path, stringsAsFactors = FALSE)

required_score_columns <- c("survey", "site", "rarity_score", "is_top_decile")
missing_score_columns <- setdiff(required_score_columns, names(scores))

if (length(missing_score_columns) > 0) {
  stop("Missing score columns: ", paste(missing_score_columns, collapse = ", "), call. = FALSE)
}

if (sum(scores$is_top_decile) != 191) {
  stop("Expected 191 top-decile Atlas 3 rarity hotspots.", call. = FALSE)
}

environment <- read.csv(
  environment_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
) |>
  rownames_to_column(var = "survey") |>
  clean_lulc_columns()

analysis_data <- scores |>
  select(survey, site, rarity_score, rarity_percentile, rarity_rank, is_top_decile) |>
  left_join(environment, by = "survey") |>
  mutate(
    hotspot_group = factor(
      hotspot_label[as.character(.data$is_top_decile)],
      levels = unname(hotspot_label)
    ),
    dominant_landcover = factor(
      dominant_landcover(pick(all_of(names(landcover_labels)))),
      levels = names(landcover_palette)
    )
  )

if (anyNA(analysis_data$tmean_year)) {
  missing_environment <- analysis_data$survey[is.na(analysis_data$tmean_year)]
  stop(
    "Some Atlas 3 rarity cells are missing environmental data: ",
    paste(head(missing_environment, 10), collapse = ", "),
    call. = FALSE
  )
}

#### SUMMARIES ####

summary_variables <- names(variable_labels)

variable_summary <- do.call(
  rbind,
  lapply(summary_variables, function(variable) summarise_one_variable(analysis_data, variable))
) |>
  mutate(variable_label = unname(variable_labels[.data$variable])) |>
  relocate("variable_label", .after = "variable")

variable_differences <- do.call(
  rbind,
  lapply(summary_variables, function(variable) variable_difference(analysis_data, variable))
) |>
  arrange(desc(abs(.data$standardized_mean_difference)))

dominant_landcover_summary <- analysis_data |>
  count(.data$hotspot_group, .data$dominant_landcover, name = "n_cells") |>
  group_by(.data$hotspot_group) |>
  mutate(group_total = sum(.data$n_cells), percent_cells = 100 * .data$n_cells / .data$group_total) |>
  ungroup() |>
  arrange(.data$hotspot_group, desc(.data$percent_cells))

forest_threshold_summary <- analysis_data |>
  group_by(.data$hotspot_group) |>
  summarise(
    n_cells = n(),
    mean_forest_percent = 100 * mean(.data$perc_forest, na.rm = TRUE),
    median_forest_percent = 100 * median(.data$perc_forest, na.rm = TRUE),
    pct_cells_forest_ge_25 = 100 * mean(.data$perc_forest >= 0.25, na.rm = TRUE),
    pct_cells_forest_ge_50 = 100 * mean(.data$perc_forest >= 0.50, na.rm = TRUE),
    pct_cells_forest_dominant = 100 * mean(.data$dominant_landcover == "Forest", na.rm = TRUE),
    .groups = "drop"
  )

#### WRITE TABLE OUTPUTS ####

write.csv(analysis_data, analysis_data_path, row.names = FALSE, na = "")
write.csv(variable_summary, variable_summary_path, row.names = FALSE, na = "")
write.csv(variable_differences, variable_difference_path, row.names = FALSE, na = "")
write.csv(dominant_landcover_summary, dominant_landcover_path, row.names = FALSE, na = "")
write.csv(forest_threshold_summary, forest_threshold_path, row.names = FALSE, na = "")

writeLines(
  c(
    "Conservation overlay status",
    "",
    "No local protected-area/conservation-area layer was found in this repo during inspection.",
    "That means this analysis cannot determine whether the Atlas 3 rarity hotspots are already covered by protected areas.",
    "",
    "To add that check, provide a polygon layer such as WDPA, Natura 2000, Danish protected areas, or IBA coverage.",
    "The next step would be to intersect that layer with the 5 km atlas grid and summarize protected-area overlap for top-decile rarity hotspots versus other Atlas 3 cells."
  ),
  con = conservation_note_path
)

#### PNG FIGURES ####

landcover_long <- analysis_data |>
  select("hotspot_group", all_of(names(landcover_labels))) |>
  pivot_longer(
    cols = all_of(names(landcover_labels)),
    names_to = "landcover",
    values_to = "proportion"
  ) |>
  mutate(
    landcover_label = factor(unname(landcover_labels[.data$landcover]), levels = names(landcover_palette)),
    percent = 100 * .data$proportion
  )

landcover_plot <- landcover_long |>
  group_by(.data$hotspot_group, .data$landcover_label) |>
  summarise(
    mean_percent = mean(.data$percent, na.rm = TRUE),
    median_percent = median(.data$percent, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(x = .data$landcover_label, y = .data$mean_percent, fill = .data$landcover_label)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.2) +
  facet_wrap(~hotspot_group, nrow = 1) +
  scale_fill_manual(values = landcover_palette, guide = "none") +
  labs(
    title = "Atlas 3 rarity hotspots are compared with other Atlas 3 cells",
    subtitle = "Bars show mean percentage cover per 5 km cell",
    x = NULL,
    y = "Mean cover (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.major.x = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = landcover_png_path,
  plot = landcover_plot,
  width = 10,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)

difference_plot <- variable_differences |>
  mutate(
    variable_label = factor(.data$variable_label, levels = rev(.data$variable_label)),
    direction = ifelse(.data$standardized_mean_difference >= 0, "Higher in hotspots", "Lower in hotspots")
  ) |>
  ggplot(aes(x = .data$standardized_mean_difference, y = .data$variable_label, fill = .data$direction)) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.35) +
  geom_col(width = 0.72) +
  scale_fill_manual(
    values = c("Higher in hotspots" = "#2166ac", "Lower in hotspots" = "#b2182b"),
    name = NULL
  ) +
  labs(
    title = "Environmental differences between Atlas 3 rarity hotspots and other cells",
    subtitle = "Standardized mean differences; larger absolute values are more distinctive",
    x = "Standardized mean difference",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = difference_png_path,
  plot = difference_plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

grid_sf <- st_read(grid_path, quiet = TRUE)

if (!"kvadratkod" %in% names(grid_sf)) {
  stop("Expected kvadratkod in atlas-grid shapefile.", call. = FALSE)
}

map_data <- grid_sf |>
  left_join(analysis_data, by = c("kvadratkod" = "site"))

dominant_map <- ggplot() +
  geom_sf(data = grid_sf, fill = "grey92", color = "white", linewidth = 0.03) +
  geom_sf(
    data = map_data |> filter(.data$is_top_decile),
    aes(fill = .data$dominant_landcover),
    color = "grey20",
    linewidth = 0.08
  ) +
  scale_fill_manual(values = landcover_palette, name = "Dominant\nland cover", drop = FALSE) +
  labs(
    title = "Dominant land cover in Atlas 3 rarity hotspots",
    subtitle = "Only top-decile rarity hotspots are coloured; all other grid cells are grey",
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey25")
  )

ggsave(
  filename = dominant_map_png_path,
  plot = dominant_map,
  width = 7.5,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

#### CONSOLE SUMMARY ####

message("Wrote analysis data: ", analysis_data_path)
message("Wrote variable summary: ", variable_summary_path)
message("Wrote variable differences: ", variable_difference_path)
message("Wrote dominant land-cover summary: ", dominant_landcover_path)
message("Wrote forest threshold summary: ", forest_threshold_path)
message("Wrote conservation overlay note: ", conservation_note_path)
message("Wrote land-cover PNG: ", landcover_png_path)
message("Wrote difference PNG: ", difference_png_path)
message("Wrote dominant land-cover map PNG: ", dominant_map_png_path)
message("Top environmental differences:")
print(
  variable_differences |>
    select(variable_label, hotspot_mean, other_mean, mean_difference, standardized_mean_difference) |>
    head(8)
)
message("Forest threshold summary:")
print(forest_threshold_summary)
