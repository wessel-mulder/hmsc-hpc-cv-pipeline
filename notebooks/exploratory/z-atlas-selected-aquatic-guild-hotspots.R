# Hotspot maps for selected aquatic foraging guilds across Danish atlas projects.
#
# Purpose:
#   Make quick maps showing where the full raw atlas dataset contains the
#   highest grid-cell richness of selected aquatic foraging guilds in Atlas 1,
#   Atlas 2, and Atlas 3.
#
# Important analysis choices:
#   - This uses the full raw atlas occurrence matrix: 201 species, 2,165 grid
#     cells in each atlas period, and no HMSC model filtering.
#   - Species groups come from the local `foraging_guild_consensus` trait field.
#   - Group richness is the number of species from the group recorded as present
#     in each 5 km atlas grid cell.
#   - Hotspots are defined separately for each group and atlas period as cells
#     in the upper 10% of non-zero group-richness values. Because richness is an
#     integer count, the retained hotspot percentage can be slightly above 10%.
#   - Plot output is PNG only.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(tidyr)
})

#### PATHS ####

occurrence_path <- file.path("Data", "Y_occurrences", "Y_occurrences.csv")
design_path <- file.path("Data", "data", "1_preprocessing", "design", "studyDesign.csv")
traits_path <- file.path(
  "Data",
  "data",
  "1_preprocessing",
  "Tr_aits",
  "traits-guild_migration.csv"
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
  "atlas-selected-aquatic-guild-hotspots"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species_path <- file.path(out_dir, "selected-aquatic-guilds-species-used.csv")
site_scores_path <- file.path(out_dir, "selected-aquatic-guilds-richness-by-site.csv")
summary_path <- file.path(out_dir, "selected-aquatic-guilds-hotspot-summary.csv")
trend_path <- file.path(out_dir, "selected-aquatic-guilds-trend-summary.csv")
map_png_path <- file.path(out_dir, "selected-aquatic-guilds-hotspot-map.png")
aquatic_map_png_path <- file.path(out_dir, "aquatic-pursuers-hotspot-map.png")
trend_png_path <- file.path(out_dir, "selected-aquatic-guilds-trends.png")

#### CONSTANTS ####

atlas_labels <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

target_guilds <- c("Dabbling ducks", "Scolopacids", "Aquatic pursuers")
expected_species <- 201
expected_sites_per_atlas <- 2165
expected_surveys <- expected_sites_per_atlas * length(atlas_labels)

#### HELPERS ####

read_raw_occurrence_matrix <- function(path) {
  y <- read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  if (nrow(y) != expected_surveys) {
    stop(
      "Expected ",
      expected_surveys,
      " raw atlas survey rows, but found ",
      nrow(y),
      ".",
      call. = FALSE
    )
  }

  if (ncol(y) != expected_species) {
    stop(
      "Expected ",
      expected_species,
      " species columns, but found ",
      ncol(y),
      ".",
      call. = FALSE
    )
  }

  y_matrix <- as.matrix(y)

  if (!is.numeric(y_matrix)) {
    stop("Expected numeric occurrence values.", call. = FALSE)
  }

  if (anyNA(y_matrix)) {
    stop("The raw occurrence matrix contains NA values.", call. = FALSE)
  }

  if (!all(y_matrix %in% c(0, 1))) {
    stop("The raw occurrence matrix should contain only 0/1 values.", call. = FALSE)
  }

  y_matrix
}

read_aligned_design <- function(path, survey_ids) {
  design <- read.csv(path, row.names = 5, stringsAsFactors = FALSE)
  missing_surveys <- setdiff(survey_ids, rownames(design))

  if (length(missing_surveys) > 0) {
    stop(
      "Some occurrence surveys are absent from the design table: ",
      paste(head(missing_surveys, 10), collapse = ", "),
      call. = FALSE
    )
  }

  design <- design[survey_ids, , drop = FALSE]

  if (!identical(rownames(design), survey_ids)) {
    stop("Design rows could not be aligned to occurrence rows.", call. = FALSE)
  }

  atlas_counts <- table(as.character(design$atlas))

  if (!all(atlas_counts[names(atlas_labels)] == expected_sites_per_atlas)) {
    stop("Expected 2,165 grid cells in each atlas period.", call. = FALSE)
  }

  design
}

read_group_species <- function(path, occurrence_species) {
  traits <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(traits)[1] <- "row_id"

  group_species <- traits |>
    filter(
      .data$foraging_guild_consensus %in% target_guilds,
      .data$latin_DOF_underscores %in% occurrence_species
    ) |>
    transmute(
      group = .data$foraging_guild_consensus,
      species = .data$latin_DOF_underscores
    ) |>
    arrange(.data$group, .data$species)

  missing_groups <- setdiff(target_guilds, unique(group_species$group))

  if (length(missing_groups) > 0) {
    stop(
      "No occurrence species found for guild(s): ",
      paste(missing_groups, collapse = ", "),
      call. = FALSE
    )
  }

  group_species
}

calculate_group_richness <- function(y_matrix, design, group_species) {
  do.call(
    rbind,
    lapply(target_guilds, function(group_name) {
      species <- group_species |>
        filter(.data$group == group_name) |>
        pull(.data$species)

      richness <- rowSums(y_matrix[, species, drop = FALSE])

      data.frame(
        survey = rownames(y_matrix),
        site = as.character(design$site),
        atlas = as.character(design$atlas),
        period = unname(atlas_labels[as.character(design$atlas)]),
        group = group_name,
        species_pool_size = length(species),
        group_richness = as.integer(richness),
        stringsAsFactors = FALSE
      )
    })
  ) |>
    mutate(
      period = factor(.data$period, levels = unname(atlas_labels)),
      group = factor(.data$group, levels = target_guilds)
    )
}

add_hotspot_flags <- function(site_scores) {
  hotspot_summary <- site_scores |>
    group_by(.data$group, .data$atlas, .data$period) |>
    summarise(
      n_cells = n(),
      species_pool_size = first(.data$species_pool_size),
      occupied_cells = sum(.data$group_richness > 0),
      max_richness = max(.data$group_richness),
      mean_richness = mean(.data$group_richness),
      hotspot_threshold = unname(
        quantile(.data$group_richness[.data$group_richness > 0], 0.9)
      ),
      .groups = "drop"
    ) |>
    mutate(
      hotspot_threshold = ceiling(.data$hotspot_threshold)
    )

  site_scores <- site_scores |>
    left_join(
      hotspot_summary |>
        select("group", "atlas", "hotspot_threshold"),
      by = c("group", "atlas")
    ) |>
    mutate(
      is_hotspot = .data$group_richness >= .data$hotspot_threshold &
        .data$group_richness > 0
    )

  hotspot_summary <- site_scores |>
    group_by(.data$group, .data$atlas, .data$period) |>
    summarise(
      n_cells = n(),
      species_pool_size = first(.data$species_pool_size),
      occupied_cells = sum(.data$group_richness > 0),
      hotspot_threshold = first(.data$hotspot_threshold),
      hotspot_cells = sum(.data$is_hotspot),
      hotspot_percent_of_all_cells = 100 * mean(.data$is_hotspot),
      max_richness = max(.data$group_richness),
      mean_richness = mean(.data$group_richness),
      .groups = "drop"
    )

  list(site_scores = site_scores, hotspot_summary = hotspot_summary)
}

#### LOAD AND ALIGN INPUTS ####

y_matrix <- read_raw_occurrence_matrix(occurrence_path)
design <- read_aligned_design(design_path, rownames(y_matrix))
group_species <- read_group_species(traits_path, colnames(y_matrix))

#### CALCULATE GROUP RICHNESS HOTSPOTS ####

site_scores <- calculate_group_richness(y_matrix, design, group_species)
hotspot_result <- add_hotspot_flags(site_scores)
site_scores <- hotspot_result$site_scores
hotspot_summary <- hotspot_result$hotspot_summary

trend_summary <- hotspot_summary |>
  mutate(
    period_year = case_when(
      .data$period == "1970s" ~ 1970,
      .data$period == "1990s" ~ 1990,
      .data$period == "2010s" ~ 2010
    ),
    occupied_percent_of_cells = 100 * .data$occupied_cells / .data$n_cells,
    total_group_occurrences = .data$mean_richness * .data$n_cells,
    mean_richness_occupied_cells = .data$total_group_occurrences /
      .data$occupied_cells
  ) |>
  arrange(.data$group, .data$period_year)

if (!all(table(site_scores$group, site_scores$atlas) == expected_sites_per_atlas)) {
  stop("Each group-atlas combination should contain 2,165 grid cells.", call. = FALSE)
}

#### JOIN TO GRID ####

grid_sf <- st_read(grid_path, quiet = TRUE)

if (!"kvadratkod" %in% names(grid_sf)) {
  stop("Expected a kvadratkod column in the atlas-grid shapefile.", call. = FALSE)
}

if (nrow(grid_sf) != expected_sites_per_atlas ||
    length(unique(grid_sf$kvadratkod)) != expected_sites_per_atlas) {
  stop("Expected 2,165 unique polygons in the atlas-grid shapefile.", call. = FALSE)
}

map_data <- do.call(
  rbind,
  lapply(seq_len(nrow(hotspot_summary)), function(i) {
    this_group <- as.character(hotspot_summary$group[[i]])
    this_atlas <- as.character(hotspot_summary$atlas[[i]])

    atlas_values <- site_scores |>
      filter(
        as.character(.data$group) == this_group,
        .data$atlas == this_atlas
      )

    joined <- grid_sf |>
      left_join(atlas_values, by = c("kvadratkod" = "site"))

    if (anyNA(joined$group_richness)) {
      stop(
        "Some grid cells failed to join for ",
        this_group,
        " Atlas ",
        this_atlas,
        ".",
        call. = FALSE
      )
    }

    joined
  })
)

#### WRITE TABLE OUTPUTS ####

write.csv(group_species, species_path, row.names = FALSE, na = "")
write.csv(site_scores, site_scores_path, row.names = FALSE, na = "")
write.csv(hotspot_summary, summary_path, row.names = FALSE, na = "")
write.csv(trend_summary, trend_path, row.names = FALSE, na = "")

#### PNG HOTSPOT MAP ####

map_data$period <- factor(map_data$period, levels = unname(atlas_labels))
map_data$group <- factor(map_data$group, levels = target_guilds)

max_richness <- max(map_data$group_richness)
richness_breaks <- seq(0, max_richness, by = 1)

hotspot_map <- ggplot(map_data) +
  geom_sf(
    aes(fill = .data$group_richness),
    color = "grey75",
    linewidth = 0.025
  ) +
  geom_sf(
    data = map_data |> filter(.data$is_hotspot),
    fill = NA,
    color = "black",
    linewidth = 0.09
  ) +
  facet_grid(group ~ period) +
  scale_fill_gradientn(
    colors = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    limits = c(0, max_richness),
    breaks = richness_breaks,
    name = "Species\nrichness"
  ) +
  labs(
    title = "Selected aquatic foraging-guild hotspots across atlas projects",
    subtitle = "Raw atlas observations; black outlines mark upper-decile occupied cells within each group and atlas",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_blank(),
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey25"),
    strip.text = element_text(face = "bold", size = 11)
  )

ggsave(
  filename = map_png_path,
  plot = hotspot_map,
  width = 10.8,
  height = 9.6,
  units = "in",
  dpi = 300,
  bg = "white"
)

aquatic_pursuer_map <- ggplot(
  map_data |> filter(.data$group == "Aquatic pursuers")
) +
  geom_sf(
    aes(fill = .data$group_richness),
    color = "grey75",
    linewidth = 0.025
  ) +
  geom_sf(
    data = map_data |>
      filter(.data$group == "Aquatic pursuers", .data$is_hotspot),
    fill = NA,
    color = "black",
    linewidth = 0.09
  ) +
  facet_wrap(~period, nrow = 1) +
  scale_fill_gradientn(
    colors = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    limits = c(0, max_richness),
    breaks = richness_breaks,
    name = "Species\nrichness"
  ) +
  labs(
    title = "Aquatic pursuer hotspots across atlas projects",
    subtitle = "Raw atlas observations; black outlines mark upper-decile occupied cells within each atlas",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text = element_blank(),
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey25"),
    strip.text = element_text(face = "bold", size = 11)
  )

ggsave(
  filename = aquatic_map_png_path,
  plot = aquatic_pursuer_map,
  width = 10.8,
  height = 3.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

#### PNG TREND PLOT ####

trend_long <- trend_summary |>
  select(
    "group",
    "period",
    "period_year",
    occupied_percent_of_cells = "occupied_percent_of_cells",
    hotspot_percent_of_all_cells = "hotspot_percent_of_all_cells",
    mean_richness = "mean_richness"
  ) |>
  pivot_longer(
    cols = c(
      "occupied_percent_of_cells",
      "hotspot_percent_of_all_cells",
      "mean_richness"
    ),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = recode(
      .data$metric,
      occupied_percent_of_cells = "Cells occupied (%)",
      hotspot_percent_of_all_cells = "Hotspot cells (%)",
      mean_richness = "Mean species richness per cell"
    ),
    metric = factor(
      .data$metric,
      levels = c(
        "Cells occupied (%)",
        "Mean species richness per cell",
        "Hotspot cells (%)"
      )
    )
  )

trend_palette <- c(
  "Dabbling ducks" = "#2166ac",
  "Scolopacids" = "#b2182b",
  "Aquatic pursuers" = "#1b9e77"
)

trend_plot <- ggplot(
  trend_long,
  aes(
    x = .data$period,
    y = .data$value,
    group = .data$group,
    color = .data$group
  )
) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.5) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_color_manual(values = trend_palette, name = NULL) +
  labs(
    title = "Guild-level atlas trends for selected aquatic foraging guilds",
    subtitle = "Raw atlas observations summarised across all 2,165 grid cells per period",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey25"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = trend_png_path,
  plot = trend_plot,
  width = 10.5,
  height = 4,
  units = "in",
  dpi = 300,
  bg = "white"
)

#### CONSOLE SUMMARY ####

message("Species used:")
print(group_species)
message("Hotspot summary:")
print(hotspot_summary)
message("Wrote species list: ", species_path)
message("Wrote site scores: ", site_scores_path)
message("Wrote hotspot summary: ", summary_path)
message("Wrote trend summary: ", trend_path)
message("Wrote PNG hotspot map: ", map_png_path)
message("Wrote PNG aquatic pursuer hotspot map: ", aquatic_map_png_path)
message("Wrote PNG trend plot: ", trend_png_path)
