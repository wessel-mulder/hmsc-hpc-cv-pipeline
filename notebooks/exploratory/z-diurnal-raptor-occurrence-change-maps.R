# Diurnal-raptor occurrence-change maps from the 1990s to the 2010s.
#
# Purpose:
#   Map how all diurnal raptor species in the full Danish atlas occurrence
#   matrix changed between Atlas 2 (1990s) and Atlas 3 (2010s).
#
# Important analysis choices:
#   - This script uses the full raw atlas dataset, not the filtered HMSC model
#     input. The raw matrix has 2,165 grid cells in each atlas period and 201
#     species columns.
#   - "Diurnal raptors" are defined from the local taxonomy table as species
#     in the orders Accipitriformes and Falconiformes. Owls are excluded.
#   - Occurrence change is binary at the 5 km grid-cell level:
#       Gain         = absent in the 1990s, present in the 2010s
#       Loss         = present in the 1990s, absent in the 2010s
#       Present both = present in both atlas periods
#       Absent both  = absent in both atlas periods
#   - Only the 1990s-to-2010s comparison is calculated.
#   - Plot output is PNG only.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(sf)
})

source(file.path("support_scripts", "data_helpers.R"))

#### PATHS ####

occurrence_path <- file.path("Data", "Y_occurrences", "Y_occurrences.csv")
design_path <- file.path("Data", "data", "1_preprocessing", "design", "studyDesign.csv")
taxonomy_path <- file.path("Data", "data", "1_preprocessing", "Taxonomy", "taxonomy.csv")
sti_path <- file.path("Data", "sti", "STI_Devictor.csv")
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
  "diurnal-raptor-occurrence-change-maps"
)

species_map_dir <- file.path(out_dir, "species-change-maps")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(species_map_dir, recursive = TRUE, showWarnings = FALSE)

species_path <- file.path(out_dir, "diurnal-raptor-species-used.csv")
change_data_path <- file.path(out_dir, "diurnal-raptor-occurrence-change-by-site.csv")
change_summary_path <- file.path(out_dir, "diurnal-raptor-occurrence-change-summary.csv")
richness_change_path <- file.path(out_dir, "diurnal-raptor-richness-change-by-site.csv")
richness_change_png_path <- file.path(out_dir, "diurnal-raptor-richness-change-map.png")
overview_png_path <- file.path(out_dir, "diurnal-raptor-species-change-overview.png")

#### CONSTANTS ####

start_atlas <- "2"
end_atlas <- "3"
comparison_label <- "1990s to 2010s"

raptor_orders <- c("Accipitriformes", "Falconiformes")

change_levels <- c("Loss", "Gain", "Present both", "Absent both")

change_palette <- c(
  "Loss" = "#b2182b",
  "Gain" = "#2166ac",
  "Present both" = "#525252",
  "Absent both" = "#f0f0f0"
)

#### HELPERS ####

make_species_slug <- function(species) {
  tolower(gsub("[^A-Za-z0-9]+", "-", species))
}

pretty_species_name <- function(species) {
  gsub("_", " ", species, fixed = TRUE)
}

read_raw_occurrence_matrix <- function(path) {
  y <- read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  if (nrow(y) != 6495) {
    stop("Expected 6,495 raw atlas survey rows, but found ", nrow(y), ".", call. = FALSE)
  }

  if (ncol(y) != 201) {
    stop("Expected 201 raw species columns, but found ", ncol(y), ".", call. = FALSE)
  }

  y_matrix <- as.matrix(y)

  if (!is.numeric(y_matrix)) {
    stop("Expected numeric 0/1 values in the raw occurrence matrix.", call. = FALSE)
  }

  if (anyNA(y_matrix)) {
    stop("Raw occurrence matrix contains NA values.", call. = FALSE)
  }

  if (!all(y_matrix %in% c(0, 1))) {
    stop("Raw occurrence matrix should contain only 0/1 values.", call. = FALSE)
  }

  y_matrix
}

read_design <- function(path, survey_ids) {
  design <- read.csv(path, row.names = 5, stringsAsFactors = FALSE)

  missing_design_rows <- setdiff(survey_ids, rownames(design))

  if (length(missing_design_rows) > 0) {
    stop(
      "Some occurrence survey IDs are missing from the design file: ",
      paste(head(missing_design_rows, 10), collapse = ", "),
      call. = FALSE
    )
  }

  design <- design[survey_ids, , drop = FALSE]

  if (!identical(rownames(design), survey_ids)) {
    stop("Design rows could not be aligned to occurrence rows.", call. = FALSE)
  }

  if (!all(table(design$atlas) == 2165)) {
    stop("Expected 2,165 full-atlas grid cells in each atlas period.", call. = FALSE)
  }

  design
}

read_raptor_species <- function(path, occurrence_species) {
  taxonomy <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(taxonomy)[1] <- "row_id"

  raptors <- taxonomy |>
    filter(
      .data$Order2_AVONET %in% raptor_orders,
      .data$latin_DOF_underscores %in% occurrence_species
    ) |>
    select(
      species = "latin_DOF_underscores",
      avonet_species = "Species2_AVONET",
      family = "Family2_AVONET",
      order = "Order2_AVONET",
      genus = "genus"
    ) |>
    arrange(.data$order, .data$family, .data$species)

  if (nrow(raptors) != 14) {
    stop(
      "Expected 14 diurnal raptor species in the full atlas matrix, but found ",
      nrow(raptors),
      ".",
      call. = FALSE
    )
  }

  raptors
}

attach_common_names <- function(raptors, path) {
  sti <- standardise_sti_species_names(
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  )

  raptors |>
    left_join(
      sti |>
        select(species = "SPECIES", common_name = "ENGL NAME"),
      by = "species"
    ) |>
    mutate(
      common_name = ifelse(
        is.na(.data$common_name),
        pretty_species_name(.data$species),
        .data$common_name
      ),
      plot_label = paste0(.data$common_name, "\n", pretty_species_name(.data$species))
    )
}

make_species_change_frame <- function(y_matrix, design, species, species_label) {
  species_data <- data.frame(
    site = as.character(design$site),
    atlas = as.character(design$atlas),
    occurrence = as.integer(y_matrix[, species]),
    stringsAsFactors = FALSE
  )

  start_df <- species_data |>
    filter(.data$atlas == start_atlas) |>
    select(site, start_occurrence = "occurrence")

  end_df <- species_data |>
    filter(.data$atlas == end_atlas) |>
    select(site, end_occurrence = "occurrence")

  change_df <- start_df |>
    inner_join(end_df, by = "site") |>
    mutate(
      species = species,
      species_label = species_label,
      comparison = comparison_label,
      change_class = case_when(
        .data$start_occurrence == 1 & .data$end_occurrence == 0 ~ "Loss",
        .data$start_occurrence == 0 & .data$end_occurrence == 1 ~ "Gain",
        .data$start_occurrence == 1 & .data$end_occurrence == 1 ~ "Present both",
        TRUE ~ "Absent both"
      )
    )

  if (nrow(change_df) != 2165) {
    stop(
      "Expected 2,165 grid-cell comparisons for ",
      species,
      ", but found ",
      nrow(change_df),
      ".",
      call. = FALSE
    )
  }

  change_df
}

summarise_species_change <- function(all_change_data) {
  all_change_data |>
    group_by(.data$species, .data$species_label) |>
    summarise(
      occupied_1990s = sum(.data$start_occurrence),
      occupied_2010s = sum(.data$end_occurrence),
      gains = sum(.data$change_class == "Gain"),
      losses = sum(.data$change_class == "Loss"),
      persistent_presence = sum(.data$change_class == "Present both"),
      net_change = .data$occupied_2010s - .data$occupied_1990s,
      total_turnover = .data$gains + .data$losses,
      percent_change_from_1990s = ifelse(
        .data$occupied_1990s == 0,
        NA_real_,
        100 * .data$net_change / .data$occupied_1990s
      ),
      .groups = "drop"
    ) |>
    arrange(desc(.data$net_change))
}

plot_species_change <- function(grid_sf, change_df, species_label, output_path) {
  plot_sf <- grid_sf |>
    left_join(change_df, by = c("kvadratkod" = "site")) |>
    mutate(change_class = factor(.data$change_class, levels = change_levels))

  if (anyNA(plot_sf$change_class)) {
    stop("Grid join failed for ", species_label, ".", call. = FALSE)
  }

  p <- ggplot(plot_sf) +
    geom_sf(aes(fill = .data$change_class), color = "grey75", linewidth = 0.03) +
    scale_fill_manual(
      values = change_palette,
      breaks = change_levels,
      drop = FALSE,
      name = NULL
    ) +
    labs(
      title = species_label,
      subtitle = "Occurrence change from the 1990s to the 2010s",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text = element_blank(),
      legend.position = "bottom",
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey25")
    )

  ggsave(
    filename = output_path,
    plot = p,
    width = 7.5,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

make_richness_change <- function(y_matrix, design, raptor_species) {
  richness_data <- data.frame(
    site = as.character(design$site),
    atlas = as.character(design$atlas),
    raptor_richness = rowSums(y_matrix[, raptor_species, drop = FALSE]),
    stringsAsFactors = FALSE
  )

  start_df <- richness_data |>
    filter(.data$atlas == start_atlas) |>
    select(site, richness_1990s = "raptor_richness")

  end_df <- richness_data |>
    filter(.data$atlas == end_atlas) |>
    select(site, richness_2010s = "raptor_richness")

  start_df |>
    inner_join(end_df, by = "site") |>
    mutate(delta_raptor_richness = .data$richness_2010s - .data$richness_1990s)
}

plot_richness_change <- function(grid_sf, richness_change, output_path) {
  plot_sf <- grid_sf |>
    left_join(richness_change, by = c("kvadratkod" = "site"))

  delta_limit <- max(abs(plot_sf$delta_raptor_richness), na.rm = TRUE)

  p <- ggplot(plot_sf) +
    geom_sf(aes(fill = .data$delta_raptor_richness), color = "grey70", linewidth = 0.03) +
    scale_fill_gradient2(
      low = "#b2182b",
      mid = "#f7f7f7",
      high = "#2166ac",
      midpoint = 0,
      limits = c(-delta_limit, delta_limit),
      breaks = seq(-delta_limit, delta_limit, by = 1),
      name = "Change in\ndiurnal raptors"
    ) +
    labs(
      title = "Change in diurnal-raptor species richness",
      subtitle = "Full atlas grid, 1990s to 2010s; positive values indicate more species",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text = element_blank(),
      legend.position = "bottom",
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey25")
    )

  ggsave(
    filename = output_path,
    plot = p,
    width = 7.5,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

plot_species_overview <- function(grid_sf, all_change_data, output_path) {
  overview_data <- all_change_data |>
    mutate(
      species_label = factor(.data$species_label, levels = unique(.data$species_label)),
      change_class = factor(.data$change_class, levels = change_levels)
    )

  plot_sf <- grid_sf |>
    left_join(overview_data, by = c("kvadratkod" = "site"))

  p <- ggplot(plot_sf) +
    geom_sf(aes(fill = .data$change_class), color = NA) +
    facet_wrap(~species_label, ncol = 4) +
    scale_fill_manual(
      values = change_palette,
      breaks = change_levels,
      drop = FALSE,
      name = NULL
    ) +
    labs(
      title = "Diurnal-raptor occurrence change from the 1990s to the 2010s",
      subtitle = "One map per species in the full atlas occurrence matrix",
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 8) +
    theme(
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey25"),
      strip.text = element_text(face = "italic", size = 6.5)
    )

  ggsave(
    filename = output_path,
    plot = p,
    width = 12,
    height = 11,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

#### LOAD AND CHECK DATA ####

y_matrix <- read_raw_occurrence_matrix(occurrence_path)
design <- read_design(design_path, rownames(y_matrix))
grid_sf <- st_read(grid_path, quiet = TRUE)

if (!"kvadratkod" %in% names(grid_sf)) {
  stop("Expected a kvadratkod column in the atlas grid shapefile.", call. = FALSE)
}

if (nrow(grid_sf) != 2165 || length(unique(grid_sf$kvadratkod)) != 2165) {
  stop("Expected 2,165 unique grid cells in the atlas-grid shapefile.", call. = FALSE)
}

raptors <- read_raptor_species(taxonomy_path, colnames(y_matrix)) |>
  attach_common_names(sti_path)

raptor_species <- raptors$species
write.csv(raptors, species_path, row.names = FALSE, na = "")

#### SPECIES-LEVEL CHANGE MAPS ####

all_change_data <- do.call(
  rbind,
  lapply(seq_len(nrow(raptors)), function(i) {
    species <- raptors$species[[i]]
    species_label <- raptors$plot_label[[i]]

    change_df <- make_species_change_frame(
      y_matrix = y_matrix,
      design = design,
      species = species,
      species_label = species_label
    )

    species_png_path <- file.path(
      species_map_dir,
      paste0(make_species_slug(species), "-occurrence-change-1990s-to-2010s.png")
    )

    plot_species_change(
      grid_sf = grid_sf,
      change_df = change_df,
      species_label = species_label,
      output_path = species_png_path
    )

    change_df
  })
)

change_summary <- summarise_species_change(all_change_data)

write.csv(all_change_data, change_data_path, row.names = FALSE, na = "")
write.csv(change_summary, change_summary_path, row.names = FALSE, na = "")

plot_species_overview(
  grid_sf = grid_sf,
  all_change_data = all_change_data,
  output_path = overview_png_path
)

#### AGGREGATE RAPTOR-RICHNESS CHANGE ####

richness_change <- make_richness_change(
  y_matrix = y_matrix,
  design = design,
  raptor_species = raptor_species
)

write.csv(richness_change, richness_change_path, row.names = FALSE, na = "")

plot_richness_change(
  grid_sf = grid_sf,
  richness_change = richness_change,
  output_path = richness_change_png_path
)

#### CONSOLE SUMMARY ####

message("Diurnal raptor species used: ", length(raptor_species))
message("Comparison: ", comparison_label)
message("Wrote species list: ", species_path)
message("Wrote site-level change data: ", change_data_path)
message("Wrote change summary: ", change_summary_path)
message("Wrote richness change data: ", richness_change_path)
message("Wrote overview PNG: ", overview_png_path)
message("Wrote richness-change PNG: ", richness_change_png_path)
message("Wrote species-level PNGs to: ", species_map_dir)
message("Species change summary:")
print(
  change_summary |>
    select(
      species,
      common_label = species_label,
      occupied_1990s,
      occupied_2010s,
      gains,
      losses,
      net_change,
      percent_change_from_1990s
    )
)
