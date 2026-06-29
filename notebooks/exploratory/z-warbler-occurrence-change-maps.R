# Warbler occurrence-change maps through Danish atlas time.
#
# Purpose:
#   Make quick maps showing how all warbler species in the full raw atlas
#   occurrence matrix changed across Atlas 1, Atlas 2, and Atlas 3.
#
# Important analysis choices:
#   - This script uses the full raw atlas dataset, not the filtered HMSC model
#     input. The raw matrix has 2,165 grid cells in each of the three atlas
#     periods and 201 species columns.
#   - "Warblers" are defined from the local taxonomy table as species in the
#     warbler families present in the atlas matrix: Acrocephalidae,
#     Locustellidae, Phylloscopidae, and Sylviidae.
#   - Occurrence change is mapped as binary cell-level change:
#       Gain         = absent in the start atlas, present in the end atlas
#       Loss         = present in the start atlas, absent in the end atlas
#       Present both = present in both atlases
#       Absent both  = absent in both atlases
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
  "warbler-occurrence-change-maps"
)

species_map_dir <- file.path(out_dir, "species-change-maps")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(species_map_dir, recursive = TRUE, showWarnings = FALSE)

warbler_species_path <- file.path(out_dir, "warbler-species-used.csv")
occupancy_summary_path <- file.path(out_dir, "warbler-occupancy-change-summary.csv")
site_richness_path <- file.path(out_dir, "warbler-richness-by-site-atlas.csv")
richness_change_png_path <- file.path(out_dir, "warbler-richness-change-map.png")
total_change_overview_png_path <- file.path(out_dir, "warbler-species-total-change-overview.png")

#### CONSTANTS ####

atlas_labels <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

atlas_comparisons <- data.frame(
  comparison = c("1970s to 1990s", "1990s to 2010s", "1970s to 2010s"),
  start_atlas = c("1", "2", "1"),
  end_atlas = c("2", "3", "3"),
  stringsAsFactors = FALSE
)

comparison_levels <- atlas_comparisons$comparison

warbler_families <- c(
  "Acrocephalidae",
  "Locustellidae",
  "Phylloscopidae",
  "Sylviidae"
)

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

read_warbler_species <- function(taxonomy_path, occurrence_species) {
  taxonomy <- read.csv(taxonomy_path, check.names = FALSE, stringsAsFactors = FALSE)
  names(taxonomy)[1] <- "row_id"

  warblers <- taxonomy |>
    filter(
      .data$Family2_AVONET %in% warbler_families,
      .data$latin_DOF_underscores %in% occurrence_species
    ) |>
    select(
      species = "latin_DOF_underscores",
      avonet_species = "Species2_AVONET",
      family = "Family2_AVONET",
      genus = "genus"
    ) |>
    arrange(.data$family, .data$species)

  if (nrow(warblers) == 0) {
    stop("No warbler species were found in the raw occurrence matrix.", call. = FALSE)
  }

  warblers
}

attach_common_names <- function(warblers, sti_path) {
  sti <- standardise_sti_species_names(read.csv(sti_path, check.names = FALSE, stringsAsFactors = FALSE))

  # The local STI helper already maps Sylvia communis/curruca to Curruca. Add
  # the same current-taxonomy label for Barred Warbler, which appears in the
  # atlas taxonomy as Curruca nisoria but in STI as Sylvia nisoria.
  sti$SPECIES <- ifelse(sti$SPECIES == "Sylvia_nisoria", "Curruca_nisoria", sti$SPECIES)

  warblers |>
    left_join(
      sti |>
        select(species = "SPECIES", common_name = "ENGL NAME"),
      by = "species"
    ) |>
    mutate(
      common_name = ifelse(is.na(.data$common_name), pretty_species_name(.data$species), .data$common_name),
      plot_label = paste0(.data$common_name, "\n", pretty_species_name(.data$species))
    )
}

site_atlas_matrix <- function(y_matrix, design, species) {
  data.frame(
    survey = rownames(y_matrix),
    site = as.character(design$site),
    atlas = as.character(design$atlas),
    occurrence = as.integer(y_matrix[, species]),
    stringsAsFactors = FALSE
  )
}

make_change_frame <- function(site_atlas_df, species, species_label) {
  do.call(
    rbind,
    lapply(seq_len(nrow(atlas_comparisons)), function(i) {
      comparison_label <- atlas_comparisons$comparison[[i]]
      start_atlas <- atlas_comparisons$start_atlas[[i]]
      end_atlas <- atlas_comparisons$end_atlas[[i]]

      start_df <- site_atlas_df |>
        filter(.data$atlas == start_atlas) |>
        select(site, start_occurrence = "occurrence")

      end_df <- site_atlas_df |>
        filter(.data$atlas == end_atlas) |>
        select(site, end_occurrence = "occurrence")

      change_df <- start_df |>
        inner_join(end_df, by = "site") |>
        mutate(
          species = species,
          species_label = species_label,
          comparison = comparison_label,
          start_atlas = start_atlas,
          end_atlas = end_atlas,
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
          " in ",
          comparison_label,
          ", but found ",
          nrow(change_df),
          ".",
          call. = FALSE
        )
      }

      change_df
    })
  )
}

summarise_change_frame <- function(change_df) {
  change_df |>
    group_by(.data$species, .data$species_label, .data$comparison, .data$change_class) |>
    summarise(n_cells = n(), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = "change_class",
      values_from = "n_cells",
      values_fill = 0
    ) |>
    mutate(
      net_change = .data$Gain - .data$Loss,
      turnover = .data$Gain + .data$Loss
    )
}

plot_species_change_map <- function(grid_sf, change_df, species, species_label, output_path) {
  plot_sf <- grid_sf |>
    left_join(change_df, by = c("kvadratkod" = "site")) |>
    mutate(
      comparison = factor(.data$comparison, levels = comparison_levels),
      change_class = factor(.data$change_class, levels = change_levels)
    )

  missing_change_rows <- plot_sf |>
    sf::st_drop_geometry() |>
    filter(is.na(.data$change_class))

  if (nrow(missing_change_rows) > 0) {
    stop("Grid join failed for ", species, ".", call. = FALSE)
  }

  p <- ggplot(plot_sf) +
    geom_sf(aes(fill = .data$change_class), color = "grey75", linewidth = 0.03) +
    facet_wrap(~comparison, nrow = 1) +
    scale_fill_manual(
      values = change_palette,
      breaks = change_levels,
      drop = FALSE,
      name = NULL
    ) +
    labs(
      title = species_label,
      subtitle = "Full atlas grid; binary occurrence change by cell",
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
      plot.subtitle = element_text(color = "grey25"),
      strip.text = element_text(face = "bold")
    )

  ggsave(
    filename = output_path,
    plot = p,
    width = 11,
    height = 4.6,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

build_warbler_richness <- function(y_matrix, design, warbler_species) {
  data.frame(
    survey = rownames(y_matrix),
    site = as.character(design$site),
    atlas = as.character(design$atlas),
    period = unname(atlas_labels[as.character(design$atlas)]),
    warbler_richness = rowSums(y_matrix[, warbler_species, drop = FALSE]),
    stringsAsFactors = FALSE
  )
}

build_richness_change_frame <- function(warbler_richness) {
  do.call(
    rbind,
    lapply(seq_len(nrow(atlas_comparisons)), function(i) {
      comparison_label <- atlas_comparisons$comparison[[i]]
      start_atlas <- atlas_comparisons$start_atlas[[i]]
      end_atlas <- atlas_comparisons$end_atlas[[i]]

      start_df <- warbler_richness |>
        filter(.data$atlas == start_atlas) |>
        select(site, start_warbler_richness = "warbler_richness")

      end_df <- warbler_richness |>
        filter(.data$atlas == end_atlas) |>
        select(site, end_warbler_richness = "warbler_richness")

      start_df |>
        inner_join(end_df, by = "site") |>
        mutate(
          comparison = comparison_label,
          delta_warbler_richness = .data$end_warbler_richness - .data$start_warbler_richness
        )
    })
  )
}

plot_warbler_richness_change_map <- function(grid_sf, richness_change, output_path) {
  plot_sf <- grid_sf |>
    left_join(richness_change, by = c("kvadratkod" = "site")) |>
    mutate(comparison = factor(.data$comparison, levels = comparison_levels))

  delta_limit <- max(abs(plot_sf$delta_warbler_richness), na.rm = TRUE)

  p <- ggplot(plot_sf) +
    geom_sf(aes(fill = .data$delta_warbler_richness), color = "grey70", linewidth = 0.03) +
    facet_wrap(~comparison, nrow = 1) +
    scale_fill_gradient2(
      low = "#b2182b",
      mid = "#f7f7f7",
      high = "#2166ac",
      midpoint = 0,
      limits = c(-delta_limit, delta_limit),
      breaks = seq(-delta_limit, delta_limit, by = 1),
      name = "Change in\nwarbler spp."
    ) +
    labs(
      title = "Change in warbler species richness",
      subtitle = "Full atlas grid; positive values mean more warbler species in the later atlas",
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
      plot.subtitle = element_text(color = "grey25"),
      strip.text = element_text(face = "bold")
    )

  ggsave(
    filename = output_path,
    plot = p,
    width = 11,
    height = 4.6,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

plot_total_change_overview <- function(grid_sf, all_change_frames, output_path) {
  total_change <- all_change_frames |>
    filter(.data$comparison == "1970s to 2010s") |>
    mutate(
      species_label = factor(.data$species_label, levels = unique(.data$species_label)),
      change_class = factor(.data$change_class, levels = change_levels)
    )

  plot_sf <- grid_sf |>
    left_join(total_change, by = c("kvadratkod" = "site"))

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
      title = "Warbler occurrence change from the 1970s to the 2010s",
      subtitle = "One small map per warbler species in the full atlas matrix",
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
    height = 13,
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

missing_grid_sites <- setdiff(unique(as.character(design$site)), grid_sf$kvadratkod)

if (length(missing_grid_sites) > 0) {
  stop(
    "Some full-atlas design sites were not found in the grid shapefile: ",
    paste(head(missing_grid_sites, 10), collapse = ", "),
    call. = FALSE
  )
}

warblers <- read_warbler_species(taxonomy_path, colnames(y_matrix)) |>
  attach_common_names(sti_path)

warbler_species <- warblers$species

if (length(warbler_species) != 17) {
  stop("Expected 17 warbler species in the full atlas matrix, but found ", length(warbler_species), ".", call. = FALSE)
}

write.csv(warblers, warbler_species_path, row.names = FALSE, na = "")

#### AGGREGATE WARBLER RICHNESS CHANGE ####

warbler_richness <- build_warbler_richness(
  y_matrix = y_matrix,
  design = design,
  warbler_species = warbler_species
)

write.csv(warbler_richness, site_richness_path, row.names = FALSE, na = "")

richness_change <- build_richness_change_frame(warbler_richness)

plot_warbler_richness_change_map(
  grid_sf = grid_sf,
  richness_change = richness_change,
  output_path = richness_change_png_path
)

#### SPECIES-LEVEL CHANGE MAPS ####

all_change_frames <- do.call(
  rbind,
  lapply(seq_len(nrow(warblers)), function(i) {
    species <- warblers$species[[i]]
    species_label <- warblers$plot_label[[i]]

    site_atlas_df <- site_atlas_matrix(
      y_matrix = y_matrix,
      design = design,
      species = species
    )

    change_df <- make_change_frame(
      site_atlas_df = site_atlas_df,
      species = species,
      species_label = species_label
    )

    species_png_path <- file.path(
      species_map_dir,
      paste0(make_species_slug(species), "-occurrence-change-map.png")
    )

    plot_species_change_map(
      grid_sf = grid_sf,
      change_df = change_df,
      species = species,
      species_label = species_label,
      output_path = species_png_path
    )

    change_df
  })
)

occupancy_summary <- summarise_change_frame(all_change_frames)
write.csv(occupancy_summary, occupancy_summary_path, row.names = FALSE, na = "")

plot_total_change_overview(
  grid_sf = grid_sf,
  all_change_frames = all_change_frames,
  output_path = total_change_overview_png_path
)

#### CONSOLE SUMMARY ####

message("Warbler species used: ", length(warbler_species))
message("Wrote species list: ", warbler_species_path)
message("Wrote site richness table: ", site_richness_path)
message("Wrote occupancy change summary: ", occupancy_summary_path)
message("Wrote richness change PNG: ", richness_change_png_path)
message("Wrote total-change overview PNG: ", total_change_overview_png_path)
message("Wrote species-level PNGs to: ", species_map_dir)
message("Species included:")
print(warblers |> select(species, common_name, family))
