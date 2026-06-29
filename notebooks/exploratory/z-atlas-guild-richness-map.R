# Foraging-guild richness in the 1990s and 2010s at shared Good-coverage sites.
#
# Purpose:
#   For every atlas grid cell, count how many distinct foraging guilds are
#   represented by the bird species recorded as present. Map that guild richness
#   for the 1990s and 2010s, restricted to the same Atlas 2 Good-coverage sites.
#
# Important analysis choices:
#   - This script uses the full raw atlas occurrence matrix: 201 species and
#     2,165 grid cells in each of the three atlas periods.
#   - Atlas 2 Good coverage is read from the available Atlas 1/2 effort workbook.
#     There is no independent Atlas 3 effort table in this repo, so the project's
#     established convention is followed: all available Atlas 3 cells are
#     treated as Good. The retained set is therefore Atlas 2 Good sites that
#     also occur in Atlas 3, giving identical sites in both mapped periods.
#   - Guild membership comes from `foraging_guild_consensus` in the local traits
#     table. All 201 occurrence species have one non-missing guild assignment.
#   - A guild contributes one unit to guild richness when at least one species
#     belonging to that guild is present in the cell. Multiple species from the
#     same guild do not increase guild richness further.
#   - Plot output is PNG only.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(sf)
})

source(file.path("support_scripts", "data_helpers.R"))

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
effort_data_dir <- "Data"

out_dir <- file.path(
  "notebooks",
  "exploratory",
  "outputs",
  "atlas-guild-richness"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

site_guild_richness_path <- file.path(
  out_dir,
  "atlas2-atlas3-good-shared-guild-richness-by-site.csv"
)
atlas_summary_path <- file.path(
  out_dir,
  "atlas2-atlas3-good-shared-guild-richness-summary.csv"
)
guild_species_path <- file.path(out_dir, "foraging-guild-species-membership.csv")
guild_presence_path <- file.path(
  out_dir,
  "atlas2-atlas3-good-shared-guild-presence-by-site.csv"
)
paired_change_path <- file.path(
  out_dir,
  "atlas2-atlas3-good-shared-guild-richness-paired-change.csv"
)
map_png_path <- file.path(
  out_dir,
  "atlas2-atlas3-good-shared-guild-richness-map.png"
)
violin_png_path <- file.path(
  out_dir,
  "atlas2-atlas3-good-shared-guild-richness-violin.png"
)

#### CONSTANTS ####

atlas_labels <- c(
  "2" = "1990s",
  "3" = "2010s"
)

period_levels <- unname(atlas_labels)
expected_species <- 201
expected_sites_per_atlas <- 2165
expected_guilds <- 34
expected_shared_good_sites <- 1465

#### HELPERS ####

read_occurrence_matrix <- function(path) {
  y <- read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  y_matrix <- as.matrix(y)

  if (ncol(y_matrix) != expected_species) {
    stop(
      "Expected ",
      expected_species,
      " species columns, but found ",
      ncol(y_matrix),
      ".",
      call. = FALSE
    )
  }

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

read_aligned_traits <- function(path, species_names) {
  traits <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(traits)[1] <- "row_id"

  missing_species <- setdiff(species_names, traits$latin_DOF_underscores)

  if (length(missing_species) > 0) {
    stop(
      "Some occurrence species are absent from the trait table: ",
      paste(missing_species, collapse = ", "),
      call. = FALSE
    )
  }

  traits <- traits[
    match(species_names, traits$latin_DOF_underscores),
    ,
    drop = FALSE
  ]

  if (!identical(traits$latin_DOF_underscores, species_names)) {
    stop("Trait rows could not be aligned to occurrence species.", call. = FALSE)
  }

  missing_guild <- is.na(traits$foraging_guild_consensus) |
    traits$foraging_guild_consensus == ""

  if (any(missing_guild)) {
    stop(
      "Some species have no foraging-guild assignment: ",
      paste(traits$latin_DOF_underscores[missing_guild], collapse = ", "),
      call. = FALSE
    )
  }

  traits
}

make_species_by_guild_matrix <- function(traits) {
  guild_names <- sort(unique(traits$foraging_guild_consensus))

  incidence <- vapply(
    guild_names,
    function(guild) as.integer(traits$foraging_guild_consensus == guild),
    integer(nrow(traits))
  )

  rownames(incidence) <- traits$latin_DOF_underscores
  colnames(incidence) <- guild_names
  incidence
}

summarise_atlas <- function(site_data) {
  site_data |>
    group_by(.data$atlas, .data$period) |>
    summarise(
      n_cells = n(),
      mean_guild_richness = mean(.data$guild_richness),
      median_guild_richness = median(.data$guild_richness),
      sd_guild_richness = sd(.data$guild_richness),
      q25_guild_richness = unname(quantile(.data$guild_richness, 0.25)),
      q75_guild_richness = unname(quantile(.data$guild_richness, 0.75)),
      min_guild_richness = min(.data$guild_richness),
      max_guild_richness = max(.data$guild_richness),
      mean_species_richness = mean(.data$species_richness),
      .groups = "drop"
    ) |>
    arrange(.data$atlas)
}

#### LOAD AND ALIGN INPUTS ####

y_matrix <- read_occurrence_matrix(occurrence_path)
design <- read_aligned_design(design_path, rownames(y_matrix))
traits <- read_aligned_traits(traits_path, colnames(y_matrix))

species_by_guild <- make_species_by_guild_matrix(traits)
guild_names <- colnames(species_by_guild)

if (length(guild_names) != expected_guilds) {
  stop(
    "Expected ",
    expected_guilds,
    " foraging guilds, but found ",
    length(guild_names),
    ".",
    call. = FALSE
  )
}

#### CALCULATE GUILD RICHNESS ####

# Matrix multiplication gives the number of present species from each guild in
# every survey cell. Converting values greater than zero to TRUE records whether
# each guild is represented, irrespective of how many member species are present.
guild_species_counts <- y_matrix %*% species_by_guild
guild_presence <- guild_species_counts > 0
guild_richness <- rowSums(guild_presence)
species_richness <- rowSums(y_matrix)

site_guild_richness <- data.frame(
  survey = rownames(y_matrix),
  site = as.character(design$site),
  atlas = as.character(design$atlas),
  period = factor(
    unname(atlas_labels[as.character(design$atlas)]),
    levels = period_levels
  ),
  guild_richness = as.integer(guild_richness),
  species_richness = as.integer(species_richness),
  stringsAsFactors = FALSE
)

#### RETAIN SHARED GOOD-COVERAGE SITES ####

effort <- read_effort_coverage(effort_data_dir)

good_atlas2_sites <- effort |>
  filter(
    grepl("_2$", .data$survey),
    !is.na(.data$coverage),
    .data$coverage == "Good"
  ) |>
  transmute(site = sub("_2$", "", .data$survey)) |>
  distinct(.data$site) |>
  pull(.data$site)

atlas3_sites <- site_guild_richness |>
  filter(.data$atlas == "3") |>
  distinct(.data$site) |>
  pull(.data$site)

shared_good_sites <- intersect(good_atlas2_sites, atlas3_sites)

if (length(shared_good_sites) != expected_shared_good_sites) {
  stop(
    "Expected ",
    expected_shared_good_sites,
    " shared Atlas 2 Good / Atlas 3 available sites, but found ",
    length(shared_good_sites),
    ".",
    call. = FALSE
  )
}

site_guild_richness <- site_guild_richness |>
  filter(
    .data$atlas %in% names(atlas_labels),
    .data$site %in% shared_good_sites
  ) |>
  mutate(
    period = factor(
      unname(atlas_labels[.data$atlas]),
      levels = period_levels
    )
  )

retained_counts <- table(site_guild_richness$atlas)

if (!all(retained_counts[names(atlas_labels)] == expected_shared_good_sites)) {
  stop("The retained Atlas 2 and Atlas 3 site counts are not identical.", call. = FALSE)
}

if (any(site_guild_richness$guild_richness > site_guild_richness$species_richness)) {
  stop("Guild richness cannot exceed species richness.", call. = FALSE)
}

if (max(site_guild_richness$guild_richness) > length(guild_names)) {
  stop("Guild richness exceeds the total number of guilds.", call. = FALSE)
}

guild_presence_data <- data.frame(
  survey = rownames(guild_presence),
  site = as.character(design$site),
  atlas = as.character(design$atlas),
  period = unname(atlas_labels[as.character(design$atlas)]),
  guild_presence,
  check.names = FALSE,
  stringsAsFactors = FALSE
) |>
  filter(
    .data$atlas %in% names(atlas_labels),
    .data$site %in% shared_good_sites
  ) |>
  mutate(period = unname(atlas_labels[.data$atlas]))

guild_species_membership <- traits |>
  transmute(
    species = .data$latin_DOF_underscores,
    foraging_guild = .data$foraging_guild_consensus
  ) |>
  arrange(.data$foraging_guild, .data$species)

atlas_summary <- summarise_atlas(site_guild_richness)

paired_change <- site_guild_richness |>
  select("site", "period", "guild_richness") |>
  tidyr::pivot_wider(
    names_from = "period",
    values_from = "guild_richness",
    names_prefix = "guild_richness_"
  ) |>
  mutate(
    delta_2010s_minus_1990s = .data$guild_richness_2010s - .data$guild_richness_1990s
  ) |>
  arrange(desc(.data$delta_2010s_minus_1990s))

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
  lapply(names(atlas_labels), function(atlas_id) {
    atlas_values <- site_guild_richness |>
      filter(.data$atlas == atlas_id)

    joined <- grid_sf |>
      filter(.data$kvadratkod %in% shared_good_sites) |>
      left_join(atlas_values, by = c("kvadratkod" = "site"))

    if (anyNA(joined$guild_richness)) {
      stop("Some Atlas ", atlas_id, " grid cells have no guild-richness value.", call. = FALSE)
    }

    joined
  })
)

#### WRITE TABLE OUTPUTS ####

write.csv(site_guild_richness, site_guild_richness_path, row.names = FALSE, na = "")
write.csv(atlas_summary, atlas_summary_path, row.names = FALSE, na = "")
write.csv(guild_species_membership, guild_species_path, row.names = FALSE, na = "")
write.csv(guild_presence_data, guild_presence_path, row.names = FALSE, na = "")
write.csv(paired_change, paired_change_path, row.names = FALSE, na = "")

#### PNG MAP ####

guild_richness_breaks <- seq(
  floor(min(site_guild_richness$guild_richness) / 2) * 2,
  ceiling(max(site_guild_richness$guild_richness) / 2) * 2,
  by = 2
)

guild_richness_map <- ggplot(map_data) +
  geom_sf(
    aes(fill = .data$guild_richness),
    color = "grey65",
    linewidth = 0.025
  ) +
  facet_wrap(~period, nrow = 1) +
  scale_fill_viridis_c(
    option = "C",
    limits = range(site_guild_richness$guild_richness),
    breaks = guild_richness_breaks,
    name = "Number of\nguilds"
  ) +
  labs(
    title = "Foraging-guild richness at shared Good-coverage sites",
    subtitle = paste0(
      "Atlas 2 Good sites retained in both periods (n = ",
      expected_shared_good_sites,
      "); Atlas 3 has no independent effort classification"
    ),
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
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave(
  filename = map_png_path,
  plot = guild_richness_map,
  width = 11,
  height = 4.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

#### PNG VIOLIN PLOT ####

violin_summary <- atlas_summary |>
  mutate(
    label = sprintf(
      "median %s\nmean %.1f",
      median_guild_richness,
      mean_guild_richness
    ),
    label_y = max(site_guild_richness$guild_richness) + 1.2
  )

guild_richness_violin <- ggplot(
  site_guild_richness,
  aes(x = .data$period, y = .data$guild_richness, fill = .data$period)
) +
  geom_violin(
    width = 0.78,
    trim = FALSE,
    scale = "width",
    color = "grey25",
    linewidth = 0.35
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    color = "grey15",
    linewidth = 0.4
  ) +
  geom_text(
    data = violin_summary,
    aes(x = .data$period, y = .data$label_y, label = .data$label),
    inherit.aes = FALSE,
    size = 3.5,
    lineheight = 0.95
  ) +
  scale_fill_manual(
    values = c("1990s" = "#4daf4a", "2010s" = "#377eb8"),
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = seq(0, max(site_guild_richness$guild_richness) + 2, by = 2),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(
    title = "Distribution of foraging-guild richness",
    subtitle = paste0(
      "The same ",
      expected_shared_good_sites,
      " grid cells are shown in each atlas period"
    ),
    x = NULL,
    y = "Number of guilds represented"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = violin_png_path,
  plot = guild_richness_violin,
  width = 6.5,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

#### CONSOLE SUMMARY ####

message("Foraging guilds represented in the trait table: ", length(guild_names))
message("Shared Atlas 2 Good / Atlas 3 available sites: ", length(shared_good_sites))
message("Wrote site guild richness: ", site_guild_richness_path)
message("Wrote atlas summary: ", atlas_summary_path)
message("Wrote species-to-guild membership: ", guild_species_path)
message("Wrote cell-by-guild presence table: ", guild_presence_path)
message("Wrote paired site change: ", paired_change_path)
message("Wrote PNG map: ", map_png_path)
message("Wrote PNG violin plot: ", violin_png_path)
message("Atlas guild-richness summary:")
print(atlas_summary)
