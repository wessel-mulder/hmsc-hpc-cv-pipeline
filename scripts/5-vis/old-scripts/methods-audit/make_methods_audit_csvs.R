rm(list = ls())

# Build methods audit CSVs from the Chapter 1 preprocessing files and current
# fitted model outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(readxl)
  library(sf)
  library(stringr)
  library(tidyr)
})

model_pattern <- "2026-07-27"
min_occs <- 5
coveragelevels <- "good_and_average"

base_dir <- "HmscOutputs"
output_dir <- file.path("manuscript-outputs", "tables", "methods-audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

period_lookup <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

status_flags <- c("not_breeding", "introduced", "rare", "domestic")

atlas_from_folder <- function(model_folder) {
  str_match(model_folder, "Atlas([0-9]+)")[, 2]
}

survey_atlas <- function(survey) {
  str_match(survey, "_([123])$")[, 2]
}

species_underscore <- function(x) {
  x |>
    str_replace_all(" ", "_") |>
    str_replace_all("\\.", "_")
}

model_folders <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE) |>
  keep(~ str_detect(.x, fixed(model_pattern))) |>
  sort()

if (length(model_folders) == 0) {
  stop("No model folders found for pattern: ", model_pattern, call. = FALSE)
}

read_model_output <- function(model_folder) {
  atlas <- atlas_from_folder(model_folder)
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
  model <- fit_env$fitted_model$posteriors

  list(
    atlas = atlas,
    folder = model_folder,
    Y = as.data.frame(model$Y, check.names = FALSE),
    TrData = as.data.frame(model$TrData, check.names = FALSE),
    XData = as.data.frame(model$XData, check.names = FALSE),
    studyDesign = as.data.frame(model$studyDesign, check.names = FALSE)
  )
}

model_outputs <- map(model_folders, read_model_output)
names(model_outputs) <- map_chr(model_outputs, "atlas")

model_output_for_atlas <- function(atlas) {
  matches <- keep(model_outputs, ~ identical(.x$atlas, as.character(atlas)))
  if (length(matches) != 1) {
    stop("Expected one model output for atlas ", atlas, ", found ", length(matches), call. = FALSE)
  }
  matches[[1]]
}

# -------------------------------------------------------------------------
# 1. Modelled species, traits, and prevalence

model_species_traits <- model_outputs[[1]]$TrData |>
  tibble::rownames_to_column("species") |>
  mutate(across(where(is.factor), as.character))

prevalence_long <- map_dfr(model_outputs, function(model) {
  tibble(
    species = colnames(model$Y),
    atlas = model$atlas,
    period = period_lookup[[model$atlas]],
    n_grid_cells = nrow(model$Y),
    n_breeding_records = as.integer(colSums(model$Y, na.rm = TRUE)),
    prevalence = n_breeding_records / n_grid_cells
  )
})

prevalence_wide <- prevalence_long |>
  select(
    species,
    atlas,
    n_grid_cells,
    n_breeding_records,
    prevalence
  ) |>
  pivot_wider(
    names_from = atlas,
    values_from = c(n_grid_cells, n_breeding_records, prevalence),
    names_glue = "atlas{atlas}_{.value}"
  )

trait_lookup_source <- read_csv(
  file.path("data-new", "derived", "traits", "birdTraits-merged-subset.csv"),
  show_col_types = FALSE
) |>
  mutate(
    species_primary = species_underscore(.data$latin_DOF),
    species_ioc_underscore = species_underscore(.data$species_IOC)
  )

trait_lookup <- bind_rows(
  trait_lookup_source |>
    mutate(species = .data$species_primary),
  trait_lookup_source |>
    filter(!is.na(.data$species_ioc_underscore), .data$species_ioc_underscore != .data$species_primary) |>
    mutate(species = .data$species_ioc_underscore)
) |>
  select(
    species,
    artnr_DOF,
    artnavn_DOF,
    latin_name_dof = latin_DOF,
    english_DOF,
    species_IOC,
    order_IOC,
    family_IOC,
    genus_IOC,
    Habitat_AVONET,
    Habitat.Density_AVONET,
    Trophic.Level_AVONET,
    Trophic.Niche_AVONET,
    Primary.Lifestyle_AVONET,
    trend_DOF,
    trend_a1a2_DOF,
    tend_a2a3_DOF,
    trend_a1a2a3_DOF
  ) |>
  distinct(.data$species, .keep_all = TRUE)

modelled_species_traits_prevalence <- model_species_traits |>
  left_join(trait_lookup, by = "species") |>
  left_join(prevalence_wide, by = "species") |>
  relocate(
    species,
    latin_name_dof,
    english_DOF,
    order_IOC,
    family_IOC,
    Migration_a3_DOF,
    foraging_guild_consensus,
    species_thermal_index
  ) |>
  select(
    species,
    latin_name_dof,
    english_DOF,
    Migration_a3_DOF,
    foraging_guild_consensus,
    species_thermal_index,
    atlas1_n_grid_cells,
    atlas2_n_grid_cells,
    atlas3_n_grid_cells,
    atlas1_n_breeding_records,
    atlas2_n_breeding_records,
    atlas3_n_breeding_records,
    atlas1_prevalence,
    atlas2_prevalence,
    atlas3_prevalence
  ) |>
  arrange(.data$latin_name_dof, .data$species)

write_csv(
  modelled_species_traits_prevalence,
  file.path(output_dir, "modelled_species_traits_prevalence.csv")
)

# -------------------------------------------------------------------------
# Shared preprocessing helpers for grid and species summaries

get_effort <- function(dir) {
  atlas1 <- read_excel(file.path(dir, "effort", "dækning atlas1+2.xlsx"), sheet = "kvadrater atlas1")
  atlas2 <- read_excel(file.path(dir, "effort", "dækning atlas1+2.xlsx"), sheet = "kvadrater atlas2")

  atlas2 <- atlas2 |>
    select(-UNDS_KODE_) |>
    rename(UNDS_KODE_ = UNDS_KOD_1)

  key <- read_excel(file.path(dir, "effort", "dækning atlas1+2.xlsx"), sheet = "dækning") |>
    mutate(
      coverage = case_when(
        .data$UNDSKODE == 0 ~ "Not understood",
        .data$UNDSKODE == 1 ~ "Good",
        .data$UNDSKODE == 2 ~ "Average",
        .data$UNDSKODE == 3 ~ "Bad",
        TRUE ~ NA_character_
      )
    )

  bind_rows(
    atlas1 |> mutate(atlas = 1),
    atlas2 |> mutate(atlas = 2)
  ) |>
    left_join(key |> rename(UNDS_KODE_ = UNDSKODE), by = "UNDS_KODE_") |>
    rename(site = KVADRATKOD) |>
    mutate(
      atlas = as.character(.data$atlas),
      survey = paste(.data$site, .data$atlas, sep = "_"),
      coverage = factor(.data$coverage, levels = c("Good", "Average", "Bad", "Not understood"))
    ) |>
    select(survey, site, atlas, coverage)
}

prepare_environment <- function() {
  lulc_lookup <- c(
    "0" = "ocean",
    "11" = "urban",
    "22" = "cropland",
    "33" = "pasture",
    "44" = "forest",
    "55" = "grass_shrub",
    "66" = "other",
    "77" = "water"
  )

  read.csv(
    file.path("data-new", "derived", "environmental", "X_Environmental.csv"),
    row.names = 1,
    check.names = FALSE
  ) |>
    as_tibble(rownames = "survey") |>
    rename_with(
      .cols = matches("^LULC_"),
      .fn = function(col) {
        lulc_code <- str_extract(col, "\\d+")
        paste0("perc_", lulc_lookup[lulc_code])
      }
    ) |>
    mutate(
      perc_fresh_saltwater = .data$perc_water + .data$perc_ocean
    ) |>
    select(-perc_other, -perc_ocean, -perc_water) |>
    mutate(
      dominant = {
        perc_data <- pick(starts_with("perc_"))
        factor(names(perc_data)[max.col(perc_data, ties.method = "first")])
      }
    )
}

grid_thresholds <- st_read(
  file.path(
    "data-new",
    "derived",
    "atlas-grid",
    "grids-ocean-thresholds",
    "grids_ocean_thresholds.shp"
  ),
  quiet = TRUE
) |>
  st_drop_geometry()

land_sites <- grid_thresholds |>
  filter(.data$pct_lnd >= 25) |>
  pull(.data$kvdrtkd) |>
  unique() |>
  sort()

effort <- get_effort(file.path("data-new", "raw")) |>
  bind_rows(
    tibble(
      survey = paste(land_sites, "3", sep = "_"),
      site = land_sites,
      atlas = "3",
      coverage = factor("Good", levels = c("Good", "Average", "Bad", "Not understood"))
    )
  )

if (coveragelevels == "good_and_average") {
  surveys2keep <- effort$survey[effort$coverage %in% c("Good", "Average")]
} else if (coveragelevels == "good") {
  surveys2keep <- effort$survey[effort$coverage %in% "Good"]
} else {
  stop("Unknown coveragelevels: ", coveragelevels, call. = FALSE)
}

X_environment <- prepare_environment()

X_land <- X_environment |>
  filter(str_remove(.data$survey, "_[123]$") %in% land_sites)

X_complete_land <- X_land |>
  drop_na()

X_final <- X_complete_land |>
  filter(.data$survey %in% surveys2keep)

Y_raw <- read.csv(
  file.path("data-new", "derived", "occurrences", "Y_occurrences.csv"),
  row.names = 1,
  check.names = FALSE
)

Y_after_grid_filters <- Y_raw[rownames(Y_raw) %in% X_final$survey, , drop = FALSE]

# -------------------------------------------------------------------------
# 2. Grid-cell processing summary

grid_stage_counts <- map_dfr(names(period_lookup), function(atlas) {
  atlas_id <- as.character(atlas)
  model_output <- model_output_for_atlas(atlas_id)
  all_surveys <- paste(grid_thresholds$kvdrtkd, atlas_id, sep = "_")
  land_surveys <- paste(land_sites, atlas_id, sep = "_")
  effort_surveys <- intersect(land_surveys, surveys2keep)
  final_surveys <- intersect(effort_surveys, X_complete_land$survey)
  model_sites <- model_output$studyDesign$site |> as.character() |> unique()

  tibble(
    atlas = atlas_id,
    period = period_lookup[[atlas_id]],
    stage_order = 1:5,
    stage = c(
      "all_grid_cells",
      "after_land_cover_threshold_pct_lnd_ge_25",
      "after_effort_filter_good_or_average",
      "after_missing_covariate_filter_complete_cases",
      "final_model_output_grid_cells"
    ),
    n_grid_cells = c(
      length(unique(str_remove(all_surveys, "_[123]$"))),
      length(unique(str_remove(land_surveys, "_[123]$"))),
      length(unique(str_remove(effort_surveys, "_[123]$"))),
      length(unique(str_remove(final_surveys, "_[123]$"))),
      length(model_sites)
    ),
    source_file = c(
      "data-new/derived/atlas-grid/grids-ocean-thresholds/grids_ocean_thresholds.shp",
      "data-new/derived/atlas-grid/grids-ocean-thresholds/grids_ocean_thresholds.shp",
      "data-new/raw/effort/dækning atlas1+2.xlsx; S0_Data_Definitions.R",
      "data-new/derived/environmental/X_Environmental.csv; S0_Data_Definitions.R",
      file.path(
        base_dir,
        model_output$folder,
        "Models/Fitted/HPC_samples_0250_thin_100_chains_4.Rdata"
      )
    ),
    note = c(
      "All grid polygons in the land-threshold shapefile.",
      "Cells retained where pct_lnd >= 25, matching S0_Data_Definitions.R.",
      if (identical(atlas_id, "3")) {
        "Atlas 3 effort is assumed Good for all land-threshold cells in S0_Data_Definitions.R."
      } else {
        "Cells with Good or Average effort coverage, restricted to land-threshold cells."
      },
      "Complete-case environmental covariates after the land and effort filters. S0 applies drop_na before the effort filter; this row is ordered for the requested audit table.",
      "Unique site random-effect levels in the fitted model output."
    )
  ) |>
    mutate(
      n_removed_from_previous = lag(.data$n_grid_cells) - .data$n_grid_cells,
      n_removed_from_start = first(.data$n_grid_cells) - .data$n_grid_cells
    )
})

write_csv(
  grid_stage_counts,
  file.path(output_dir, "grid_cell_processing_summary_long.csv")
)

grid_stage_value <- function(stage_vec, value_vec, stage_name) {
  value_vec[stage_vec == stage_name][1]
}

grid_cell_processing_wide <- grid_stage_counts |>
  group_by(.data$atlas, .data$period) |>
  summarise(
    total_starting_grid_cells = grid_stage_value(.data$stage, .data$n_grid_cells, "all_grid_cells"),
    filtered_area = grid_stage_value(.data$stage, .data$n_removed_from_previous, "after_land_cover_threshold_pct_lnd_ge_25"),
    filtered_incomplete_environmental_variables = grid_stage_value(.data$stage, .data$n_removed_from_previous, "after_missing_covariate_filter_complete_cases"),
    filtered_effort = grid_stage_value(.data$stage, .data$n_removed_from_previous, "after_effort_filter_good_or_average"),
    grid_cells_retained = grid_stage_value(.data$stage, .data$n_grid_cells, "final_model_output_grid_cells"),
    .groups = "drop"
  )

write_csv(
  grid_cell_processing_wide,
  file.path(output_dir, "grid_cell_processing_summary.csv")
)

# -------------------------------------------------------------------------
# 3. Species-processing summary

traits_all <- read_excel(
  file.path("data-new", "raw", "traits", "birdTraits-merged-all.xlsx"),
  .name_repair = "unique_quiet"
) |>
  mutate(
    species = species_underscore(.data$latin_DOF),
    across(all_of(status_flags), ~ replace_na(.x, 0)),
    removed_by_status_flag = rowSums(across(all_of(status_flags), ~ .x == 1)) > 0
  )

count_status_flags <- function(df) {
  tibble(
    removed_not_breeding = sum(df$not_breeding == 1, na.rm = TRUE),
    removed_introduced = sum(df$introduced == 1, na.rm = TRUE),
    removed_rare = sum(df$rare == 1, na.rm = TRUE),
    removed_domestic = sum(df$domestic == 1, na.rm = TRUE)
  )
}

species_status_stages <- map_dfr(names(period_lookup), function(atlas) {
  atlas_id <- as.character(atlas)
  model_output <- model_output_for_atlas(atlas_id)
  occ_col <- paste0("occ_a", atlas_id, "_DOF")
  recorded <- traits_all |>
    filter(!is.na(.data[[occ_col]]), .data[[occ_col]] > 0)
  flagged <- recorded |>
    filter(.data$removed_by_status_flag)
  retained <- recorded |>
    filter(!.data$removed_by_status_flag)
  occurrence_species <- Y_after_grid_filters |>
    as.data.frame(check.names = FALSE) |>
    select(where(~ TRUE))
  occurrence_recorded_n <- sum(colSums(
    occurrence_species[survey_atlas(rownames(occurrence_species)) == atlas_id, , drop = FALSE],
    na.rm = TRUE
  ) > 0)

  bind_rows(
    tibble(
      atlas = atlas_id,
      period = period_lookup[[atlas_id]],
      stage_order = 1,
      stage = "recorded_in_atlas_trait_table_occ_gt_0",
      n_species = nrow(recorded),
      n_removed_from_previous = NA_integer_,
      source_file = "data-new/raw/traits/birdTraits-merged-all.xlsx",
      note = paste0("Recorded species are rows with ", occ_col, " > 0.")
    ),
    tibble(
      atlas = atlas_id,
      period = period_lookup[[atlas_id]],
      stage_order = 2,
      stage = "removed_by_status_flags",
      n_species = nrow(flagged),
      n_removed_from_previous = nrow(flagged),
      source_file = "scripts/0-preprocessing/traits/script_used.R; data-new/raw/traits/birdTraits-merged-all.xlsx",
      note = "Status flags follow script_used.R: not_breeding, introduced, rare, or domestic equals 1."
    ) |>
      bind_cols(count_status_flags(flagged)),
    tibble(
      atlas = atlas_id,
      period = period_lookup[[atlas_id]],
      stage_order = 3,
      stage = "retained_after_status_flags",
      n_species = nrow(retained),
      n_removed_from_previous = nrow(recorded) - nrow(retained),
      source_file = "scripts/0-preprocessing/traits/script_used.R; data-new/raw/traits/birdTraits-merged-all.xlsx",
      note = "Species retained after removing not_breeding, introduced, rare, and domestic rows."
    ),
    tibble(
      atlas = atlas_id,
      period = period_lookup[[atlas_id]],
      stage_order = 4,
      stage = "recorded_after_grid_covariate_effort_filters_before_min_records",
      n_species = occurrence_recorded_n,
      n_removed_from_previous = NA_integer_,
      source_file = "data-new/derived/occurrences/Y_occurrences.csv; S0_Data_Definitions.R",
      note = "Species with at least one breeding record in the occurrence matrix after grid, covariate, and effort filters. Not a direct rowwise continuation of the trait-table status filter because the local occurrence matrix has a slightly different species set."
    ),
    tibble(
      atlas = atlas_id,
      period = period_lookup[[atlas_id]],
      stage_order = 6,
      stage = "final_modelled_species",
      n_species = length(colnames(model_output$Y)),
      n_removed_from_previous = NA_integer_,
      source_file = file.path(
        base_dir,
        model_output$folder,
        "Models/Fitted/HPC_samples_0250_thin_100_chains_4.Rdata"
      ),
      note = "Species columns in the fitted model output."
    )
  )
})

Y_iterative <- Y_after_grid_filters
min_record_removals <- map_dfr(names(period_lookup), function(atlas) {
  atlas_id <- as.character(atlas)
  atlas_rows <- survey_atlas(rownames(Y_iterative)) == atlas_id
  Y_sub <- Y_iterative[atlas_rows, , drop = FALSE]
  records_before <- colSums(Y_sub, na.rm = TRUE)
  to_filter <- names(records_before[records_before < min_occs])
  n_before <- ncol(Y_iterative)

  if (length(to_filter) > 0) {
    Y_iterative <<- Y_iterative[, !(colnames(Y_iterative) %in% to_filter), drop = FALSE]
  }

  tibble(
    atlas = atlas_id,
    period = period_lookup[[atlas_id]],
    stage_order = 5,
    stage = "removed_due_to_insufficient_breeding_records_lt5_iterative",
    n_species = length(to_filter),
    n_removed_from_previous = length(to_filter),
    n_species_before_step = n_before,
    n_species_after_step = ncol(Y_iterative),
    source_file = "S0_Data_Definitions.R; data-new/derived/occurrences/Y_occurrences.csv",
    note = paste0(
      "Matches S0_Data_Definitions.R iterative filter: remove species with < ",
      min_occs,
      " records in this atlas from all atlases before moving to the next atlas."
    )
  )
})

species_processing_summary <- species_status_stages |>
  bind_rows(min_record_removals) |>
  arrange(.data$atlas, .data$stage_order, .data$stage) |>
  mutate(
    removed_not_breeding = replace_na(.data$removed_not_breeding, NA_integer_),
    removed_introduced = replace_na(.data$removed_introduced, NA_integer_),
    removed_rare = replace_na(.data$removed_rare, NA_integer_),
    removed_domestic = replace_na(.data$removed_domestic, NA_integer_)
  )

write_csv(
  species_processing_summary,
  file.path(output_dir, "species_processing_summary_long.csv")
)

recorded_any_atlas <- traits_all |>
  filter(if_any(all_of(paste0("occ_a", names(period_lookup), "_DOF")), ~ !is.na(.x) & .x > 0))

after_uncertain_breeding <- recorded_any_atlas |>
  filter(.data$not_breeding != 1)

introduced_domestic_filter <- (after_uncertain_breeding$introduced == 1) |
  (after_uncertain_breeding$domestic == 1)

after_introduced_domestic <- after_uncertain_breeding |>
  filter(!introduced_domestic_filter)

model_species_sets <- map(model_outputs, ~ colnames(.x$Y))
if (!all(map_lgl(model_species_sets, ~ setequal(.x, model_species_sets[[1]])))) {
  stop("Modelled species sets differ among atlas outputs.", call. = FALSE)
}

n_final_modelled_species <- length(model_species_sets[[1]])

species_processing_wide <- tibble(
  total_species_across_all_atlases_before_filtering = nrow(recorded_any_atlas),
  filtered_uncertain_breeding = nrow(recorded_any_atlas) - nrow(after_uncertain_breeding),
  filtered_introduced_domestic = sum(introduced_domestic_filter, na.rm = TRUE),
  filtered_rarity_min_records = nrow(after_introduced_domestic) - n_final_modelled_species,
  species_retained = n_final_modelled_species
)

write_csv(
  species_processing_wide,
  file.path(output_dir, "species_processing_summary.csv")
)

message("Wrote: ", file.path(output_dir, "modelled_species_traits_prevalence.csv"))
message("Wrote: ", file.path(output_dir, "grid_cell_processing_summary.csv"))
message("Wrote: ", file.path(output_dir, "grid_cell_processing_summary_long.csv"))
message("Wrote: ", file.path(output_dir, "species_processing_summary.csv"))
message("Wrote: ", file.path(output_dir, "species_processing_summary_long.csv"))
