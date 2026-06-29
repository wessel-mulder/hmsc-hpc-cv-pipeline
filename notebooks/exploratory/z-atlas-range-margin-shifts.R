# Range-margin shifts across adjacent Danish bird atlas projects.
#
# Purpose:
#   Recreate the logic of the British atlas range-margin analysis for the
#   Danish atlas data. For species with southerly start-period distributions,
#   this script measures movement of the northern margin. For species with
#   northerly start-period distributions, it measures movement of the southern
#   margin. The regression intercept estimates directional margin movement for
#   species with no overall change in occupied grid-cell count.
#
# Important analysis choices:
#   - Adjacent contrasts only: 1970s to 1990s and 1990s to 2010s.
#   - Species universe: the 157 HMSC retained species from the Good/Average
#     model-input RData file, excluding diurnal raptors in Accipitriformes and
#     Falconiformes. Owls are retained.
#   - Site universe: the full raw atlas grid, 2,165 5 km cells in each atlas.
#   - Marginal records: 40 occupied 25 km2 cells, chosen to area-match the
#     paper's ten 100 km2 cells.
#   - Main model excludes species with fewer than 40 occupied cells in either
#     atlas and species occupying more than 2,000 cells in either atlas.
#   - A sensitivity model includes the near-ubiquitous species.
#   - Plot output is PNG only.

rm(list = ls())

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(ggplot2)
  library(sf)
})

#### PATHS ####

occurrence_path <- file.path("Data", "Y_occurrences", "Y_occurrences.csv")
design_path <- file.path("Data", "data", "1_preprocessing", "design", "studyDesign.csv")
hmsc_input_path <- file.path(
  "Data",
  "preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData"
)
taxonomy_path <- file.path("Data", "data", "1_preprocessing", "Taxonomy", "taxonomy.csv")
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
  "atlas-range-margin-shifts"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species_data_path <- file.path(out_dir, "atlas-range-margin-shifts-no-diurnal-raptors-species.csv")
excluded_raptors_path <- file.path(
  out_dir,
  "atlas-range-margin-shifts-excluded-diurnal-raptors.csv"
)
exclusions_path <- file.path(out_dir, "atlas-range-margin-shifts-no-diurnal-raptors-exclusions.csv")
main_model_path <- file.path(
  out_dir,
  "atlas-range-margin-shifts-no-diurnal-raptors-main-model-summary.csv"
)
sensitivity_model_path <- file.path(
  out_dir,
  "atlas-range-margin-shifts-no-diurnal-raptors-sensitivity-model-summary.csv"
)
plot_data_path <- file.path(out_dir, "atlas-range-margin-shifts-no-diurnal-raptors-main-plot-data.csv")
main_png_path <- file.path(out_dir, "atlas-range-margin-shifts-no-diurnal-raptors-main.png")

#### CONSTANTS ####

atlas_labels <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

atlas_comparisons <- data.frame(
  comparison = c("1970s to 1990s", "1990s to 2010s"),
  start_atlas = c("1", "2"),
  end_atlas = c("2", "3"),
  stringsAsFactors = FALSE
)

expected_raw_surveys <- 6495
expected_raw_species <- 201
expected_hmsc_species <- 157
expected_cells_per_atlas <- 2165
margin_cells <- 40
ubiquitous_threshold <- 2000
excluded_raptor_orders <- c("Accipitriformes", "Falconiformes")

distribution_levels <- c("southerly", "northerly")

distribution_labels <- c(
  "southerly" = "Southerly species\nnorthern margin",
  "northerly" = "Northerly species\nsouthern margin"
)

variant_labels <- c(
  "main" = "Main model",
  "omit_1_strongest_decline" = "Omit strongest decline",
  "omit_2_strongest_declines" = "Omit two strongest declines",
  "plus_near_ubiquitous" = "Include near-ubiquitous species"
)

#### HELPERS ####

read_raw_occurrence_matrix <- function(path) {
  y <- read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  if (nrow(y) != expected_raw_surveys) {
    stop("Expected 6,495 raw atlas survey rows, but found ", nrow(y), ".", call. = FALSE)
  }

  if (ncol(y) != expected_raw_species) {
    stop("Expected 201 raw species columns, but found ", ncol(y), ".", call. = FALSE)
  }

  y_matrix <- as.matrix(y)

  if (!is.numeric(y_matrix)) {
    stop("Expected numeric 0/1 values in the raw occurrence matrix.", call. = FALSE)
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
      "Some raw occurrence surveys are missing from the design table: ",
      paste(head(missing_surveys, 10), collapse = ", "),
      call. = FALSE
    )
  }

  design <- design[survey_ids, , drop = FALSE]

  if (!identical(rownames(design), survey_ids)) {
    stop("Design rows could not be aligned to raw occurrence rows.", call. = FALSE)
  }

  atlas_counts <- table(as.character(design$atlas))

  if (!all(atlas_counts[names(atlas_labels)] == expected_cells_per_atlas)) {
    stop("Expected 2,165 full-grid cells in each atlas period.", call. = FALSE)
  }

  design
}

read_hmsc_species <- function(path) {
  rdata <- new.env(parent = emptyenv())
  load(path, envir = rdata)

  if (!exists("Y", envir = rdata)) {
    stop("Expected object `Y` in ", path, call. = FALSE)
  }

  species <- colnames(rdata$Y)

  if (length(species) != expected_hmsc_species) {
    stop(
      "Expected ",
      expected_hmsc_species,
      " HMSC retained species, but found ",
      length(species),
      ".",
      call. = FALSE
    )
  }

  species
}

read_excluded_diurnal_raptors <- function(path, hmsc_species) {
  taxonomy <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(taxonomy)[1] <- "row_id"

  excluded <- taxonomy |>
    filter(
      .data$latin_DOF_underscores %in% hmsc_species,
      .data$Order2_AVONET %in% excluded_raptor_orders
    ) |>
    transmute(
      species = .data$latin_DOF_underscores,
      family = .data$Family2_AVONET,
      order = .data$Order2_AVONET,
      exclusion_reason = "diurnal raptor"
    ) |>
    arrange(.data$order, .data$family, .data$species)

  if (nrow(excluded) == 0) {
    warning("No diurnal raptors were found in the HMSC species list.")
  }

  excluded
}

read_site_northings <- function(path, expected_sites) {
  grid_sf <- st_read(path, quiet = TRUE)

  if (!"kvadratkod" %in% names(grid_sf)) {
    stop("Expected a `kvadratkod` column in the atlas grid shapefile.", call. = FALSE)
  }

  if (nrow(grid_sf) != expected_cells_per_atlas ||
      length(unique(grid_sf$kvadratkod)) != expected_cells_per_atlas) {
    stop("Expected 2,165 unique polygons in the atlas grid shapefile.", call. = FALSE)
  }

  missing_sites <- setdiff(expected_sites, grid_sf$kvadratkod)

  if (length(missing_sites) > 0) {
    stop(
      "Some raw atlas sites are missing from the grid shapefile: ",
      paste(head(missing_sites, 10), collapse = ", "),
      call. = FALSE
    )
  }

  centroids <- st_coordinates(st_centroid(st_geometry(grid_sf)))

  data.frame(
    site = grid_sf$kvadratkod,
    northing_m = centroids[, "Y"],
    easting_m = centroids[, "X"],
    stringsAsFactors = FALSE
  )
}

survey_ids <- function(sites, atlas) {
  paste(sites, atlas, sep = "_")
}

mean_margin_northing <- function(northings, distribution_class, n_margin_cells) {
  if (length(northings) < n_margin_cells) {
    return(NA_real_)
  }

  if (distribution_class == "southerly") {
    return(mean(sort(northings, decreasing = TRUE)[seq_len(n_margin_cells)]))
  }

  if (distribution_class == "northerly") {
    return(mean(sort(northings, decreasing = FALSE)[seq_len(n_margin_cells)]))
  }

  NA_real_
}

exclusion_reason <- function(start_occupied, end_occupied, near_ubiquitous, distribution_class) {
  reasons <- character(0)

  if (is.na(distribution_class)) {
    reasons <- c(reasons, "start-period centroid equals national mean northing")
  }

  if (start_occupied < margin_cells) {
    reasons <- c(reasons, paste0("start occupied cells < ", margin_cells))
  }

  if (end_occupied < margin_cells) {
    reasons <- c(reasons, paste0("end occupied cells < ", margin_cells))
  }

  if (near_ubiquitous) {
    reasons <- c(
      reasons,
      paste0("occupied cells > ", ubiquitous_threshold, " in at least one atlas")
    )
  }

  if (length(reasons) == 0) {
    return("")
  }

  paste(reasons, collapse = "; ")
}

calculate_species_contrast <- function(
    species,
    comparison,
    y_matrix,
    sites,
    site_northings,
    national_mean_northing) {
  start_atlas <- comparison$start_atlas
  end_atlas <- comparison$end_atlas
  start_surveys <- survey_ids(sites, start_atlas)
  end_surveys <- survey_ids(sites, end_atlas)

  start_occurrence <- as.integer(y_matrix[start_surveys, species])
  end_occurrence <- as.integer(y_matrix[end_surveys, species])

  start_occupied_sites <- sites[start_occurrence == 1]
  end_occupied_sites <- sites[end_occurrence == 1]
  start_occupied <- length(start_occupied_sites)
  end_occupied <- length(end_occupied_sites)

  start_northings <- site_northings$northing_m[
    match(start_occupied_sites, site_northings$site)
  ]
  end_northings <- site_northings$northing_m[
    match(end_occupied_sites, site_northings$site)
  ]

  start_mean_northing <- if (start_occupied > 0) mean(start_northings) else NA_real_
  distribution_class <- ifelse(
    is.na(start_mean_northing),
    NA_character_,
    ifelse(
      start_mean_northing < national_mean_northing,
      "southerly",
      ifelse(start_mean_northing > national_mean_northing, "northerly", NA_character_)
    )
  )

  start_margin_northing <- mean_margin_northing(
    start_northings,
    distribution_class,
    margin_cells
  )
  end_margin_northing <- mean_margin_northing(
    end_northings,
    distribution_class,
    margin_cells
  )

  near_ubiquitous <- start_occupied > ubiquitous_threshold ||
    end_occupied > ubiquitous_threshold
  margin_eligible <- !is.na(distribution_class) &&
    start_occupied >= margin_cells &&
    end_occupied >= margin_cells
  main_included <- margin_eligible && !near_ubiquitous

  data.frame(
    comparison = comparison$comparison,
    start_atlas = start_atlas,
    end_atlas = end_atlas,
    start_period = unname(atlas_labels[start_atlas]),
    end_period = unname(atlas_labels[end_atlas]),
    species = species,
    n_shared_cells = length(sites),
    start_occupied_cells = start_occupied,
    end_occupied_cells = end_occupied,
    status_change_log10 = ifelse(
      start_occupied > 0 && end_occupied > 0,
      log10(end_occupied / start_occupied),
      NA_real_
    ),
    start_mean_northing_m = start_mean_northing,
    national_mean_northing_m = national_mean_northing,
    distribution_class = distribution_class,
    margin_type = ifelse(
      distribution_class == "southerly",
      "northern margin",
      ifelse(distribution_class == "northerly", "southern margin", NA_character_)
    ),
    start_margin_northing_m = start_margin_northing,
    end_margin_northing_m = end_margin_northing,
    margin_shift_km = (end_margin_northing - start_margin_northing) / 1000,
    near_ubiquitous = near_ubiquitous,
    margin_eligible = margin_eligible,
    main_included = main_included,
    exclusion_reason = exclusion_reason(
      start_occupied,
      end_occupied,
      near_ubiquitous,
      distribution_class
    ),
    stringsAsFactors = FALSE
  )
}

fit_margin_model <- function(data, contrast, distribution_class, variant_id, variant_label) {
  model_data <- data[
    data$comparison == contrast &
      data$distribution_class == distribution_class,
    ,
    drop = FALSE
  ]

  model_data <- model_data[
    is.finite(model_data$status_change_log10) &
      is.finite(model_data$margin_shift_km),
    ,
    drop = FALSE
  ]

  empty_result <- data.frame(
    comparison = contrast,
    distribution_class = distribution_class,
    distribution_label = unname(distribution_labels[distribution_class]),
    model_variant = variant_id,
    model_variant_label = variant_label,
    n_species = nrow(model_data),
    intercept_km = NA_real_,
    intercept_se = NA_real_,
    intercept_p_value = NA_real_,
    slope_km_per_log10_ratio = NA_real_,
    slope_se = NA_real_,
    slope_p_value = NA_real_,
    r_squared = NA_real_,
    adjusted_r_squared = NA_real_,
    strongest_decline_species = ifelse(
      nrow(model_data) > 0,
      model_data$species[which.min(model_data$status_change_log10)],
      NA_character_
    ),
    stringsAsFactors = FALSE
  )

  if (nrow(model_data) < 3 ||
      length(unique(model_data$status_change_log10)) < 2) {
    return(empty_result)
  }

  model <- lm(margin_shift_km ~ status_change_log10, data = model_data)
  coefficients <- broom::tidy(model)
  model_fit <- broom::glance(model)

  data.frame(
    comparison = contrast,
    distribution_class = distribution_class,
    distribution_label = unname(distribution_labels[distribution_class]),
    model_variant = variant_id,
    model_variant_label = variant_label,
    n_species = nrow(model_data),
    intercept_km = coefficients$estimate[coefficients$term == "(Intercept)"],
    intercept_se = coefficients$std.error[coefficients$term == "(Intercept)"],
    intercept_p_value = coefficients$p.value[coefficients$term == "(Intercept)"],
    slope_km_per_log10_ratio =
      coefficients$estimate[coefficients$term == "status_change_log10"],
    slope_se = coefficients$std.error[coefficients$term == "status_change_log10"],
    slope_p_value = coefficients$p.value[coefficients$term == "status_change_log10"],
    r_squared = model_fit$r.squared,
    adjusted_r_squared = model_fit$adj.r.squared,
    strongest_decline_species =
      model_data$species[which.min(model_data$status_change_log10)],
    stringsAsFactors = FALSE
  )
}

make_sensitivity_data <- function(species_data, contrast, distribution_class, variant_id) {
  if (variant_id == "plus_near_ubiquitous") {
    data <- species_data[
      species_data$comparison == contrast &
        species_data$distribution_class == distribution_class &
        species_data$margin_eligible,
      ,
      drop = FALSE
    ]
  } else {
    data <- species_data[
      species_data$comparison == contrast &
        species_data$distribution_class == distribution_class &
        species_data$main_included,
      ,
      drop = FALSE
    ]
  }

  data <- data[
    is.finite(data$status_change_log10) &
      is.finite(data$margin_shift_km),
    ,
    drop = FALSE
  ]

  if (variant_id == "omit_1_strongest_decline" && nrow(data) >= 2) {
    data <- data[order(data$status_change_log10), , drop = FALSE]
    data <- data[-1, , drop = FALSE]
  }

  if (variant_id == "omit_2_strongest_declines" && nrow(data) >= 3) {
    data <- data[order(data$status_change_log10), , drop = FALSE]
    data <- data[-seq_len(2), , drop = FALSE]
  }

  data
}

format_p_value <- function(p) {
  ifelse(is.na(p), "p = NA", ifelse(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 2))))
}

#### LOAD AND CHECK INPUTS ####

y_matrix <- read_raw_occurrence_matrix(occurrence_path)
design <- read_aligned_design(design_path, rownames(y_matrix))
hmsc_species <- read_hmsc_species(hmsc_input_path)
excluded_diurnal_raptors <- read_excluded_diurnal_raptors(
  taxonomy_path,
  hmsc_species
)
analysis_species <- setdiff(hmsc_species, excluded_diurnal_raptors$species)

missing_hmsc_species <- setdiff(hmsc_species, colnames(y_matrix))

if (length(missing_hmsc_species) > 0) {
  stop(
    "Some HMSC retained species are absent from the raw atlas matrix: ",
    paste(missing_hmsc_species, collapse = ", "),
    call. = FALSE
  )
}

if (length(analysis_species) + nrow(excluded_diurnal_raptors) != length(hmsc_species)) {
  stop("The analysis species count does not reconcile with excluded raptors.", call. = FALSE)
}

sites_by_atlas <- split(as.character(design$site), as.character(design$atlas))
full_shared_sites <- Reduce(intersect, sites_by_atlas[names(atlas_labels)])

if (length(full_shared_sites) != expected_cells_per_atlas) {
  stop(
    "Expected ",
    expected_cells_per_atlas,
    " sites shared by all raw atlas periods, but found ",
    length(full_shared_sites),
    ".",
    call. = FALSE
  )
}

site_northings <- read_site_northings(grid_path, full_shared_sites)
site_northings <- site_northings[
  match(full_shared_sites, site_northings$site),
  ,
  drop = FALSE
]

if (!identical(site_northings$site, full_shared_sites)) {
  stop("Site northings could not be aligned to shared atlas sites.", call. = FALSE)
}

national_mean_northing <- mean(site_northings$northing_m)

#### CALCULATE SPECIES-LEVEL MARGIN SHIFTS ####

species_margin_data <- do.call(
  rbind,
  lapply(seq_len(nrow(atlas_comparisons)), function(i) {
    comparison <- atlas_comparisons[i, ]
    contrast_sites <- intersect(
      sites_by_atlas[[comparison$start_atlas]],
      sites_by_atlas[[comparison$end_atlas]]
    )

    if (length(contrast_sites) != expected_cells_per_atlas) {
      stop(
        "Expected ",
        expected_cells_per_atlas,
        " cells for ",
        comparison$comparison,
        ", but found ",
        length(contrast_sites),
        ".",
        call. = FALSE
      )
    }

    contrast_sites <- full_shared_sites[full_shared_sites %in% contrast_sites]

    do.call(
      rbind,
      lapply(analysis_species, calculate_species_contrast,
             comparison = comparison,
             y_matrix = y_matrix,
             sites = contrast_sites,
             site_northings = site_northings,
             national_mean_northing = national_mean_northing)
    )
  })
) |>
  mutate(
    comparison = factor(.data$comparison, levels = atlas_comparisons$comparison),
    distribution_class = factor(.data$distribution_class, levels = distribution_levels),
    distribution_label = unname(distribution_labels[as.character(.data$distribution_class)])
  )

if (any(species_margin_data$margin_eligible & !is.finite(species_margin_data$status_change_log10))) {
  stop("Eligible species include a zero occupancy count, making log10 ratio invalid.", call. = FALSE)
}

if (any(
  species_margin_data$main_included &
    (species_margin_data$start_occupied_cells < margin_cells |
       species_margin_data$end_occupied_cells < margin_cells)
)) {
  stop("A main-analysis species has fewer than 40 occupied cells.", call. = FALSE)
}

#### FIT MAIN AND SENSITIVITY MODELS ####

model_keys <- expand.grid(
  comparison = atlas_comparisons$comparison,
  distribution_class = distribution_levels,
  stringsAsFactors = FALSE
)

sensitivity_models <- do.call(
  rbind,
  lapply(seq_len(nrow(model_keys)), function(i) {
    contrast <- model_keys$comparison[[i]]
    distribution_class <- model_keys$distribution_class[[i]]

    do.call(
      rbind,
      lapply(names(variant_labels), function(variant_id) {
        model_data <- make_sensitivity_data(
          species_margin_data,
          contrast,
          distribution_class,
          variant_id
        )

        fit_margin_model(
          model_data,
          contrast,
          distribution_class,
          variant_id,
          unname(variant_labels[variant_id])
        )
      })
    )
  })
)

main_models <- sensitivity_models |>
  filter(.data$model_variant == "main")

plot_data <- species_margin_data |>
  filter(.data$main_included) |>
  mutate(
    comparison = factor(.data$comparison, levels = atlas_comparisons$comparison),
    distribution_label = factor(
      .data$distribution_label,
      levels = unname(distribution_labels)
    )
  )

model_annotations <- main_models |>
  mutate(
    comparison = factor(.data$comparison, levels = atlas_comparisons$comparison),
    distribution_label = factor(
      .data$distribution_label,
      levels = unname(distribution_labels)
    ),
    label = paste0(
      "intercept = ",
      round(.data$intercept_km, 1),
      " km\n",
      format_p_value(.data$intercept_p_value),
      "\nn = ",
      .data$n_species,
      "; R2 = ",
      round(.data$r_squared, 2)
    )
  )

#### WRITE TABLE OUTPUTS ####

write.csv(species_margin_data, species_data_path, row.names = FALSE, na = "")
write.csv(excluded_diurnal_raptors, excluded_raptors_path, row.names = FALSE, na = "")
write.csv(
  species_margin_data |> filter(!.data$main_included),
  exclusions_path,
  row.names = FALSE,
  na = ""
)
write.csv(main_models, main_model_path, row.names = FALSE, na = "")
write.csv(sensitivity_models, sensitivity_model_path, row.names = FALSE, na = "")
write.csv(plot_data, plot_data_path, row.names = FALSE, na = "")

#### PNG FIGURE ####

margin_plot <- ggplot(
  plot_data,
  aes(
    x = .data$status_change_log10,
    y = .data$margin_shift_km
  )
) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_point(alpha = 0.72, size = 1.7, color = "#2b2b2b") +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "#2166ac",
    fill = "#9ecae1",
    linewidth = 0.85
  ) +
  geom_text(
    data = model_annotations,
    aes(label = .data$label),
    x = -Inf,
    y = Inf,
    hjust = -0.05,
    vjust = 1.08,
    size = 3.1,
    lineheight = 0.95,
    inherit.aes = FALSE
  ) +
  facet_grid(distribution_label ~ comparison) +
  labs(
    title = "Range-margin shifts across adjacent Danish atlas projects",
    subtitle = paste0(
      "HMSC retained species excluding diurnal raptors; margins use the ",
      margin_cells,
      " most marginal occupied 5 km cells"
    ),
    x = "Change in status: log10(occupied cells in later atlas / earlier atlas)",
    y = "Margin shift (km; positive = northward)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey25"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = main_png_path,
  plot = margin_plot,
  width = 10.5,
  height = 7.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

#### CONSOLE SUMMARY ####

message("Raw occurrence matrix: ", nrow(y_matrix), " rows x ", ncol(y_matrix), " species")
message("HMSC retained species: ", length(hmsc_species))
message("Excluded diurnal raptors: ", nrow(excluded_diurnal_raptors))
message("Analysis species after raptor exclusion: ", length(analysis_species))
message("Full raw shared atlas cells: ", length(full_shared_sites))
message("Main-analysis species counts:")
print(
  plot_data |>
    count(.data$comparison, .data$distribution_class, name = "n_species")
)
message("Main model summary:")
print(main_models)
message("Wrote species margin data: ", species_data_path)
message("Wrote excluded diurnal raptors: ", excluded_raptors_path)
message("Wrote exclusions: ", exclusions_path)
message("Wrote main model summary: ", main_model_path)
message("Wrote sensitivity model summary: ", sensitivity_model_path)
message("Wrote main plot data: ", plot_data_path)
message("Wrote PNG figure: ", main_png_path)
