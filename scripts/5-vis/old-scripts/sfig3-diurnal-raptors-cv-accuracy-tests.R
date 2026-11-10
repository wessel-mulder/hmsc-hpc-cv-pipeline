#!/usr/bin/env Rscript

# Panel-wise permutation tests for Supplementary Figure 3.
#
# Tests whether the observed median MFCV accuracy difference between diurnal
# raptors and other species is larger than expected under random guild labels
# within each metric-period panel.
#
# Output:
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig3-diurnal-raptors-cv-accuracy-permutation-tests.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
})

base_dir <- "HmscOutputs"
model_pattern <- "2026-07-27"
samples_slug <- "samples_0250_thin_100_chains_4"
nfolds <- 10
evaluation <- "MFCV"
guild_column <- "foraging_guild_consensusDiurnal raptors"
n_permutations <- 99999
random_seed <- 42

output_dir <- file.path("manuscript-outputs", "tables", "supplement")
output_path <- file.path(
  output_dir,
  paste0(model_pattern, "-sfig3-diurnal-raptors-cv-accuracy-permutation-tests.csv")
)

period_lookup <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)
period_levels <- unname(period_lookup)

metric_lookup <- tribble(
  ~metric, ~metric_label,
  "AUC", "AUC",
  "TjurR2", "Tjur R2",
  "RMSE", "RMSE"
)

comparison_levels <- c("Diurnal raptors", "Other species")

extract_atlas <- function(model_folder) {
  atlas <- str_match(basename(model_folder), "Atlas([0-9]+)")[, 2]

  if (is.na(atlas) || !atlas %in% names(period_lookup)) {
    stop("Could not extract atlas number from model folder: ", model_folder, call. = FALSE)
  }

  atlas
}

read_model_accuracy <- function(model_folder) {
  atlas <- extract_atlas(model_folder)
  period <- period_lookup[[atlas]]

  fit_path <- file.path(
    base_dir,
    model_folder,
    "Models",
    "Fitted",
    paste0("HPC_", samples_slug, ".Rdata")
  )
  accuracy_path <- file.path(
    base_dir,
    model_folder,
    "Models",
    "Fitted",
    paste0("MF_", samples_slug, "_nfolds_", nfolds, ".rdata")
  )

  if (!file.exists(fit_path)) {
    stop("Missing fitted model: ", fit_path, call. = FALSE)
  }
  if (!file.exists(accuracy_path)) {
    stop("Missing model-fit metrics: ", accuracy_path, call. = FALSE)
  }

  fit_env <- new.env(parent = emptyenv())
  load(fit_path, envir = fit_env)
  model <- fit_env$fitted_model$posteriors

  if (is.null(model$Tr) || !guild_column %in% colnames(model$Tr)) {
    stop("Missing trait column: ", guild_column, call. = FALSE)
  }

  accuracy_env <- new.env(parent = emptyenv())
  load(accuracy_path, envir = accuracy_env)

  if (!exists(evaluation, envir = accuracy_env, inherits = FALSE)) {
    stop("Metric file does not contain ", evaluation, ": ", accuracy_path, call. = FALSE)
  }

  metric_set <- accuracy_env[[evaluation]]
  species <- rownames(model$Tr)
  comparison <- if_else(
    model$Tr[, guild_column] == 1,
    "Diurnal raptors",
    "Other species"
  )

  map_dfr(metric_lookup$metric, function(metric_name) {
    values <- metric_set[[metric_name]]

    if (is.null(values)) {
      stop("Missing ", evaluation, " metric ", metric_name, " in ", model_folder, call. = FALSE)
    }
    if (length(values) != length(species)) {
      stop(
        "Metric/species length mismatch for ",
        model_folder,
        " / ",
        metric_name,
        ": ",
        length(values),
        " values vs ",
        length(species),
        " species",
        call. = FALSE
      )
    }

    tibble(
      model_folder = model_folder,
      atlas = atlas,
      period = period,
      species = species,
      comparison = comparison,
      metric = metric_name,
      value = as.numeric(values)
    )
  })
}

probability_superiority <- function(raptor_values, other_values) {
  comparison_matrix <- outer(raptor_values, other_values, FUN = "-")
  mean(comparison_matrix > 0, na.rm = TRUE) +
    0.5 * mean(comparison_matrix == 0, na.rm = TRUE)
}

panel_permutation_test <- function(data, n_perm = n_permutations) {
  complete <- data |>
    filter(!is.na(.data$value), !is.na(.data$comparison))

  raptor_values <- complete |>
    filter(.data$comparison == "Diurnal raptors") |>
    pull(.data$value)
  other_values <- complete |>
    filter(.data$comparison == "Other species") |>
    pull(.data$value)

  n_raptors <- length(raptor_values)
  n_other <- length(other_values)
  all_values <- c(raptor_values, other_values)
  n_total <- length(all_values)

  if (n_raptors == 0 || n_other == 0) {
    stop("Each panel must contain both diurnal raptors and other species.", call. = FALSE)
  }

  observed_difference <- median(raptor_values) - median(other_values)
  null_differences <- replicate(n_perm, {
    pseudo_raptors <- sample.int(n_total, size = n_raptors, replace = FALSE)
    median(all_values[pseudo_raptors]) - median(all_values[-pseudo_raptors])
  })

  p_two_sided <- (sum(abs(null_differences) >= abs(observed_difference)) + 1) /
    (n_perm + 1)

  superiority <- probability_superiority(raptor_values, other_values)

  tibble(
    n_raptors = n_raptors,
    n_other_species = n_other,
    median_raptors = median(raptor_values),
    median_other_species = median(other_values),
    median_difference_raptors_minus_other = observed_difference,
    probability_raptor_greater_than_other = superiority,
    cliffs_delta = 2 * superiority - 1,
    n_permutations = n_perm,
    p_permutation_two_sided = p_two_sided
  )
}

set.seed(random_seed)

model_folders <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE) |>
  keep(~ str_detect(.x, fixed(model_pattern))) |>
  keep(~ str_detect(.x, "Atlas[0-9]+")) |>
  sort()

if (length(model_folders) == 0) {
  stop("No model folders found for pattern: ", model_pattern, call. = FALSE)
}

accuracy_df <- map_dfr(model_folders, read_model_accuracy) |>
  left_join(metric_lookup, by = "metric") |>
  mutate(
    period = factor(.data$period, levels = period_levels),
    metric_label = factor(.data$metric_label, levels = metric_lookup$metric_label),
    comparison = factor(.data$comparison, levels = comparison_levels)
  )

test_results <- accuracy_df |>
  group_by(.data$metric, .data$metric_label, .data$atlas, .data$period) |>
  group_modify(~ panel_permutation_test(.x)) |>
  ungroup() |>
  arrange(.data$metric_label, .data$period) |>
  mutate(
    p_fdr_bh = p.adjust(.data$p_permutation_two_sided, method = "BH"),
    random_seed = random_seed
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(test_results, output_path, row.names = FALSE)

print(test_results, n = Inf)
message("Wrote: ", output_path)
