#!/usr/bin/env Rscript

# Supplementary figure: temporal changes in cross-validated accuracy for
# diurnal raptors and four common guilds with increasing prevalence.
#
# Outputs:
#   manuscript-outputs/figures/supplement/<model-pattern>-sfig4-guild-cv-accuracy-temporal-change.png
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig4-guild-cv-accuracy-period-summary.csv
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig4-guild-cv-accuracy-temporal-change-summary.csv
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig4-guild-cv-accuracy-temporal-change-species-values.csv
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig4-guild-cv-accuracy-temporal-change-raptor-comparisons.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

base_dir <- "HmscOutputs"
model_pattern <- "2026-07-27"
samples_slug <- "samples_0250_thin_100_chains_4"
nfolds <- 10
evaluation <- "MFCV"
trait_table_path <- file.path(
  "manuscript-outputs",
  "tables",
  "methods-audit",
  "modelled_species_traits_prevalence.csv"
)
n_permutations <- 99999
random_seed <- 42

selected_guilds <- c(
  "Diurnal raptors",
  "Tit-like birds",
  "Passerine seedeaters",
  "Low flycatching feeders",
  "Foliage gleaners"
)

output_figure_dir <- file.path("manuscript-outputs", "figures", "supplement")
output_table_dir <- file.path("manuscript-outputs", "tables", "supplement")

output_figure_path <- file.path(
  output_figure_dir,
  paste0(model_pattern, "-sfig4-guild-cv-accuracy-temporal-change.png")
)
period_summary_path <- file.path(
  output_table_dir,
  paste0(model_pattern, "-sfig4-guild-cv-accuracy-period-summary.csv")
)
change_summary_path <- file.path(
  output_table_dir,
  paste0(model_pattern, "-sfig4-guild-cv-accuracy-temporal-change-summary.csv")
)
species_change_path <- file.path(
  output_table_dir,
  paste0(model_pattern, "-sfig4-guild-cv-accuracy-temporal-change-species-values.csv")
)
raptor_comparison_path <- file.path(
  output_table_dir,
  paste0(model_pattern, "-sfig4-guild-cv-accuracy-temporal-change-raptor-comparisons.csv")
)

period_lookup <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)
period_levels <- unname(period_lookup)

metric_lookup <- tribble(
  ~metric, ~metric_label, ~improvement_direction,
  "AUC", "AUC", 1,
  "TjurR2", "Tjur R2", 1,
  "RMSE", "RMSE", -1
)

contrast_levels <- c(
  "1990s - 1970s",
  "2010s - 1990s",
  "2010s - 1970s"
)

guild_cols <- c(
  "Diurnal raptors" = "#b23a48",
  "Tit-like birds" = "#26847d",
  "Passerine seedeaters" = "#b9832d",
  "Low flycatching feeders" = "#4b6fb3",
  "Foliage gleaners" = "#4f8f3a"
)

extract_atlas <- function(model_folder) {
  atlas <- str_match(basename(model_folder), "Atlas([0-9]+)")[, 2]

  if (is.na(atlas) || !atlas %in% names(period_lookup)) {
    stop("Could not extract atlas number from model folder: ", model_folder, call. = FALSE)
  }

  atlas
}

load_species_guild_lookup <- function(path) {
  if (!file.exists(path)) {
    stop("Missing species trait/prevalence table: ", path, call. = FALSE)
  }

  traits <- read.csv(path, check.names = FALSE)
  required_columns <- c("species", "foraging_guild_consensus")
  missing_columns <- setdiff(required_columns, colnames(traits))

  if (length(missing_columns) > 0) {
    stop("Missing columns in trait/prevalence table: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  if (anyDuplicated(traits$species) > 0) {
    stop("Species names are not unique in trait/prevalence table: ", path, call. = FALSE)
  }

  traits |>
    transmute(
      species = .data$species,
      guild = .data$foraging_guild_consensus
    )
}

species_guild_lookup <- load_species_guild_lookup(trait_table_path)

assign_selected_guild <- function(species) {
  missing_species <- setdiff(species, species_guild_lookup$species)

  if (length(missing_species) > 0) {
    stop(
      "Species missing from trait/prevalence table: ",
      paste(missing_species, collapse = ", "),
      call. = FALSE
    )
  }

  guild <- species_guild_lookup$guild[match(species, species_guild_lookup$species)]
  if_else(guild %in% selected_guilds, guild, NA_character_)
}

read_model_accuracy <- function(model_folder) {
  atlas <- extract_atlas(model_folder)
  period <- period_lookup[[atlas]]
  period_key <- paste0("atlas", atlas)

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

  if (is.null(model$Tr)) {
    stop("Missing trait matrix in fitted model: ", fit_path, call. = FALSE)
  }

  accuracy_env <- new.env(parent = emptyenv())
  load(accuracy_path, envir = accuracy_env)

  if (!exists(evaluation, envir = accuracy_env, inherits = FALSE)) {
    stop("Metric file does not contain ", evaluation, ": ", accuracy_path, call. = FALSE)
  }

  metric_set <- accuracy_env[[evaluation]]
  species <- rownames(model$Tr)
  guild <- assign_selected_guild(species)

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
      period_key = period_key,
      period = period,
      species = species,
      guild = guild,
      metric = metric_name,
      value = as.numeric(values)
    )
  })
}

probability_superiority <- function(target_values, comparison_values) {
  comparison_matrix <- outer(target_values, comparison_values, FUN = "-")
  mean(comparison_matrix > 0, na.rm = TRUE) +
    0.5 * mean(comparison_matrix == 0, na.rm = TRUE)
}

change_permutation_test <- function(data, comparison_guild, n_perm = n_permutations) {
  complete <- data |>
    filter(!is.na(.data$change), .data$guild %in% c("Diurnal raptors", comparison_guild))

  raptor_values <- complete |>
    filter(.data$guild == "Diurnal raptors") |>
    pull(.data$change)
  comparison_values <- complete |>
    filter(.data$guild == comparison_guild) |>
    pull(.data$change)

  n_raptors <- length(raptor_values)
  n_comparison <- length(comparison_values)
  all_values <- c(raptor_values, comparison_values)
  n_total <- length(all_values)

  if (n_raptors == 0 || n_comparison == 0) {
    stop("Each test must contain diurnal raptors and the comparison guild.", call. = FALSE)
  }

  observed_difference <- median(raptor_values) - median(comparison_values)
  null_differences <- replicate(n_perm, {
    pseudo_raptors <- sample.int(n_total, size = n_raptors, replace = FALSE)
    median(all_values[pseudo_raptors]) - median(all_values[-pseudo_raptors])
  })

  p_two_sided <- (sum(abs(null_differences) >= abs(observed_difference)) + 1) /
    (n_perm + 1)

  superiority <- probability_superiority(raptor_values, comparison_values)

  tibble(
    comparison_guild = comparison_guild,
    n_raptors = n_raptors,
    n_comparison_guild = n_comparison,
    median_change_raptors = median(raptor_values),
    median_change_comparison_guild = median(comparison_values),
    median_change_difference_raptors_minus_comparison = observed_difference,
    probability_raptor_change_greater_than_comparison = superiority,
    cliffs_delta = 2 * superiority - 1,
    n_permutations = n_perm,
    p_permutation_two_sided = p_two_sided
  )
}

model_folders <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE) |>
  keep(~ str_detect(.x, fixed(model_pattern))) |>
  keep(~ str_detect(.x, "Atlas[0-9]+")) |>
  sort()

if (length(model_folders) == 0) {
  stop("No model folders found for pattern: ", model_pattern, call. = FALSE)
}

accuracy_df <- map_dfr(model_folders, read_model_accuracy) |>
  filter(.data$guild %in% selected_guilds) |>
  left_join(metric_lookup, by = "metric") |>
  mutate(
    guild = factor(.data$guild, levels = selected_guilds),
    period = factor(.data$period, levels = period_levels),
    metric_label = factor(.data$metric_label, levels = metric_lookup$metric_label)
  )

period_summary <- accuracy_df |>
  group_by(.data$guild, .data$metric, .data$metric_label, .data$period) |>
  summarise(
    n_species = sum(!is.na(.data$value)),
    mean_metric = mean(.data$value, na.rm = TRUE),
    median_metric = median(.data$value, na.rm = TRUE),
    q25_metric = quantile(.data$value, 0.25, na.rm = TRUE),
    q75_metric = quantile(.data$value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

change_df <- accuracy_df |>
  select(
    "species",
    "guild",
    "metric",
    "metric_label",
    "improvement_direction",
    "period_key",
    "value"
  ) |>
  pivot_wider(names_from = "period_key", values_from = "value") |>
  mutate(
    `1990s - 1970s` = .data$atlas2 - .data$atlas1,
    `2010s - 1990s` = .data$atlas3 - .data$atlas2,
    `2010s - 1970s` = .data$atlas3 - .data$atlas1
  ) |>
  pivot_longer(
    cols = all_of(contrast_levels),
    names_to = "contrast",
    values_to = "change"
  ) |>
  mutate(
    contrast = factor(.data$contrast, levels = contrast_levels),
    improved = case_when(
      .data$improvement_direction > 0 ~ .data$change > 0,
      .data$improvement_direction < 0 ~ .data$change < 0,
      TRUE ~ NA
    )
  )

change_summary <- change_df |>
  group_by(.data$guild, .data$metric, .data$metric_label, .data$contrast) |>
  summarise(
    n_species = sum(!is.na(.data$change)),
    mean_change = mean(.data$change, na.rm = TRUE),
    median_change = median(.data$change, na.rm = TRUE),
    q25_change = quantile(.data$change, 0.25, na.rm = TRUE),
    q75_change = quantile(.data$change, 0.75, na.rm = TRUE),
    prop_species_with_improved_accuracy = mean(.data$improved, na.rm = TRUE),
    .groups = "drop"
  )

set.seed(random_seed)

comparison_guilds <- setdiff(selected_guilds, "Diurnal raptors")

raptor_comparisons <- change_df |>
  group_by(.data$metric, .data$metric_label, .data$contrast) |>
  group_modify(~ map_dfr(comparison_guilds, \(guild_name) {
    change_permutation_test(.x, comparison_guild = guild_name)
  })) |>
  ungroup() |>
  mutate(
    p_fdr_bh = p.adjust(.data$p_permutation_two_sided, method = "BH"),
    random_seed = random_seed
  ) |>
  arrange(.data$metric_label, .data$contrast, .data$comparison_guild)

guild_labels <- accuracy_df |>
  distinct(.data$guild, .data$species) |>
  count(.data$guild, name = "n_species") |>
  mutate(label = paste0(str_wrap(as.character(.data$guild), width = 12), "\n(n = ", .data$n_species, ")")) |>
  select("guild", "label") |>
  deframe()

summary_points <- change_summary |>
  mutate(
    guild = factor(.data$guild, levels = selected_guilds),
    metric_label = factor(.data$metric_label, levels = metric_lookup$metric_label),
    contrast = factor(.data$contrast, levels = contrast_levels)
  )

temporal_change_plot <- ggplot(
  change_df,
  aes(x = .data$guild, y = .data$change)
) +
  geom_hline(yintercept = 0, linewidth = 0.28, colour = "grey45") +
  geom_boxplot(
    aes(fill = .data$guild),
    width = 0.54,
    outlier.shape = NA,
    alpha = 0.62,
    colour = "grey20",
    linewidth = 0.30,
    na.rm = TRUE
  ) +
  geom_point(
    aes(fill = .data$guild),
    position = position_jitter(width = 0.09, height = 0, seed = 42),
    shape = 21,
    size = 1.35,
    stroke = 0.14,
    colour = "grey18",
    alpha = 0.68,
    na.rm = TRUE
  ) +
  geom_point(
    data = summary_points,
    aes(x = .data$guild, y = .data$median_change),
    inherit.aes = FALSE,
    shape = 23,
    size = 2.25,
    stroke = 0.35,
    fill = "white",
    colour = "black"
  ) +
  facet_grid(rows = vars(.data$metric_label), cols = vars(.data$contrast), scales = "free_y") +
  scale_x_discrete(labels = guild_labels) +
  scale_fill_manual(values = guild_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  labs(
    x = NULL,
    y = "Change in 10-fold cross-validated metric value"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(0.62, "lines"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text = element_text(face = "bold", colour = "black", size = 9),
    axis.text.x = element_text(colour = "black", size = 7, lineheight = 0.92),
    axis.text.y = element_text(colour = "black", size = 8),
    axis.title.y = element_text(size = 9),
    legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA)
  )

dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  output_figure_path,
  temporal_change_plot,
  width = 8.8,
  height = 7.4,
  units = "in",
  dpi = 300,
  bg = "white"
)

write.csv(period_summary, period_summary_path, row.names = FALSE)
write.csv(change_summary, change_summary_path, row.names = FALSE)
write.csv(change_df, species_change_path, row.names = FALSE)
write.csv(raptor_comparisons, raptor_comparison_path, row.names = FALSE)

print(change_summary, n = Inf)
print(raptor_comparisons, n = Inf)
message("Wrote: ", output_figure_path)
message("Wrote: ", period_summary_path)
message("Wrote: ", change_summary_path)
message("Wrote: ", species_change_path)
message("Wrote: ", raptor_comparison_path)
