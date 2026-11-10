#!/usr/bin/env Rscript

# Supplementary diagnostic: species-level trajectories in full-model Tjur R2
# and atlas prevalence. All diurnal raptors are highlighted over non-raptor
# species that increased in prevalence from the 1970s to the 2010s.
#
# Outputs:
#   manuscript-outputs/figures/supplement/<model-pattern>-sfig-raptor-species-full-tjur-prevalence-trends.png
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig-raptor-species-full-tjur-prevalence-trends-values.csv
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig-raptor-species-full-tjur-prevalence-trends-summary.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(patchwork)
})

base_dir <- "HmscOutputs"
model_pattern <- "2026-07-27"
samples_slug <- "samples_0250_thin_100_chains_4"

trait_table_path <- file.path(
  "manuscript-outputs",
  "tables",
  "methods-audit",
  "modelled_species_traits_prevalence.csv"
)

output_figure_dir <- file.path("manuscript-outputs", "figures", "supplement")
output_table_dir <- file.path("manuscript-outputs", "tables", "supplement")

output_figure_path <- file.path(
  output_figure_dir,
  paste0(model_pattern, "-sfig-raptor-species-full-tjur-prevalence-trends.png")
)
trend_values_path <- file.path(
  output_table_dir,
  paste0(model_pattern, "-sfig-raptor-species-full-tjur-prevalence-trends-values.csv")
)
trend_summary_path <- file.path(
  output_table_dir,
  paste0(model_pattern, "-sfig-raptor-species-full-tjur-prevalence-trends-summary.csv")
)

period_lookup <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)
period_levels <- unname(period_lookup)

highlight_colours <- c(
  "Other increasing species" = "grey62",
  "Diurnal raptors" = "#b23a48"
)

extract_atlas <- function(model_folder) {
  atlas <- str_match(basename(model_folder), "Atlas([0-9]+)")[, 2]

  if (is.na(atlas) || !atlas %in% names(period_lookup)) {
    stop("Could not extract atlas number from model folder: ", model_folder, call. = FALSE)
  }

  atlas
}

load_species_order <- function(model_folder) {
  fit_path <- file.path(
    base_dir,
    model_folder,
    "Models",
    "Fitted",
    paste0("HPC_", samples_slug, ".Rdata")
  )

  if (!file.exists(fit_path)) {
    stop("Missing fitted model: ", fit_path, call. = FALSE)
  }

  fit_env <- new.env(parent = emptyenv())
  load(fit_path, envir = fit_env)

  if (!exists("fitted_model", envir = fit_env, inherits = FALSE)) {
    stop("Fitted model file does not contain fitted_model: ", fit_path, call. = FALSE)
  }

  species <- rownames(fit_env$fitted_model$posteriors$Tr)
  if (is.null(species)) {
    stop("Could not recover species order from fitted model: ", fit_path, call. = FALSE)
  }

  species
}

read_full_model_tjur_r2 <- function(model_folder) {
  atlas <- extract_atlas(model_folder)
  period <- period_lookup[[atlas]]
  period_key <- paste0("atlas", atlas)
  species <- load_species_order(model_folder)

  metric_path <- file.path(
    base_dir,
    model_folder,
    "Models",
    "Fitted",
    paste0("MF_preCV_", samples_slug, ".Rdata")
  )

  if (!file.exists(metric_path)) {
    stop("Missing full-model metric file: ", metric_path, call. = FALSE)
  }

  metric_env <- new.env(parent = emptyenv())
  load(metric_path, envir = metric_env)

  if (!exists("MF", envir = metric_env, inherits = FALSE)) {
    stop("Metric file does not contain MF: ", metric_path, call. = FALSE)
  }

  values <- metric_env$MF$TjurR2

  if (is.null(values)) {
    stop("Missing full-model TjurR2 in: ", metric_path, call. = FALSE)
  }
  if (length(values) != length(species)) {
    stop(
      "Metric/species length mismatch for ",
      model_folder,
      ": ",
      length(values),
      " values vs ",
      length(species),
      " species",
      call. = FALSE
    )
  }

  tibble(
    species = species,
    period = period,
    period_key = period_key,
    period_index = as.integer(atlas),
    full_model_tjur_r2 = as.numeric(values)
  )
}

if (!file.exists(trait_table_path)) {
  stop("Missing species trait/prevalence table: ", trait_table_path, call. = FALSE)
}

species_traits <- read.csv(trait_table_path, check.names = FALSE)
required_trait_columns <- c(
  "species",
  "english_DOF",
  "foraging_guild_consensus",
  "atlas1_prevalence",
  "atlas2_prevalence",
  "atlas3_prevalence"
)
missing_trait_columns <- setdiff(required_trait_columns, colnames(species_traits))
if (length(missing_trait_columns) > 0) {
  stop("Missing trait/prevalence columns: ", paste(missing_trait_columns, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(species_traits$species) > 0) {
  stop("Species names are not unique in: ", trait_table_path, call. = FALSE)
}

prevalence_long <- species_traits |>
  select(
    "species",
    "english_DOF",
    "foraging_guild_consensus",
    "atlas1_prevalence",
    "atlas2_prevalence",
    "atlas3_prevalence"
  ) |>
  mutate(
    prevalence_change_2010s_minus_1970s = .data$atlas3_prevalence - .data$atlas1_prevalence,
    prevalence_increased_2010s_minus_1970s = .data$prevalence_change_2010s_minus_1970s > 0
  ) |>
  pivot_longer(
    cols = starts_with("atlas") & ends_with("prevalence"),
    names_to = "period_key",
    values_to = "prevalence"
  ) |>
  mutate(
    period = period_lookup[str_remove(str_remove(.data$period_key, "_prevalence$"), "^atlas")],
    period_key = str_remove(.data$period_key, "_prevalence$"),
    period_index = as.integer(str_remove(.data$period_key, "^atlas"))
  )

model_folders <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE) |>
  keep(~ str_detect(.x, fixed(model_pattern))) |>
  keep(~ str_detect(.x, "Atlas[0-9]+")) |>
  sort()

if (length(model_folders) == 0) {
  stop("No model folders found for pattern: ", model_pattern, call. = FALSE)
}

trend_values_all <- map_dfr(model_folders, read_full_model_tjur_r2) |>
  left_join(
    prevalence_long,
    by = c("species", "period", "period_key", "period_index")
  ) |>
  mutate(
    period = factor(.data$period, levels = period_levels),
    is_raptor = .data$foraging_guild_consensus == "Diurnal raptors"
  )

if (any(is.na(trend_values_all$prevalence))) {
  missing_species <- trend_values_all |>
    filter(is.na(.data$prevalence)) |>
    distinct(.data$species) |>
    pull(.data$species)
  stop("Missing prevalence values for: ", paste(missing_species, collapse = ", "), call. = FALSE)
}

trend_values <- trend_values_all |>
  filter(.data$is_raptor | .data$prevalence_increased_2010s_minus_1970s) |>
  mutate(
    highlight = factor(
      if_else(.data$is_raptor, "Diurnal raptors", "Other increasing species"),
      levels = names(highlight_colours)
    )
  )

wide_values <- trend_values |>
  select(
    "species",
    "english_DOF",
    "foraging_guild_consensus",
    "is_raptor",
    "prevalence_change_2010s_minus_1970s",
    "period_key",
    "full_model_tjur_r2",
    "prevalence"
  ) |>
  pivot_wider(
    names_from = "period_key",
    values_from = c("full_model_tjur_r2", "prevalence")
  )

species_trends <- trend_values |>
  select(
    "species",
    "english_DOF",
    "foraging_guild_consensus",
    "is_raptor",
    "period_index",
    "full_model_tjur_r2",
    "prevalence"
  ) |>
  pivot_longer(
    cols = c("full_model_tjur_r2", "prevalence"),
    names_to = "measure",
    values_to = "value"
  ) |>
  group_by(
    .data$species,
    .data$english_DOF,
    .data$foraging_guild_consensus,
    .data$is_raptor,
    .data$measure
  ) |>
  summarise(
    slope_per_atlas_step = coef(lm(value ~ period_index))[["period_index"]],
    value_1970s = value[period_index == 1][[1]],
    value_1990s = value[period_index == 2][[1]],
    value_2010s = value[period_index == 3][[1]],
    change_2010s_minus_1970s = value_2010s - value_1970s,
    increased = change_2010s_minus_1970s > 0,
    .groups = "drop"
  )

compare_raptor_trends <- function(data) {
  raptor_values <- data$change_2010s_minus_1970s[data$is_raptor]
  other_values <- data$change_2010s_minus_1970s[!data$is_raptor]
  slope_model <- lm(slope_per_atlas_step ~ is_raptor, data = data)
  slope_coefs <- summary(slope_model)$coefficients

  tibble(
    n_raptors = sum(data$is_raptor),
    n_other_increasing_species = sum(!data$is_raptor),
    n_raptors_with_increasing_measure = sum(data$increased[data$is_raptor]),
    n_other_species_with_increasing_measure = sum(data$increased[!data$is_raptor]),
    median_change_raptors = median(raptor_values, na.rm = TRUE),
    median_change_other_increasing_species = median(other_values, na.rm = TRUE),
    median_change_difference_raptors_minus_other = median(raptor_values, na.rm = TRUE) -
      median(other_values, na.rm = TRUE),
    mean_slope_raptors = mean(data$slope_per_atlas_step[data$is_raptor], na.rm = TRUE),
    mean_slope_other_increasing_species = mean(data$slope_per_atlas_step[!data$is_raptor], na.rm = TRUE),
    lm_slope_difference_raptors_minus_other = slope_coefs["is_raptorTRUE", "Estimate"],
    lm_slope_difference_p_value = slope_coefs["is_raptorTRUE", "Pr(>|t|)"],
    wilcoxon_change_p_value = wilcox.test(
      change_2010s_minus_1970s ~ is_raptor,
      data = data,
      exact = FALSE
    )$p.value
  )
}

trend_summary <- species_trends |>
  group_by(.data$measure) |>
  group_modify(~ compare_raptor_trends(.x)) |>
  ungroup() |>
  mutate(
    measure = recode(
      .data$measure,
      "full_model_tjur_r2" = "Full-model Tjur R2",
      "prevalence" = "Prevalence"
    )
  )

base_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.title = element_text(colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA)
  )

make_trend_plot <- function(data, value_col, y_label, y_scale) {
  ggplot(data, aes(x = .data$period, y = .data[[value_col]], group = .data$species)) +
    geom_line(
      data = filter(data, !.data$is_raptor),
      aes(colour = "Other increasing species"),
      linewidth = 0.28,
      alpha = 0.28,
      lineend = "round"
    ) +
    geom_point(
      data = filter(data, !.data$is_raptor),
      aes(colour = "Other increasing species"),
      size = 0.65,
      alpha = 0.25
    ) +
    geom_line(
      data = filter(data, .data$is_raptor),
      aes(colour = "Diurnal raptors"),
      linewidth = 0.9,
      alpha = 0.95,
      lineend = "round"
    ) +
    geom_point(
      data = filter(data, .data$is_raptor),
      aes(colour = "Diurnal raptors"),
      size = 1.8,
      alpha = 0.98
    ) +
    scale_colour_manual(
      values = highlight_colours,
      breaks = names(highlight_colours)
    ) +
    y_scale +
    labs(x = NULL, y = y_label) +
    base_theme
}

tjur_plot <- make_trend_plot(
  trend_values,
  value_col = "full_model_tjur_r2",
  y_label = "Full-model Tjur R2",
  y_scale = scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.08)))
)

prevalence_plot <- make_trend_plot(
  trend_values,
  value_col = "prevalence",
  y_label = "Prevalence in retained atlas grid cells",
  y_scale = scale_y_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0.02, 0.04))
  )
)

species_trend_plot <- tjur_plot + prevalence_plot +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  output_figure_path,
  species_trend_plot,
  width = 9.2,
  height = 4.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

write.csv(trend_values, trend_values_path, row.names = FALSE)
write.csv(trend_summary, trend_summary_path, row.names = FALSE)

print(trend_summary)
message("Wrote: ", output_figure_path)
message("Wrote: ", trend_values_path)
message("Wrote: ", trend_summary_path)
