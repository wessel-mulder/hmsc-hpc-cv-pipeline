#!/usr/bin/env Rscript

# Supplementary figure: cross-validated accuracy for diurnal raptors vs other species.
#
# Output:
#   manuscript-outputs/figures/supplement/<model-pattern>-sfig3-diurnal-raptors-cv-accuracy.png

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
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

output_dir <- file.path("manuscript-outputs", "figures", "supplement")
output_path <- file.path(
  output_dir,
  paste0(model_pattern, "-sfig3-diurnal-raptors-cv-accuracy.png")
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
comparison_cols <- c(
  "Diurnal raptors" = "#b23a48",
  "Other species" = "grey72"
)

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

comparison_labels <- accuracy_df |>
  distinct(.data$comparison, .data$species) |>
  count(.data$comparison, name = "n_species") |>
  mutate(label = paste0(.data$comparison, "\n(n = ", .data$n_species, ")")) |>
  select("comparison", "label") |>
  deframe()

summary_df <- accuracy_df |>
  group_by(.data$period, .data$metric_label, .data$comparison) |>
  summarise(
    n_species = sum(!is.na(.data$value)),
    median = median(.data$value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    accuracy_df |>
      group_by(.data$period, .data$metric_label) |>
      summarise(label_y = max(.data$value, na.rm = TRUE), .groups = "drop"),
    by = c("period", "metric_label")
  ) |>
  mutate(median_label = paste0("median = ", sprintf("%.2f", .data$median)))

diurnal_raptor_accuracy_plot <- ggplot(
  accuracy_df,
  aes(x = .data$comparison, y = .data$value)
) +
  geom_violin(
    aes(fill = .data$comparison),
    width = 0.82,
    trim = FALSE,
    alpha = 0.50,
    colour = NA,
    na.rm = TRUE
  ) +
  geom_point(
    aes(fill = .data$comparison),
    position = position_jitter(width = 0.055, height = 0, seed = 42),
    shape = 21,
    size = 0.85,
    stroke = 0.12,
    colour = "grey25",
    alpha = 0.55,
    na.rm = TRUE
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    fill = "white",
    colour = "grey20",
    linewidth = 0.28,
    na.rm = TRUE
  ) +
  geom_point(
    data = summary_df,
    aes(x = .data$comparison, y = .data$median),
    inherit.aes = FALSE,
    shape = 23,
    size = 2.4,
    stroke = 0.35,
    fill = "white",
    colour = "black"
  ) +
  geom_text(
    data = summary_df,
    aes(x = .data$comparison, y = .data$label_y, label = .data$median_label),
    inherit.aes = FALSE,
    vjust = -0.55,
    size = 2.45,
    colour = "grey15"
  ) +
  facet_grid(rows = vars(.data$metric_label), cols = vars(.data$period), scales = "free_y") +
  scale_x_discrete(labels = comparison_labels) +
  scale_fill_manual(values = comparison_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.20))) +
  labs(
    x = NULL,
    y = "10-fold cross-validated metric value"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(0.55, "lines"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text = element_text(face = "bold", colour = "black", size = 10),
    axis.text.x = element_text(colour = "black", size = 8),
    axis.text.y = element_text(colour = "black", size = 8),
    axis.title.y = element_text(size = 9),
    legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA)
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  output_path,
  diurnal_raptor_accuracy_plot,
  width = 7.6,
  height = 7.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Wrote: ", output_path)
