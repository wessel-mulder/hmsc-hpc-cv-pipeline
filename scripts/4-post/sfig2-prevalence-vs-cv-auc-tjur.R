#!/usr/bin/env Rscript

# Supplementary figure: species prevalence vs cross-validated AUC and Tjur R2.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(scales)
  library(stringr)
  library(tibble)
})

base_dir <- "HmscOutputs"
model_pattern <- "2026-07-27"
samples_slug <- "samples_0250_thin_100_chains_4"
nfolds <- 10
output_dir <- file.path("manuscript-outputs", "figures", "supplement")
output_path <- file.path(output_dir, paste0(model_pattern, "-sfig2-prevalence-vs-cv-auc-tjur.png"))

period_lookup <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)
period_levels <- unname(period_lookup)

metric_lookup <- tibble(
  metric = c("AUC", "TjurR2"),
  metric_label = c("AUC", "Tjur R2")
)

period_cols <- c(
  "1970s" = "#4c78a8",
  "1990s" = "#f58518",
  "2010s" = "#54a24b"
)

extract_atlas <- function(model_folder) {
  atlas <- str_match(basename(model_folder), "Atlas([0-9]+)")[, 2]

  if (is.na(atlas) || !atlas %in% names(period_lookup)) {
    stop("Could not extract atlas number from model folder: ", model_folder, call. = FALSE)
  }

  atlas
}

load_fitted_model <- function(model_folder) {
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

  fit_env$fitted_model
}

read_prevalence <- function(model_folder) {
  atlas <- extract_atlas(model_folder)
  period <- period_lookup[[atlas]]
  fitted_model <- load_fitted_model(model_folder)
  y_matrix <- fitted_model$posteriors$Y

  if (is.null(y_matrix) || is.null(dim(y_matrix)) || is.null(colnames(y_matrix))) {
    stop("Could not recover fitted-model response matrix for: ", model_folder, call. = FALSE)
  }

  tibble(
    model_folder = model_folder,
    atlas = atlas,
    period = period,
    species = colnames(y_matrix),
    n_sites = nrow(y_matrix),
    n_presence = colSums(y_matrix, na.rm = TRUE),
    prevalence = colMeans(y_matrix, na.rm = TRUE)
  )
}

get_species_order <- function(model_folder) {
  fitted_model <- load_fitted_model(model_folder)
  species <- rownames(fitted_model$posteriors$Tr)

  if (is.null(species)) {
    stop("Could not recover species order from fitted model: ", model_folder, call. = FALSE)
  }

  species
}

read_cv_accuracy <- function(model_folder) {
  atlas <- extract_atlas(model_folder)
  period <- period_lookup[[atlas]]
  accuracy_path <- file.path(
    base_dir,
    model_folder,
    "Models",
    "Fitted",
    paste0("MF_", samples_slug, "_nfolds_", nfolds, ".rdata")
  )

  if (!file.exists(accuracy_path)) {
    stop("Missing model-fit metrics: ", accuracy_path, call. = FALSE)
  }

  metric_env <- new.env(parent = emptyenv())
  load(accuracy_path, envir = metric_env)

  if (!exists("MFCV", envir = metric_env, inherits = FALSE)) {
    stop("Metric file does not contain MFCV: ", accuracy_path, call. = FALSE)
  }

  species <- get_species_order(model_folder)
  map_dfr(metric_lookup$metric, function(metric_name) {
    values <- metric_env$MFCV[[metric_name]]

    if (is.null(values)) {
      stop("Missing MFCV metric ", metric_name, " in ", model_folder, call. = FALSE)
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
      metric = metric_name,
      value = as.numeric(values)
    )
  })
}

polynomial_lm_label <- function(data) {
  complete <- data |>
    filter(!is.na(.data$prevalence), !is.na(.data$value))

  if (nrow(complete) < 4 || length(unique(complete$prevalence)) < 3) {
    return("quadratic lm\nR2 = NA")
  }

  model <- lm(value ~ poly(prevalence, 2, raw = TRUE), data = complete)

  paste0(
    "quadratic lm\nR2 = ",
    sprintf("%.2f", summary(model)$r.squared)
  )
}

model_folders <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE) |>
  keep(~ str_detect(.x, fixed(model_pattern))) |>
  keep(~ str_detect(.x, "Atlas[0-9]+")) |>
  sort()

if (length(model_folders) == 0) {
  stop("No model folders found for pattern: ", model_pattern, call. = FALSE)
}

plot_df <- left_join(
  map_dfr(model_folders, read_cv_accuracy),
  map_dfr(model_folders, read_prevalence),
  by = c("model_folder", "atlas", "period", "species")
) |>
  left_join(metric_lookup, by = "metric") |>
  mutate(
    period = factor(.data$period, levels = period_levels),
    metric_label = factor(.data$metric_label, levels = metric_lookup$metric_label)
  )

annotation_df <- plot_df |>
  group_by(.data$metric_label) |>
  summarise(label = polynomial_lm_label(pick(everything())), .groups = "drop")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

prevalence_accuracy_plot <- ggplot(
  plot_df,
  aes(x = .data$prevalence, y = .data$value)
) +
  geom_point(aes(colour = .data$period), size = 1.35, alpha = 0.58) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ poly(x, 2, raw = TRUE),
    se = TRUE,
    colour = "grey15",
    fill = "grey75",
    linewidth = 0.55,
    alpha = 0.35
  ) +
  geom_text(
    data = annotation_df,
    aes(label = .data$label),
    x = Inf,
    y = Inf,
    hjust = 1.08,
    vjust = 1.18,
    size = 3.1,
    colour = "grey20",
    inherit.aes = FALSE
  ) +
  facet_wrap(vars(.data$metric_label), ncol = 2, scales = "free_y") +
  scale_x_continuous(
    labels = label_percent(accuracy = 1),
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),
    expand = expansion(mult = c(0.015, 0.045))
  ) +
  scale_colour_manual(values = period_cols, name = "Atlas period") +
  labs(
    x = "Species prevalence in retained atlas grid cells",
    y = "10-fold cross-validated metric value"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.25),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text = element_text(face = "bold", colour = "black"),
    plot.title = element_text(face = "bold", size = 12),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", colour = NA)
  )

ggsave(
  output_path,
  prevalence_accuracy_plot,
  width = 8.2,
  height = 4.6,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Wrote: ", output_path)
