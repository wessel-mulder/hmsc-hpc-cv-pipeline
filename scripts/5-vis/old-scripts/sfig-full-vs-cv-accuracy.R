#!/usr/bin/env Rscript

# Supplementary figure: full-model and cross-validated model accuracy.
#
# Usage:
#   Rscript scripts/4-data-visualisation/sfig-full-vs-cv-accuracy.R
#   Rscript scripts/4-data-visualisation/sfig-full-vs-cv-accuracy.R --model-pattern 2026-07-27
#
# Inputs:
#   HmscOutputs/<model>/Models/Fitted/MF_samples_0250_thin_100_chains_4_nfolds_10.rdata
#   HmscOutputs/<model>/Models/Fitted/MF_preCV_samples_0250_thin_100_chains_4.Rdata
#   HmscOutputs/<model>/Models/Fitted/HPC_samples_0250_thin_100_chains_4.Rdata
#
# Output:
#   manuscript-outputs/figures/supplement/<model-pattern>-sfig-full-vs-cv-accuracy.png

required_packages <- c("dplyr", "ggplot2", "purrr", "stringr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install required package(s) before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(tibble)
})

default_options <- list(
  model_pattern = "2026-07-27",
  model_dir = "HmscOutputs",
  output_dir = file.path("manuscript-outputs", "figures", "supplement"),
  output_slug = "sfig-full-vs-cv-accuracy",
  samples_slug = "samples_0250_thin_100_chains_4",
  nfolds = 10,
  width = 7.6,
  height = 7.2,
  dpi = 300
)

usage_text <- paste(
  "Usage:",
  "  Rscript scripts/4-data-visualisation/sfig-full-vs-cv-accuracy.R [options]",
  "",
  "Options:",
  "  --model-pattern VALUE   Model-folder date/pattern to include.",
  "  --model-dir PATH        Directory containing HMSC output folders.",
  "  --output-dir PATH       Directory for the PNG output.",
  "  --output-slug VALUE     Filename suffix without extension.",
  "  --samples-slug VALUE    HMSC sample-file slug.",
  "  --nfolds VALUE          Number of cross-validation folds.",
  "  --width VALUE           Figure width in inches.",
  "  --height VALUE          Figure height in inches.",
  "  --dpi VALUE             Output PNG resolution.",
  "  --help                  Show this message.",
  sep = "\n"
)

parse_cli_options <- function(args, defaults = default_options) {
  options <- defaults
  option_names <- names(defaults)

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]

    if (key %in% c("-h", "--help")) {
      cat(usage_text, "\n")
      quit(save = "no", status = 0)
    }

    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }

    if (grepl("=", key, fixed = TRUE)) {
      key_value <- strsplit(sub("^--", "", key), "=", fixed = TRUE)[[1]]
      option_name <- key_value[[1]]
      option_value <- paste(key_value[-1], collapse = "=")
    } else {
      option_name <- sub("^--", "", key)
      if (i == length(args)) {
        stop("Missing value for option: --", option_name, call. = FALSE)
      }
      i <- i + 1
      option_value <- args[[i]]
    }

    option_name <- gsub("-", "_", option_name, fixed = TRUE)
    if (!option_name %in% option_names) {
      stop("Unknown option: --", gsub("_", "-", option_name), call. = FALSE)
    }

    options[[option_name]] <- option_value
    i <- i + 1
  }

  options$nfolds <- as.integer(options$nfolds)
  options$width <- as.numeric(options$width)
  options$height <- as.numeric(options$height)
  options$dpi <- as.integer(options$dpi)

  if (is.na(options$nfolds) || options$nfolds < 2) {
    stop("--nfolds must be an integer greater than one.", call. = FALSE)
  }
  if (any(is.na(c(options$width, options$height))) || options$width <= 0 || options$height <= 0) {
    stop("--width and --height must be positive numbers.", call. = FALSE)
  }
  if (is.na(options$dpi) || options$dpi <= 0) {
    stop("--dpi must be a positive integer.", call. = FALSE)
  }

  options
}

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)

  repeat {
    has_project_file <- file.exists(file.path(current, "hmsc-hpc-cv-pipeline.Rproj"))
    has_project_dirs <- dir.exists(file.path(current, "HmscOutputs")) &&
      dir.exists(file.path(current, "scripts"))

    if (has_project_file || has_project_dirs) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop(
        "Could not find the project root. Run this script from inside the repository.",
        call. = FALSE
      )
    }
    current <- parent
  }
}

is_absolute_path <- function(path) {
  grepl("^(/|~|[A-Za-z]:[/\\\\])", path)
}

resolve_project_path <- function(path, project_root) {
  if (is_absolute_path(path)) {
    return(path)
  }

  file.path(project_root, path)
}

period_lookup <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

metric_lookup <- tribble(
  ~metric, ~metric_label,
  "AUC", "AUC",
  "TjurR2", "Tjur R2",
  "RMSE", "RMSE"
)

max_or_na <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  max(x, na.rm = TRUE)
}

quantile_or_na <- function(x, probability) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  as.numeric(quantile(x, probability, na.rm = TRUE))
}

extract_atlas <- function(model_folder) {
  atlas <- str_match(basename(model_folder), "Atlas([0-9]+)")[, 2]

  if (is.na(atlas) || !atlas %in% names(period_lookup)) {
    stop("Could not extract atlas number from model folder: ", model_folder, call. = FALSE)
  }

  atlas
}

find_model_folders <- function(model_dir, model_pattern) {
  list.dirs(model_dir, recursive = FALSE, full.names = TRUE) |>
    keep(~ str_detect(basename(.x), fixed(model_pattern))) |>
    keep(~ str_detect(basename(.x), "Atlas[0-9]+")) |>
    sort()
}

get_species_order <- function(model_folder, samples_slug) {
  fit_path <- file.path(
    model_folder,
    "Models",
    "Fitted",
    paste0("HPC_", samples_slug, ".Rdata")
  )

  if (!file.exists(fit_path)) {
    return(NULL)
  }

  fit_env <- new.env(parent = emptyenv())
  load(fit_path, envir = fit_env)

  if (!exists("fitted_model", envir = fit_env, inherits = FALSE)) {
    return(NULL)
  }

  species <- rownames(fit_env$fitted_model$posteriors$Tr)
  if (is.null(species)) {
    return(NULL)
  }

  species
}

check_full_model_metrics <- function(model_folder, metric_env, metric_names, samples_slug) {
  precv_path <- file.path(
    model_folder,
    "Models",
    "Fitted",
    paste0("MF_preCV_", samples_slug, ".Rdata")
  )

  if (!file.exists(precv_path)) {
    warning("No pre-CV full-model metric file found for ", basename(model_folder), call. = FALSE)
    return(invisible(FALSE))
  }

  precv_env <- new.env(parent = emptyenv())
  load(precv_path, envir = precv_env)

  if (!exists("MF", envir = precv_env, inherits = FALSE)) {
    warning("Pre-CV metric file does not contain MF for ", basename(model_folder), call. = FALSE)
    return(invisible(FALSE))
  }

  mismatch <- metric_names |>
    map_lgl(~ !isTRUE(all.equal(
      metric_env$MF[[.x]],
      precv_env$MF[[.x]],
      check.attributes = FALSE
    )))

  if (any(mismatch)) {
    warning(
      "MF metrics in the CV file differ from MF_preCV metrics for ",
      basename(model_folder),
      ": ",
      paste(metric_names[mismatch], collapse = ", "),
      call. = FALSE
    )
    return(invisible(FALSE))
  }

  invisible(TRUE)
}

read_model_accuracy <- function(model_folder, samples_slug, nfolds, evaluation_lookup) {
  atlas <- extract_atlas(model_folder)
  period <- period_lookup[[atlas]]

  accuracy_path <- file.path(
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

  required_objects <- c("MF", "MFCV")
  missing_objects <- required_objects[
    !vapply(required_objects, exists, logical(1), envir = metric_env, inherits = FALSE)
  ]
  if (length(missing_objects) > 0) {
    stop(
      "Metric file is missing expected object(s): ",
      paste(missing_objects, collapse = ", "),
      call. = FALSE
    )
  }

  check_full_model_metrics(
    model_folder = model_folder,
    metric_env = metric_env,
    metric_names = metric_lookup$metric,
    samples_slug = samples_slug
  )

  species <- get_species_order(model_folder, samples_slug)
  if (is.null(species)) {
    n_species <- length(metric_env$MF[[metric_lookup$metric[[1]]]])
    species <- paste0("species_", seq_len(n_species))
  }

  map_dfr(evaluation_lookup$evaluation, function(evaluation_name) {
    metric_set <- metric_env[[evaluation_name]]

    map_dfr(metric_lookup$metric, function(metric_name) {
      values <- metric_set[[metric_name]]

      if (is.null(values)) {
        stop(
          "Missing metric ",
          metric_name,
          " for ",
          evaluation_name,
          " in ",
          basename(model_folder),
          call. = FALSE
        )
      }

      if (length(values) != length(species)) {
        stop(
          "Metric/species length mismatch for ",
          basename(model_folder),
          " / ",
          evaluation_name,
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
        model_folder = basename(model_folder),
        atlas = atlas,
        period = period,
        species = species,
        evaluation = evaluation_name,
        metric = metric_name,
        value = as.numeric(values)
      )
    })
  })
}

prepare_accuracy_data <- function(model_folders, samples_slug, nfolds, evaluation_lookup) {
  map_dfr(
    model_folders,
    read_model_accuracy,
    samples_slug = samples_slug,
    nfolds = nfolds,
    evaluation_lookup = evaluation_lookup
  ) |>
    left_join(metric_lookup, by = "metric") |>
    left_join(evaluation_lookup, by = "evaluation") |>
    mutate(
      period = factor(.data$period, levels = unname(period_lookup)),
      metric_label = factor(.data$metric_label, levels = metric_lookup$metric_label),
      evaluation_label = factor(
        .data$evaluation_label,
        levels = evaluation_lookup$evaluation_label
      )
    )
}

summarise_accuracy <- function(accuracy_df) {
  accuracy_df |>
    group_by(.data$period, .data$metric_label, .data$evaluation_label) |>
    summarise(
      n_species = sum(!is.na(.data$value)),
      median = median(.data$value, na.rm = TRUE),
      q25 = quantile_or_na(.data$value, 0.25),
      q75 = quantile_or_na(.data$value, 0.75),
      .groups = "drop"
    ) |>
    left_join(
      accuracy_df |>
        group_by(.data$period, .data$metric_label) |>
        summarise(label_y = max_or_na(.data$value), .groups = "drop"),
      by = c("period", "metric_label")
    ) |>
    mutate(median_label = paste0("median = ", sprintf("%.2f", .data$median)))
}

plot_accuracy <- function(accuracy_df, summary_df, evaluation_lookup) {
  ggplot(
    accuracy_df,
    aes(x = .data$evaluation_label, y = .data$value)
  ) +
    geom_violin(
      aes(fill = .data$evaluation_label),
      width = 0.82,
      trim = FALSE,
      alpha = 0.46,
      colour = NA,
      na.rm = TRUE
    ) +
    geom_point(
      aes(fill = .data$evaluation_label),
      position = position_jitter(width = 0.055, height = 0, seed = 42),
      shape = 21,
      size = 0.85,
      stroke = 0.12,
      colour = "grey25",
      alpha = 0.48,
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
      aes(x = .data$evaluation_label, y = .data$median),
      inherit.aes = FALSE,
      shape = 23,
      size = 2.4,
      stroke = 0.35,
      fill = "white",
      colour = "black"
    ) +
    geom_text(
      data = summary_df,
      aes(x = .data$evaluation_label, y = .data$label_y, label = .data$median_label),
      inherit.aes = FALSE,
      vjust = -0.55,
      size = 2.45,
      colour = "grey15"
    ) +
    facet_grid(rows = vars(.data$metric_label), cols = vars(.data$period), scales = "free_y") +
    scale_x_discrete(labels = c("Full model" = "Full\nmodel", "10-fold CV" = "10-fold\nCV")) +
    scale_fill_manual(
      values = setNames(c("grey72", "#2a9d8f"), evaluation_lookup$evaluation_label),
      name = NULL
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.20))) +
    labs(
      x = NULL,
      y = "Metric value"
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
}

main <- function() {
  options <- parse_cli_options(commandArgs(trailingOnly = TRUE))
  project_root <- find_project_root()
  model_dir <- resolve_project_path(options$model_dir, project_root)
  output_dir <- resolve_project_path(options$output_dir, project_root)

  if (!dir.exists(model_dir)) {
    stop("Model directory does not exist: ", model_dir, call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  evaluation_lookup <- tribble(
    ~evaluation, ~evaluation_label,
    "MF", "Full model",
    "MFCV", paste0(options$nfolds, "-fold CV")
  )

  model_folders <- find_model_folders(model_dir, options$model_pattern)
  if (length(model_folders) == 0) {
    stop("No model folders found for pattern: ", options$model_pattern, call. = FALSE)
  }

  message("Using model folders: ", paste(basename(model_folders), collapse = ", "))

  accuracy_df <- prepare_accuracy_data(
    model_folders = model_folders,
    samples_slug = options$samples_slug,
    nfolds = options$nfolds,
    evaluation_lookup = evaluation_lookup
  )
  summary_df <- summarise_accuracy(accuracy_df)

  message("Median model-performance summaries:")
  print(summary_df, n = Inf)

  accuracy_plot <- plot_accuracy(
    accuracy_df = accuracy_df,
    summary_df = summary_df,
    evaluation_lookup = evaluation_lookup
  )

  safe_model_pattern <- gsub("[^A-Za-z0-9._-]+", "-", options$model_pattern)
  output_path <- file.path(
    output_dir,
    paste0(safe_model_pattern, "-", options$output_slug, ".png")
  )

  ggsave(
    output_path,
    accuracy_plot,
    width = options$width,
    height = options$height,
    units = "in",
    dpi = options$dpi,
    bg = "white"
  )

  message("Wrote: ", output_path)
  invisible(output_path)
}

main()
