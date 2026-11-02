source(file.path("support_scripts", "project_paths.R"))

hmsc_fitted_file <- "HPC_samples_0250_thin_100_chains_4.Rdata"
hmsc_model_fit_file <- "MF_samples_0250_thin_100_chains_4_nfolds_10.rdata"

atlas_numbers <- function(model_folders) {
  as.numeric(stringr::str_extract(basename(model_folders), "(?<=Atlas)\\d+"))
}

name_by_atlas <- function(x, model_folders) {
  names(x) <- atlas_numbers(model_folders)
  x
}

figure_model_folders <- function(pattern, base_dir = "HmscOutputs") {
  folders <- find_model_folders(base_dir = base_dir, pattern = pattern)
  if (length(folders) == 0) {
    stop("No model folders found for pattern: ", pattern, call. = FALSE)
  }
  folders
}

load_hmsc_posteriors <- function(model_folders, base_dir = "HmscOutputs",
                                 fitted_file = hmsc_fitted_file) {
  models <- lapply(model_folders, function(model_folder) {
    path <- file.path(base_dir, model_folder, "Models", "Fitted", fitted_file)
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    if (!exists("fitted_model", envir = env, inherits = FALSE)) {
      stop("Expected `fitted_model` in ", path, call. = FALSE)
    }
    env$fitted_model$posteriors
  })
  name_by_atlas(models, model_folders)
}

hmsc_study_design <- function(model) {
  model$studyDesign |>
    tibble::rownames_to_column(var = "survey") |>
    dplyr::left_join(
      model$ranLevels$site$s |> tibble::rownames_to_column(var = "site"),
      by = "site"
    )
}

load_hmsc_study_designs <- function(models) {
  lapply(models, hmsc_study_design)
}

load_hmsc_model_fit <- function(model_folders, base_dir = "HmscOutputs",
                                model_fit_file = hmsc_model_fit_file) {
  fits <- lapply(model_folders, function(model_folder) {
    path <- file.path(base_dir, model_folder, "Models", "Fitted", model_fit_file)
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    if (!exists("MF", envir = env, inherits = FALSE) ||
        !exists("MFCV", envir = env, inherits = FALSE)) {
      stop("Expected `MF` and `MFCV` in ", path, call. = FALSE)
    }
    list(MF = env$MF, MFCV = env$MFCV)
  })
  name_by_atlas(fits, model_folders)
}

load_or_compute_site_predictions <- function(models, model_folders,
                                             base_dir = "HmscOutputs",
                                             prediction_file = file.path("Results", "Preds", "predicted-values.rdata"),
                                             cache = TRUE) {
  predictions <- Map(function(model, folder) {
    target_path <- file.path(base_dir, folder, prediction_file)

    if (file.exists(target_path)) {
      message("Loading existing fitted-site predictions for: ", folder)
      return(readRDS(target_path))
    }

    message("Computing fitted-site predictions for: ", folder)
    pred_temp <- pcomputePredictedValues(model)
    pred_y <- rowMeans(pred_temp, dims = 2)

    if (cache) {
      dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(pred_y, target_path)
    }

    pred_y
  }, models, model_folders)

  name_by_atlas(predictions, model_folders)
}

load_gradient_predictions <- function(model_folders, base_dir = "HmscOutputs",
                                      fitted_file = hmsc_fitted_file) {
  predictions <- lapply(model_folders, function(model_folder) {
    path <- file.path(
      base_dir,
      model_folder,
      "Results",
      "Preds",
      sprintf("Preds_%s_%s", model_folder, fitted_file)
    )
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    if (!exists("Preds", envir = env, inherits = FALSE)) {
      stop("Expected `Preds` in ", path, call. = FALSE)
    }
    env$Preds
  })
  name_by_atlas(predictions, model_folders)
}

load_empirical_responses <- function(models) {
  lapply(models, function(model) model$Y)
}

community_weighted_means <- function(site_predictions, models) {
  Map(function(pred_y, model) {
    site_richness <- rowSums(pred_y)
    (pred_y %*% model$Tr) / matrix(rep(site_richness, model$nt), ncol = model$nt)
  }, site_predictions, models)
}

predicted_richness_frames <- function(site_predictions, designs) {
  purrr::map2(site_predictions, designs, function(pred_y, design) {
    pred_y |>
      as.data.frame() |>
      dplyr::mutate(richness = rowSums(pred_y)) |>
      tibble::rownames_to_column(var = "survey") |>
      dplyr::select(survey, richness) |>
      dplyr::left_join(design, by = "survey")
  })
}

load_vp_estimates <- function(model_folders, base_dir = "HmscOutputs",
                              scaled = FALSE) {
  vps <- lapply(model_folders, function(model_folder) {
    path <- file.path(
      base_dir,
      model_folder,
      "Results",
      sprintf("%sparameter_estimates_VP_.csv", model_folder)
    )
    df <- read.csv(path) |>
      tibble::column_to_rownames(var = "X")

    if (isTRUE(scaled)) {
      r2_values <- as.numeric(df["TjurR2", ])
      df <- df |>
        dplyr::filter(rownames(df) != "TjurR2") |>
        purrr::map2_dfc(r2_values, ~ .x * .y) |>
        dplyr::mutate(Factor = rownames(df)[rownames(df) != "TjurR2"]) |>
        tibble::column_to_rownames(var = "Factor")
    }

    df
  })
  name_by_atlas(vps, model_folders)
}

read_parameter_effects <- function(pattern, effect = "Beta",
                                   base_dir = "HmscOutputs") {
  effect_type <- ifelse(tolower(effect) == "gamma", "Gamma", "Beta")
  id_col_name <- ifelse(effect_type == "Gamma", "Traits", "Species")
  rename_to <- ifelse(effect_type == "Gamma", "traits", "species")
  model_folders <- figure_model_folders(pattern = pattern, base_dir = base_dir)

  outputs <- lapply(model_folders, function(model_folder) {
    path <- file.path(
      base_dir,
      model_folder,
      "Results",
      sprintf("%sparameter_estimates_%s_.xlsx", model_folder, effect_type)
    )

    long_effects <- readxl::read_excel(path, sheet = "Posterior mean") |>
      tidyr::pivot_longer(
        cols = -dplyr::all_of(id_col_name),
        names_to = "variable",
        values_to = "effect_size"
      )

    long_sigs <- readxl::read_excel(path, sheet = "Pr(x>0)") |>
      tidyr::pivot_longer(
        cols = -dplyr::all_of(id_col_name),
        names_to = "variable",
        values_to = "sig_val"
      )

    long_effects |>
      dplyr::left_join(long_sigs, by = c(id_col_name, "variable")) |>
      dplyr::mutate(effect_size = ifelse(sig_val > 0.05 & sig_val < 0.95, NA, effect_size)) |>
      dplyr::rename(!!rename_to := dplyr::all_of(id_col_name)) |>
      dplyr::select(-sig_val) |>
      dplyr::mutate(atlas = as.numeric(sub(".*Atlas([0-9]+).*", "\\1", model_folder)))
  })

  dplyr::bind_rows(outputs)
}
