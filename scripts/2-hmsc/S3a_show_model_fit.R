############################################################
# Calculate and visualise in-sample model fit before CV
#
# This script is intended to be run after:
#   S2_Model_Fitting_HPC_Version.R
#   S2b_HPC_output_merger.R
#
# It loads the merged Hmsc posterior, calculates fitted values and model-fit
# metrics, and creates one PNG per metric. No cross-validation output is required.
#
# Optional command-line arguments:
#   1: folder-name pattern (default: "2026-07-27")
#   2: samples             (default: 250)
#   3: thin                (default: 100)
#   4: chains              (default: 4)
#   5: parallel workers    (default: 4)
#
# Example:
#   Rscript S3a_show_model_fit.R 2026-07-27 250 100 4 4
############################################################

print("loading packages")
.libPaths(c("~/Rlibs", .libPaths()))

if (!requireNamespace("Hmsc", quietly = TRUE)) {
  stop("Package 'Hmsc' is required.")
}
if (!requireNamespace("cli", quietly = TRUE)) {
  stop("Package 'cli' is required.")
}

args <- commandArgs(trailingOnly = TRUE)
pattern2match <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "2026-11-11"
samples <- if (length(args) >= 2) as.integer(args[2]) else 250L
thin <- if (length(args) >= 3) as.integer(args[3]) else 100L
nChains <- if (length(args) >= 4) as.integer(args[4]) else 4L
nParallel <- if (length(args) >= 5) as.integer(args[5]) else 4L

numeric_args <- c(samples = samples, thin = thin, chains = nChains,
                  parallel_workers = nParallel)
if (anyNA(numeric_args) || any(numeric_args < 1)) {
  stop("Samples, thin, chains, and parallel workers must be positive integers.")
}

# Metric-specific display settings adapted from S5_show_model_fit.R.
metric_settings <- function(type, values) {
  finite_values <- values[is.finite(values)]

  if (type == "RMSE") {
    upper <- if (length(finite_values)) max(finite_values) else 1
    upper <- max(upper, .Machine$double.eps)
    return(list(xlim = c(0, upper), reference = 0))
  }

  if (type %in% c("AUC", "O.AUC")) {
    return(list(xlim = c(0, 1), reference = 0.5))
  }

  # R2-like metrics can be negative, although their usual upper limit is one.
  lower <- if (length(finite_values)) min(0, finite_values) else 0
  upper <- if (length(finite_values)) max(1, finite_values) else 1
  list(xlim = c(lower, upper), reference = 0)
}

plot_model_fit <- function(type, values, modelname, samples, thin, nChains) {
  values <- as.numeric(values)
  ok <- is.finite(values)
  settings <- metric_settings(type, values)
  metric_mean <- mean(values, na.rm = TRUE)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(1, 2), mar = c(5, 5, 5, 1))

  if (!any(ok)) {
    plot.new()
    title(main = sprintf("%s: no finite values", type))
    plot.new()
    return(invisible(NULL))
  }

  species_index <- seq_along(values)
  plot(values, species_index,
       xlim = settings$xlim,
       ylim = c(0.5, length(values) + 0.5),
       pch = 19, cex = 0.55,
       xlab = type, ylab = "Species index",
       main = sprintf(
         "%s\n%s\nthin = %i, samples = %i, chains = %i\nmean = %.4f",
         type, modelname, thin, samples, nChains, metric_mean
       ))
  abline(v = settings$reference, lty = 2, col = "grey45")
  abline(v = metric_mean, lwd = 2, col = "#0072B2")

  hist(values[ok],
       breaks = "FD",
       xlim = settings$xlim,
       col = "grey80", border = "white",
       xlab = type, main = sprintf("Distribution of %s", type))
  abline(v = settings$reference, lty = 2, col = "grey45")
  abline(v = metric_mean, lwd = 2, col = "#0072B2")
  legend("topright",
         legend = c(sprintf("Mean = %.4f", metric_mean), "Reference"),
         col = c("#0072B2", "grey45"),
         lwd = c(2, 1), lty = c(1, 2), bty = "n")

  invisible(NULL)
}

matching_folders <- list.dirs("HmscOutputs", recursive = FALSE, full.names = FALSE)
matching_folders <- matching_folders[
  grepl(pattern2match, basename(matching_folders), fixed = TRUE)
]

if (!length(matching_folders)) {
  stop(sprintf(
    "No folders directly below 'HmscOutputs' matched pattern '%s'.",
    pattern2match
  ))
}

for (models_description in matching_folders) {
  localDir <- file.path("HmscOutputs", models_description)
  ModelDir <- file.path(localDir, "Models", "Fitted")
  ResultDir <- file.path(localDir, "Results")

  input_file <- file.path(
    ModelDir,
    sprintf(
      "HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",
      samples, thin, nChains
    )
  )
  fit_file <- file.path(
    ModelDir,
    sprintf(
      "MF_preCV_samples_%.4d_thin_%.2d_chains_%.1d.Rdata",
      samples, thin, nChains
    )
  )


  cli::cli_h1("Model {models_description}")

  if (!file.exists(input_file)) {
    cli::cli_alert_warning("Merged posterior not found: {.file {input_file}}")
    next
  }

  if (!dir.exists(ResultDir)) {
    dir.create(ResultDir, recursive = TRUE)
  }

  if (file.exists(fit_file)) {
    cli::cli_alert_info("Loading cached pre-CV model fit: {.file {fit_file}}")
    fit_environment <- new.env(parent = emptyenv())
    load(fit_file, envir = fit_environment)
    if (!exists("MF", envir = fit_environment, inherits = FALSE)) {
      stop(sprintf("Cached file does not contain an object named MF: %s", fit_file))
    }
    MF <- fit_environment$MF
  } else {
    cli::cli_progress_step("Loading merged posterior")
    model_environment <- new.env(parent = emptyenv())
    load(input_file, envir = model_environment)

    if (!exists("fitted_model", envir = model_environment, inherits = FALSE) ||
        is.null(model_environment$fitted_model$posteriors)) {
      stop(sprintf(
        "Expected fitted_model$posteriors was not found in: %s",
        input_file
      ))
    }
    hM <- model_environment$fitted_model$posteriors
    rm(model_environment)
    cli::cli_progress_done()

    cli::cli_progress_step("Computing in-sample predicted values")
    preds <- Hmsc::computePredictedValues(hM, nParallel = nParallel)
    cli::cli_progress_done()

    cli::cli_progress_step("Evaluating in-sample model fit")
    MF <- Hmsc::evaluateModelFit(hM, predY = preds)
    cli::cli_progress_done()

    save(MF, file = fit_file)
    rm(hM, preds)
    cli::cli_alert_success("Saved pre-CV model fit: {.file {fit_file}}")
  }

  if (!is.list(MF) || !length(MF)) {
    stop(sprintf("Model-fit result is empty for model '%s'.", models_description))
  }

  for (type in names(MF)) {
    png_file <- file.path(
      ResultDir,
      sprintf("%s_model_fit_preCV_samples_%.4d_thin_%.2d_chains_%.1d_%s.png",
              models_description, samples, thin, nChains, make.names(type))
    )
    grDevices::png(png_file, width = 1800, height = 1050, res = 150)
    tryCatch(
      {
        cli::cli_progress_step("Plotting {type}")
        plot_model_fit(type, MF[[type]], models_description,
                       samples, thin, nChains)
        cli::cli_progress_done()
      },
      finally = grDevices::dev.off()
    )
    cli::cli_alert_success("Saved model-fit figure: {.file {png_file}}")
  }
}
