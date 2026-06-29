# Predicted environmental-gradient richness direction
#
# This exploratory script reads the saved HMSC environmental-gradient prediction
# objects from the "CoverageGoodAverage" atlas models. For each environmental
# covariate, it sums species occurrence probabilities at each gradient point to
# get predicted richness, then compares the high end of the gradient with the
# low end for each posterior draw.
#
# The main question is whether the direction of the richness response changes
# among atlas periods. The script reports both the median direction and a simple
# posterior-support direction. A direction is treated as supported when at least
# 95% of posterior draws agree that the high-minus-low richness delta is above
# or below zero.
#
# No plots are produced here; the outputs are CSV summaries.

remove(list = ls())

model_dirs <- c(
  "HmscOutputs/2026-03-13_06-58-56_Atlas1_MinOccs5_CoverageGoodAverage",
  "HmscOutputs/2026-03-13_06-58-56_Atlas2_MinOccs5_CoverageGoodAverage",
  "HmscOutputs/2026-03-13_06-58-56_Atlas3_MinOccs5_CoverageGoodAverage"
)

out_dir <- "notebooks/exploratory/outputs/predicted-gradient-richness-direction"
run_date <- format(Sys.Date())

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(
  out_dir,
  sprintf("%s-richness-gradient-direction-summary.csv", run_date)
)
comparison_path <- file.path(
  out_dir,
  sprintf("%s-richness-gradient-direction-comparison.csv", run_date)
)

extract_atlas <- function(model_dir) {
  model_name <- basename(model_dir)
  atlas <- regmatches(model_name, regexpr("Atlas[0-9]+", model_name))

  if (length(atlas) == 0 || is.na(atlas)) {
    stop("Could not identify atlas number from model folder: ", model_dir, call. = FALSE)
  }

  atlas
}

find_gradient_prediction_file <- function(model_dir) {
  pred_dir <- file.path(model_dir, "Results", "Preds")
  pred_file <- list.files(
    pred_dir,
    pattern = "^Preds_.*[.]Rdata$",
    full.names = TRUE
  )

  if (length(pred_file) != 1L) {
    stop(
      "Expected one environmental-gradient Preds file in ",
      pred_dir,
      " but found ",
      length(pred_file),
      call. = FALSE
    )
  }

  pred_file
}

summarise_prediction <- function(pred_list, gradient, focal_var) {
  # HMSC stores one matrix per posterior draw. Rows are gradient points and
  # columns are species, so rowSums gives expected richness per gradient point.
  richness_by_draw <- do.call(cbind, lapply(pred_list, rowSums))
  focal_values <- gradient$XDataNew[[focal_var]]

  if (is.null(focal_values)) {
    stop("The focal variable is missing from the saved gradient: ", focal_var, call. = FALSE)
  }

  low_idx <- 1L
  high_idx <- length(focal_values)

  delta <- richness_by_draw[high_idx, ] - richness_by_draw[low_idx, ]
  slope <- apply(
    richness_by_draw,
    2,
    function(y) unname(coef(lm(y ~ focal_values))[2])
  )

  data.frame(
    gradient_min = focal_values[low_idx],
    gradient_max = focal_values[high_idx],
    n_gradient_points = length(focal_values),
    n_posterior_draws = ncol(richness_by_draw),
    n_species = ncol(pred_list[[1]]),
    richness_low_median = median(richness_by_draw[low_idx, ]),
    richness_high_median = median(richness_by_draw[high_idx, ]),
    delta_median = median(delta),
    delta_q025 = unname(quantile(delta, 0.025)),
    delta_q975 = unname(quantile(delta, 0.975)),
    pr_delta_positive = mean(delta > 0),
    slope_median = median(slope),
    slope_pr_positive = mean(slope > 0)
  )
}

summarise_model_dir <- function(model_dir) {
  model_name <- basename(model_dir)
  atlas <- extract_atlas(model_dir)
  pred_file <- find_gradient_prediction_file(model_dir)

  pred_env <- new.env(parent = emptyenv())
  load(pred_file, envir = pred_env)

  if (!exists("Preds", envir = pred_env, inherits = FALSE)) {
    stop("The prediction file does not contain an object named Preds: ", pred_file, call. = FALSE)
  }

  Preds <- pred_env$Preds
  rows <- list()

  for (covariate in names(Preds)) {
    # "total" follows the default constructGradient behaviour used in
    # S6_Making_Predictions.R. "marginal" uses the Gradient2 object, where
    # non-focal covariates were held fixed by non.focalVariables = 1.
    total <- summarise_prediction(Preds[[covariate]]$predY, Preds[[covariate]]$Gradient, covariate)
    marginal <- summarise_prediction(Preds[[covariate]]$predY2, Preds[[covariate]]$Gradient2, covariate)

    rows[[length(rows) + 1L]] <- cbind(
      model_name = model_name,
      atlas = atlas,
      covariate = covariate,
      effect = "total",
      total
    )
    rows[[length(rows) + 1L]] <- cbind(
      model_name = model_name,
      atlas = atlas,
      covariate = covariate,
      effect = "marginal",
      marginal
    )
  }

  do.call(rbind, rows)
}

add_direction_labels <- function(summary_df) {
  summary_df$direction <- ifelse(
    summary_df$delta_median > 0,
    "positive",
    ifelse(summary_df$delta_median < 0, "negative", "flat")
  )
  summary_df$supported_direction <- ifelse(
    summary_df$pr_delta_positive >= 0.95,
    "positive",
    ifelse(summary_df$pr_delta_positive <= 0.05, "negative", "uncertain")
  )

  summary_df
}

make_comparison_table <- function(summary_df) {
  keys <- unique(summary_df[, c("effect", "covariate")])
  comparison_rows <- list()

  for (i in seq_len(nrow(keys))) {
    subset_df <- summary_df[
      summary_df$effect == keys$effect[i] &
        summary_df$covariate == keys$covariate[i],
    ]

    median_dirs <- setNames(subset_df$direction, subset_df$atlas)
    supported_dirs <- setNames(subset_df$supported_direction, subset_df$atlas)
    deltas <- setNames(round(subset_df$delta_median, 3), subset_df$atlas)

    supported_non_uncertain <- supported_dirs[supported_dirs != "uncertain"]

    comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
      effect = keys$effect[i],
      covariate = keys$covariate[i],
      Atlas1_median_direction = unname(median_dirs["Atlas1"]),
      Atlas2_median_direction = unname(median_dirs["Atlas2"]),
      Atlas3_median_direction = unname(median_dirs["Atlas3"]),
      Atlas1_supported_direction = unname(supported_dirs["Atlas1"]),
      Atlas2_supported_direction = unname(supported_dirs["Atlas2"]),
      Atlas3_supported_direction = unname(supported_dirs["Atlas3"]),
      Atlas1_delta = unname(deltas["Atlas1"]),
      Atlas2_delta = unname(deltas["Atlas2"]),
      Atlas3_delta = unname(deltas["Atlas3"]),
      median_direction_change = length(unique(median_dirs)) > 1L,
      supported_direction_change = length(unique(supported_non_uncertain)) > 1L
    )
  }

  do.call(rbind, comparison_rows)
}

summary_df <- do.call(rbind, lapply(model_dirs, summarise_model_dir))
summary_df <- add_direction_labels(summary_df)
summary_df <- summary_df[order(summary_df$effect, summary_df$covariate, summary_df$atlas), ]

comparison_df <- make_comparison_table(summary_df)

write.csv(summary_df, summary_path, row.names = FALSE)
write.csv(comparison_df, comparison_path, row.names = FALSE)

message("Saved summary: ", summary_path)
message("Saved comparison: ", comparison_path)
print(comparison_df, row.names = FALSE)
