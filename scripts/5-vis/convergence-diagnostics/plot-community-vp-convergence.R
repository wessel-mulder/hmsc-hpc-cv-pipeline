# Run from the repository root:
# Rscript scripts/5-vis/convergence-diagnostics/plot-community-vp-convergence.R
suppressPackageStartupMessages(library(ggplot2))
stopifnot(requireNamespace("posterior", quietly = TRUE))
diagnostic_helpers <- new.env()
source("scripts/2-hmsc/model-diagnostics/helpers.R", local = diagnostic_helpers)
source("scripts/4-post/compare-scaled-vp-changes.R")

results_dir <- "scripts/4-post/compare-scaled-vp-changes"
output_dir <- "scripts/5-vis/convergence-diagnostics"
model_dirs <- c("Environment only" = results_dir, "Spatial" = file.path(results_dir, "spatial"))
atlas_labels <- c("1970s", "1990s", "2010s")
variable_labels <- c(
  tmean_breeding = "Temperature", prec_breeding = "Precipitation",
  hh = "Land-use heterogeneity", perc_urban = "Urban cover",
  perc_cropland = "Cropland cover", perc_pasture = "Pasture cover",
  perc_forest_deciduous = "Deciduous forest cover",
  perc_forest_evergreen = "Evergreen forest cover",
  perc_forest_mixed_other = "Mixed/other forest cover",
  perc_grass_shrub = "Grass/shrubland cover",
  "(Intercept)" = "Mean intercept (probit)"
)
trace_variables <- c("tmean_breeding", "prec_breeding", "perc_urban", "(Intercept)")
metrics <- c(
  rhat = sprintf("Rank-normalised R-hat  (< %s)", vp_diagnostic_thresholds$rhat),
  ess_bulk = "Bulk ESS  (≥ 400)",
  ess_tail = "Tail ESS  (≥ 400)",
  mcse_sd_ratio = "MCSE of mean / posterior SD  (≤ 5%)"
)
flag_columns <- c("flag_rhat", "flag_bulk", "flag_tail", "flag_mcse")
diagnostic_rows <- list()
trace_rows <- list()
intercept_rows <- list()

read_mean_intercept <- function(result) {
  file <- result$provenance$input
  stopifnot(identical(unname(tools::md5sum(file)), result$provenance$input_md5))
  fitted <- new.env(parent = emptyenv())
  load(file, envir = fitted)
  model <- fitted$fitted_model$posteriors
  intercept <- match("(Intercept)", model$covNames)
  stopifnot(!is.na(intercept), setequal(model$spNames, result$species),
    length(model$postList) == 4L, all(lengths(model$postList) == 250L))
  vapply(model$postList, function(chain) {
    vapply(chain, function(draw) mean(draw$Beta[intercept, ]), numeric(1))
  }, numeric(250L))
}

for (model in names(model_dirs)) {
  comparison <- readRDS(file.path(model_dirs[[model]], "comparison.rds"))
  stopifnot(isTRUE(comparison$provenance$settings$spatial) == (model == "Spatial"),
    identical(comparison$provenance$diagnostic_thresholds, vp_diagnostic_thresholds))
  diagnostic_helpers$diag_thresholds <- comparison$provenance$diagnostic_thresholds
  for (atlas in seq_along(atlas_labels)) {
    result <- readRDS(file.path(model_dirs[[model]], sprintf("atlas%d.rds", atlas)))
    stopifnot(result$atlas == atlas, identical(dim(result$community_percent), c(250L, 4L, 10L)),
      identical(result$provenance$input_md5, unname(comparison$provenance$model_md5[atlas])))
    summary <- result$community_summary
    stopifnot(setequal(summary$variable, setdiff(names(variable_labels), "(Intercept)")),
      !anyDuplicated(summary$variable))
    key <- paste(model, atlas)
    intercept_draws <- read_mean_intercept(result)
    intercept_summary <- diagnostic_helpers$diagnose_matrix(intercept_draws)
    intercept_rows[[key]] <- cbind(data.frame(model = model, atlas = atlas_labels[atlas]),
      intercept_summary[c("rhat", "ess_bulk", "ess_tail", "mcse_sd_ratio", "status")])
    summary_columns <- c(names(metrics), flag_columns)
    summary <- rbind(summary[c("variable", summary_columns)],
      cbind(variable = "(Intercept)", intercept_summary[summary_columns]))
    diagnostic_rows[[key]] <- do.call(rbind, lapply(seq_along(metrics), function(i) {
      values <- summary[[names(metrics)[i]]]
      flag <- summary[[flag_columns[i]]]
      label <- switch(names(metrics)[i], rhat = sprintf("%.3f", values),
        mcse_sd_ratio = sprintf("%.1f%%", values * 100), sprintf("%.0f", values))
      data.frame(model = model, atlas = atlas_labels[atlas], variable = summary$variable,
        metric = unname(metrics[i]), label = ifelse(is.finite(values), label, "NA"),
        status = ifelse(is.na(flag) | !is.finite(values), "Unavailable",
          ifelse(flag, "Flagged", "Pass")))
    }))
    vp_variables <- setdiff(trace_variables, "(Intercept)")
    draws <- result$community_percent[, , vp_variables, drop = FALSE]
    trace_rows[[key]] <- data.frame(
      expand.grid(iteration = seq_len(dim(draws)[1]), chain = seq_len(dim(draws)[2]),
        variable = trace_variables, stringsAsFactors = FALSE),
      value = c(as.vector(draws), as.vector(intercept_draws)), atlas = atlas_labels[atlas], model = model
    )
  }
}

diagnostics <- do.call(rbind, diagnostic_rows)
diagnostics$metric <- factor(diagnostics$metric, levels = unname(metrics))
diagnostics$variable <- factor(unname(variable_labels[diagnostics$variable]),
  levels = rev(unname(variable_labels)))
column_labels <- c(paste("Env. only", atlas_labels, sep = "\n"),
  paste("Spatial", atlas_labels, sep = "\n"))
diagnostics$column <- factor(paste(ifelse(diagnostics$model == "Spatial", "Spatial", "Env. only"),
  diagnostics$atlas, sep = "\n"), levels = column_labels)
stopifnot(nrow(diagnostics) == 264L, !anyNA(diagnostics$variable), !anyNA(diagnostics$column))

heatmap <- ggplot(diagnostics, aes(column, variable, fill = status)) +
  geom_tile(colour = "white", linewidth = .7) +
  geom_text(aes(label = label), size = 3.4, colour = "#202020") +
  geom_vline(xintercept = 3.5, colour = "white", linewidth = 1.4) +
  geom_hline(yintercept = 1.5, colour = "grey65", linewidth = .35) +
  facet_wrap(~metric, ncol = 2) +
  scale_fill_manual(values = c(Pass = "#E6EFE8", Flagged = "#F5B9A5", Unavailable = "#E5E5E5")) +
  scale_x_discrete(position = "top", expand = expansion(add = .02)) +
  scale_y_discrete(expand = expansion(add = .02)) +
  labs(title = "Convergence of community scaled VP and mean intercept",
    subtitle = "2 November spatial and 6 November environment-only fits",
    x = NULL, y = NULL, fill = NULL,
    caption = paste("Each cell uses four chains of 250 retained draws. The intercept is the mean baseline coefficient across 155 species per draw.",
      "The MCSE panel shows the numerical error of the posterior mean as a percentage of the posterior standard deviation.",
      sep = "\n")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
    axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 10),
    panel.spacing = grid::unit(1.4, "lines"),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 18),
    plot.caption = element_text(hjust = 0),
    plot.margin = margin(12, 15, 12, 12))

traces <- do.call(rbind, trace_rows)
traces$chain <- factor(traces$chain)
traces$model <- factor(traces$model, levels = names(model_dirs))
trace_names <- variable_labels
trace_names["(Intercept)"] <- "Mean intercept"
trace_labels <- unlist(lapply(unname(trace_names[trace_variables]),
  function(label) paste(label, atlas_labels, sep = "\n")), use.names = FALSE)
traces$panel <- factor(paste(unname(trace_names[traces$variable]), traces$atlas, sep = "\n"),
  levels = trace_labels)
stopifnot(nrow(traces) == 24000L, all(is.finite(traces$value)))
trace_plot <- ggplot(traces, aes(iteration, value, colour = chain)) +
  geom_line(linewidth = .28, alpha = .8) +
  facet_grid(panel ~ model, scales = "free_y") +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#009E73", "#9C5691")) +
  labs(title = "Saved chains for community scaled VP and mean intercept",
    subtitle = "Temperature, precipitation, urban cover and the baseline coefficient",
    x = "Retained draw within chain", y = "Community mean", colour = "Chain",
    caption = paste("Predictor rows show scaled VP (%). Intercept rows show coefficients on the probit scale.",
      "The vertical scale is shared across models within each row.", sep = "\n")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
    panel.spacing.y = grid::unit(.6, "lines"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 18),
    plot.caption = element_text(hjust = 0),
    plot.margin = margin(12, 15, 12, 12))

ggsave(file.path(output_dir, "community-vp-diagnostic-flags.png"), heatmap,
  width = 14, height = 10.8, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "community-vp-chain-traces.png"), trace_plot,
  width = 11, height = 16, dpi = 300, bg = "white")
print(with(subset(diagnostics, status == "Flagged"), table(model, atlas, metric)))
print(do.call(rbind, intercept_rows), row.names = FALSE)
