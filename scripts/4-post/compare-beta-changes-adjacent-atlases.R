# Compare species-environment Beta posteriors between adjacent bird atlases.
# Run from the repository root:
#   Rscript scripts/4-post/compare-beta-changes-adjacent-atlases.R
#
# Betas are converted from each atlas's within-atlas standardized predictor
# scale back to the original predictor units before adjacent-period contrasts.

suppressPackageStartupMessages({
  library(ggplot2)
})

set.seed(20261103)

# Plot selection control: set to NULL to show every estimate. Otherwise omit
# each plotted estimate whose posterior support for either direction is at or
# below this threshold.
plot_min_direction_support <- 0.85

model_root <- file.path("data-new", "2-hmsc", "derived", "models")
input_file <- file.path(
  "data-new", "2-hmsc", "derived",
  "preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData"
)
output_dir <- file.path("scripts", "4-post", "compare-beta-changes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Full posterior draws are local derived data; tables and PNGs remain alongside the script.
draw_dir <- file.path("data-new", "4-post", "derived", "compare-beta-changes")
dir.create(draw_dir, recursive = TRUE, showWarnings = FALSE)

model_stamp <- "2026-11-02_06-55-21"
model_file <- function(atlas) {
  file.path(
    model_root,
    sprintf("%s_Atlas%d_MinOccs5_CoverageGoodAverage", model_stamp, atlas),
    "Models", "Fitted", "HPC_samples_0250_thin_100_chains_4.Rdata"
  )
}

if (!file.exists(input_file)) stop("Missing preprocessed model input: ", input_file)
model_files <- vapply(1:3, model_file, character(1))
if (any(!file.exists(model_files))) {
  stop("Missing November 2 fitted model(s): ", paste(model_files[!file.exists(model_files)], collapse = ", "))
}

input_env <- new.env(parent = emptyenv())
load(input_file, envir = input_env)
if (!all(c("X", "Y", "Design") %in% ls(input_env))) {
  stop("Preprocessed input must contain X, Y and Design.")
}
raw_X <- input_env$X

extract_model <- function(path, atlas) {
  model_env <- new.env(parent = emptyenv())
  load(path, envir = model_env)
  if (!"fitted_model" %in% ls(model_env)) stop("No fitted_model object in: ", path)
  m <- model_env$fitted_model$posteriors

  if (length(m$postList) != 4L || any(lengths(m$postList) != 250L)) {
    stop("Expected four chains with 250 retained draws in Atlas ", atlas, ".")
  }
  if (is.null(m$covNames) || is.null(m$spNames)) stop("Missing coefficient/species names in Atlas ", atlas, ".")
  if (!all(rownames(m$X) %in% rownames(raw_X))) {
    stop("Fitted sites do not all occur in the bundled preprocessing input for Atlas ", atlas, ".")
  }
  if (!identical(colnames(m$XData), setdiff(m$covNames, "(Intercept)"))) {
    stop("XData and Beta predictor names/order differ in Atlas ", atlas, ".")
  }

  predictors <- setdiff(m$covNames, "(Intercept)")
  raw_atlas <- raw_X[rownames(m$X), predictors, drop = FALSE]
  model_atlas <- m$XData[, predictors, drop = FALSE]
  if (!identical(rownames(raw_atlas), rownames(model_atlas))) stop("Site order mismatch in Atlas ", atlas, ".")

  # The saved Hmsc XScalePar describes the already standardized matrix. Verify
  # that the bundled raw inputs reproduce the exact model XData scaling.
  raw_mean <- vapply(raw_atlas, mean, numeric(1))
  raw_sd <- vapply(raw_atlas, sd, numeric(1))
  if (any(!is.finite(raw_sd) | raw_sd <= 0)) stop("Invalid raw predictor SD in Atlas ", atlas, ".")
  restandardized <- scale(raw_atlas, center = raw_mean, scale = raw_sd)
  max_abs_error <- max(abs(restandardized - model_atlas))
  if (!is.finite(max_abs_error) || max_abs_error > 1e-8) {
    stop("Raw predictor values do not reproduce fitted model scaling in Atlas ", atlas,
         " (max abs error = ", signif(max_abs_error, 4), ").")
  }

  # Store all retained draws, preserving coefficient and species order.
  beta_chains <- lapply(m$postList, function(chain) {
    do.call(rbind, lapply(chain, function(draw) {
      b <- draw$Beta
      if (!identical(dim(b), c(length(m$covNames), length(m$spNames)))) {
        stop("Unexpected Beta dimensions in Atlas ", atlas, ".")
      }
      as.vector(b)
    }))
  })
  parameter_names <- as.vector(outer(m$covNames, m$spNames, paste, sep = "::"))
  beta_chains <- lapply(beta_chains, function(x) {
    colnames(x) <- parameter_names
    x
  })
  # Classical per-parameter Gelman-Rubin R-hat, vectorized over all Betas.
  # Avoids the high overhead of constructing thousands of separate coda chains.
  chain_array <- simplify2array(beta_chains)
  n_iter <- dim(chain_array)[1]
  chain_means <- apply(chain_array, c(2, 3), mean)
  chain_vars <- apply(chain_array, c(2, 3), var)
  within_var <- rowMeans(chain_vars)
  between_var <- n_iter * apply(chain_means, 1, var)
  var_hat <- ((n_iter - 1) / n_iter) * within_var + between_var / n_iter
  rhat <- sqrt(var_hat / within_var)
  names(rhat) <- parameter_names

  # Separate scaling already occurred before Hmsc was called: the fitted
  # matrices should have mean 0 and SD 1, while raw means/SDs come from X.
  scale_audit <- data.frame(
    atlas = atlas,
    predictor = predictors,
    raw_mean = unname(raw_mean[predictors]),
    raw_sd = unname(raw_sd[predictors]),
    fitted_XData_mean = vapply(model_atlas, mean, numeric(1)),
    fitted_XData_sd = vapply(model_atlas, sd, numeric(1)),
    max_standardization_abs_error = max_abs_error,
    stringsAsFactors = FALSE
  )
  diagnostics <- data.frame(
    atlas = atlas,
    parameter = names(rhat),
    rhat = unname(rhat),
    stringsAsFactors = FALSE
  )
  rhat_matrix <- matrix(rhat, nrow = length(m$covNames), ncol = length(m$spNames),
                        dimnames = list(m$covNames, m$spNames))
  eligible <- rhat_matrix[predictors, , drop = FALSE] < 1.1

  # Convert each draw to an original-unit slope: beta_std = beta_raw * SD.
  beta_array <- array(NA_real_, dim = c(nrow(do.call(rbind, beta_chains)), length(m$covNames), length(m$spNames)),
                      dimnames = list(NULL, m$covNames, m$spNames))
  for (chain_i in seq_along(beta_chains)) {
    chain_start <- (chain_i - 1L) * nrow(beta_chains[[chain_i]]) + 1L
    chain_end <- chain_start + nrow(beta_chains[[chain_i]]) - 1L
    beta_array[chain_start:chain_end, , ] <- array(beta_chains[[chain_i]],
                                                   dim = c(nrow(beta_chains[[chain_i]]),
                                                           length(m$covNames), length(m$spNames)))
  }
  for (p in predictors) beta_array[, p, ] <- beta_array[, p, , drop = FALSE] / raw_sd[[p]]

  list(
    beta = beta_array,
    predictors = predictors,
    species = m$spNames,
    eligible = eligible,
    rhat = rhat_matrix,
    scale_audit = scale_audit,
    diagnostics = diagnostics,
    file = path
  )
}

models <- lapply(1:3, function(a) extract_model(model_files[[a]], a))
names(models) <- paste0("Atlas", 1:3)

species <- models[[1]]$species
predictors <- models[[1]]$predictors
for (a in 2:3) {
  if (!identical(models[[a]]$species, species)) stop("Species lists/order differ between atlas models.")
  if (!identical(models[[a]]$predictors, predictors)) stop("Predictor lists/order differ between atlas models.")
}

write.csv(do.call(rbind, lapply(models, `[[`, "scale_audit")),
          file.path(output_dir, "beta-scaling-audit.csv"), row.names = FALSE)
diagnostics <- do.call(rbind, lapply(models, `[[`, "diagnostics"))
write.csv(diagnostics, file.path(output_dir, "beta-chain-diagnostics.csv"), row.names = FALSE)
diagnostics$predictor <- sub("::.*$", "", diagnostics$parameter)
diagnostics$coefficient_type <- ifelse(diagnostics$predictor == "(Intercept)", "intercept", "predictor")
diagnostic_summary <- do.call(rbind, lapply(
  split(diagnostics, interaction(diagnostics$atlas, diagnostics$coefficient_type, diagnostics$predictor, drop = TRUE)),
  function(z) data.frame(
    atlas = z$atlas[[1]], coefficient_type = z$coefficient_type[[1]],
    predictor = z$predictor[[1]], n_parameters = nrow(z),
    median_rhat = median(z$rhat, na.rm = TRUE),
    q95_rhat = unname(quantile(z$rhat, 0.95, na.rm = TRUE)),
    max_rhat = max(z$rhat, na.rm = TRUE),
    n_rhat_ge_1_1 = sum(z$rhat >= 1.1, na.rm = TRUE)
  )
))
write.csv(diagnostic_summary,
          file.path(output_dir, "beta-chain-diagnostic-summary.csv"), row.names = FALSE)

summarise_draws <- function(x) {
  c(
    median = median(x),
    q025 = unname(quantile(x, 0.025)),
    q975 = unname(quantile(x, 0.975)),
    probability_positive = mean(x > 0),
    probability_negative = mean(x < 0)
  )
}

atlas_summary <- list()
for (a in 1:3) {
  beta <- models[[a]]$beta
  for (p in predictors) {
    keep <- models[[a]]$eligible[p, ]
    vals <- beta[, p, keep, drop = FALSE][, 1, ]
    q <- t(apply(vals, 2, summarise_draws))
    atlas_summary[[length(atlas_summary) + 1L]] <- data.frame(
      atlas = a, species = species[keep], predictor = p, q,
      row.names = NULL, check.names = FALSE
    )
  }
}
atlas_summary <- do.call(rbind, atlas_summary)
write.csv(atlas_summary, file.path(output_dir, "beta-atlas-posterior-summary-original-units.csv"), row.names = FALSE)

n_paired_draws <- 10000L
change_summaries <- list()
change_probability <- list()
paired_draws <- list()
exclusion_audit <- list()
species_index <- data.frame(species_index = seq_along(species), species = species)
write.csv(species_index, file.path(output_dir, "species-model-order.csv"), row.names = FALSE)

for (pair in list(c(1L, 2L), c(2L, 3L))) {
  early <- pair[[1]]
  late <- pair[[2]]
  beta_early <- models[[early]]$beta
  beta_late <- models[[late]]$beta

  # Independent fits: sample independently from each atlas posterior. A shared
  # sampled-draw index is used across coefficients to preserve draw dependence.
  idx_early <- sample.int(dim(beta_early)[1], n_paired_draws, replace = TRUE)
  idx_late <- sample.int(dim(beta_late)[1], n_paired_draws, replace = TRUE)
  delta <- beta_late[idx_late, , , drop = FALSE] - beta_early[idx_early, , , drop = FALSE]
  key <- paste0("Atlas", early, "_to_Atlas", late)
  eligible_pair <- models[[early]]$eligible & models[[late]]$eligible
  pair_draws <- list()

  for (p in predictors) {
    keep <- eligible_pair[p, ]
    excluded <- which(!keep)
    if (length(excluded)) {
      exclusion_audit[[length(exclusion_audit) + 1L]] <- data.frame(
        atlas_pair = key, species = species[excluded], predictor = p,
        rhat_earlier = models[[early]]$rhat[p, excluded],
        rhat_later = models[[late]]$rhat[p, excluded],
        failed_earlier_model = models[[early]]$rhat[p, excluded] >= 1.1,
        failed_later_model = models[[late]]$rhat[p, excluded] >= 1.1,
        stringsAsFactors = FALSE
      )
    }
    vals <- delta[, p, keep, drop = FALSE][, 1, , drop = FALSE]
    vals <- matrix(vals, nrow = n_paired_draws, ncol = sum(keep),
                   dimnames = list(NULL, species[keep]))
    pair_draws[[p]] <- vals
    q <- t(apply(vals, 2, summarise_draws))
    change_summaries[[length(change_summaries) + 1L]] <- data.frame(
      atlas_pair = key, species = species[keep], predictor = p,
      n_paired_posterior_draws = n_paired_draws, q,
      row.names = NULL, check.names = FALSE
    )
    change_probability[[length(change_probability) + 1L]] <- data.frame(
      atlas_pair = key, species_index = which(keep), species = species[keep],
      predictor = p, probability_delta_positive = colMeans(vals > 0),
      row.names = NULL
    )
  }
  paired_draws[[key]] <- pair_draws
}

change_summary <- do.call(rbind, change_summaries)
change_probability <- do.call(rbind, change_probability)
write.csv(change_summary, file.path(output_dir, "beta-adjacent-atlas-change-posterior-summary-original-units.csv"), row.names = FALSE)
write.csv(change_probability, file.path(output_dir, "beta-adjacent-atlas-change-probability-positive.csv"), row.names = FALSE)
saveRDS(paired_draws, file.path(draw_dir, "beta-adjacent-atlas-paired-posterior-draws.rds"), compress = "xz")
exclusion_audit <- do.call(rbind, exclusion_audit)
write.csv(exclusion_audit, file.path(output_dir, "beta-rhat-exclusion-audit.csv"), row.names = FALSE)

# Summarize the average slope change across the included species within each
# posterior draw. This is an overview statistic; opposing species changes can cancel.
aggregate_shift <- list()
for (key in names(paired_draws)) {
  for (p in predictors) {
    vals <- paired_draws[[key]][[p]]
    avg_delta <- rowMeans(vals)
    q <- summarise_draws(avg_delta)
    aggregate_shift[[length(aggregate_shift) + 1L]] <- data.frame(
      atlas_pair = key, predictor = p, n_species = ncol(vals),
      median_average_change = q[["median"]], q025 = q[["q025"]], q975 = q[["q975"]],
      probability_average_change_positive = q[["probability_positive"]],
      probability_average_change_negative = q[["probability_negative"]]
    )
  }
}
aggregate_shift <- do.call(rbind, aggregate_shift)
write.csv(aggregate_shift, file.path(output_dir, "beta-change-average-across-species.csv"), row.names = FALSE)

# Select named species with the strongest directional posterior evidence for
# an easier species-level view. This is a display selection, not a new test.
directional_examples <- do.call(rbind, lapply(
  split(change_summary, interaction(change_summary$atlas_pair, change_summary$predictor, drop = TRUE)),
  function(z) {
    z <- z[order(z$probability_positive), ]
    rbind(head(z, 3L), tail(z, 3L))
  }
))
write.csv(directional_examples, file.path(output_dir, "beta-change-directional-examples.csv"), row.names = FALSE)

predictor_labels <- c(
  tmean_breeding = "Breeding-season mean temperature (per °C)",
  prec_breeding = "Breeding-season precipitation (per mm)",
  hh = "Landscape diversity (per unit)",
  perc_cropland = "Cropland (per percentage point)",
  perc_urban = "Urban cover (per percentage point)",
  perc_pasture = "Pasture (per percentage point)",
  perc_grass_shrub = "Grass/shrub cover (per percentage point)",
  perc_forest_deciduous = "Deciduous forest (per percentage point)",
  perc_forest_evergreen = "Evergreen forest (per percentage point)",
  perc_forest_mixed_other = "Mixed/other forest (per percentage point)"
)
# Manual display filter. This affects figures only; saved summaries and draws
# remain complete. Individual species-predictor estimates are filtered, and
# predictors remain if at least one eligible estimate passes in either pair.
change_probability$direction_support <- pmax(
  change_probability$probability_delta_positive,
  1 - change_probability$probability_delta_positive
)
supported_predictors <- if (is.null(plot_min_direction_support)) predictors else {
  unique(change_probability$predictor[
    change_probability$direction_support > plot_min_direction_support
  ])
}
plot_predictors <- predictors[predictors %in% supported_predictors]
if (!length(plot_predictors)) stop("No Beta predictors meet plot_min_direction_support.")
plot_suffix <- if (is.null(plot_min_direction_support)) "all" else
  paste0("support-over-", round(plot_min_direction_support * 100))
change_probability$predictor_label <- predictor_labels[change_probability$predictor]
change_probability$predictor_label <- factor(
  change_probability$predictor_label,
  levels = unname(predictor_labels[predictors])
)

for (key in unique(change_probability$atlas_pair)) {
  plot_data <- change_probability[change_probability$atlas_pair == key &
                                    change_probability$predictor %in% plot_predictors, ]
  if (!is.null(plot_min_direction_support)) {
    plot_data <- plot_data[plot_data$direction_support > plot_min_direction_support, ]
  }
  all_tiles <- expand.grid(species_index = seq_along(species), predictor = plot_predictors,
                           KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  all_tiles <- merge(all_tiles, plot_data[, c("species_index", "predictor", "probability_delta_positive")],
                     by = c("species_index", "predictor"), all.x = TRUE, sort = FALSE)
  all_tiles$predictor_label <- factor(predictor_labels[all_tiles$predictor],
                                      levels = unname(predictor_labels[plot_predictors]))
  observed_tiles <- all_tiles[!is.na(all_tiles$probability_delta_positive), ]
  p <- ggplot() +
    geom_tile(data = all_tiles, aes(x = species_index, y = predictor_label), fill = "grey85") +
    geom_tile(data = observed_tiles,
              aes(x = species_index, y = predictor_label, fill = probability_delta_positive)) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0.5,
      limits = c(0, 1), breaks = c(0, 0.5, 1),
      name = "Posterior\nPr(change > 0)"
    ) +
    labs(
      title = paste("Species–environment Beta change:", gsub("_", " ", key)),
      subtitle = "Original predictor units; only species-predictor Betas with R-hat < 1.1 in both models are included",
      x = "Species in model order (see species-model-order.csv)", y = NULL,
      caption = paste0("Grey cells are excluded by R-hat or have direction support at or below ",
                       if (is.null(plot_min_direction_support)) "no threshold" else paste0(round(plot_min_direction_support * 100), "%"),
                       ". A positive change means the later slope is higher (more positive or less negative)." )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  ggsave(file.path(output_dir, paste0("beta-change-direction-probability-", plot_suffix, "-", tolower(key), ".png")),
         p, width = 13, height = 6.5, dpi = 300, bg = "white")
}

aggregate_shift_plot <- aggregate_shift[aggregate_shift$predictor %in% plot_predictors, ]
aggregate_shift_plot$predictor_label <- factor(
  predictor_labels[aggregate_shift_plot$predictor], levels = rev(unname(predictor_labels[plot_predictors]))
)
aggregate_shift_plot$atlas_pair_label <- factor(
  aggregate_shift_plot$atlas_pair,
  levels = c("Atlas1_to_Atlas2", "Atlas2_to_Atlas3"),
  labels = c("Atlas 1 to 2", "Atlas 2 to 3")
)
aggregate_shift_plot$direction <- ifelse(aggregate_shift_plot$probability_average_change_positive >= 0.5,
                                    "More positive", "More negative")
aggregate_shift_plot$direction_support <- pmax(
  aggregate_shift_plot$probability_average_change_positive,
  aggregate_shift_plot$probability_average_change_negative
)
if (!is.null(plot_min_direction_support)) {
  aggregate_shift_plot <- aggregate_shift_plot[
    aggregate_shift_plot$direction_support > plot_min_direction_support, ]
}
aggregate_shift_plot$support_band <- ifelse(
  aggregate_shift_plot$direction_support > 0.95, ">95%",
  ifelse(aggregate_shift_plot$direction_support > 0.90, ">90–95%",
         ifelse(aggregate_shift_plot$direction_support > 0.85, ">85–90%", "≤85%"))
)
write.csv(aggregate_shift, file.path(output_dir, "beta-change-average-across-species.csv"), row.names = FALSE)
aggregate_shift_plot$direction_support_label <- paste(aggregate_shift_plot$direction, aggregate_shift_plot$support_band, sep = " | ")
heatmap_levels <- c(
  "More positive | ≤85%", "More positive | >85–90%", "More positive | >90–95%", "More positive | >95%",
  "More negative | ≤85%", "More negative | >85–90%", "More negative | >90–95%", "More negative | >95%"
)
aggregate_shift_plot$direction_support_label <- factor(aggregate_shift_plot$direction_support_label, levels = heatmap_levels)
aggregate_shift_plot$atlas_pair_label <- factor(
  aggregate_shift_plot$atlas_pair,
  levels = c("Atlas1_to_Atlas2", "Atlas2_to_Atlas3"),
  labels = c("Atlas 1 → 2", "Atlas 2 → 3")
)
aggregate_shift_plot$predictor_label <- factor(
  predictor_labels[aggregate_shift_plot$predictor], levels = rev(unname(predictor_labels[plot_predictors]))
)
support_levels <- c("≤85%", ">85–90%", ">90–95%", ">95%")
aggregate_shift_plot$support_band <- factor(aggregate_shift_plot$support_band, levels = support_levels)
community_bubble <- ggplot(aggregate_shift_plot,
                           aes(x = atlas_pair_label, y = predictor_label,
                               color = direction, size = support_band)) +
  geom_point(shape = 16) +
  scale_color_manual(values = c("More positive" = "#20853F", "More negative" = "#B2182B"),
                     limits = c("More positive", "More negative"), name = "Direction") +
  scale_size_manual(values = c("≤85%" = 4.5, ">85–90%" = 7,
                               ">90–95%" = 9.5, ">95%" = 12),
                    limits = support_levels, breaks = support_levels, drop = FALSE,
                    name = "Posterior support") +
  labs(
    title = "Community-wide direction of species–environment change",
    subtitle = paste0("Circle color shows direction; size shows support. Only circles above ",
                      if (is.null(plot_min_direction_support)) "no threshold" else paste0(round(plot_min_direction_support * 100), "%"), " support are shown."),
    x = NULL, y = NULL,
    caption = "Size is posterior support for the average species-level Beta change.\nOnly contrasts with R-hat < 1.1 in both atlases are included; counts are in the CSV."
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 9), plot.title = element_text(face = "bold"),
        legend.position = "right")
ggsave(file.path(output_dir, paste0("beta-change-community-wide-circles-", plot_suffix, ".png")),
       community_bubble, width = 12, height = 9, dpi = 300, bg = "white")

directional_examples_plot <- directional_examples[directional_examples$predictor %in% plot_predictors, ]
directional_examples_plot$direction_support <- pmax(
  directional_examples_plot$probability_positive,
  1 - directional_examples_plot$probability_positive
)
if (!is.null(plot_min_direction_support)) {
  directional_examples_plot <- directional_examples_plot[
    directional_examples_plot$direction_support > plot_min_direction_support, ]
}
directional_examples_plot$predictor_label <- predictor_labels[directional_examples_plot$predictor]
directional_examples_plot$facet_label <- paste0(
  directional_examples_plot$predictor_label, " | ",
  ifelse(directional_examples_plot$atlas_pair == "Atlas1_to_Atlas2", "Atlas 1 to 2", "Atlas 2 to 3")
)
facet_levels <- unique(directional_examples_plot$facet_label)
species_levels <- unlist(lapply(facet_levels, function(f) {
  z <- directional_examples_plot[directional_examples_plot$facet_label == f, ]
  z <- z[order(z$median), ]
  paste(z$species, f, sep = "||| ")
}), use.names = FALSE)
directional_examples_plot$species_facet <- factor(
  paste(directional_examples_plot$species, directional_examples_plot$facet_label, sep = "||| "),
  levels = species_levels
)
examples_plot <- ggplot(directional_examples_plot, aes(y = species_facet)) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_segment(aes(x = q025, xend = q975, yend = species_facet), linewidth = 0.6) +
  geom_point(aes(x = median, color = probability_positive), size = 1.6) +
  scale_color_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0.5,
    limits = c(0, 1), name = "Posterior\nPr(change > 0)"
  ) +
  scale_y_discrete(labels = function(x) sub("\\|\\|\\|.*$", "", x)) +
  facet_wrap(~facet_label, scales = "free", ncol = 2) +
  labs(
    title = "Species with the clearest directional Beta changes",
    subtitle = paste0("Up to three highest and three lowest posterior probabilities per predictor and atlas pair; estimates must exceed ",
                      if (is.null(plot_min_direction_support)) "no threshold" else paste0(round(plot_min_direction_support * 100), "% support")),
    x = "Change in slope (original predictor units)", y = NULL,
    caption = "Intervals show 95% credible intervals. Selection ranks posterior direction evidence; selected intervals can still include zero."
  ) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold", size = 8), plot.title = element_text(face = "bold"))
ggsave(file.path(output_dir, paste0("beta-change-directional-species-examples-", plot_suffix, ".png")),
       examples_plot, width = 14, height = 17, dpi = 300, bg = "white")

writeLines(c(
  "Species-environment Beta posterior contrasts across adjacent atlases",
  "Model set: November 2, 2026; Atlas1/2/3 MinOccs5 CoverageGoodAverage.",
  paste0("Full paired posterior draws (local only): ", file.path(draw_dir, "beta-adjacent-atlas-paired-posterior-draws.rds")),
  "Each model contains four chains x 250 retained draws; all retained draws were concatenated within atlas.",
  "The preprocessed raw environmental matrix reproduces fitted XData after per-atlas standardization (max absolute error recorded in beta-scaling-audit.csv).",
  "Each Beta draw was divided by the SD of its predictor among that atlas's exact fitted sites, yielding slopes in original predictor units.",
  "Adjacent contrasts are Atlas2 minus Atlas1 and Atlas3 minus Atlas2.",
  paste0("For each contrast, ", n_paired_draws, " posterior draw pairs were sampled independently from the two fitted atlas posteriors (seed 20261103)."),
  paste0("PNG display filter: ", if (is.null(plot_min_direction_support)) "all estimates" else paste0("individual estimates with direction support above ", round(plot_min_direction_support * 100), "%"), ". Full posterior tables and draw files are not filtered."),
  "These are posterior contrasts between independently fitted models, not paired experimental tests.",
  "Evidence is reported continuously using posterior medians, 95% credible intervals and sign probabilities; no binary significance cutoff or multiplicity adjustment is applied.",
  paste0("Only Betas with R-hat < 1.1 are included. Eligible species-predictor contrasts: ",
         paste(vapply(list(c(1L, 2L), c(2L, 3L)), function(pair) {
           sum(models[[pair[[1]]]]$eligible & models[[pair[[2]]]]$eligible)
         }, integer(1)), collapse = ", "), " of ", length(predictors) * length(species), " for Atlas1-2 and Atlas2-3."),
  "beta-rhat-exclusion-audit.csv lists excluded species-predictor contrasts and the atlas model(s) that failed the threshold.",
  paste0("beta-change-community-wide-circles-", plot_suffix, ".png shows average slope changes across eligible species; only dots above the display support threshold are plotted."),
  paste0("beta-change-directional-species-examples-", plot_suffix, ".png names the three most positive and three most negative directional posterior contrasts per retained predictor and atlas pair."),
  paste0("Non-intercept Beta parameters with R-hat >= 1.1: ", paste(vapply(1:3, function(a) {
    sum(diagnostics$atlas == a & !grepl("^\\(Intercept\\)::", diagnostics$parameter) & diagnostics$rhat >= 1.1, na.rm = TRUE)
  }, integer(1)), collapse = ", "), " for Atlases 1, 2 and 3, respectively (out of 1,550 per atlas)."),
  "Treat contrasts involving poorly mixing parameters as exploratory; improve sampling/diagnostics before strong inferential claims.",
  "Posterior uncertainty is conditional on the fitted model structures, predictor set, priors and training data."
), file.path(output_dir, "README.txt"))

message("Saved adjacent-atlas Beta contrasts and PNG figures under: ", output_dir)
message("Maximum Beta R-hat across atlas models: ", signif(max(diagnostics$rhat, na.rm = TRUE), 4))
