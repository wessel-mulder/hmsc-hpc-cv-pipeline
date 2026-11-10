# Compare trait-environment Gamma posteriors between adjacent bird atlases.
# Run from the repository root:
#   Rscript scripts/4-post/compare-gamma-changes-adjacent-atlases.R
#
# Gamma draws are converted to original environmental units. Continuous STI
# terms are also converted from the model's standardized trait scale to raw STI
# units; factor terms remain 0/1 contrasts against their common reference level.

suppressPackageStartupMessages(library(ggplot2))
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
output_dir <- file.path("scripts", "4-post", "compare-gamma-changes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Full posterior draws are local derived data; tables and PNGs remain alongside the script.
draw_dir <- file.path("data-new", "4-post", "derived", "compare-gamma-changes")
dir.create(draw_dir, recursive = TRUE, showWarnings = FALSE)
model_stamp <- "2026-11-02_06-55-21"
model_file <- function(atlas) file.path(
  model_root,
  sprintf("%s_Atlas%d_MinOccs5_CoverageGoodAverage", model_stamp, atlas),
  "Models", "Fitted", "HPC_samples_0250_thin_100_chains_4.Rdata"
)

if (!file.exists(input_file)) stop("Missing preprocessed model input: ", input_file)
model_files <- vapply(1:3, model_file, character(1))
if (any(!file.exists(model_files))) {
  stop("Missing November 2 fitted model(s): ", paste(model_files[!file.exists(model_files)], collapse = ", "))
}
input_env <- new.env(parent = emptyenv())
load(input_file, envir = input_env)
if (!all(c("X", "Y", "Tr") %in% ls(input_env))) stop("Preprocessed input must contain X, Y and Tr.")
raw_X <- input_env$X

extract_model <- function(path, atlas) {
  model_env <- new.env(parent = emptyenv())
  load(path, envir = model_env)
  if (!"fitted_model" %in% ls(model_env)) stop("No fitted_model object in: ", path)
  m <- model_env$fitted_model$posteriors
  if (length(m$postList) != 4L || any(lengths(m$postList) != 250L)) {
    stop("Expected four chains with 250 retained draws in Atlas ", atlas, ".")
  }
  if (is.null(m$covNames) || is.null(m$trNames)) stop("Missing coefficient or trait-term names in Atlas ", atlas, ".")
  if (!all(rownames(m$X) %in% rownames(raw_X))) stop("Fitted sites absent from raw input in Atlas ", atlas, ".")
  if (!identical(colnames(m$XData), setdiff(m$covNames, "(Intercept)"))) {
    stop("Environmental predictor names/order differ in Atlas ", atlas, ".")
  }
  if (!identical(colnames(m$TrScaled), m$trNames)) stop("Trait-term names/order differ in Atlas ", atlas, ".")

  predictors <- setdiff(m$covNames, "(Intercept)")
  trait_terms <- setdiff(m$trNames, "(Intercept)")
  raw_atlas <- raw_X[rownames(m$X), predictors, drop = FALSE]
  model_atlas <- m$XData[, predictors, drop = FALSE]
  raw_mean <- vapply(raw_atlas, mean, numeric(1))
  raw_sd <- vapply(raw_atlas, sd, numeric(1))
  if (any(!is.finite(raw_sd) | raw_sd <= 0)) stop("Invalid raw environmental SD in Atlas ", atlas, ".")
  max_x_error <- max(abs(scale(raw_atlas, center = raw_mean, scale = raw_sd) - model_atlas))
  if (!is.finite(max_x_error) || max_x_error > 1e-8) {
    stop("Raw environmental input does not reproduce fitted XData in Atlas ", atlas,
         " (max abs error = ", signif(max_x_error, 4), ").")
  }

  # Confirm the numeric trait was scaled in the model and retain its raw-unit SD.
  sti_term <- "species_thermal_index"
  if (!sti_term %in% m$trNames || !sti_term %in% names(m$TrData)) stop("STI term missing in Atlas ", atlas, ".")
  sti_sd <- sd(m$TrData[[sti_term]], na.rm = TRUE)
  sti_term_i <- match(sti_term, m$trNames)
  if (!is.finite(sti_sd) || abs(sti_sd - m$TrScalePar[2, sti_term_i]) > 1e-8) {
    stop("Stored STI scale does not match raw STI SD in Atlas ", atlas, ".")
  }
  sti_z <- as.numeric(scale(m$TrData[[sti_term]]))
  sti_z_error <- max(abs(sti_z - m$TrScaled[, sti_term]), na.rm = TRUE)
  if (!is.finite(sti_z_error) || sti_z_error > 1e-8) stop("STI scaling mismatch in Atlas ", atlas, ".")

  gamma_chains <- lapply(m$postList, function(chain) {
    do.call(rbind, lapply(chain, function(draw) {
      g <- draw$Gamma
      if (!identical(dim(g), c(length(m$covNames), length(m$trNames)))) {
        stop("Unexpected Gamma dimensions in Atlas ", atlas, ".")
      }
      as.vector(g)
    }))
  })
  parameter_names <- as.vector(outer(m$covNames, m$trNames, paste, sep = "::"))
  gamma_chains <- lapply(gamma_chains, function(x) {
    colnames(x) <- parameter_names
    x
  })

  # Classical per-parameter Gelman-Rubin R-hat, vectorized over Gamma terms.
  chain_array <- simplify2array(gamma_chains)
  n_iter <- dim(chain_array)[1]
  chain_means <- apply(chain_array, c(2, 3), mean)
  chain_vars <- apply(chain_array, c(2, 3), var)
  within_var <- rowMeans(chain_vars)
  between_var <- n_iter * apply(chain_means, 1, var)
  var_hat <- ((n_iter - 1) / n_iter) * within_var + between_var / n_iter
  rhat <- sqrt(var_hat / within_var)
  names(rhat) <- parameter_names
  rhat_matrix <- matrix(rhat, nrow = length(m$covNames), ncol = length(m$trNames),
                        dimnames = list(m$covNames, m$trNames))

  # Convert Gamma: per raw environmental unit, and for STI per raw trait unit.
  trait_unit_sd <- rep(1, length(m$trNames))
  names(trait_unit_sd) <- m$trNames
  trait_unit_sd[[sti_term]] <- sti_sd
  gamma_array <- array(NA_real_,
                       dim = c(sum(lengths(m$postList)), length(m$covNames), length(m$trNames)),
                       dimnames = list(NULL, m$covNames, m$trNames))
  for (chain_i in seq_along(gamma_chains)) {
    chain_start <- (chain_i - 1L) * nrow(gamma_chains[[chain_i]]) + 1L
    chain_end <- chain_start + nrow(gamma_chains[[chain_i]]) - 1L
    gamma_array[chain_start:chain_end, , ] <- array(
      gamma_chains[[chain_i]],
      dim = c(nrow(gamma_chains[[chain_i]]), length(m$covNames), length(m$trNames))
    )
  }
  for (p in predictors) {
    for (tr in trait_terms) gamma_array[, p, tr] <- gamma_array[, p, tr] / (raw_sd[[p]] * trait_unit_sd[[tr]])
  }

  migration_ref <- levels(m$TrData$Migration_a3_DOF)[[1]]
  guild_ref <- levels(m$TrData$foraging_guild_consensus)[[1]]
  term_group <- ifelse(startsWith(trait_terms, "Migration_a3_DOF"), "Migration strategy",
                       ifelse(startsWith(trait_terms, "foraging_guild_consensus"), "Foraging guild",
                              "Species thermal index"))
  term_label <- ifelse(term_group == "Migration strategy",
                       sub("^Migration_a3_DOF", "", trait_terms),
                       ifelse(term_group == "Foraging guild",
                              sub("^foraging_guild_consensus", "", trait_terms), "STI"))
  trait_key <- data.frame(
    trait_term = trait_terms, trait_group = term_group, trait_label = term_label,
    trait_unit = ifelse(trait_terms == sti_term, "per STI unit", "factor contrast vs reference"),
    stringsAsFactors = FALSE
  )
  diagnostics <- data.frame(atlas = atlas, parameter = names(rhat), rhat = unname(rhat), stringsAsFactors = FALSE)
  scale_audit <- expand.grid(atlas = atlas, predictor = predictors, trait_term = trait_terms,
                             stringsAsFactors = FALSE)
  scale_audit$environmental_raw_sd <- raw_sd[scale_audit$predictor]
  scale_audit$trait_raw_unit_sd <- trait_unit_sd[scale_audit$trait_term]
  scale_audit$max_environmental_standardization_error <- max_x_error
  scale_audit$sti_standardization_error <- sti_z_error
  scale_audit$migration_reference <- migration_ref
  scale_audit$foraging_guild_reference <- guild_ref

  eligible <- rhat_matrix[predictors, trait_terms, drop = FALSE] < 1.1
  list(gamma = gamma_array, predictors = predictors, trait_terms = trait_terms,
       eligible = eligible, rhat = rhat_matrix, diagnostics = diagnostics,
       scale_audit = scale_audit, trait_key = trait_key,
       TrScaled = m$TrScaled, file = path)
}

models <- lapply(1:3, function(a) extract_model(model_files[[a]], a))
names(models) <- paste0("Atlas", 1:3)
predictors <- models[[1]]$predictors
trait_terms <- models[[1]]$trait_terms
for (a in 2:3) {
  if (!identical(models[[a]]$predictors, predictors)) stop("Environmental predictor lists/order differ.")
  if (!identical(models[[a]]$trait_terms, trait_terms)) stop("Trait terms/reference coding differ between atlases.")
  if (!identical(models[[a]]$TrScaled, models[[1]]$TrScaled)) stop("Trait design matrices differ between atlas models.")
}
trait_key <- models[[1]]$trait_key
write.csv(trait_key, file.path(output_dir, "gamma-trait-term-key.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(models, `[[`, "scale_audit")),
          file.path(output_dir, "gamma-scaling-audit.csv"), row.names = FALSE)
diagnostics <- do.call(rbind, lapply(models, `[[`, "diagnostics"))
write.csv(diagnostics, file.path(output_dir, "gamma-chain-diagnostics.csv"), row.names = FALSE)

diagnostics$predictor <- sub("::.*$", "", diagnostics$parameter)
diagnostics$trait_term <- sub("^.*::", "", diagnostics$parameter)
diagnostics$trait_group <- trait_key$trait_group[match(diagnostics$trait_term, trait_key$trait_term)]
diagnostics$coefficient_type <- ifelse(diagnostics$predictor == "(Intercept)" |
                                         diagnostics$trait_term == "(Intercept)", "intercept", "trait-predictor")
diag_groups <- split(diagnostics, interaction(diagnostics$atlas, diagnostics$coefficient_type,
                                              diagnostics$predictor, diagnostics$trait_group, drop = TRUE))
diagnostic_summary <- do.call(rbind, lapply(diag_groups, function(z) data.frame(
  atlas = z$atlas[[1]], coefficient_type = z$coefficient_type[[1]],
  predictor = z$predictor[[1]], trait_group = z$trait_group[[1]], n_parameters = nrow(z),
  median_rhat = median(z$rhat, na.rm = TRUE), q95_rhat = unname(quantile(z$rhat, 0.95, na.rm = TRUE)),
  max_rhat = max(z$rhat, na.rm = TRUE), n_rhat_ge_1_1 = sum(z$rhat >= 1.1, na.rm = TRUE)
)))
write.csv(diagnostic_summary, file.path(output_dir, "gamma-chain-diagnostic-summary.csv"), row.names = FALSE)

summarise_draws <- function(x) c(
  median = median(x), q025 = unname(quantile(x, 0.025)), q975 = unname(quantile(x, 0.975)),
  probability_positive = mean(x > 0), probability_negative = mean(x < 0)
)

atlas_summary <- list()
for (a in 1:3) {
  g <- models[[a]]$gamma
  for (p in predictors) for (tr in trait_terms) {
    if (!models[[a]]$eligible[p, tr]) next
    q <- summarise_draws(g[, p, tr])
    atlas_summary[[length(atlas_summary) + 1L]] <- data.frame(
      atlas = a, predictor = p, trait_key[match(tr, trait_key$trait_term), ], q,
      row.names = NULL, check.names = FALSE
    )
  }
}
atlas_summary <- do.call(rbind, atlas_summary)
write.csv(atlas_summary, file.path(output_dir, "gamma-atlas-posterior-summary-original-units.csv"), row.names = FALSE)

n_paired_draws <- 10000L
pair_order <- list(c(1L, 2L), c(2L, 3L))
pair_draws <- list()
change_summary <- list()
change_probability <- list()
exclusion_audit <- list()
for (pair in pair_order) {
  early <- pair[[1]]
  late <- pair[[2]]
  early_gamma <- models[[early]]$gamma
  late_gamma <- models[[late]]$gamma
  i_early <- sample.int(dim(early_gamma)[1], n_paired_draws, replace = TRUE)
  i_late <- sample.int(dim(late_gamma)[1], n_paired_draws, replace = TRUE)
  delta <- late_gamma[i_late, , , drop = FALSE] - early_gamma[i_early, , , drop = FALSE]
  pair_name <- paste0("Atlas", early, "_to_Atlas", late)
  eligible_pair <- models[[early]]$eligible & models[[late]]$eligible
  pair_draws[[pair_name]] <- list()
  for (p in predictors) {
    keep <- eligible_pair[p, ]
    dropped <- which(!keep)
    if (length(dropped)) exclusion_audit[[length(exclusion_audit) + 1L]] <- data.frame(
      atlas_pair = pair_name, predictor = p,
      trait_key[match(trait_terms[dropped], trait_key$trait_term), ],
      rhat_earlier = models[[early]]$rhat[p, trait_terms[dropped]],
      rhat_later = models[[late]]$rhat[p, trait_terms[dropped]],
      stringsAsFactors = FALSE
    )
    vals <- matrix(delta[, p, keep, drop = FALSE][, 1, , drop = FALSE],
                   nrow = n_paired_draws, ncol = sum(keep),
                   dimnames = list(NULL, trait_terms[keep]))
    pair_draws[[pair_name]][[p]] <- vals
    for (tr in trait_terms[keep]) {
      q <- summarise_draws(vals[, tr])
      term_info <- trait_key[match(tr, trait_key$trait_term), ]
      change_summary[[length(change_summary) + 1L]] <- data.frame(
        atlas_pair = pair_name, predictor = p, term_info,
        n_paired_posterior_draws = n_paired_draws, q,
        row.names = NULL, check.names = FALSE
      )
      change_probability[[length(change_probability) + 1L]] <- data.frame(
        atlas_pair = pair_name, predictor = p, term_info,
        probability_delta_positive = mean(vals[, tr] > 0),
        probability_delta_negative = mean(vals[, tr] < 0),
        row.names = NULL, check.names = FALSE
      )
    }
  }
}
change_summary <- do.call(rbind, change_summary)
change_probability <- do.call(rbind, change_probability)
write.csv(change_summary, file.path(output_dir, "gamma-adjacent-atlas-change-posterior-summary-original-units.csv"), row.names = FALSE)
write.csv(change_probability, file.path(output_dir, "gamma-adjacent-atlas-change-direction-probability.csv"), row.names = FALSE)
saveRDS(pair_draws, file.path(draw_dir, "gamma-adjacent-atlas-paired-posterior-draws.rds"), compress = "xz")
if (length(exclusion_audit)) {
  exclusion_audit <- do.call(rbind, exclusion_audit)
} else {
  exclusion_audit <- data.frame(atlas_pair = character(), predictor = character(),
                                trait_term = character(), rhat_earlier = numeric(), rhat_later = numeric())
}
write.csv(exclusion_audit, file.path(output_dir, "gamma-rhat-exclusion-audit.csv"), row.names = FALSE)

predictor_labels <- c(
  tmean_breeding = "Temperature (per °C)", prec_breeding = "Precipitation (per mm)",
  hh = "Landscape diversity (per unit)", perc_cropland = "Cropland (per percentage point)",
  perc_urban = "Urban cover (per percentage point)", perc_pasture = "Pasture (per percentage point)",
  perc_grass_shrub = "Grass/shrub cover (per percentage point)",
  perc_forest_deciduous = "Deciduous forest (per percentage point)",
  perc_forest_evergreen = "Evergreen forest (per percentage point)",
  perc_forest_mixed_other = "Mixed/other forest (per percentage point)"
)
change_probability$direction <- ifelse(change_probability$probability_delta_positive >= 0.5,
                                        "More positive", "More negative")
change_probability$direction_support <- pmax(change_probability$probability_delta_positive,
                                              change_probability$probability_delta_negative)
change_probability$support_band <- factor(
  ifelse(change_probability$direction_support > 0.95, ">95%",
         ifelse(change_probability$direction_support > 0.90, ">90–95%",
                ifelse(change_probability$direction_support > 0.85, ">85–90%", "≤85%"))),
  levels = c("≤85%", ">85–90%", ">90–95%", ">95%")
)
# Manual display filter. It changes the PNGs only; posterior tables and draws
# remain complete. The filter is applied to individual predictor-trait-atlas
# estimates, so no low-support dots remain in the figures.
if (is.null(plot_min_direction_support)) {
  supported_combinations <- unique(change_probability[c("predictor", "trait_group")])
  plot_suffix <- "all"
} else {
  supported_combinations <- unique(change_probability[
    change_probability$direction_support > plot_min_direction_support,
    c("predictor", "trait_group")
  ])
  plot_suffix <- paste0("support-over-", round(plot_min_direction_support * 100))
}
change_probability$combination_key <- paste(change_probability$predictor, change_probability$trait_group, sep = "||| ")
supported_keys <- paste(supported_combinations$predictor, supported_combinations$trait_group, sep = "||| ")
change_probability_plot <- change_probability[change_probability$combination_key %in% supported_keys, ]
if (!is.null(plot_min_direction_support)) {
  change_probability_plot <- change_probability[
    change_probability$direction_support > plot_min_direction_support, ]
}
if (!nrow(change_probability_plot)) stop("No Gamma predictor-by-trait-group combinations meet plot_min_direction_support.")
change_probability$atlas_pair_label <- factor(
  change_probability$atlas_pair,
  levels = c("Atlas1_to_Atlas2", "Atlas2_to_Atlas3"),
  labels = c("Atlas 1 → 2", "Atlas 2 → 3")
)
change_probability$predictor_label <- factor(
  predictor_labels[change_probability$predictor], levels = rev(unname(predictor_labels[predictors]))
)
change_probability$trait_group <- factor(change_probability$trait_group,
                                          levels = c("Species thermal index", "Migration strategy", "Foraging guild"))
change_probability$trait_label <- factor(change_probability$trait_label,
                                          levels = unique(trait_key$trait_label))

trait_plot_specs <- list(
  "Species thermal index" = list(width = 11, height = 7, angle = 0, label_size = 9),
  "Migration strategy" = list(width = 13, height = 8, angle = 0, label_size = 8),
  "Foraging guild" = list(width = 22, height = 10, angle = 90, label_size = 6)
)
for (group_name in names(trait_plot_specs)) {
  spec <- trait_plot_specs[[group_name]]
  plot_data <- change_probability_plot[change_probability_plot$trait_group == group_name, ]
  if (!nrow(plot_data)) next
  plot_data$trait_label <- droplevels(factor(plot_data$trait_label, levels = unique(trait_key$trait_label)))
  plot_data$predictor_label <- factor(predictor_labels[plot_data$predictor],
                                     levels = rev(unname(predictor_labels[predictors])))
  plot_data$atlas_pair_label <- factor(
    plot_data$atlas_pair,
    levels = c("Atlas1_to_Atlas2", "Atlas2_to_Atlas3"),
    labels = c("Atlas 1 → 2", "Atlas 2 → 3")
  )
  group_title <- switch(group_name,
                        "Species thermal index" = "Species thermal index",
                        "Migration strategy" = "Migration strategy",
                        "Foraging guild" = "Foraging guild")
  trait_bubble <- ggplot(plot_data,
                         aes(x = trait_label, y = predictor_label,
                             color = direction, size = support_band)) +
    geom_point() +
    scale_color_manual(values = c("More positive" = "#20853F", "More negative" = "#B2182B"),
                       limits = c("More positive", "More negative"), name = "Direction") +
    scale_size_manual(values = c("≤85%" = 4, ">85–90%" = 6.5,
                                 ">90–95%" = 9, ">95%" = 11.5),
                      limits = c("≤85%", ">85–90%", ">90–95%", ">95%"),
                      breaks = c("≤85%", ">85–90%", ">90–95%", ">95%"), drop = FALSE,
                      name = "Posterior support") +
    facet_grid(. ~ atlas_pair_label, scales = "free_x", space = "free_x") +
    scale_x_discrete(labels = function(x) ifelse(group_name == "Foraging guild" & nchar(x) > 22,
                                                  paste0(substr(x, 1, 19), "…"), x)) +
    labs(
      title = paste("Change in trait-modified environmental associations:", group_title),
      subtitle = paste0("Green/red shows change direction; circle size shows support. Only estimates above ",
                        if (is.null(plot_min_direction_support)) "no threshold" else paste0(round(plot_min_direction_support * 100), "%"), "."),
      x = "Trait term", y = NULL,
      caption = "Factor Gamma terms are relative to migration category 3 or foraging guild Tit-like birds. STI Gamma is per raw STI unit. Only R-hat < 1.1 terms in both atlases are shown."
    ) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = spec$angle, hjust = if (spec$angle == 90) 1 else 0.5,
                                     vjust = if (spec$angle == 90) 0.5 else 1, size = spec$label_size),
          axis.text.y = element_text(size = 8), strip.text = element_text(face = "bold", size = 10),
          plot.title = element_text(face = "bold", size = 14), legend.position = "right")
  filename <- switch(group_name,
                     "Species thermal index" = paste0("gamma-change-species-thermal-index-circles-", plot_suffix, ".png"),
                     "Migration strategy" = paste0("gamma-change-migration-strategy-circles-", plot_suffix, ".png"),
                     "Foraging guild" = paste0("gamma-change-foraging-guild-circles-", plot_suffix, ".png"))
  ggsave(file.path(output_dir, filename), trait_bubble,
         width = spec$width, height = spec$height, dpi = 300, bg = "white")
}

passing_counts <- vapply(pair_order, function(pair) sum(models[[pair[[1]]]]$eligible & models[[pair[[2]]]]$eligible), integer(1))
writeLines(c(
  "Trait-environment Gamma posterior contrasts across adjacent atlases",
  "Model set: November 2, 2026; Atlas1/2/3 MinOccs5 CoverageGoodAverage.",
  paste0("Full paired posterior draws (local only): ", file.path(draw_dir, "gamma-adjacent-atlas-paired-posterior-draws.rds")),
  "All retained Gamma draws from four chains x 250 draws were retained within atlas; posterior differences pair draws independently across fits.",
  "Environmental predictors were converted to original units using SDs from the exact fitted sites; scaling reproduces fitted XData to numerical tolerance.",
  "Trait design matrices and term coding are identical across atlases. Factor Gamma terms are 0/1 contrasts against Migration category 3 and Foraging guild Tit-like birds.",
  "Species thermal index is standardized in the fitted trait matrix; its Gamma terms were divided by the common raw STI SD and are reported per STI unit.",
  "Atlas pairs are Atlas 1 to 2 and Atlas 2 to 3. Only Gamma parameters with R-hat < 1.1 in both models are included.",
  paste0("PNG display filter: ", if (is.null(plot_min_direction_support)) "all estimates" else paste0("individual estimates with direction support above ", round(plot_min_direction_support * 100), "%"), ". Full posterior tables are not filtered."),
  paste0("Eligible predictor-by-trait contrasts: ", paste(passing_counts, collapse = ", "), " of ", length(predictors) * length(trait_terms), " for Atlas1-2 and Atlas2-3."),
  "The PNG circle plots show individual Gamma terms; they do not average across trait categories.",
  "The trait-level Gamma parameters describe modeled moderation of species-environment slopes; they are not causal effects."
), file.path(output_dir, "README.txt"))
message("Saved adjacent-atlas Gamma contrasts and circle plot under: ", output_dir)
message("Maximum non-intercept Gamma R-hat: ", signif(max(diagnostics$rhat[!startsWith(diagnostics$parameter, "(Intercept)::") &
                                                                  !endsWith(diagnostics$parameter, "::(Intercept)")], na.rm = TRUE), 4))
