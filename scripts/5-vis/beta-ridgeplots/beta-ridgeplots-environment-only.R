# Run from project root:
# Rscript scripts/5-vis/beta-ridgeplots/beta-ridgeplots-environment-only.R
# Betas retain the fitted atlas-specific predictor scaling: effect per one
# within-atlas predictor SD, on the probit linear-predictor scale.
# No conversion to original predictor units and no significance filtering.

required <- c("ggplot2", "ggridges", "posterior")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(ggplot2); library(ggridges)})
output_dir <- "scripts/5-vis/beta-ridgeplots"
model_root <- "data-new/2-hmsc/derived/models"
model_stamp <- "2026-11-06_05-11-30"
atlas_labels <- c("1970s", "1990s", "2010s")
variable_labels <- c(
  tmean_breeding = "Temperature", prec_breeding = "Precipitation",
  hh = "Land-use heterogeneity", perc_urban = "Urban cover", perc_cropland = "Cropland cover",
  perc_pasture = "Pasture cover", perc_forest_deciduous = "Deciduous forest cover",
  perc_forest_evergreen = "Evergreen forest cover", perc_forest_mixed_other = "Mixed/other forest cover",
  perc_grass_shrub = "Grass/shrubland cover")
variable_colours <- c(
  "Temperature" = "firebrick3", "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3", "Urban cover" = "snow3",
  "Cropland cover" = "goldenrod1", "Pasture cover" = "darkorange",
  "Deciduous forest cover" = "#A6D96A", "Evergreen forest cover" = "#006400",
  "Mixed/other forest cover" = "#87CEEB", "Grass/shrubland cover" = "#B39DDB")
input_env <- new.env(parent = emptyenv())
load("data-new/2-hmsc/derived/preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData",
     envir = input_env)
beta_draws <- list(); mean_rows <- list(); audit_rows <- list(); diagnostic_rows <- list()
reference_species <- NULL
for (atlas in 1:3) {
  folder <- sprintf("%s_Atlas%d_MinOccs5_CoverageGoodAverage_EnvironmentOnly", model_stamp, atlas)
  e <- new.env(parent = emptyenv())
  load(file.path(model_root, folder, "Models", "Fitted", "HPC_samples_0250_thin_100_chains_4.Rdata"), e)
  m <- e$fitted_model$posteriors
  predictors <- setdiff(m$covNames, "(Intercept)")
  stopifnot(m$nr == 0L, m$ns == 155L, setequal(predictors, names(variable_labels)),
            length(m$postList) == 4L, all(lengths(m$postList) == 250L))
  if (is.null(reference_species)) reference_species <- m$spNames
  stopifnot(setequal(m$spNames, reference_species), !anyDuplicated(m$spNames),
            all(rownames(m$XData) %in% rownames(input_env$X)),
            identical(rownames(m$XData), rownames(m$X)))
  pred_ix <- match(names(variable_labels), m$covNames)
  species_ix <- match(reference_species, m$spNames)
  raw <- input_env$X[rownames(m$XData), names(variable_labels), drop = FALSE]
  fitted_x <- as.matrix(m$X[, pred_ix, drop = FALSE])
  scaled_raw <- scale(raw)
  max_error <- max(abs(scaled_raw - fitted_x))
  stopifnot(is.finite(max_error), max_error < 1e-8,
            max(abs(colMeans(fitted_x))) < 1e-8,
            max(abs(apply(fitted_x, 2, sd) - 1)) < 1e-8)
  audit_rows[[atlas]] <- data.frame(atlas = atlas, variable = names(variable_labels),
    n_cells = nrow(raw), raw_mean = colMeans(raw), raw_sd = vapply(raw, sd, numeric(1)),
    fitted_mean = colMeans(fitted_x), fitted_sd = apply(fitted_x, 2, sd), max_scaling_error = max_error)
  draws <- array(NA_real_, c(250L, 4L, 10L, 155L),
    dimnames = list(NULL, NULL, names(variable_labels), reference_species))
  for (chain in 1:4) for (iteration in 1:250) {
    b <- m$postList[[chain]][[iteration]]$Beta
    stopifnot(identical(dim(b), c(11L, 155L)), all(is.finite(b)))
    draws[iteration, chain, , ] <- b[pred_ix, species_ix, drop = FALSE]
  }
  beta_draws[[atlas]] <- draws
  means <- apply(draws, c(3, 4), mean)
  # Verify selected species/predictor means by direct original-chain extraction.
  for (v in c(1L, 2L, 10L)) for (s in c(1L, 78L, 155L)) {
    direct <- unlist(lapply(m$postList, function(ch) vapply(ch,
      function(d) d$Beta[pred_ix[v], species_ix[s]], numeric(1))))
    stopifnot(abs(mean(direct) - means[v, s]) < 1e-12)
  }
  mean_rows[[atlas]] <- data.frame(atlas = atlas,
    variable = rep(names(variable_labels), times = 155L),
    species = rep(reference_species, each = 10L), beta_mean = as.vector(means))
  diagnostic_rows[[atlas]] <- do.call(rbind, lapply(1:10, function(v) {
    rh <- vapply(1:155, function(s) posterior::rhat(draws[, , v, s]), numeric(1))
    stopifnot(all(is.finite(rh)))
    data.frame(atlas = atlas, variable = names(variable_labels)[v],
      median_rhat = median(rh), max_rhat = max(rh), n_rhat_ge_1_1 = sum(rh >= 1.1))
  }))
  message("Atlas ", atlas, ": 155 species, ten predictors, 1,000 draws; scaling error ", signif(max_error, 3))
}
species_beta_means <- do.call(rbind, mean_rows)
predictor_scale_audit <- do.call(rbind, audit_rows)
beta_diagnostic_summary <- do.call(rbind, diagnostic_rows)
stopifnot(nrow(species_beta_means) == 4650L,
          all(table(species_beta_means$atlas, species_beta_means$variable) == 155L))

make_ridges <- function(mode = c("species-means", "posterior-mixture")) {
  mode <- match.arg(mode)
  get_values <- function(v, a) {
    if (mode == "species-means") {
      species_beta_means$beta_mean[species_beta_means$atlas == a &
                                    species_beta_means$variable == names(variable_labels)[v]]
    } else as.vector(beta_draws[[a]][, , v, ])
  }
  bandwidths <- vapply(1:10, function(v) {
    bw.nrd0(unlist(lapply(1:3, function(a) get_values(v, a)), use.names = FALSE))
  }, numeric(1))
  stopifnot(all(is.finite(bandwidths)), all(bandwidths > 0))
  max_beta <- max(vapply(1:10, function(v) max(abs(unlist(lapply(1:3,
    function(a) get_values(v, a)), use.names = FALSE))), numeric(1)))
  radius <- ceiling((max_beta + 4 * max(bandwidths)) / 0.25) * 0.25
  density_rows <- lapply(1:10, function(v) {
    d <- lapply(1:3, function(a) {
      x <- get_values(v, a)
      stopifnot(length(x) == if (mode == "species-means") 155L else 155000L)
      den <- density(x, bw = bandwidths[v], from = -radius, to = radius, n = 2048)
      area <- sum(diff(den$x) * (head(den$y, -1) + tail(den$y, -1)) / 2)
      stopifnot(abs(area - 1) < .01)
      data.frame(variable = unname(variable_labels[v]), atlas = atlas_labels[a],
                 x = den$x, density = den$y, bandwidth = bandwidths[v])
    })
    d <- do.call(rbind, d)
    # One common height multiplier for all three atlas ridges of this variable.
    d$height <- d$density / max(d$density)
    d
  })
  ridge_data <- do.call(rbind, density_rows)
  ridge_data$variable <- factor(ridge_data$variable, levels = unname(variable_labels))
  ridge_data$atlas <- factor(ridge_data$atlas, levels = rev(atlas_labels))
  is_mean <- mode == "species-means"
  p <- ggplot(ridge_data, aes(x = x, y = atlas, height = height,
                              fill = variable, alpha = atlas, group = atlas)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = .4, linetype = "dashed") +
    ggridges::geom_ridgeline(stat = "identity", scale = .8, colour = NA, linewidth = 0,
                              min_height = 0) +
    facet_grid(variable ~ ., switch = "y") +
    scale_fill_manual(values = variable_colours) +
    scale_alpha_manual(values = c("1970s" = .45, "1990s" = .7, "2010s" = 1)) +
    scale_x_continuous(breaks = seq(-1, 1, .25), expand = expansion(mult = c(0, 0))) +
    coord_cartesian(xlim = c(-1, 1), expand = FALSE) +
    scale_y_discrete(expand = expansion(add = c(.15, 1))) +
    labs(title = "Environment-only models: species–environment effect sizes",
      subtitle = if (is_mean) "Each ridge: 155 species' posterior mean Betas; three atlas ridges per variable"
                 else "Each ridge: all 1,000 Beta draws for each of 155 species, equally weighted",
      x = "Beta: change in probit linear predictor per 1 within-atlas SD of the environmental predictor",
      y = NULL,
      caption = paste0(
        "November 6 models; all ten predictors; no intercept or significance filtering. Dashed line: zero effect.\n",
        if (is_mean) "Ridges describe differences among species; they do not show uncertainty in each species' coefficient.\n"
        else "Ridges mix differences among species with posterior uncertainty; they are not uncertainty intervals for a community mean.\n",
        "One bandwidth per variable shared across atlases. Display cropped to [−1, +1]; densities use all values. No outlines.\n",
        "Atlas-specific SDs differ: temporal shifts can reflect changing predictor variability as well as changing raw-unit slopes.")) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none", panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      strip.placement = "outside", strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 10),
      axis.text.y = element_text(size = 9), panel.spacing.y = grid::unit(.2, "lines"),
      plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0, size = 9),
      plot.margin = margin(12, 12, 12, 12))
  list(plot = p, data = ridge_data, radius = radius, bandwidths = bandwidths)
}
ridge_results <- lapply(c("species-means", "posterior-mixture"), make_ridges)
names(ridge_results) <- c("species-means", "posterior-mixture")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (mode in names(ridge_results)) {
  path <- file.path(output_dir, paste0("beta-ridgeplots-environment-only-", mode, ".png"))
  ggsave(path, ridge_results[[mode]]$plot, width = 13, height = 15, dpi = 300, bg = "white")
  message("Saved: ", path)
}
print(beta_diagnostic_summary, row.names = FALSE)
print(predictor_scale_audit[, c("atlas", "variable", "raw_sd")], row.names = FALSE)
