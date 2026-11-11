# Run from the repository root:
# Rscript scripts/5-vis/vp-distributions/vp-geom-sankey-bump-environment-only.R
# Compare species-mean Tjur R2-scaled VP for November spatial and environment-only fits.
# Original model inputs are paired within atlas; no common-cell refit or RSFA.
required_packages <- c("ggplot2", "ggsankey", "dplyr", "patchwork")
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})

model_dir <- "data-new/2-hmsc/derived/models"
output_dir <- "scripts/5-vis/vp-distributions"
atlas_labels <- c("1970s", "1990s", "2010s")
variable_labels <- c(
  tmean_breeding = "Temperature", prec_breeding = "Precipitation",
  hh = "Land-use heterogeneity", perc_urban = "Urban (% coverage)",
  perc_cropland = "Cropland (% coverage)", perc_pasture = "Pasture (% coverage)",
  perc_forest_deciduous = "Deciduous forest (% coverage)",
  perc_forest_evergreen = "Evergreen forest (% coverage)",
  perc_forest_mixed_other = "Mixed/other forest (% coverage)",
  perc_grass_shrub = "Grass/Shrubland (% coverage)"
)
variable_colours <- c(
  "Temperature" = "firebrick3", "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3", "Urban (% coverage)" = "snow3",
  "Cropland (% coverage)" = "goldenrod1", "Pasture (% coverage)" = "darkorange",
  "Deciduous forest (% coverage)" = "#A6D96A", "Evergreen forest (% coverage)" = "#006400",
  "Mixed/other forest (% coverage)" = "#87CEEB", "Grass/Shrubland (% coverage)" = "#B39DDB"
)
model_folder <- function(atlas, environment_only) {
  if (environment_only) sprintf("2026-11-06_05-11-30_Atlas%d_MinOccs5_CoverageGoodAverage_EnvironmentOnly", atlas)
  else sprintf("2026-11-02_06-55-21_Atlas%d_MinOccs5_CoverageGoodAverage", atlas)
}
load_fit <- function(folder) {
  e <- new.env(parent = emptyenv())
  load(file.path(model_dir, folder, "Models", "Fitted", "HPC_samples_0250_thin_100_chains_4.Rdata"), e)
  e$fitted_model$posteriors
}
read_vp <- function(folder, environment_only) {
  raw <- read.csv(file.path(model_dir, folder, "Results", paste0(folder, "parameter_estimates_VP_.csv")),
                   row.names = 1, check.names = FALSE)
  expected <- c(names(variable_labels), if (!environment_only) "Random: site", "TjurR2")
  if (!setequal(rownames(raw), expected) || ncol(raw) != 155L || anyDuplicated(colnames(raw))) {
    stop("Unexpected VP rows/species in ", folder)
  }
  if (!all(vapply(raw, is.numeric, logical(1))) || any(!is.finite(as.matrix(raw)))) {
    stop("Invalid VP values: ", folder)
  }
  parts <- as.matrix(raw[setdiff(rownames(raw), "TjurR2"), , drop = FALSE])
  if (any(parts < 0) || max(abs(colSums(parts) - 1)) > 1e-8) stop("Invalid relative VP partition.")
  raw
}

rows <- list(); totals <- list(); reference_species <- NULL
for (atlas in 1:3) {
  folders <- c(spatial = model_folder(atlas, FALSE), environment_only = model_folder(atlas, TRUE))
  old <- load_fit(folders[["spatial"]]); new <- load_fit(folders[["environment_only"]])
  for (field in c("Y", "X", "XData", "Tr", "studyDesign", "spNames", "covNames")) {
    if (!isTRUE(all.equal(do.call("$", list(old, field)), do.call("$", list(new, field))))) {
      stop("Inputs differ within Atlas ", atlas, ": ", field)
    }
  }
  stopifnot(old$nr == 1L, new$nr == 0L)
  for (type in names(folders)) {
    env_only <- type == "environment_only"
    raw <- read_vp(folders[[type]], env_only)
    if (is.null(reference_species)) reference_species <- colnames(raw)
    if (!setequal(colnames(raw), reference_species)) stop("Species sets differ.")
    raw <- raw[, reference_species, drop = FALSE]
    fit <- if (env_only) new else old
    stopifnot(identical(colnames(raw), fit$spNames))
    tjur <- as.numeric(raw["TjurR2", ])
    scaled <- sweep(as.matrix(raw[names(variable_labels), , drop = FALSE]), 2, tjur, "*")
    if (any(scaled < 0)) stop("Negative ribbon widths; inspect Tjur R2 rather than clipping.")
    means <- rowMeans(scaled)
    expected_total <- if (env_only) mean(tjur) else mean(tjur * (1 - as.numeric(raw["Random: site", ])))
    stopifnot(abs(sum(means) - expected_total) < 1e-10)
    key <- paste(atlas, type, sep = "::")
    rows[[key]] <- data.frame(model_type = type, atlas = atlas_labels[[atlas]],
                              variable = unname(variable_labels), mean_VP = unname(means))
    totals[[key]] <- data.frame(atlas = atlas, model_type = type, n_species = length(tjur),
                                environmental_scaled_vp = sum(means), mean_tjur = mean(tjur))
  }
  rm(old, new); gc(verbose = FALSE)
}
plot_data <- do.call(rbind, rows)
total_summary <- do.call(rbind, totals)
plot_data$atlas <- factor(plot_data$atlas, levels = atlas_labels)
plot_data$variable <- factor(plot_data$variable, levels = unname(variable_labels))
stopifnot(nrow(plot_data) == 60L)
y_max <- max(0.1, ceiling(max(total_summary$environmental_scaled_vp) / 0.025) * 0.025)
percent_labels <- function(x) paste0(sub("\\.0$", "", sprintf("%.1f", x * 100)), "%")
make_plot <- function(type, title, subtitle) {
  dat <- plot_data[plot_data$model_type == type, ]
  ggplot(dat, aes(x = atlas, node = variable, fill = variable, value = mean_VP)) +
    ggsankey::geom_sankey_bump(space = 0, type = "alluvial", alpha = 0.8, colour = NA, linewidth = 0) +
    scale_y_continuous(limits = c(0, y_max), breaks = seq(0, y_max, 0.025),
                       minor_breaks = seq(0, y_max, 0.0125), labels = percent_labels) +
    scale_fill_manual(values = variable_colours, drop = FALSE) +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Tjur R²-scaled VP (%)", fill = "Variable") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_line(colour = "grey90", linewidth = 0.25),
          plot.background = element_rect(fill = "grey99", colour = NA),
          axis.text.x = element_text(face = "bold", size = 12),
          plot.title = element_text(face = "bold"), legend.position = "right")
}
spatial_plot <- make_plot("spatial", "Environment + spatial effects", "Models: 2 November 2026")
environment_plot <- make_plot("environment_only", "Environment only", "Models: 6 November 2026")
caption <- paste0("155 matched species; mean across species of (VP × Tjur R²); all ten environmental predictors.\n",
  "Spatial-model environmental shares are not renormalised; ribbon heights reflect both VP allocation and model fit.")
environment_standalone <- environment_plot + labs(caption = paste0(
  "155 matched species; mean across species of (VP × Tjur R²); all ten environmental predictors.\n",
  "Ribbons sum to mean Tjur R²; heights reflect both VP allocation and model fit.")) +
  theme(plot.caption = element_text(hjust = 0, size = 10))
comparison <- (spatial_plot + environment_plot) + plot_layout(guides = "collect") +
  plot_annotation(title = "Variable importance: spatial versus environment-only models", caption = caption,
    theme = theme(plot.title = element_text(face = "bold", size = 17),
                  plot.caption = element_text(hjust = 0, size = 10),
                  plot.background = element_rect(fill = "grey99", colour = NA)))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "vp-geom-sankey-bump-environment-only.png"), environment_standalone,
       width = 11, height = 6.5, dpi = 300, bg = "grey99")
ggsave(file.path(output_dir, "vp-geom-sankey-bump-model-comparison.png"), comparison,
       width = 17, height = 7, dpi = 300, bg = "grey99")
message("Saved outline-free environment-only and model-comparison PNGs; common upper axis = ", y_max * 100, "%.")
print(total_summary, row.names = FALSE)
