# Descriptive distributions across species, not posterior distributions.
# Run from the project root: Rscript scripts/5-vis/vp-distributions/vp-distributions.R
# November 2 models; same species, ten environmental predictors, no spatial panel.
# Relative VP retains the spatial component in its original denominator.
# Scaled VP = exported relative VP * exported species-level Tjur R2,
# exactly as in vp-geom-sankey-bump.R. No significance tests are performed.

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2.")
suppressPackageStartupMessages(library(ggplot2))

model_stamp <- "2026-11-02_06-55-21"
model_root <- file.path("data-new", "2-hmsc", "derived", "models")
output_dir <- file.path("scripts", "5-vis", "vp-distributions")
atlas_labels <- c("1970s", "1990s", "2010s")
atlas_colours <- c("1970s" = "#0072B2", "1990s" = "#D55E00", "2010s" = "#009E73")
variable_labels <- c(
  tmean_breeding = "Temperature",
  prec_breeding = "Precipitation",
  hh = "Land-use heterogeneity",
  perc_urban = "Urban cover",
  perc_cropland = "Cropland cover",
  perc_pasture = "Pasture cover",
  perc_forest_deciduous = "Deciduous forest cover",
  perc_forest_evergreen = "Evergreen forest cover",
  perc_forest_mixed_other = "Mixed/other forest cover",
  perc_grass_shrub = "Grass/shrubland cover"
)

vp_tables <- lapply(1:3, function(atlas) {
  model <- sprintf("%s_Atlas%d_MinOccs5_CoverageGoodAverage", model_stamp, atlas)
  path <- file.path(model_root, model, "Results", paste0(model, "parameter_estimates_VP_.csv"))
  if (!file.exists(path)) stop("Missing November VP export: ", path)
  raw <- read.csv(path, row.names = 1, check.names = FALSE)
  predictors <- setdiff(rownames(raw), c("TjurR2", "Random: site"))
  if (!setequal(predictors, names(variable_labels)) ||
      !all(c("TjurR2", "Random: site") %in% rownames(raw))) {
    stop("Unexpected VP rows; update predictor labels before plotting: ", path)
  }
  if (!all(vapply(raw, is.numeric, logical(1))) || any(!is.finite(as.matrix(raw)))) {
    stop("Non-numeric or missing VP values: ", path)
  }
  if (any(as.matrix(raw) < 0) || any(as.matrix(raw) > 1)) {
    stop("VP or Tjur R2 outside [0, 1]: ", path)
  }
  if (max(abs(colSums(raw[setdiff(rownames(raw), "TjurR2"), ]) - 1)) > 1e-8) {
    stop("Relative VP components do not sum to one: ", path)
  }
  raw
})
species <- colnames(vp_tables[[1]])
if (length(species) != 155L || anyDuplicated(species) ||
    !all(vapply(vp_tables, function(x) setequal(colnames(x), species), logical(1)))) {
  stop("Expected the same 155 unique species in all three atlases.")
}

species_vp <- do.call(rbind, lapply(1:3, function(atlas) {
  raw <- vp_tables[[atlas]][, species, drop = FALSE]
  relative <- as.matrix(raw[names(variable_labels), , drop = FALSE])
  scaled <- sweep(relative, 2, as.numeric(raw["TjurR2", ]), "*")
  do.call(rbind, lapply(c("relative", "scaled"), function(scale_name) {
    values <- if (scale_name == "relative") relative else scaled
    data.frame(
      atlas = atlas_labels[[atlas]],
      variable = rep(names(variable_labels), times = length(species)),
      species = rep(species, each = length(variable_labels)),
      scale = scale_name,
      vp_percent = as.vector(values) * 100
    )
  }))
}))
species_vp$atlas <- factor(species_vp$atlas, levels = atlas_labels)
species_vp$variable <- factor(species_vp$variable, levels = names(variable_labels),
                              labels = unname(variable_labels))
stopifnot(nrow(species_vp) == 3L * 10L * 155L * 2L,
          all(table(species_vp$atlas, species_vp$variable, species_vp$scale) == 155L))

percent_labels <- function(x) paste0(format(x, trim = TRUE, scientific = FALSE), "%")
plot_theme <- theme_minimal(base_size = 12) + theme(
  legend.position = "top",
  panel.grid.minor = element_blank(),
  strip.text = element_text(face = "bold", size = 11),
  panel.spacing = grid::unit(1.2, "lines"),
  plot.title = element_text(face = "bold"),
  plot.caption = element_text(hjust = 0, size = 10)
)

# Gaussian KDE reflected at both boundaries (0 and 100 percentage points).
# A single pooled-species bandwidth and grid per predictor/scale ensure that
# smoothing differences cannot be caused by atlas-specific bandwidth choices.
# Evaluate out to four bandwidths beyond the largest value, including all species.
make_density <- function(dat) {
  do.call(rbind, lapply(split(dat, dat$variable, drop = TRUE), function(panel) {
    bandwidth <- stats::bw.nrd0(panel$vp_percent)
    if (!is.finite(bandwidth) || bandwidth <= 0) stop("Invalid KDE bandwidth.")
    grid_x <- seq(0, min(100, max(panel$vp_percent) + 4 * bandwidth), length.out = 2048)
    do.call(rbind, lapply(split(panel, panel$atlas, drop = TRUE), function(group) {
      x <- group$vp_percent
      kernel <- function(centres) {
        rowMeans(stats::dnorm(outer(grid_x, centres, "-"), sd = bandwidth))
      }
      density_y <- kernel(x) + kernel(-x) + kernel(200 - x)
      mass <- sum(diff(grid_x) * (head(density_y, -1) + tail(density_y, -1)) / 2)
      if (abs(mass - 1) > 0.01) stop("KDE numerical integral differs from one: ", mass)
      data.frame(variable = group$variable[[1]], atlas = group$atlas[[1]],
                 vp_percent = grid_x, density = density_y, bandwidth = bandwidth)
    }))
  }))
}

for (scale_name in c("relative", "scaled")) {
  dat <- species_vp[species_vp$scale == scale_name, ]
  density_data <- make_density(dat)
  axis_title <- if (scale_name == "relative") {
    "Relative VP (% of model-explained variation)"
  } else {
    "Tjur R²-scaled VP (%)"
  }
  scale_title <- if (scale_name == "relative") "Relative VP" else "Tjur R²-scaled VP"
  caption <- paste0(
    "November 2 models; 155 matched species per atlas; one exported estimate per species and predictor.\n",
    "Spatial component excluded from panels; original VP denominator retained. Descriptive estimates, without posterior uncertainty."
  )
  density_plot <- ggplot(density_data, aes(vp_percent, density, colour = atlas)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~variable, ncol = 2, scales = "free") +
    scale_colour_manual(values = atlas_colours, drop = FALSE) +
    scale_x_continuous(labels = percent_labels, expand = expansion(mult = c(0, 0.03))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = paste0(scale_title, ": distribution across species"),
         subtitle = "Boundary-corrected densities; common atlas bandwidth within each panel; panel axes vary",
         x = axis_title, y = "Density (per percentage point)", colour = "Atlas", caption = caption) +
    plot_theme

  ecdf_plot <- ggplot(dat, aes(vp_percent, colour = atlas)) +
    stat_ecdf(geom = "step", linewidth = 0.8, pad = TRUE) +
    facet_wrap(~variable, ncol = 2, scales = "free_x") +
    scale_colour_manual(values = atlas_colours, drop = FALSE) +
    scale_x_continuous(limits = c(0, NA), labels = percent_labels,
                       expand = expansion(mult = c(0, 0.03))) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = function(x) percent_labels(x * 100)) +
    labs(title = paste0(scale_title, ": cumulative distribution across species"),
         subtitle = "Height = percentage of species with VP at or below the horizontal-axis value; panel horizontal axes vary",
         x = axis_title, y = "Species at or below VP (%)", colour = "Atlas", caption = caption) +
    plot_theme

  for (plot_type in c("density", "ecdf")) {
    path <- file.path(output_dir, paste0("vp-distributions-", scale_name, "-", plot_type, ".png"))
    ggsave(path, if (plot_type == "density") density_plot else ecdf_plot,
           width = 12, height = 15, dpi = 300, bg = "white")
    message("Saved: ", path)
  }
}
message("Validated 155 matched species, all ten predictors, both VP scales; no significance tests performed.")
