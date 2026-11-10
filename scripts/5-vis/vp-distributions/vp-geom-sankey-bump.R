# Compare variable importance across atlases using the VP panel from
# old-scripts/fig1-vp-effect-sizes.R. Run from the project root with:
# Rscript scripts/5-vis/vp-distributions/vp-geom-sankey-bump.R
# Ribbon widths are species-mean VP fractions multiplied by species Tjur R2.
# As in the original figure, the site random effect is excluded.

required_packages <- c("ggplot2", "ggsankey", "dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required package(s): ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}
suppressPackageStartupMessages({
  library(ggplot2)
  # This ggsankey version calls rowwise() from the search path.
  library(dplyr)
})

model_pattern <- "2026-11-02_06-55-21"
scaled <- TRUE
model_dir <- file.path("data-new", "2-hmsc", "derived", "models")
figure_path <- file.path("scripts", "5-vis", "vp-distributions", "vp-geom-sankey-bump.png")
atlas_labels <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")
variable_labels <- c(
  tmean_breeding = "Temperature",
  prec_breeding = "Precipitation",
  hh = "Land-use heterogeneity",
  perc_urban = "Urban (% coverage)",
  perc_cropland = "Cropland (% coverage)",
  perc_pasture = "Pasture (% coverage)",
  perc_forest_deciduous = "Deciduous forest (% coverage)",
  perc_forest_evergreen = "Evergreen forest (% coverage)",
  perc_forest_mixed_other = "Mixed/other forest (% coverage)",
  perc_grass_shrub = "Grass/Shrubland (% coverage)"
)
variable_colours <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban (% coverage)" = "snow3",
  "Cropland (% coverage)" = "goldenrod1",
  "Pasture (% coverage)" = "darkorange",
  "Deciduous forest (% coverage)" = "#A6D96A",
  "Evergreen forest (% coverage)" = "#006400",
  "Mixed/other forest (% coverage)" = "#87CEEB",
  "Grass/Shrubland (% coverage)" = "#B39DDB"
)

model_folders <- list.dirs(model_dir, recursive = FALSE, full.names = FALSE)
model_folders <- model_folders[grepl(model_pattern, model_folders, fixed = TRUE)]
atlas_numbers <- sub(".*Atlas([0-9]+).*", "\\1", model_folders)
if (length(model_folders) != 3L || anyDuplicated(atlas_numbers) ||
    !setequal(atlas_numbers, names(atlas_labels))) {
  stop("Expected exactly one model for each of Atlases 1–3 matching ",
       model_pattern, " under ", model_dir, call. = FALSE)
}
model_folders <- model_folders[match(names(atlas_labels), atlas_numbers)]

species_names <- NULL
plot_data <- do.call(rbind, lapply(seq_along(model_folders), function(i) {
  model <- model_folders[[i]]
  vp_path <- file.path(model_dir, model, "Results",
                       paste0(model, "parameter_estimates_VP_.csv"))
  raw <- read.csv(vp_path, check.names = FALSE, row.names = 1)
  required_rows <- c(names(variable_labels), "TjurR2")
  if (!all(required_rows %in% rownames(raw))) {
    stop("Missing VP rows in ", vp_path, ": ",
         paste(setdiff(required_rows, rownames(raw)), collapse = ", "))
  }
  # Fail rather than silently omit environmental predictors added to the model.
  environmental_rows <- setdiff(rownames(raw), c("TjurR2", "Random: site"))
  if (!setequal(environmental_rows, names(variable_labels))) {
    stop("Update variable labels/colours for unrecognised VP rows in ", vp_path,
         ": ", paste(setdiff(environmental_rows, names(variable_labels)), collapse = ", "))
  }
  if (is.null(species_names)) species_names <<- colnames(raw)
  if (!setequal(species_names, colnames(raw))) {
    stop("Species differ between atlases; reconcile species before comparing VP.")
  }
  values <- as.matrix(raw[names(variable_labels), species_names, drop = FALSE])
  tjur_r2 <- as.numeric(raw["TjurR2", species_names])
  if (!is.numeric(values) || any(!is.finite(values)) ||
      any(!is.finite(tjur_r2))) {
    stop("Non-numeric or missing VP/Tjur R2 values in ", vp_path)
  }
  if (isTRUE(scaled)) values <- sweep(values, 2, tjur_r2, "*")
  if (any(values < 0)) stop("Negative ribbon widths in ", vp_path)
  data.frame(
    atlas = unname(atlas_labels[[i]]),
    variable = unname(variable_labels),
    mean_VP = rowMeans(values)
  )
}))
plot_data$atlas <- factor(plot_data$atlas, levels = unname(atlas_labels))
plot_data$variable <- factor(plot_data$variable, levels = unname(variable_labels))

percent_labels <- function(x) {
  paste0(sub("\\.0$", "", sprintf("%.1f", x * 100)), "%")
}
# Retain the original 0–10% axis, extending it if the plotted totals require it.
y_max <- max(0.1, ceiling(max(tapply(plot_data$mean_VP, plot_data$atlas, sum)) /
                           0.025) * 0.025)
vp_plot <- ggplot(plot_data, aes(x = atlas, node = variable,
                                fill = variable, value = mean_VP)) +
  ggsankey::geom_sankey_bump(space = 0, type = "alluvial", alpha = 0.8,
                            colour = "black", linewidth = 0.25) +
  scale_y_continuous(limits = c(0, y_max), breaks = seq(0, y_max, 0.025),
                     minor_breaks = seq(0, y_max, 0.0125), labels = percent_labels) +
  scale_fill_manual(values = variable_colours) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = NULL, fill = "Variable", y = if (isTRUE(scaled)) {
    "Absolute variance explained (%)"
  } else {
    "Relative variance explained (%)"
  }) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_line(colour = "grey90", linewidth = 0.25),
    plot.background = element_rect(fill = "grey99", colour = NA),
    axis.text.x = element_text(face = "bold", size = 12)
  )

dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)
ggsave(figure_path, vp_plot, width = 11, height = 6, dpi = 300, bg = "grey99")
message("Saved figure: ", figure_path)
