# Run from the repository root: Rscript scripts/5-vis/spatial-latents/plot-environment-correlations.R
# Descriptive Pearson correlations between all ten model predictors, using
# exactly the same grid cells in each atlas. Only PNG output is generated.
suppressPackageStartupMessages(library(ggplot2))

model_root <- "data-new/2-hmsc/derived/models"
model_stamp <- "2026-11-02_06-55-21"
output_dir <- "scripts/5-vis/spatial-latents"
variable_labels <- c(
  tmean_breeding = "Temperature",
  prec_breeding = "Precipitation",
  hh = "Land-use heterogeneity",
  perc_cropland = "Cropland",
  perc_urban = "Urban",
  perc_pasture = "Pasture",
  perc_grass_shrub = "Grass/shrubland",
  perc_forest_deciduous = "Deciduous forest",
  perc_forest_evergreen = "Evergreen forest",
  perc_forest_mixed_other = "Mixed/other forest"
)
variables <- names(variable_labels)

models <- lapply(1:3, function(atlas) {
  path <- file.path(model_root,
    sprintf("%s_Atlas%d_MinOccs5_CoverageGoodAverage", model_stamp, atlas),
    "Models", "Unfitted", "unfitted_models.RData")
  e <- new.env(parent = emptyenv())
  load(path, e)
  m <- e$models[[1]]
  stopifnot(setequal(names(m$XData), variables),
    identical(rownames(m$XData), rownames(m$studyDesign)),
    !anyDuplicated(as.character(m$studyDesign$site)))
  list(x = m$XData[, variables, drop = FALSE],
    sites = as.character(m$studyDesign$site))
})
common_sites <- sort(Reduce(intersect, lapply(models, `[[`, "sites")))
stopifnot(length(common_sites) > 2L)

correlation_matrices <- lapply(models, function(m) {
  ix <- match(common_sites, m$sites)
  stopifnot(!anyNA(ix), identical(m$sites[ix], common_sites))
  x <- as.matrix(m$x[ix, , drop = FALSE])
  stopifnot(all(is.finite(x)), all(apply(x, 2, sd) > 0))
  cor(x, method = "pearson", use = "everything")
})
atlas_labels <- c("Atlas 1 (1970s)", "Atlas 2 (1990s)", "Atlas 3 (2010s)")
plot_data <- do.call(rbind, lapply(1:3, function(a) {
  d <- as.data.frame(as.table(correlation_matrices[[a]]), stringsAsFactors = FALSE)
  names(d) <- c("row", "column", "r")
  d$atlas <- atlas_labels[a]
  d
}))
plot_data$row <- factor(plot_data$row, levels = rev(variables))
plot_data$column <- factor(plot_data$column, levels = variables)
plot_data$atlas <- factor(plot_data$atlas, levels = atlas_labels)
plot_data$label <- sprintf("%.2f", plot_data$r)
plot_data$label[abs(plot_data$r) < 0.005] <- "0.00"

p <- ggplot(plot_data, aes(column, row, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = label, colour = abs(r) > 0.65), size = 3.0) +
  scale_colour_manual(values = c("FALSE" = "#222222", "TRUE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B35806",
    midpoint = 0, limits = c(-1, 1), breaks = seq(-1, 1, 0.5), name = "Pearson r") +
  scale_x_discrete(labels = variable_labels, expand = c(0, 0)) +
  scale_y_discrete(labels = variable_labels, expand = c(0, 0)) +
  facet_wrap(~atlas, nrow = 1) +
  coord_fixed() +
  labs(title = "Correlations between environmental predictors",
    subtitle = sprintf("The same %s grid cells in all three atlases; all ten model predictors",
      format(length(common_sites), big.mark = ",")),
    x = NULL, y = NULL,
    caption = paste("Descriptive correlations across cells; no significance tests. Cover variables are proportions of the same landscape.",
      "Inputs: November 2026 model-definition objects. Predictor standardisation does not change Pearson correlations.", sep = "\n")) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, colour = "#222222"),
    axis.text.y = element_text(colour = "#222222"),
    strip.text = element_text(face = "bold", size = 14, margin = margin(b = 12)),
    panel.spacing = grid::unit(1.0, "lines"),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(margin = margin(b = 16)),
    plot.caption = element_text(hjust = 0, margin = margin(t = 15)),
    legend.position = "bottom", legend.key.width = grid::unit(1.2, "cm"),
    plot.margin = margin(15, 20, 15, 15))

output_path <- file.path(output_dir, "environment-predictor-correlations-atlases-1-2-3.png")
ggsave(output_path, p, width = 21, height = 9.2, units = "in", dpi = 300, bg = "white")
message("Saved: ", output_path, " (", length(common_sites), " common cells)")
for (a in 1:3) {
  r <- correlation_matrices[[a]]
  ij <- which(upper.tri(r), arr.ind = TRUE)
  strongest <- ij[order(abs(r[ij]), decreasing = TRUE)[1:5], , drop = FALSE]
  cat("\n", atlas_labels[a], ": strongest five absolute correlations\n", sep = "")
  print(data.frame(variable_1 = rownames(r)[strongest[, 1]],
    variable_2 = colnames(r)[strongest[, 2]], r = round(r[strongest], 3)), row.names = FALSE)
}
