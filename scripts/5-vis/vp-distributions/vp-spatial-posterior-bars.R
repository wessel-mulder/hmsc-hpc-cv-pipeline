# Run from the repository root:
# Rscript scripts/5-vis/vp-distributions/vp-spatial-posterior-bars.R
suppressPackageStartupMessages(library(ggplot2))
analysis_dir <- "scripts/4-post/compare-scaled-vp-changes/spatial"
comparison_file <- file.path(analysis_dir, "comparison.rds")
if (!file.exists(comparison_file)) {
  stop("Run Rscript scripts/4-post/compare-scaled-vp-changes.R --spatial first.")
}
comparison <- readRDS(comparison_file)
stopifnot(isTRUE(comparison$provenance$settings$spatial),
  identical(comparison$provenance$settings$model_prefix, "2026-11-02_06-55-21"))
atlas_labels <- c("1970s", "1990s", "2010s")
variable_labels <- c(
  tmean_breeding = "Temperature", prec_breeding = "Precipitation",
  hh = "Land-use\nheterogeneity", perc_urban = "Urban\ncover",
  perc_cropland = "Cropland\ncover", perc_pasture = "Pasture\ncover",
  perc_forest_deciduous = "Deciduous\nforest cover",
  perc_forest_evergreen = "Evergreen\nforest cover",
  perc_forest_mixed_other = "Mixed/other\nforest cover",
  perc_grass_shrub = "Grass/shrubland\ncover"
)
variables <- names(variable_labels)
atlases <- lapply(1:3, function(atlas) {
  result <- readRDS(file.path(analysis_dir, paste0("atlas", atlas, ".rds")))
  stopifnot(result$atlas == atlas, length(result$species) == 155L,
    setequal(result$community_summary$variable, variables),
    identical(result$provenance$input_md5, unname(comparison$provenance$model_md5[atlas])))
  result$community_summary$mean[match(variables, result$community_summary$variable)]
})
bars <- expand.grid(variable = variables, atlas = atlas_labels, stringsAsFactors = FALSE)
bars$mean <- unlist(atlases, use.names = FALSE)
offsets <- c(-.25, 0, .25)
bars$x <- match(bars$variable, variables) + offsets[match(bars$atlas, atlas_labels)]
bars$atlas <- factor(bars$atlas, levels = atlas_labels)

pairs <- list(c(1L, 2L), c(2L, 3L), c(1L, 3L))
pair_labels <- vapply(pairs, function(pair) paste(atlas_labels[pair], collapse = " to "), character(1))
lines <- comparison$community_summary
stopifnot(nrow(lines) == 30L, setequal(lines$variable, variables),
  setequal(lines$pair, pair_labels), !anyDuplicated(lines[c("variable", "pair")]))
lines$direction <- factor(ifelse(lines$joint_lower > 0, "Increase",
  ifelse(lines$joint_upper < 0, "Decrease", "Uncertain")), levels = c("Increase", "Decrease", "Uncertain"))
lines$diagnostics <- ifelse(lines$source_diagnostic_status == "pass", "Passed", "Flagged")
pair_index <- match(lines$pair, pair_labels)
centres <- match(lines$variable, variables)
earlier <- c(1L, 2L, 1L)[pair_index]
later <- c(2L, 3L, 3L)[pair_index]
lines$start <- centres + offsets[earlier]
lines$end <- centres + offsets[later]
bar_max <- tapply(bars$mean, bars$variable, max)
gap <- .075 * max(bars$mean)
lines$y <- unname(bar_max[lines$variable]) + gap * pair_index
lines$tick <- lines$y - .26 * gap
bar_means <- matrix(bars$mean, nrow = length(variables))
bar_differences <- bar_means[cbind(centres, later)] - bar_means[cbind(centres, earlier)]
stopifnot(max(abs(lines$mean - bar_differences)) < 1e-10)

caption <- paste(
  "Bars show posterior mean Tjur R²-scaled VP across 155 species in the 2 November spatial models.",
  "Lines connect atlas pairs. Green and red mark joint 95% intervals entirely above and below zero. Grey intervals include zero.",
  "The joint intervals cover 30 environmental comparisons. The spatial component is included in the partition.",
  sep = "\n"
)
if (any(lines$diagnostics == "Flagged")) {
  caption <- paste(caption, "Dashed lines mark comparisons with flagged posterior diagnostics.", sep = "\n")
}
p <- ggplot(bars, aes(x, mean, fill = atlas)) +
  geom_col(width = .205) +
  geom_segment(data = lines, aes(x = start, xend = end, y = y, yend = y,
    colour = direction, linetype = diagnostics), inherit.aes = FALSE, linewidth = .8) +
  geom_segment(data = lines, aes(x = start, xend = start, y = y, yend = tick,
    colour = direction), inherit.aes = FALSE, linewidth = .65, show.legend = FALSE) +
  geom_segment(data = lines, aes(x = end, xend = end, y = y, yend = tick,
    colour = direction), inherit.aes = FALSE, linewidth = .65, show.legend = FALSE) +
  scale_fill_manual(values = c("1970s" = "#C7DDEB", "1990s" = "#77A8C8", "2010s" = "#32668A")) +
  scale_colour_manual(values = c(Increase = "#20853F", Decrease = "#B2182B", Uncertain = "#999999")) +
  scale_linetype_manual(values = c(Passed = "solid", Flagged = "dashed"), guide = "none") +
  scale_x_continuous(breaks = seq_along(variables), labels = unname(variable_labels),
    expand = expansion(add = .5)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, .045))) +
  labs(title = "Environmental importance across atlases",
    subtitle = "Community means and posterior changes in the 2 November spatial models",
    x = NULL, y = "Community mean Tjur R²-scaled VP (%)", fill = "Atlas",
    colour = "Direction of change", caption = caption) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", legend.justification = "left",
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10, colour = "grey20"),
    plot.title = element_text(face = "bold", size = 18),
    plot.caption = element_text(hjust = 0, margin = margin(t = 12)),
    plot.margin = margin(14, 18, 14, 14))

output <- "scripts/5-vis/vp-distributions/vp-spatial-posterior-bars.png"
ggsave(output, p, width = 15, height = 6.5, dpi = 300, bg = "white")
print(table(lines$direction, lines$diagnostics))
message("Saved: ", normalizePath(output))
