# Run from the repository root:
# Rscript scripts/5-vis/spatial-latents/climate-spatial-residualisation.R
# Exploratory first-stage spatial detrending, not a refit of the bird models.
# Objects retained in memory: climate_decomposition, smooth_summary,
# smooth_fits, variogram_summary, pair_counts, common_sites, neighbour_graph.

required <- c("mgcv", "sf", "spdep", "ggplot2", "cowplot")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(mgcv); library(sf); library(ggplot2); library(cowplot)})
set.seed(20261105)
output_dir <- file.path("scripts", "5-vis", "spatial-latents")
prefix <- "climate-spatial-residualisation"
input_path <- file.path("data-new", "2-hmsc", "derived",
                        "preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData")
e <- new.env(parent = emptyenv()); load(input_path, envir = e)
X <- e$X; design <- e$Design
stopifnot(identical(rownames(X), rownames(design)), !anyDuplicated(rownames(X)))
variables <- c("tmean_breeding", "prec_breeding")
variable_labels <- c(tmean_breeding = "Temperature", prec_breeding = "Precipitation")
unit_labels <- c(tmean_breeding = "Temperature (°C)", prec_breeding = "Precipitation (mm/season)")
atlas_labels <- c("1970s", "1990s", "2010s")
basis_dimensions <- c(50L, 100L, 200L)
main_k <- 100L
sites_by_atlas <- split(as.character(design$site), as.character(design$atlas))
stopifnot(setequal(names(sites_by_atlas), as.character(1:3)),
          identical(as.integer(lengths(sites_by_atlas)), c(1747L, 1821L, 1875L)),
          !any(vapply(sites_by_atlas, anyDuplicated, integer(1))))
common_sites <- sort(Reduce(intersect, sites_by_atlas))
stopifnot(length(common_sites) == 1697L)
row_indices <- lapply(1:3, function(a) match(paste0(common_sites, "_", a), rownames(X)))
stopifnot(!anyNA(unlist(row_indices)))
coords <- as.matrix(design[row_indices[[1]], c("lon", "lat")])
stopifnot(all(is.finite(coords)), !anyDuplicated(data.frame(coords)))
for (a in 1:3) {
  stopifnot(identical(as.character(design$site[row_indices[[a]]]), common_sites),
            max(abs(coords - as.matrix(design[row_indices[[a]], c("lon", "lat")]))) < 1e-6,
            all(is.finite(as.matrix(X[row_indices[[a]], variables]))))
}
# S0 projects these coordinates to EPSG:23032 (metres).
map_points <- st_transform(st_as_sf(data.frame(site = common_sites, x = coords[, 1], y = coords[, 2]),
                                    coords = c("x", "y"), crs = 23032), 4326)
lonlat <- st_coordinates(map_points)
spatial_data <- data.frame(east_km = (coords[, 1] - mean(coords[, 1])) / 1000,
                           north_km = (coords[, 2] - mean(coords[, 2])) / 1000)
# Union of directed four-nearest-neighbour links, then row standardise.
neighbour_graph <- suppressWarnings(spdep::knn2nb(spdep::knearneigh(coords, k = 4), sym = TRUE))
weights <- spdep::nb2listw(neighbour_graph, style = "W")
graph_components <- spdep::n.comp.nb(neighbour_graph)$nc
moran_i <- function(value) {
  z <- value - mean(value)
  sum(z * spdep::lag.listw(weights, z)) / sum(z^2) * length(z) / spdep::Szero(weights)
}
stopifnot(abs(spdep::Szero(weights) - length(common_sites)) < 1e-8)

summary_rows <- list(); decomposition_rows <- list(); smooth_fits <- list()
index <- 0L
for (a in 1:3) for (variable in variables) for (k in basis_dimensions) {
  index <- index + 1L
  dat <- spatial_data
  dat$value <- X[row_indices[[a]], variable]
  fit <- mgcv::gam(value ~ s(east_km, north_km, bs = "tp", k = k),
                    data = dat, family = gaussian(), method = "REML")
  if (!isTRUE(fit$converged) || !identical(fit$outer.info$conv, "full convergence")) {
    stop("GAM failed to converge: Atlas ", a, ", ", variable, ", k=", k)
  }
  fitted_trend <- as.numeric(predict(fit, type = "response"))
  residual <- dat$value - fitted_trend
  stopifnot(all(is.finite(fitted_trend)), all(is.finite(residual)),
            max(abs(dat$value - fitted_trend - residual)) < 1e-10,
            abs(mean(residual)) < 1e-7)
  kc <- mgcv::k.check(fit, subsample = nrow(dat), n.rep = 200)
  fs <- summary(fit)
  original_sd <- sd(dat$value); residual_sd <- sd(residual)
  # Check the descriptive Moran calculation against the package implementation.
  for (v in list(dat$value, residual)) {
    direct_i <- spdep::moran(v, weights, length(v), spdep::Szero(weights))$I
    stopifnot(abs(moran_i(v) - direct_i) < 1e-12)
  }
  summary_rows[[index]] <- data.frame(
    atlas = a, variable = variable, k = k, n_cells = nrow(dat),
    original_mean = mean(dat$value), original_sd = original_sd,
    residual_mean = mean(residual), residual_sd = residual_sd,
    residual_sd_percent = 100 * residual_sd / original_sd,
    residual_variance_percent = 100 * var(residual) / var(dat$value),
    original_min = min(dat$value), original_max = max(dat$value),
    residual_min = min(residual), residual_max = max(residual),
    edf = unname(fs$s.table[1, "edf"]), basis_max_edf = k - 1,
    k_index = unname(kc[1, "k-index"]),
    original_moran = moran_i(dat$value), residual_moran = moran_i(residual),
    gam_converged = fit$converged, outer_convergence = fit$outer.info$conv)
  decomposition_rows[[index]] <- data.frame(
    site = common_sites, atlas = a, variable = variable, k = k,
    lon = lonlat[, 1], lat = lonlat[, 2],
    original = dat$value, trend = fitted_trend, residual = residual)
  smooth_fits[[paste(a, variable, k, sep = "::")]] <- fit
  message("Atlas ", a, "; ", variable, "; k=", k,
          ": residual SD ", round(100 * residual_sd / original_sd, 1),
          "% of original; Moran I ", round(moran_i(residual), 3), ".")
}
smooth_summary <- do.call(rbind, summary_rows)
climate_decomposition <- do.call(rbind, decomposition_rows)
stopifnot(nrow(smooth_summary) == 18L,
          nrow(climate_decomposition) == 18L * 1697L,
          all(smooth_summary$residual_variance_percent >= 0 &
                smooth_summary$residual_variance_percent <= 100 + 1e-6))

# Reuse unordered cell pairs and bins for all fields. All directions included.
pair_indices <- which(lower.tri(matrix(FALSE, nrow(coords), nrow(coords))), arr.ind = TRUE)
pair_distance_km <- as.vector(dist(coords)) / 1000
stopifnot(nrow(pair_indices) == choose(1697, 2),
          max(abs(pair_distance_km[c(1, 100, nrow(pair_indices))] -
                    sqrt(rowSums((coords[pair_indices[c(1, 100, nrow(pair_indices)), 1], , drop = FALSE] -
                                    coords[pair_indices[c(1, 100, nrow(pair_indices)), 2], , drop = FALSE])^2)) / 1000)) < 1e-9)
pair_bins <- cut(pair_distance_km, breaks = seq(0, 200, 10), right = TRUE, include.lowest = FALSE)
pair_counts <- data.frame(bin = 1:20, distance_km = seq(5, 195, 10),
                          n_pairs = as.integer(table(pair_bins)))
variogram_summary <- do.call(rbind, lapply(1:3, function(a) {
  do.call(rbind, lapply(variables, function(variable) {
    # Use explicit indexing rather than non-standard evaluation of loop names.
    dat <- climate_decomposition[climate_decomposition$atlas == a &
      climate_decomposition$variable == variable & climate_decomposition$k == main_k, ]
    stopifnot(identical(dat$site, common_sites))
    do.call(rbind, lapply(c("original", "residual"), function(field) {
      v <- dat[[field]]
      semivariance <- (v[pair_indices[, 1]] - v[pair_indices[, 2]])^2 / 2
      means <- as.numeric(tapply(semivariance, pair_bins, mean))
      data.frame(atlas = a, variable = variable, field = field, pair_counts,
                 semivariance = means, normalised_semivariance = means / var(v))
    }))
  }))
}))

mainland_bbox <- c(xmin = 7.95, xmax = 13, ymin = 54.4, ymax = 57.85)
bornholm_bbox <- c(xmin = 14.65, xmax = 15.2, ymin = 54.95, ymax = 55.32)
inside <- function(dat, b) with(dat, lon >= b[["xmin"]] & lon <= b[["xmax"]] &
                                 lat >= b[["ymin"]] & lat <= b[["ymax"]])
stopifnot(all(inside(climate_decomposition, mainland_bbox) | inside(climate_decomposition, bornholm_bbox)))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_png <- function(plot, suffix, width, height) {
  path <- file.path(output_dir, paste0(prefix, "-", suffix, ".png"))
  ggsave(path, plot, width = width, height = height, dpi = 300, bg = "white")
  message("Saved: ", path)
}
for (variable in variables) {
  z <- climate_decomposition[climate_decomposition$variable == variable & climate_decomposition$k == main_k, ]
  absolute_limits <- range(z$original, z$trend)
  residual_limit <- max(abs(z$residual))
  scale_map <- function(field) {
    if (field == "residual") scale_fill_gradient2(low = "#2C7FB8", mid = "#FFF7E6", high = "#D95F0E",
      midpoint = 0, limits = c(-residual_limit, residual_limit), name = paste("Residual", unit_labels[[variable]]))
    else scale_fill_gradientn(colours = c("#2C7FB8", "#FFF7E6", "#D95F0E"), limits = absolute_limits,
                               name = unit_labels[[variable]])
  }
  map <- function(d, field, bbox, border = FALSE, size = 1.05) {
    d$value <- d[[field]]
    ggplot(d, aes(lon, lat, fill = value)) +
      geom_point(shape = 21, colour = "white", stroke = 0.10, size = size) + scale_map(field) +
      coord_quickmap(xlim = bbox[c("xmin", "xmax")], ylim = bbox[c("ymin", "ymax")], expand = FALSE) +
      theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = if (border) element_rect(fill = NA, colour = "grey35", linewidth = 0.5) else element_blank())
  }
  fields <- c("original", "trend", "residual")
  panels <- unlist(lapply(1:3, function(a) lapply(fields, function(field) {
    d <- z[z$atlas == a, ]
    main <- map(d, field, mainland_bbox)
    inset <- map(d, field, bornholm_bbox, border = TRUE, size = 1.3)
    body <- ggdraw(main) + draw_plot(inset, x = 0.79, y = 0.68, width = 0.18, height = 0.22)
    field_label <- c(original = "Original", trend = "Fitted spatial trend", residual = "Residual")[[field]]
    strip <- ggdraw() + draw_label(paste0(atlas_labels[[a]], " | ", field_label), size = 11, fontface = "bold") +
      theme(plot.background = element_rect(fill = "grey92", colour = "white"))
    plot_grid(strip, body, ncol = 1, rel_heights = c(0.07, 1))
  })), recursive = FALSE)
  legends <- lapply(c("original", "residual"), function(field) get_legend(
    map(z, field, mainland_bbox) + theme(legend.position = "bottom") +
      guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = grid::unit(11, "cm")))))
  header <- ggdraw() + draw_label(paste0(variable_labels[[variable]], ": original climate, spatial trend and residual"),
                                  x = 0.01, hjust = 0, fontface = "bold", size = 18)
  footer <- ggdraw() + draw_label(
    "Same 1,697 cells in each atlas; separate REML thin-plate smooths (k = 100). Trend includes the atlas mean.\nOriginal/trend share a scale; residuals use a separate zero-centred scale. Bornholm insets; no bird models refitted.",
    x = 0.01, hjust = 0, size = 10)
  composite <- plot_grid(header, plot_grid(plotlist = panels, ncol = 3),
                         plot_grid(plotlist = legends, nrow = 1), footer, ncol = 1,
                         rel_heights = c(0.06, 1, 0.075, 0.065))
  save_png(composite, paste0(variable, "-maps"), 15, 17)
}

atlas_colours <- c("1970s" = "#0072B2", "1990s" = "#D55E00", "2010s" = "#009E73")
sensitivity_data <- rbind(
  data.frame(smooth_summary, metric = "Residual SD (% of original)", value = smooth_summary$residual_sd_percent),
  data.frame(smooth_summary, metric = "Residual variance (% of original)", value = smooth_summary$residual_variance_percent),
  data.frame(smooth_summary, metric = "Residual Moran's I", value = smooth_summary$residual_moran))
sensitivity_data$climate <- factor(variable_labels[sensitivity_data$variable], levels = unname(variable_labels))
sensitivity_data$atlas_label <- factor(atlas_labels[sensitivity_data$atlas], levels = atlas_labels)
sensitivity_data$metric <- factor(sensitivity_data$metric, levels = unique(sensitivity_data$metric))
p <- ggplot(sensitivity_data, aes(edf, value, colour = atlas_label, group = atlas_label)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  geom_text(data = sensitivity_data[sensitivity_data$atlas_label == "1970s", ],
            aes(label = k), vjust = -0.8, size = 3, show.legend = FALSE) +
  facet_grid(metric ~ climate, scales = "free_y") + scale_colour_manual(values = atlas_colours) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(title = "Sensitivity to spatial smooth flexibility", subtitle = "Labels on blue curves show k (50, 100, 200); horizontal axis is fitted effective degrees of freedom",
       x = "Smooth effective degrees of freedom", y = NULL, colour = "Atlas",
       caption = "REML selects smoothing separately for every fit. More residual removal is not automatically a better adjustment.\nMoran's I uses a common symmetrised four-nearest-neighbour graph with row-standardised weights.") +
  theme_minimal(base_size = 12) + theme(legend.position = "top", panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"), plot.caption = element_text(hjust = 0))
save_png(p, "sensitivity", 12, 11)

for (variable in variables) {
  z <- variogram_summary[variogram_summary$variable == variable, ]
  z$atlas_label <- factor(atlas_labels[z$atlas], levels = atlas_labels)
  z$field <- factor(z$field, levels = c("original", "residual"), labels = c("Original", "Residual"))
  base_plot <- function(response, y_label) {
    ggplot(z, aes(distance_km, .data[[response]], colour = field, group = field)) +
      geom_line(linewidth = 0.7, na.rm = TRUE) +
      geom_point(aes(shape = n_pairs < 100), size = 1.3, na.rm = TRUE) +
      scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1), guide = "none") +
      scale_colour_manual(values = c("Original" = "#0072B2", "Residual" = "#D55E00")) +
      facet_wrap(~atlas_label, nrow = 1) +
      labs(x = "Distance (km)", y = y_label, colour = NULL) +
      theme_minimal(base_size = 11) + theme(legend.position = "top", panel.grid.minor = element_blank())
  }
  physical_units <- if (variable == "tmean_breeding") "Semivariance (°C²)" else "Semivariance ((mm/season)²)"
  physical <- base_plot("semivariance", physical_units)
  normalised <- base_plot("normalised_semivariance", "Semivariance / own field variance") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55")
  counts <- unique(z[, c("atlas_label", "distance_km", "n_pairs")])
  count_plot <- ggplot(counts, aes(distance_km, n_pairs)) + geom_col(fill = "grey55", width = 8) +
    facet_wrap(~atlas_label, nrow = 1) + labs(x = "Distance (km)", y = "Unordered cell pairs") +
    theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
  header <- ggdraw() + draw_label(paste0(variable_labels[[variable]], ": spatial structure before and after detrending"),
                                  x = 0.01, hjust = 0, fontface = "bold", size = 16)
  footer <- ggdraw() + draw_label(
    "k = 100; same cell pairs in every atlas and field. Bins (0,10], …, (190,200] km. Open points: fewer than 100 pairs.\nNormalisation uses each field's own sample variance. A value near one does not establish independence.",
    x = 0.01, hjust = 0, size = 10)
  save_png(plot_grid(header, physical, normalised, count_plot, footer, ncol = 1,
                     rel_heights = c(0.12, 1, 1, 0.8, 0.15)), paste0(variable, "-variograms"), 13, 12)
}
message("Validated: 18 converged GAMs, 1,697 aligned cells, exact decompositions, common spatial weights and pair bins.")
print(smooth_summary[smooth_summary$k == main_k,
  c("atlas", "variable", "original_sd", "residual_sd", "residual_variance_percent", "edf", "original_moran", "residual_moran")], row.names = FALSE)
