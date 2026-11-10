# Run from the project root:
# Rscript scripts/5-vis/spatial-latents/plot-spatial-latents.R
# Map Eta in the saved within-model alignment; do not align between atlases.
# These are latent-factor coordinates, not species-specific Eta %*% Lambda.

required <- c("sf", "ggplot2", "cowplot")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(sf); library(ggplot2); library(cowplot)})

model_stamp <- "2026-11-02_06-55-21"
model_root <- file.path("data-new", "2-hmsc", "derived", "models")
output_dir <- file.path("scripts", "5-vis", "spatial-latents")
atlas_labels <- c("1970s", "1990s", "2010s")
expected_cells <- c(1747L, 1821L, 1875L)
grid_path <- file.path("data-new", "0-pre", "raw", "DOF_Shapefiles_",
                       "DK5km_ED50grid_approx_kvadrkod_DOF.shp")
grid <- st_read(grid_path, quiet = TRUE)
if (anyDuplicated(grid$kvadratkod)) stop("Duplicate grid cell IDs.")
grid_centers <- suppressWarnings(st_centroid(st_transform(grid, 4326)))

extract_eta <- function(atlas) {
  model <- sprintf("%s_Atlas%d_MinOccs5_CoverageGoodAverage", model_stamp, atlas)
  path <- file.path(model_root, model, "Models", "Fitted",
                    "HPC_samples_0250_thin_100_chains_4.Rdata")
  if (!file.exists(path)) stop("Missing fitted model: ", path)
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  m <- e$fitted_model$posteriors
  if (length(m$postList) != 4L || any(lengths(m$postList) != 250L) ||
      m$nr != 1L || m$ny != expected_cells[[atlas]]) stop("Unexpected model dimensions: ", model)
  sites <- as.character(m$studyDesign$site)
  unit_ids <- levels(m$dfPi$site)
  if (anyDuplicated(sites) || !identical(as.character(m$dfPi$site), sites) ||
      !identical(as.integer(m$Pi[, 1]), match(sites, unit_ids)) ||
      !all(sites %in% grid$kvadratkod)) stop("Invalid site indexing: ", model)
  coords <- m$ranLevels[[1]]$s[sites, , drop = FALSE]
  if (ncol(coords) != 2L || any(!is.finite(as.matrix(coords)))) stop("Invalid coordinates.")
  # S0_Data_Definitions.R projects Design coordinates to EPSG:23032.
  points <- st_as_sf(data.frame(site = sites, X = coords[, 1], Y = coords[, 2]),
                     coords = c("X", "Y"), crs = 23032)
  points <- st_transform(points, 4326)
  grid_match <- match(sites, grid$kvadratkod)
  offset <- as.numeric(st_distance(points, grid_centers[grid_match, ], by_element = TRUE))
  if (any(!is.finite(offset)) || max(offset) > 1000) {
    stop("Fitted coordinates disagree with grid centroids; check CRS/site mapping.")
  }
  lonlat <- st_coordinates(points)

  # Retain chain and iteration dimensions while extracting all saved draws.
  eta <- array(NA_real_, dim = c(250L, 4L, length(sites), 5L))
  for (chain in 1:4) for (iteration in 1:250) {
    draw <- m$postList[[chain]][[iteration]]$Eta[[1]]
    if (!identical(dim(draw), c(length(unit_ids), 5L)) || any(!is.finite(draw))) {
      stop("Invalid Eta dimensions/values: ", model)
    }
    eta[iteration, chain, , ] <- draw[m$Pi[, 1], , drop = FALSE]
  }
  mean_eta <- apply(eta, c(3, 4), mean)
  sd_eta <- apply(eta, c(3, 4), sd)
  # Independent direct extraction checks for selected cells and every factor.
  for (site_i in unique(c(1L, ceiling(length(sites) / 2), length(sites)))) {
    for (factor_i in 1:5) {
      direct <- unlist(lapply(m$postList, function(chain) {
        vapply(chain, function(draw) draw$Eta[[1]][m$Pi[site_i, 1], factor_i], numeric(1))
      }), use.names = FALSE)
      stopifnot(abs(mean_eta[site_i, factor_i] - mean(direct)) < 1e-12,
                abs(sd_eta[site_i, factor_i] - sd(direct)) < 1e-12)
    }
  }
  if (any(!is.finite(mean_eta)) || any(!is.finite(sd_eta)) || any(sd_eta < 0)) {
    stop("Invalid posterior summaries.")
  }
  message("Atlas ", atlas, ": ", length(sites), " cells, five factors, 1,000 draws; ",
          "max grid-centroid offset = ", round(max(offset), 1), " m.")
  do.call(rbind, lapply(1:5, function(factor_i) data.frame(
    atlas = atlas, factor = factor_i, site = sites, lon = lonlat[, 1], lat = lonlat[, 2],
    posterior_mean = mean_eta[, factor_i], posterior_sd = sd_eta[, factor_i]
  )))
}

summaries <- do.call(rbind, lapply(1:3, extract_eta))
stopifnot(nrow(summaries) == sum(expected_cells) * 5L)
mean_limit <- max(abs(summaries$posterior_mean))
sd_limit <- max(summaries$posterior_sd)
mainland_bbox <- c(xmin = 7.95, xmax = 13.0, ymin = 54.4, ymax = 57.85)
bornholm_bbox <- c(xmin = 14.65, xmax = 15.2, ymin = 54.95, ymax = 55.32)
in_bbox <- function(d, b) with(d, lon >= b[["xmin"]] & lon <= b[["xmax"]] &
                                  lat >= b[["ymin"]] & lat <= b[["ymax"]])
if (!all(in_bbox(summaries, mainland_bbox) | in_bbox(summaries, bornholm_bbox))) {
  stop("Map extents would omit fitted cells.")
}
for (atlas in 1:3) if (!any(in_bbox(summaries[summaries$atlas == atlas, ], bornholm_bbox))) {
  stop("Missing Bornholm cells in Atlas ", atlas)
}

fill_scale <- function(statistic) {
  if (statistic == "posterior_mean") {
    scale_fill_gradient2(low = "#2C7FB8", mid = "#FFF7E6", high = "#D95F0E",
                          midpoint = 0, limits = c(-mean_limit, mean_limit),
                          name = "Posterior mean Eta")
  } else {
    scale_fill_gradient(low = "#FFF7E6", high = "#D95F0E", limits = c(0, sd_limit),
                        name = "Posterior SD of Eta")
  }
}
point_map <- function(dat, statistic, bbox, border = FALSE, point_size = 1.2) {
  dat$value <- dat[[statistic]]
  ggplot(dat, aes(lon, lat, fill = value)) +
    geom_point(shape = 21, colour = "white", stroke = 0.10, size = point_size) +
    fill_scale(statistic) +
    coord_quickmap(xlim = bbox[c("xmin", "xmax")], ylim = bbox[c("ymin", "ymax")],
                   expand = FALSE) +
    theme_void() +
    theme(legend.position = "none", plot.background = element_rect(fill = "white", colour = NA),
          panel.border = if (border) element_rect(colour = "grey35", fill = NA, linewidth = 0.55)
                         else element_blank(),
          plot.margin = margin(0, 0, 0, 0))
}
factor_panel <- function(dat, statistic, factor_i) {
  main <- point_map(dat, statistic, mainland_bbox)
  inset <- point_map(dat, statistic, bornholm_bbox, border = TRUE, point_size = 1.45)
  map <- ggdraw(main) + draw_plot(inset, x = 0.79, y = 0.68, width = 0.18, height = 0.22)
  strip <- ggdraw() + draw_label(paste("Factor", factor_i), fontface = "bold", size = 12) +
    theme(plot.background = element_rect(fill = "grey92", colour = "white"))
  plot_grid(strip, map, ncol = 1, rel_heights = c(0.07, 1))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (atlas in 1:3) for (statistic in c("posterior_mean", "posterior_sd")) {
  dat <- summaries[summaries$atlas == atlas, ]
  panels <- lapply(1:5, function(factor_i) {
    factor_panel(dat[dat$factor == factor_i, ], statistic, factor_i)
  })
  legend <- get_legend(point_map(dat, statistic, mainland_bbox) +
                         theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
                         guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                                       barwidth = grid::unit(12, "cm"))))
  panels[[6]] <- ggdraw() + draw_label("Within-atlas factors\n\nBornholm shown in each inset\n\nNo cross-atlas factor matching",
                                      size = 12, lineheight = 1.3, colour = "grey30")
  body <- plot_grid(plotlist = panels, ncol = 3)
  label <- if (statistic == "posterior_mean") "Posterior means" else "Posterior uncertainty (SD)"
  header <- ggdraw() +
    draw_label(paste0("Atlas ", atlas, " (", atlas_labels[[atlas]], ") | Spatial latent factors"),
               x = 0.02, y = 0.72, hjust = 0, fontface = "bold", size = 19) +
    draw_label(paste0(label, " | ", expected_cells[[atlas]], " fitted cells | November 2 models"),
               x = 0.02, y = 0.23, hjust = 0, size = 12)
  footer <- ggdraw() + draw_label(
    "4 chains × 250 retained draws; saved within-model alignment. Factor signs/order are arbitrary across atlases.\nSD is conditional on this alignment; these maps do not establish convergence or measure species-specific spatial effects.",
    x = 0.02, hjust = 0, size = 10, lineheight = 1.15)
  composite <- plot_grid(header, body, legend, footer, ncol = 1,
                          rel_heights = c(0.10, 1, 0.08, 0.07))
  filename <- sprintf("atlas%d-spatial-latents-%s.png", atlas, gsub("_", "-", statistic))
  ggsave(file.path(output_dir, filename), composite, width = 15, height = 12,
         dpi = 300, bg = "white")
  message("Saved: ", file.path(output_dir, filename))
}
message("Validated all 15 factors and their mean/SD maps; no fitted cells omitted.")
