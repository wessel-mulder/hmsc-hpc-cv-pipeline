# Predicted richness relative to the maximum within each atlas period.
# Uses the same cached richness inputs, map extents and inset layout as GWR Fig. 3.
library(ggplot2)
library(dplyr)
library(sf)
library(cowplot)
library(readr)
inputs <- readRDS("intermediate-outputs/gwr-output/fig10-gwr-richness-turnover/gwr/models/2026-03-13-gwr-analysis-inputs-100km.rds")
periods <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")
values <- bind_rows(lapply(names(inputs$richness_dfs), function(a) {
  inputs$richness_dfs[[a]] |> transmute(site, period = periods[[a]], richness,
    period_maximum = max(richness, na.rm = TRUE),
    relative_richness = 100 * richness / period_maximum)
}))
stopifnot(all(is.finite(values$relative_richness)), all(values$relative_richness <= 100))
grid_path <- "data-new/derived/atlas-grid/dof-5km-grid/DK5km_ED50grid_approx_kvadrkod_DOF.shp"
shape_sf <- st_read(grid_path, quiet = TRUE) |>
  st_transform(25832)

mainland_bbox <- st_bbox(
  c(xmin = 430000, xmax = 745000, ymin = 6040000, ymax = 6410000),
  crs = st_crs(25832)
)
bornholm_bbox <- st_bbox(
  c(xmin = 858000, xmax = 897000, ymin = 6107000, ymax = 6145000),
  crs = st_crs(25832)
)

mainland_width <- as.numeric(mainland_bbox["xmax"] - mainland_bbox["xmin"])
mainland_height <- as.numeric(mainland_bbox["ymax"] - mainland_bbox["ymin"])
bornholm_width <- as.numeric(bornholm_bbox["xmax"] - bornholm_bbox["xmin"])
bornholm_height <- as.numeric(bornholm_bbox["ymax"] - bornholm_bbox["ymin"])

inset_w <- bornholm_width / mainland_width
inset_h <- bornholm_height / mainland_height
inset_margin_x <- 0.2
inset_margin_y <- 0.2

map_theme <- function() {
  theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
}

facet_strip <- function(label) {
  ggdraw() +
    draw_label(
      label,
      x = 0.5,
      y = 0.52,
      hjust = 0.5,
      vjust = 0.5,
      fontface = "bold",
      size = 12,
      colour = "black"
    ) +
    theme(
      plot.background = element_rect(fill = "grey92", colour = "white"),
      plot.margin = margin(0, 0, 1.2, 0)
    )
}

period_strip <- function(label) {
  ggdraw() +
    draw_label(
      label,
      x = 0.5,
      y = 0.5,
      hjust = 0.5,
      vjust = 0.5,
      angle = 90,
      fontface = "bold",
      size = 12,
      colour = "black"
    ) +
    theme(
      plot.background = element_rect(fill = "grey92", colour = "white"),
      plot.margin = margin(0, 1.2, 0, 0)
    )
}


# Match the older richness maps in misc-figures/sfig1-richness-maps-accuracies.R.
richness_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
colour_scale <- function() scale_colour_gradientn(colours = richness_palette(100), limits = c(0, 100),
  breaks = seq(0, 100, 25), name = "Relative richness\n(% of period maximum)")
base_map <- function(points, bbox, border = FALSE, point_size = 1.75) {
  ggplot(points) + geom_sf(aes(colour = relative_richness), size = point_size) +
    colour_scale() +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE) +
    map_theme() + theme(panel.border = if (border)
      element_rect(colour = "grey30", fill = NA, linewidth = 0.5) else element_blank())
}
points <- shape_sf |> inner_join(values, by = c("kvadratkod" = "site"))
stopifnot(nrow(points) == nrow(values))
points <- suppressWarnings(st_point_on_surface(points))
map_panel <- function(period_name) {
  p <- points |> filter(period == period_name)
  main <- base_map(p, mainland_bbox)
  inset <- base_map(p, bornholm_bbox, TRUE, 2.05)
  ggdraw(main) + draw_plot(inset, x = 1 - inset_w - inset_margin_x,
    y = 1 - inset_h - inset_margin_y, width = inset_w, height = inset_h)
}
rows <- lapply(unname(periods), function(p) plot_grid(period_strip(p), map_panel(p),
  nrow = 1, rel_widths = c(0.055, 1), align = "hv", axis = "tblr"))
header <- plot_grid(NULL, facet_strip("Relative predicted species richness"),
  nrow = 1, rel_widths = c(0.055, 1))
maps <- plot_grid(header, plot_grid(plotlist = rows, ncol = 1),
  ncol = 1, rel_heights = c(0.035, 1))
legend <- get_legend(base_map(points, mainland_bbox) + theme(legend.position = "right",
  legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  guides(colour = guide_colourbar(barheight = grid::unit(2.6, "in"))))
figure <- plot_grid(maps, legend, nrow = 1, rel_widths = c(1, 0.4))
out <- "manuscript-outputs/figures/exploratory"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
png <- file.path(out, "2026-03-13-relative-predicted-richness-by-atlas.png")
ggsave(png, figure, width = 8.5, height = 14, units = "in", dpi = 300, bg = "white")
write_csv(values, file.path(out, "2026-03-13-relative-predicted-richness-by-atlas.csv"))
print(values |> group_by(period) |> summarise(n = n(), max_richness = first(period_maximum),
  min_percent = min(relative_richness), max_percent = max(relative_richness)))
message("Saved: ", png)
