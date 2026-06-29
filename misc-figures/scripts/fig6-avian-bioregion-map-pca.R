rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, cowplot, RColorBrewer, grid)

pattern <- "2026-03-13"
in_dir <- file.path("notebooks", "exploratory", "outputs", "bioregion-pca")
analysis_path <- file.path(in_dir, "models", paste0(pattern, "-avian-bioregion-pca-analysis.rds"))
out_dir <- file.path("misc-figures", "outputs", "main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis <- readRDS(analysis_path)

year_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")
reference_atlas <- "1"
cluster_levels <- paste("Cluster", 1:3)

# Re-label clusters for manuscript readability:
# old Cluster 1 -> new Cluster 1, old Cluster 3 -> new Cluster 2,
# old Cluster 2 -> new Cluster 3.
cluster_recode <- c(
  "Cluster 1" = "Cluster 1",
  "Cluster 3" = "Cluster 2",
  "Cluster 2" = "Cluster 3"
)
relabel_clusters <- function(x) {
  factor(unname(cluster_recode[as.character(x)]), levels = cluster_levels)
}

bioregion_assignments <- analysis$bioregion_assignments |>
  map(~ mutate(.x, bioregion = relabel_clusters(.data$bioregion)))
cluster_occupancy <- analysis$cluster_occupancy |>
  mutate(bioregion = relabel_clusters(.data$bioregion)) |>
  group_by(period, bioregion) |>
  summarise(n_sites = sum(.data$n_sites), .groups = "drop") |>
  group_by(period) |>
  mutate(pct_sites = .data$n_sites / sum(.data$n_sites)) |>
  ungroup() |>
  arrange(period, bioregion)
pca_scores <- analysis$pca_scores |>
  mutate(bioregion = relabel_clusters(.data$bioregion))
pca_centroids <- analysis$pca_centroids |>
  mutate(bioregion = relabel_clusters(.data$bioregion))
bioregion_palette <- c(
  "Cluster 1" = "#8ECAE6",
  "Cluster 2" = "#FDE68A",
  "Cluster 3" = "#D6A11F"
)

pca_variance <- analysis$reference_pca$sdev^2 / sum(analysis$reference_pca$sdev^2)
pc1_axis_label <- paste0("PCA axis 1 (", scales::percent(pca_variance[[1]], accuracy = 0.1), ")")
pc2_axis_label <- paste0("PCA axis 2 (", scales::percent(pca_variance[[2]], accuracy = 0.1), ")")

mainland_bbox <- c(xmin = 400000, xmax = 750000, ymin = 6000000, ymax = 6450000)
bornholm_bbox <- c(xmin = 855000, xmax = 905000, ymin = 6100000, ymax = 6160000)
mainland_width <- mainland_bbox[["xmax"]] - mainland_bbox[["xmin"]]
mainland_height <- mainland_bbox[["ymax"]] - mainland_bbox[["ymin"]]
bornholm_width <- bornholm_bbox[["xmax"]] - bornholm_bbox[["xmin"]]
bornholm_height <- bornholm_bbox[["ymax"]] - bornholm_bbox[["ymin"]]
bornholm_inset_width <- bornholm_width / mainland_width
bornholm_inset_height <- bornholm_height / mainland_height

bioregion_theme <- function(legend_position = "bottom") {
  theme_minimal(base_size = 10) +
    theme(
      legend.position = legend_position,
      legend.title = element_text(face = "bold"),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(0, 0, 0, 0)
    )
}

plot_bioregion_base <- function(df, fill_col, palette, title = NULL,
                                bbox = mainland_bbox, show_legend = FALSE,
                                border = FALSE) {
  ggplot(df) +
    geom_point(aes(x = X, y = Y, colour = .data[[fill_col]]), size = 1.25, alpha = 0.95) +
    scale_colour_manual(values = palette, drop = FALSE, name = "Cluster") +
    coord_fixed(
      xlim = c(bbox[["xmin"]], bbox[["xmax"]]),
      ylim = c(bbox[["ymin"]], bbox[["ymax"]]),
      expand = FALSE
    ) +
    bioregion_theme(if (show_legend) "bottom" else "none") +
    theme(
      panel.border = if (border) {
        element_rect(colour = "grey35", fill = NA, linewidth = 0.45)
      } else {
        element_blank()
      }
    ) +
    labs(title = title)
}

plot_bioregion_map <- function(df, fill_col, palette, title) {
  p_main <- plot_bioregion_base(df, fill_col, palette, title = title, bbox = mainland_bbox)
  p_inset <- plot_bioregion_base(df, fill_col, palette, bbox = bornholm_bbox, border = TRUE) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.7)
    )

  ggdraw(p_main) +
    draw_plot(
      p_inset,
      x = 1 - bornholm_inset_width - 0.2,
      y = 1 - bornholm_inset_height - 0.2,
      width = bornholm_inset_width,
      height = bornholm_inset_height
    )
}

bioregion_legend <- function(df, fill_col, palette) {
  plot_bioregion_base(df, fill_col, palette, show_legend = TRUE) |>
    get_legend()
}

plot_cluster_occupancy <- function(df, palette) {
  df |>
    arrange(bioregion) |>
    mutate(
      xmin = lag(cumsum(.data$pct_sites), default = 0),
      xmax = cumsum(.data$pct_sites),
      xmid = (.data$xmin + .data$xmax) / 2,
      label = paste0(
        str_replace(as.character(.data$bioregion), "Cluster ", "C"),
        " ",
        scales::percent(.data$pct_sites, accuracy = 1)
      )
    ) |>
    ggplot() +
    geom_rect(
      aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = 0.75,
        ymax = 1.25,
        fill = .data$bioregion
      ),
      colour = "white",
      linewidth = 0.35
    ) +
    geom_text(
      aes(x = .data$xmid, y = 1, label = .data$label),
      size = 2.6,
      colour = "grey15",
      fontface = "bold"
    ) +
    scale_fill_manual(values = palette, drop = FALSE, name = "Cluster") +
    scale_x_continuous(
      breaks = NULL,
      labels = NULL,
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(limits = c(0.7, 1.3), expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(0, 10, 4, 10)
    )
}

plot_bioregion_map_with_bar <- function(map_df, occupancy_df, palette, title) {
  map_plot <- plot_bioregion_map(map_df, "bioregion", palette, title)
  bar_plot <- plot_cluster_occupancy(occupancy_df, palette)

  plot_grid(map_plot, bar_plot, ncol = 1, rel_heights = c(1, 0.1), align = "v")
}

bioregion_period_maps <- imap(bioregion_assignments, function(df, atlas) {
  period <- year_lookup[[as.character(atlas)]]
  plot_bioregion_map_with_bar(
    df,
    cluster_occupancy |> filter(.data$period == .env$period),
    bioregion_palette,
    period
  )
})

bioregion_map_panel <- plot_grid(
  plotlist = bioregion_period_maps,
  nrow = 1,
  align = "h",
  rel_widths = c(1, 1, 1)
)

pca_centroids_plot <- pca_centroids |>
  mutate(period = factor(period, levels = unname(year_lookup))) |>
  arrange(bioregion, period)

bioregion_ordination_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, colour = bioregion)) +
  geom_point(alpha = 0.07, size = 0.75) +
  geom_path(
    data = pca_centroids_plot,
    aes(group = bioregion),
    linewidth = 0.8,
    arrow = arrow(length = unit(0.15, "cm")),
    alpha = 0.8
  ) +
  geom_point(data = pca_centroids_plot, aes(size = n_sites), alpha = 0.95) +
  geom_text(
    data = pca_centroids_plot,
    aes(label = period),
    size = 3,
    nudge_y = 0.015,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = bioregion_palette, drop = FALSE, name = "Cluster") +
  scale_size_continuous(name = "Sites", range = c(2, 5)) +
  labs(
    title = "Projected movement of modelled avian clusters",
    subtitle = "Atlas 2 and 3 communities projected into the reference Hellinger-PCA space",
    x = pc1_axis_label,
    y = pc2_axis_label
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

bioregion_publication_figure <- plot_grid(
  bioregion_map_panel,
  bioregion_ordination_plot +
    labs(title = NULL, subtitle = NULL) +
    theme(legend.position = "none"),
  labels = c("A", "B"),
  label_size = 16,
  label_fontface = "bold",
  ncol = 1,
  rel_heights = c(1.05, 1)
)

png_path <- file.path(out_dir, paste0(pattern, "-fig6-avian-bioregion-map-pca.png"))
pdf_path <- file.path(out_dir, paste0(pattern, "-fig6-avian-bioregion-map-pca.pdf"))

ggsave(png_path, bioregion_publication_figure, width = 13.5, height = 11, units = "in", dpi = 300, bg = "white")
ggsave(pdf_path, bioregion_publication_figure, width = 13.5, height = 11, units = "in", bg = "white")

message("Saved: ", png_path)
message("Saved: ", pdf_path)
