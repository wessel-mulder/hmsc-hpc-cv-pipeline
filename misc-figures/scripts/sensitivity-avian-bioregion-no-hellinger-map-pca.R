rm(list = ls())

# Quick sensitivity figure: avian bioregions without Hellinger transformation.
#
# This intentionally mirrors the manuscript bioregion map/PCA workflow, but uses
# raw modelled predicted occurrence probabilities as the community matrix. We
# fix k = 3 so the output is visually comparable with the Hellinger-based
# reference figure.

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, vegan, cluster, cowplot, grid, scales)

source(file.path("support_scripts", "figure_data_helpers.R"))

pattern <- "2026-03-13"
base_dir <- "HmscOutputs"
out_dir <- file.path("misc-figures", "outputs", "main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

year_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")
reference_atlas <- "1"
cluster_levels <- paste("Cluster", 1:3)

cluster_palette <- c(
  "Cluster 1" = "#8ECAE6",
  "Cluster 2" = "#FDE68A",
  "Cluster 3" = "#D6A11F"
)

model_folders <- figure_model_folders(pattern = pattern, base_dir = base_dir)
model_nums <- as.character(atlas_numbers(model_folders))
names(model_folders) <- model_nums

prediction_paths <- file.path(
  base_dir,
  model_folders,
  "Results", "Preds", "predicted-values.rdata"
)
if (any(!file.exists(prediction_paths))) {
  stop(
    "Missing prediction files:\n",
    paste(prediction_paths[!file.exists(prediction_paths)], collapse = "\n"),
    call. = FALSE
  )
}

preds <- map(prediction_paths, readRDS)
names(preds) <- names(model_folders)

site_frame_list <- read_csv(
  file.path(
    "notebooks", "exploratory", "outputs", "bioregion-pca",
    paste0(pattern, "-bioregion-assignments.csv")
  ),
  show_col_types = FALSE
) |>
  mutate(atlas = as.character(.data$atlas)) |>
  select(.data$atlas, .data$period, .data$survey, .data$site, .data$X, .data$Y) |>
  group_split(.data$atlas)
site_frames <- set_names(site_frame_list, map_chr(site_frame_list, ~unique(.x$atlas)))

align_site_frame <- function(pred, site_frame) {
  missing_from_sites <- setdiff(rownames(pred), site_frame$survey)
  if (length(missing_from_sites) > 0) {
    stop("Prediction rows missing from site frame: ", length(missing_from_sites), call. = FALSE)
  }
  site_frame |>
    slice(match(rownames(pred), .data$survey))
}

raw_communities <- preds
reference_raw <- raw_communities[[reference_atlas]]

reference_dist <- vegan::vegdist(reference_raw, method = "bray")
reference_fit <- cluster::pam(reference_dist, k = 3, diss = TRUE)

# Reorder clusters by Atlas 1 size so the two smaller clusters are C1/C2 and
# the largest is C3, matching the current manuscript convention.
reference_cluster_sizes <- sort(table(reference_fit$clustering), decreasing = FALSE)
cluster_order <- as.integer(names(reference_cluster_sizes))
cluster_recode <- setNames(seq_along(cluster_order), cluster_order)
reference_clusters <- unname(cluster_recode[as.character(reference_fit$clustering)])

medoid_indices <- reference_fit$id.med[cluster_order]
medoid_matrix <- reference_raw[medoid_indices, , drop = FALSE]
rownames(medoid_matrix) <- cluster_levels

bray_to_medoids <- function(query_matrix, medoid_matrix) {
  map_dfc(seq_len(nrow(medoid_matrix)), function(i) {
    medoid <- medoid_matrix[i, ]
    numerator <- rowSums(abs(sweep(query_matrix, 2, medoid, "-")), na.rm = TRUE)
    denominator <- rowSums(sweep(query_matrix, 2, medoid, "+"), na.rm = TRUE)
    tibble(!!cluster_levels[[i]] := ifelse(denominator == 0, 0, numerator / denominator))
  })
}

assign_to_medoids <- function(query_matrix, site_frame, medoid_matrix) {
  distances <- bray_to_medoids(query_matrix, medoid_matrix)
  distance_matrix <- as.matrix(distances)
  nearest_idx <- max.col(-distance_matrix, ties.method = "first")
  sorted_distances <- t(apply(distance_matrix, 1, sort))

  site_frame |>
    mutate(
      bioregion = factor(cluster_levels[nearest_idx], levels = cluster_levels),
      nearest_distance = distance_matrix[cbind(seq_len(nrow(distance_matrix)), nearest_idx)],
      second_distance = sorted_distances[, 2],
      assignment_margin = .data$second_distance - .data$nearest_distance
    )
}

bioregion_assignments <- imap(raw_communities, function(raw_matrix, atlas) {
  assign_to_medoids(
    raw_matrix,
    align_site_frame(raw_matrix, site_frames[[atlas]]),
    medoid_matrix
  )
})

# Confirm Atlas 1 nearest-medoid assignment reproduces the reordered PAM labels.
atlas1_validation <- tibble(
  survey = rownames(reference_raw),
  reference_cluster = factor(cluster_levels[reference_clusters], levels = cluster_levels)
) |>
  left_join(
    bioregion_assignments[[reference_atlas]] |> select(.data$survey, assigned_cluster = .data$bioregion),
    by = "survey"
  ) |>
  summarise(match_rate = mean(.data$reference_cluster == .data$assigned_cluster, na.rm = TRUE))

cluster_occupancy <- imap_dfr(bioregion_assignments, function(df, atlas) {
  df |>
    count(.data$bioregion, name = "n_sites") |>
    complete(
      bioregion = factor(cluster_levels, levels = cluster_levels),
      fill = list(n_sites = 0)
    ) |>
    mutate(
      atlas = .env$atlas,
      period = year_lookup[[as.character(.env$atlas)]],
      pct_sites = .data$n_sites / sum(.data$n_sites)
    )
})

reference_pca <- prcomp(reference_raw, center = TRUE, scale. = FALSE)
pca_variance <- reference_pca$sdev^2 / sum(reference_pca$sdev^2)
pca_scores <- imap_dfr(raw_communities, function(raw_matrix, atlas) {
  predict(reference_pca, newdata = raw_matrix)[, 1:2, drop = FALSE] |>
    as_tibble() |>
    bind_cols(
      bioregion_assignments[[atlas]] |>
        select(.data$survey, .data$site, .data$X, .data$Y, .data$bioregion)
    ) |>
    mutate(
      atlas = .env$atlas,
      period = year_lookup[[as.character(.env$atlas)]]
    )
})

pca_centroids <- pca_scores |>
  group_by(.data$period, .data$bioregion) |>
  summarise(
    PC1 = mean(.data$PC1, na.rm = TRUE),
    PC2 = mean(.data$PC2, na.rm = TRUE),
    n_sites = n(),
    .groups = "drop"
  ) |>
  mutate(period = factor(.data$period, levels = unname(year_lookup))) |>
  arrange(.data$bioregion, .data$period)

mainland_bbox <- c(xmin = 400000, xmax = 750000, ymin = 6000000, ymax = 6450000)
bornholm_bbox <- c(xmin = 855000, xmax = 905000, ymin = 6100000, ymax = 6160000)
mainland_width <- mainland_bbox[["xmax"]] - mainland_bbox[["xmin"]]
mainland_height <- mainland_bbox[["ymax"]] - mainland_bbox[["ymin"]]
bornholm_width <- bornholm_bbox[["xmax"]] - bornholm_bbox[["xmin"]]
bornholm_height <- bornholm_bbox[["ymax"]] - bornholm_bbox[["ymin"]]
bornholm_inset_width <- bornholm_width / mainland_width
bornholm_inset_height <- bornholm_height / mainland_height

map_theme <- function(legend_position = "none") {
  theme_minimal(base_size = 10) +
    theme(
      legend.position = legend_position,
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(0, 0, 0, 0)
    )
}

plot_map_base <- function(df, bbox, title = NULL, border = FALSE) {
  ggplot(df) +
    geom_point(aes(x = .data$X, y = .data$Y, colour = .data$bioregion), size = 1.25, alpha = 0.95) +
    scale_colour_manual(values = cluster_palette, drop = FALSE, name = "Cluster") +
    coord_fixed(
      xlim = c(bbox[["xmin"]], bbox[["xmax"]]),
      ylim = c(bbox[["ymin"]], bbox[["ymax"]]),
      expand = FALSE
    ) +
    map_theme() +
    labs(title = title) +
    theme(
      panel.border = if (border) {
        element_rect(colour = "grey35", fill = NA, linewidth = 0.55)
      } else {
        element_blank()
      }
    )
}

plot_map <- function(df, title) {
  p_main <- plot_map_base(df, mainland_bbox, title = title)
  p_inset <- plot_map_base(df, bornholm_bbox, border = TRUE) +
    theme_void() +
    theme(
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

plot_cluster_occupancy <- function(df) {
  df |>
    arrange(.data$bioregion) |>
    mutate(
      xmin = lag(cumsum(.data$pct_sites), default = 0),
      xmax = cumsum(.data$pct_sites),
      xmid = (.data$xmin + .data$xmax) / 2,
      label = paste0(str_replace(as.character(.data$bioregion), "Cluster ", "C"), " ", percent(.data$pct_sites, accuracy = 1))
    ) |>
    ggplot() +
    geom_rect(
      aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0.75, ymax = 1.25, fill = .data$bioregion),
      colour = "white",
      linewidth = 0.35
    ) +
    geom_text(aes(x = .data$xmid, y = 1, label = .data$label), size = 2.6, fontface = "bold") +
    scale_fill_manual(values = cluster_palette, drop = FALSE) +
    scale_x_continuous(limits = c(0, 1), breaks = NULL, expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0.7, 1.3), expand = expansion(mult = c(0, 0))) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = 9) +
    theme(plot.margin = margin(0, 10, 4, 10))
}

plot_map_with_bar <- function(df, occupancy_df, title) {
  plot_grid(
    plot_map(df, title),
    plot_cluster_occupancy(occupancy_df),
    ncol = 1,
    rel_heights = c(1, 0.1),
    align = "v"
  )
}

map_panel_list <- imap(bioregion_assignments, function(df, atlas) {
  period <- year_lookup[[as.character(atlas)]]
  plot_map_with_bar(
    df,
    cluster_occupancy |> filter(.data$period == .env$period),
    period
  )
})
map_panel <- plot_grid(plotlist = map_panel_list, nrow = 1, align = "h")

pca_plot <- ggplot(pca_scores, aes(x = .data$PC1, y = .data$PC2, colour = .data$bioregion)) +
  geom_point(alpha = 0.07, size = 0.75) +
  geom_path(
    data = pca_centroids,
    aes(group = .data$bioregion),
    linewidth = 0.8,
    arrow = arrow(length = unit(0.15, "cm")),
    alpha = 0.85
  ) +
  geom_point(data = pca_centroids, aes(size = .data$n_sites), alpha = 0.95) +
  geom_text(
    data = pca_centroids,
    aes(label = .data$period),
    size = 3,
    nudge_y = 0.015,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = cluster_palette, drop = FALSE, name = "Cluster") +
  scale_size_continuous(name = "Sites", range = c(2, 5)) +
  labs(
    x = paste0("PCA axis 1 (", percent(pca_variance[[1]], accuracy = 0.1), ")"),
    y = paste0("PCA axis 2 (", percent(pca_variance[[2]], accuracy = 0.1), ")")
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

figure <- plot_grid(
  map_panel,
  pca_plot + theme(legend.position = "none"),
  labels = c("A", "B"),
  label_size = 16,
  label_fontface = "bold",
  ncol = 1,
  rel_heights = c(1.05, 1)
)

prefix <- paste0(pattern, "-sensitivity-no-hellinger-bioregion-map-pca")
png_path <- file.path(out_dir, paste0(prefix, ".png"))
pdf_path <- file.path(out_dir, paste0(prefix, ".pdf"))
assignment_path <- file.path(out_dir, paste0(prefix, "-assignments.csv"))
diagnostic_path <- file.path(out_dir, paste0(prefix, "-diagnostics.csv"))

write_csv(bind_rows(bioregion_assignments), assignment_path)
write_csv(
  tibble(
    check = c("fixed_k", "atlas1_validation_match_rate"),
    value = c(3, atlas1_validation$match_rate)
  ),
  diagnostic_path
)

ggsave(png_path, figure, width = 13.5, height = 11, units = "in", dpi = 300, bg = "white")
ggsave(pdf_path, figure, width = 13.5, height = 11, units = "in", bg = "white")

message("Saved: ", png_path)
message("Saved: ", pdf_path)
