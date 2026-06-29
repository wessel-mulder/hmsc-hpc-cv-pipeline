rm(list = ls())

# Figure 8b concepts: environmental-dependency similarity among avian
# bioregion clusters within and across atlas periods.
#
# This script starts from the temporal dependency summary created by
# `fig8-avian-bioregion-environmental-dependency-waffle-pca.R`. Each
# period x cluster profile is a 16-feature vector:
#
#   8 environmental variables x 2 supported Beta directions
#
# Values are weighted supported dependencies, then normalized so each
# period x cluster sums to one across all positive and negative features.
# The plots below are deliberately exploratory alternatives to the waffle chart:
#
#   A. Pairwise profile similarity heatmap among all 6 period-cluster profiles.
#   B. PCoA ordination of those 6 profiles, with arrows from 1970s to 2010s.
#   C. Cross-period similarity matrix: 1970s clusters vs 2010s clusters.
#   D. Similarity lollipop comparing within-period and across-period contrasts.

# ---- Packages ----
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, patchwork, scales)

# ---- Paths ----
pattern <- "2026-03-13"
out_dir <- file.path("misc-figures", "outputs", "main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dependency_summary_path <- file.path(
  out_dir,
  paste0(pattern, "-fig8-avian-bioregion-environmental-dependency-1970s-2010s-waffle-summary.csv")
)

if (!file.exists(dependency_summary_path)) {
  stop(
    "Missing temporal dependency summary: ", dependency_summary_path,
    "\nRun misc-figures/scripts/fig8-avian-bioregion-environmental-dependency-waffle-pca.R first.",
    call. = FALSE
  )
}

# ---- Shared labels and palettes ----
cluster_levels <- paste("Cluster", 1:3)
period_levels <- c("1970s", "2010s")

cluster_palette <- c(
  "Cluster 1" = "#76B7D8",
  "Cluster 2" = "#F6E58D",
  "Cluster 3" = "#C49A00"
)

environment_order <- c(
  "Temperature",
  "Precipitation",
  "Land-use heterogeneity",
  "Urban (% coverage)",
  "Cropland (% coverage)",
  "Pasture (% coverage)",
  "Forest (% coverage)",
  "Grass/Shrubland (% coverage)"
)

feature_order <- as.vector(outer(
  c("Positive", "Negative"),
  environment_order,
  paste,
  sep = " | "
))

profile_levels <- tidyr::expand_grid(
  period = period_levels,
  bioregion = cluster_levels
) |>
  transmute(profile_id = paste(.data$period, .data$bioregion, sep = " | ")) |>
  pull(.data$profile_id)

# ---- Build dependency-profile matrix ----
dependency_summary <- read_csv(dependency_summary_path, show_col_types = FALSE) |>
  mutate(
    period = case_when(
      as.character(.data$period) %in% c("1970", "1970s") ~ "1970s",
      as.character(.data$period) %in% c("2010", "2010s") ~ "2010s",
      TRUE ~ as.character(.data$period)
    ),
    period = factor(.data$period, levels = period_levels),
    bioregion = factor(.data$bioregion, levels = cluster_levels),
    effect_direction = factor(.data$effect_direction, levels = c("Positive", "Negative")),
    variable_label = factor(.data$variable_label, levels = environment_order),
    feature = paste(.data$effect_direction, .data$variable_label, sep = " | ")
  )

profile_long <- dependency_summary |>
  group_by(.data$period, .data$bioregion) |>
  mutate(
    total_dependency = sum(.data$weighted_dependency, na.rm = TRUE),
    profile_share = if_else(
      .data$total_dependency > 0,
      .data$weighted_dependency / .data$total_dependency,
      0
    )
  ) |>
  ungroup() |>
  complete(
    period = factor(period_levels, levels = period_levels),
    bioregion = factor(cluster_levels, levels = cluster_levels),
    feature = feature_order,
    fill = list(profile_share = 0)
  ) |>
  mutate(
    profile_id = paste(.data$period, .data$bioregion, sep = " | "),
    profile_id = factor(
      .data$profile_id,
      levels = profile_levels
    ),
    feature = factor(.data$feature, levels = feature_order)
  )

profile_matrix_df <- profile_long |>
  select("profile_id", "period", "bioregion", "feature", "profile_share") |>
  distinct() |>
  pivot_wider(
    id_cols = c("profile_id", "period", "bioregion"),
    names_from = "feature",
    values_from = "profile_share",
    values_fill = 0
  ) |>
  arrange(.data$profile_id)

profile_matrix <- profile_matrix_df |>
  select(all_of(feature_order)) |>
  as.matrix()
rownames(profile_matrix) <- as.character(profile_matrix_df$profile_id)

# ---- Similarity helpers ----
bray_dissimilarity <- function(x, y) {
  denominator <- sum(x + y, na.rm = TRUE)
  if (!is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }
  sum(abs(x - y), na.rm = TRUE) / denominator
}

profile_ids <- rownames(profile_matrix)
similarity_matrix <- matrix(
  NA_real_,
  nrow = length(profile_ids),
  ncol = length(profile_ids),
  dimnames = list(profile_ids, profile_ids)
)

for (i in seq_along(profile_ids)) {
  for (j in seq_along(profile_ids)) {
    similarity_matrix[i, j] <- 1 - bray_dissimilarity(profile_matrix[i, ], profile_matrix[j, ])
  }
}

dissimilarity_matrix <- 1 - similarity_matrix

similarity_long <- as.data.frame(as.table(similarity_matrix)) |>
  as_tibble() |>
  rename(profile_a = Var1, profile_b = Var2, similarity = Freq) |>
  separate(.data$profile_a, into = c("period_a", "cluster_a"), sep = " \\| ", remove = FALSE) |>
  separate(.data$profile_b, into = c("period_b", "cluster_b"), sep = " \\| ", remove = FALSE) |>
  mutate(
    profile_a = factor(.data$profile_a, levels = profile_ids),
    profile_b = factor(.data$profile_b, levels = rev(profile_ids)),
    period_a = factor(.data$period_a, levels = period_levels),
    period_b = factor(.data$period_b, levels = period_levels),
    cluster_a = factor(.data$cluster_a, levels = cluster_levels),
    cluster_b = factor(.data$cluster_b, levels = cluster_levels),
    is_same_profile = .data$profile_a == as.character(.data$profile_b)
  )

# ---- Panel A: all-pair similarity heatmap ----
all_pair_heatmap <- similarity_long |>
  ggplot(aes(x = .data$profile_a, y = .data$profile_b, fill = .data$similarity)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = percent(.data$similarity, accuracy = 1)), size = 2.8) +
  scale_fill_viridis_c(
    option = "mako",
    direction = -1,
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    name = "Bray-Curtis\nsimilarity"
  ) +
  labs(x = NULL, y = NULL) +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right"
  )

# ---- Panel B: profile ordination through time ----
pcoa <- cmdscale(
  as.dist(dissimilarity_matrix),
  eig = TRUE,
  k = 2
)

positive_eigen_sum <- sum(pcoa$eig[pcoa$eig > 0], na.rm = TRUE)
pcoa_variance <- pcoa$eig[seq_len(2)] / positive_eigen_sum

pcoa_scores <- as.data.frame(pcoa$points) |>
  set_names(c("Axis1", "Axis2")) |>
  rownames_to_column("profile_id") |>
  separate(.data$profile_id, into = c("period", "bioregion"), sep = " \\| ", remove = FALSE) |>
  mutate(
    period = factor(.data$period, levels = period_levels),
    bioregion = factor(.data$bioregion, levels = cluster_levels)
  )

pcoa_arrow_df <- pcoa_scores |>
  select("period", "bioregion", "Axis1", "Axis2") |>
  pivot_wider(
    names_from = "period",
    values_from = c("Axis1", "Axis2")
  ) |>
  mutate(label_x = (.data$Axis1_1970s + .data$Axis1_2010s) / 2)

profile_ordination <- ggplot(pcoa_scores, aes(x = .data$Axis1, y = .data$Axis2)) +
  geom_hline(yintercept = 0, colour = "grey86", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey86", linewidth = 0.3) +
  geom_segment(
    data = pcoa_arrow_df,
    aes(
      x = .data$Axis1_1970s,
      y = .data$Axis2_1970s,
      xend = .data$Axis1_2010s,
      yend = .data$Axis2_2010s,
      colour = .data$bioregion
    ),
    inherit.aes = FALSE,
    linewidth = 0.75,
    arrow = arrow(length = unit(0.09, "in"))
  ) +
  geom_point(
    aes(colour = .data$bioregion, shape = .data$period),
    size = 4,
    stroke = 0.8
  ) +
  geom_text(
    aes(
      label = paste0(str_replace(.data$bioregion, "Cluster ", "C"), " ", .data$period)
    ),
    size = 3,
    nudge_y = 0.032,
    check_overlap = TRUE
  ) +
  scale_colour_manual(values = cluster_palette, name = "Cluster") +
  scale_shape_manual(values = c("1970s" = 16, "2010s" = 17), name = "Period") +
  labs(
    x = paste0("PCoA 1 (", percent(pcoa_variance[[1]], accuracy = 0.1), ")"),
    y = paste0("PCoA 2 (", percent(pcoa_variance[[2]], accuracy = 0.1), ")")
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(5.5, 25, 5.5, 18)
  )

# ---- Panel C: cross-period similarity matrix ----
cross_period_similarity <- similarity_long |>
  filter(.data$period_a == "1970s", .data$period_b == "2010s") |>
  mutate(
    diagonal_match = .data$cluster_a == .data$cluster_b,
    cluster_a = factor(.data$cluster_a, levels = cluster_levels),
    cluster_b = factor(.data$cluster_b, levels = rev(cluster_levels))
  )

cross_period_heatmap <- cross_period_similarity |>
  ggplot(aes(x = .data$cluster_a, y = .data$cluster_b, fill = .data$similarity)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_tile(
    data = filter(cross_period_similarity, .data$diagonal_match),
    fill = NA,
    colour = "grey10",
    linewidth = 0.85
  ) +
  geom_text(aes(label = percent(.data$similarity, accuracy = 1)), size = 3.6) +
  scale_fill_viridis_c(
    option = "mako",
    direction = -1,
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    guide = "none"
  ) +
  labs(
    x = "1970s cluster",
    y = "2010s cluster"
  ) +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank())

# ---- Panel D: contrast lollipop ----
same_cluster_temporal <- similarity_long |>
  filter(
    .data$period_a == "1970s",
    .data$period_b == "2010s",
    .data$cluster_a == .data$cluster_b
  ) |>
  transmute(
    contrast_group = "Same cluster across periods",
    contrast = paste(.data$cluster_a, "1970s-2010s"),
    cluster = .data$cluster_a,
    similarity
  )

within_period_pairs <- similarity_long |>
  filter(
    .data$period_a == .data$period_b,
    as.integer(.data$cluster_a) < as.integer(.data$cluster_b)
  ) |>
  transmute(
    contrast_group = paste(.data$period_a, "within period"),
    contrast = paste(.data$cluster_a, "vs", .data$cluster_b),
    cluster = NA,
    similarity
  )

nearest_cross_period <- cross_period_similarity |>
  group_by(.data$cluster_a) |>
  slice_max(.data$similarity, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(
    contrast_group = "Nearest 2010s profile",
    contrast = paste(.data$cluster_a, "1970s ->", .data$cluster_b, "2010s"),
    cluster = .data$cluster_a,
    similarity
  )

contrast_summary <- bind_rows(
  same_cluster_temporal,
  within_period_pairs,
  nearest_cross_period
) |>
  mutate(
    contrast_group = factor(
      .data$contrast_group,
      levels = c(
        "1970s within period",
        "2010s within period",
        "Same cluster across periods",
        "Nearest 2010s profile"
      )
    ),
    contrast = fct_reorder(.data$contrast, .data$similarity)
  )

contrast_lollipop <- contrast_summary |>
  ggplot(aes(x = .data$similarity, y = .data$contrast)) +
  geom_segment(
    aes(x = 0, xend = .data$similarity, yend = .data$contrast),
    colour = "grey78",
    linewidth = 0.5
  ) +
  geom_point(
    aes(colour = .data$cluster),
    size = 3.4
  ) +
  facet_grid(.data$contrast_group ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = cluster_palette, na.value = "grey35", guide = "none") +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0.01, 0.08))
  ) +
  labs(
    x = "Bray-Curtis similarity",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold", hjust = 0)
  )

# ---- Composite and exports ----
similarity_concepts <- (
  all_pair_heatmap | profile_ordination
) / (
  cross_period_heatmap | contrast_lollipop
) +
  plot_annotation(
    tag_levels = "A",
    title = "Environmental-dependency profile similarity among avian bioregion clusters",
    subtitle = "Profiles are proportional weighted supported VP x Beta-direction dependencies across environmental variables."
  ) &
  theme(
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold", size = 13)
  )

summary_csv_path <- file.path(
  out_dir,
  paste0(pattern, "-fig8b-avian-bioregion-environmental-dependency-similarity-summary.csv")
)
profile_csv_path <- file.path(
  out_dir,
  paste0(pattern, "-fig8b-avian-bioregion-environmental-dependency-profile-vectors.csv")
)
png_path <- file.path(
  out_dir,
  paste0(pattern, "-fig8b-avian-bioregion-environmental-dependency-similarity-concepts.png")
)
pdf_path <- file.path(
  out_dir,
  paste0(pattern, "-fig8b-avian-bioregion-environmental-dependency-similarity-concepts.pdf")
)

write_csv(similarity_long, summary_csv_path)
write_csv(profile_long, profile_csv_path)

ggsave(
  png_path,
  similarity_concepts,
  width = 14,
  height = 10.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  pdf_path,
  similarity_concepts,
  width = 14,
  height = 10.5,
  units = "in",
  bg = "white"
)

message("Saved: ", png_path)
message("Saved: ", pdf_path)
message("Saved: ", summary_csv_path)
message("Saved: ", profile_csv_path)
