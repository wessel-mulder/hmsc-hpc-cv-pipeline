library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(patchwork)

set.seed(42)

output_dir <- "notebooks/exploratory/outputs/bioregion-cluster-vp-effects/concept-driver-clouds"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

environment_order <- c(
  "Temperature",
  "Precipitation",
  "Land-use heterogeneity",
  "Urban",
  "Cropland",
  "Pasture",
  "Forest",
  "Grass/Shrubland"
)

environment_colours <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban" = "grey65",
  "Cropland" = "goldenrod1",
  "Pasture" = "darkorange",
  "Forest" = "springgreen4",
  "Grass/Shrubland" = "springgreen2"
)

clusters <- paste("Cluster", 1:3)

fake_species <- tibble(
  cluster = rep(clusters, c(34, 29, 45)),
  species = sprintf("species_%03d", seq_len(108))
) |>
  group_by(cluster) |>
  mutate(species_rank = row_number()) |>
  ungroup()

cluster_variable_profiles <- expand_grid(cluster = clusters, variable = environment_order) |>
  mutate(
    support_rate = case_when(
      cluster == "Cluster 1" & variable %in% c("Precipitation", "Forest", "Grass/Shrubland") ~ 0.50,
      cluster == "Cluster 2" & variable %in% c("Cropland", "Pasture", "Land-use heterogeneity") ~ 0.52,
      cluster == "Cluster 3" & variable %in% c("Temperature", "Urban", "Cropland") ~ 0.48,
      TRUE ~ 0.20
    ),
    positive_bias = case_when(
      variable %in% c("Forest", "Grass/Shrubland", "Precipitation") ~ 0.68,
      variable %in% c("Urban", "Cropland") ~ 0.30,
      variable == "Temperature" & cluster == "Cluster 3" ~ 0.72,
      variable == "Pasture" & cluster == "Cluster 2" ~ 0.35,
      TRUE ~ 0.50
    ),
    vp_scale = case_when(
      variable %in% c("Temperature", "Cropland", "Forest") ~ 1.5,
      variable %in% c("Precipitation", "Pasture") ~ 1.2,
      TRUE ~ 0.9
    )
  )

fake_driver_data <- fake_species |>
  crossing(variable = environment_order) |>
  left_join(cluster_variable_profiles, by = c("cluster", "variable")) |>
  mutate(
    is_supported = runif(n()) < support_rate,
    direction = if_else(runif(n()) < positive_bias, "Positive", "Negative"),
    vp = pmin(0.38, rgamma(n(), shape = 1.6, scale = 0.035) * vp_scale),
    side = if_else(direction == "Positive", "Positive association", "Negative association"),
    x_side = if_else(direction == "Positive", -1, 1),
    variable = factor(variable, levels = environment_order),
    cluster = factor(cluster, levels = clusters)
  ) |>
  filter(is_supported)

grid_cols <- 10
grid_col_spacing <- 0.085
grid_row_spacing <- 0.075

grid_driver_data <- fake_driver_data |>
  arrange(cluster, variable, direction, desc(vp), species_rank) |>
  group_by(cluster, variable, direction) |>
  mutate(
    grid_rank = row_number() - 1,
    grid_col = grid_rank %% grid_cols,
    grid_row = grid_rank %/% grid_cols,
    direction_sign = if_else(direction == "Positive", -1, 1),
    grid_x = direction_sign * (0.18 + grid_col * grid_col_spacing),
    grid_y_offset = grid_row * grid_row_spacing
  ) |>
  ungroup()

plot_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9),
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text = element_text(face = "bold"),
    axis.title = element_text(size = 9)
  )

# Concept A: the most literal reading of the sketch. Each variable gets a
# small panel, clusters are horizontal bands, and supported species are packed
# into 10-column grids on either side of the positive/negative divider.
concept_a <- ggplot(
  grid_driver_data |>
    mutate(y_position = as.numeric(cluster) + grid_y_offset),
  aes(
    x = grid_x,
    y = y_position,
    colour = variable,
    size = vp,
    shape = direction
  )
) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey35") +
  geom_point(alpha = 0.82, stroke = 0.75) +
  facet_wrap(~variable, nrow = 2) +
  scale_x_continuous(
    breaks = c(-0.55, 0, 0.55),
    labels = c("Positive", "", "Negative"),
    limits = c(-1.05, 1.05)
  ) +
  scale_y_continuous(
    breaks = seq_along(clusters),
    labels = clusters,
    trans = "reverse"
  ) +
  scale_colour_manual(values = environment_colours, guide = "none") +
  scale_shape_manual(values = c("Positive" = 16, "Negative" = 4)) +
  scale_size_continuous(range = c(1.2, 5.2), name = "Variable importance") +
  labs(
    title = "Concept A: mirrored 10-column species grids per variable",
    subtitle = "Each symbol is one supported species; dots are positive, crosses are negative, and size is species-level VP.",
    x = NULL,
    y = NULL,
    shape = "Supported effect"
  ) +
  plot_theme +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(0.28, "lines")
  )

# Concept B: keeps each cluster as one row and stacks variables within the row.
# This gives a denser but more poster-like cloud, closer to the hand sketch.
variable_offsets <- tibble(
  variable = factor(environment_order, levels = environment_order),
  variable_offset = seq(-0.38, 0.38, length.out = length(environment_order))
)

concept_b_data <- grid_driver_data |>
  mutate(cluster_index = as.numeric(cluster)) |>
  left_join(variable_offsets, by = "variable") |>
  mutate(
    y_position = cluster_index + variable_offset + grid_y_offset * 0.44,
    grid_x = direction_sign * (0.16 + grid_col * 0.052)
  )

concept_b <- ggplot(
  concept_b_data,
  aes(
    x = grid_x,
    y = y_position,
    colour = variable,
    size = vp,
    shape = direction
  )
) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "grey35") +
  geom_hline(yintercept = c(1.5, 2.5), linewidth = 0.35, colour = "grey80") +
  geom_point(alpha = 0.84, stroke = 0.75) +
  scale_y_continuous(
    breaks = seq_along(clusters),
    labels = clusters,
    trans = "reverse"
  ) +
  scale_x_continuous(
    breaks = c(-0.44, 0, 0.44),
    labels = c("Positive", "", "Negative"),
    limits = c(-0.92, 0.92)
  ) +
  scale_colour_manual(values = environment_colours, name = "Variable") +
  scale_shape_manual(values = c("Positive" = 16, "Negative" = 4)) +
  scale_size_continuous(range = c(1.1, 5.3), name = "Variable importance") +
  labs(
    title = "Concept B: compact 10-column grids inside cluster rows",
    subtitle = "Variables are separated by small vertical offsets; within each variable, species are packed into dot grids.",
    x = NULL,
    y = NULL,
    shape = "Supported effect"
  ) +
  plot_theme +
  theme(
    panel.grid.major.y = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

# Concept C: a summarized version. Instead of showing individual species, each
# bubble is the cluster-variable-direction total VP across supported species.
concept_c_data <- fake_driver_data |>
  group_by(cluster, variable, direction) |>
  summarise(
    species_n = n_distinct(species),
    total_vp = sum(vp, na.rm = TRUE),
    mean_vp = mean(vp, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    x_position = as.numeric(variable) + if_else(direction == "Positive", -0.18, 0.18),
    cluster = fct_rev(cluster)
  )

concept_c <- ggplot(
  concept_c_data,
  aes(
    x = x_position,
    y = cluster,
    colour = variable,
    size = total_vp,
    shape = direction
  )
) +
  geom_vline(
    xintercept = seq_along(environment_order) + 0,
    linewidth = 0.22,
    colour = "grey85"
  ) +
  geom_point(alpha = 0.86, stroke = 0.95) +
  scale_x_continuous(
    breaks = seq_along(environment_order),
    labels = environment_order,
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_colour_manual(values = environment_colours, guide = "none") +
  scale_shape_manual(values = c("Positive" = 16, "Negative" = 4)) +
  scale_size_continuous(range = c(2.2, 12), name = "Summed VP") +
  labs(
    title = "Concept C: compressed direction-bubble summary",
    subtitle = "Each symbol aggregates supported species within a cluster, variable, and direction.",
    x = NULL,
    y = NULL,
    shape = "Supported effect"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.major.x = element_blank()
  )

ggsave(
  file.path(output_dir, "concept-a-mirrored-variable-clouds.png"),
  concept_a,
  width = 12,
  height = 7.5,
  dpi = 320
)

ggsave(
  file.path(output_dir, "concept-b-cluster-row-clouds.png"),
  concept_b,
  width = 10.5,
  height = 6.5,
  dpi = 320
)

ggsave(
  file.path(output_dir, "concept-c-aggregate-direction-bubbles.png"),
  concept_c,
  width = 10.5,
  height = 5.8,
  dpi = 320
)

concept_sheet <- (concept_a / concept_b / concept_c) +
  plot_annotation(
    title = "Fake-data concepts for species-level environmental dependency clouds",
    subtitle = "These are layout drafts only; no real species, VP, or Beta values are shown."
  )

ggsave(
  file.path(output_dir, "concept-sheet-driver-cloud-layouts.png"),
  concept_sheet,
  width = 13,
  height = 19,
  dpi = 320
)
