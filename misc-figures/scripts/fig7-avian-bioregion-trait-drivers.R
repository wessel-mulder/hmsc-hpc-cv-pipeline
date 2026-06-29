rm(list = ls())

# Figure 7: trait drivers of avian bioregion shifts
#
# This script rebuilds the manuscript-ready trait-driver figure from saved
# exploratory outputs. It does not rerun the HMSC models, bioregion clustering,
# or trait-delta calculations. The figure focuses on cells that moved into the
# relabelled Cluster 3 by the 2010s.

# ---- Packages ----
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, patchwork, scales, grid, cowplot)

# ---- Paths and figure controls ----
pattern <- "2026-03-13"
in_dir <- file.path("notebooks", "exploratory", "outputs", "bioregion-trait-composition")
trait_delta_path <- file.path(in_dir, paste0(pattern, "-cluster2-shift-trait-deltas-size-normalized.csv"))
trait_metadata_path <- file.path("Data", "preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData")
bioregion_analysis_path <- file.path(
  "notebooks", "exploratory", "outputs", "bioregion-pca", "models",
  paste0(pattern, "-avian-bioregion-pca-analysis.rds")
)
out_dir <- file.path("misc-figures", "outputs", "main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cluster2_shift_focal_comparison <- "2010s minus 1970s"

# ---- Trait orders and palettes ----
# Thermal colours match the updated Fig. 4 CWM-STI map palette.
thermal_levels <- c("very cold", "cold", "medium", "warm", "very warm")
migration_levels <- c(
  "long-distance",
  "short-and long-distance",
  "short-distance",
  "sedentary and short-distance",
  "sedentary"
)

thermal_colours <- c("#313695", "#74ADD1", "#FEE090", "#F46D43", "#9E0142")
thermal_palette <- setNames(
  thermal_colours,
  thermal_levels
)
migration_palette <- setNames(
  colorRampPalette(c("#542788", "#f2e5ff"))(length(migration_levels)),
  migration_levels
)

# ---- Load saved analysis outputs ----
cluster2_shift_trait_deltas <- read_csv(trait_delta_path, show_col_types = FALSE)
bioregion_analysis <- readRDS(bioregion_analysis_path)

# ---- Display labels for relabelled bioregion transitions ----
# The exploratory analysis identified shifts into the old Cluster 2. In the
# manuscript figure, clusters were reordered so those target cells are Cluster 3.
transition_display_labels <- c(
  "Cluster 1 to Cluster 2" = "Cluster 1 to Cluster 3",
  "Cluster 3 to Cluster 2" = "Cluster 2 to Cluster 3"
)
transition_display_levels <- unname(transition_display_labels)
transition_axis_labels <- c(
  "Cluster 1 to Cluster 3" = "Cluster 1\nto Cluster 3",
  "Cluster 2 to Cluster 3" = "Cluster 2\nto Cluster 3"
)

load(trait_metadata_path)

# ---- Cluster relabelling and colours ----
# Reorder the saved exploratory clusters so the two smaller regions become
# Clusters 1 and 2, and the large target region becomes Cluster 3.
cluster_levels <- paste("Cluster", 1:3)
cluster_recode <- c(
  "Cluster 1" = "Cluster 1",
  "Cluster 3" = "Cluster 2",
  "Cluster 2" = "Cluster 3"
)
relabel_clusters <- function(x) {
  factor(unname(cluster_recode[as.character(x)]), levels = cluster_levels)
}
bioregion_palette <- c(
  "Cluster 1" = "#8ECAE6",
  "Cluster 2" = "#FDE68A",
  "Cluster 3" = "#D6A11F"
)

# ---- Map extents ----
# Bornholm is plotted as an inset but uses the same point size as the mainland.
mainland_bbox <- c(xmin = 400000, xmax = 750000, ymin = 6000000, ymax = 6450000)
bornholm_bbox <- c(xmin = 855000, xmax = 905000, ymin = 6100000, ymax = 6160000)
mainland_width <- mainland_bbox[["xmax"]] - mainland_bbox[["xmin"]]
mainland_height <- mainland_bbox[["ymax"]] - mainland_bbox[["ymin"]]
bornholm_width <- bornholm_bbox[["xmax"]] - bornholm_bbox[["xmin"]]
bornholm_height <- bornholm_bbox[["ymax"]] - bornholm_bbox[["ymin"]]
bornholm_inset_width <- bornholm_width / mainland_width
bornholm_inset_height <- bornholm_height / mainland_height

# ---- Bioregion map helpers ----
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

# Manual label positions for the 1970s reference map.
atlas1_label_positions <- tibble(
  map_class = factor(cluster_levels, levels = cluster_levels),
  X = c(470000, 460000, 585000),
  Y = c(6310000, 6205000, 6138000),
  label = cluster_levels
)

plot_bioregion_base <- function(df, colour_col, palette, bbox = mainland_bbox,
                                show_legend = FALSE, border = FALSE,
                                label_df = NULL, legend_title = "Cluster",
                                legend_rows = 1) {
  ggplot(df) +
    geom_point(aes(x = X, y = Y, colour = .data[[colour_col]]), size = 1.25, alpha = 0.95) +
    {
      if (!is.null(label_df)) {
        geom_label(
          data = label_df,
          aes(x = X, y = Y, label = label),
          inherit.aes = FALSE,
          linewidth = 0,
          alpha = 0.88,
          fill = "white",
          colour = "grey15",
          fontface = "bold",
          size = 3.2,
          show.legend = FALSE
        )
      }
    } +
    scale_colour_manual(
      values = palette,
      drop = FALSE,
      name = legend_title,
      guide = guide_legend(nrow = legend_rows, byrow = TRUE)
    ) +
    coord_fixed(
      xlim = c(bbox[["xmin"]], bbox[["xmax"]]),
      ylim = c(bbox[["ymin"]], bbox[["ymax"]]),
      expand = FALSE
    ) +
    bioregion_theme(if (show_legend) "bottom" else "none") +
    theme(
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_blank(),
      panel.border = if (border) {
        element_rect(colour = "grey35", fill = NA, linewidth = 0.45)
      } else {
        element_blank()
      }
    )
}

# The inset is placed in the upper-right of the map panel using relative panel
# coordinates. Its relative size is based on the Bornholm/mainland map extents.
plot_bioregion_map <- function(df, colour_col, palette, title, label_df = NULL,
                               legend_title = "Cluster", legend_rows = 1) {
  p_main <- plot_bioregion_base(
    df,
    colour_col = colour_col,
    palette = palette,
    bbox = mainland_bbox,
    show_legend = TRUE,
    label_df = label_df,
    legend_title = legend_title,
    legend_rows = legend_rows
  ) +
    labs(title = title)
  p_inset <- plot_bioregion_base(
    df,
    colour_col = colour_col,
    palette = palette,
    bbox = bornholm_bbox,
    border = TRUE,
    legend_title = legend_title,
    legend_rows = legend_rows
  ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.7)
    )

  p_main +
    inset_element(
      p_inset,
      left = 0.73,
      bottom = 0.67,
      right = 0.73 + bornholm_inset_width,
      top = 0.67 + bornholm_inset_height,
      align_to = "panel",
      on_top = TRUE
    )
}

# ---- Panel D: 1970s reference bioregion map ----
atlas1_bioregion_df <- bioregion_analysis$bioregion_assignments[["1"]] |>
  mutate(
    bioregion = relabel_clusters(.data$bioregion),
    map_class = .data$bioregion
  )

atlas1_bioregion_map <- plot_bioregion_map(
  atlas1_bioregion_df,
  colour_col = "map_class",
  palette = bioregion_palette,
  title = "1970s cluster distribution",
  label_df = atlas1_label_positions
)

# ---- Species trait metadata ----
# Thermal groups are quintiles of species thermal index, matching the STI
# grouping used elsewhere in this chapter.
thermal_groups <- function(sti) {
  cut(
    sti,
    breaks = quantile(sti, seq(0, 1, 0.2), na.rm = TRUE),
    include.lowest = TRUE,
    labels = thermal_levels
  )
}

species_trait_table <- Tr |>
  as.data.frame() |>
  rownames_to_column("species") |>
  select(
    species,
    migratory_strategy = Migration_a3_DOF,
    foraging_guild = foraging_guild_consensus,
    species_thermal_index
  ) |>
  mutate(
    thermal_group = thermal_groups(.data$species_thermal_index),
    thermal_group = factor(.data$thermal_group, levels = thermal_levels)
  )

# ---- Standalone thermal and migratory bar plot data ----
# This plot is kept as a separate export, even though the main composite now uses
# compact heatmaps for the same thermal/migratory summaries.
thermal_migration_driver_df <- cluster2_shift_trait_deltas |>
  filter(
    comparison == cluster2_shift_focal_comparison,
    trait_family %in% c("Thermal affinity", "Migratory strategy")
  ) |>
  mutate(
    trait_value = factor(trait_value, levels = c(thermal_levels, rev(migration_levels))),
    trait_colour_group = as.character(trait_value),
    transition = factor(transition_display_labels[.data$transition], levels = transition_display_levels)
  )

cluster2_shift_thermal_migration_plot <- ggplot(
  thermal_migration_driver_df,
  aes(
    x = mean_expected_richness_delta_per_species,
    y = trait_value,
    fill = trait_colour_group
  )
) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.35) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  facet_grid(trait_family ~ transition, scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = c(thermal_palette, migration_palette),
    breaks = c(rev(thermal_levels), migration_levels),
    name = NULL
  ) +
  labs(
    title = "Size-normalized thermal and migratory drivers of shifts into Cluster 3 by the 2010s",
    subtitle = "Bars show mean expected-richness delta per species in each trait group.",
    x = "Mean expected richness delta per species",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

# ---- Foraging-guild driver data ----
aquatic_wetland_guilds <- c(
  "Aquatic pursuers", "Dabbling ducks", "Dippers", "Diving ducks",
  "Grazing waterfowl", "Gulls", "Kingfishers", "Marsh warblers",
  "Plovers", "Rails", "Reedlings", "Scolopacids", "Terns", "Wading birds"
  )

# Order guild rows by their average size-normalized delta across the two focal
# transitions, so declining guilds sit at the bottom and increasing guilds at
# the top of the heatmap.
foraging_driver_order <- cluster2_shift_trait_deltas |>
  filter(
    comparison == cluster2_shift_focal_comparison,
    trait_family == "Foraging guild"
  ) |>
  group_by(trait_value) |>
  summarise(mean_delta = mean(mean_expected_richness_delta_per_species, na.rm = TRUE), .groups = "drop") |>
  arrange(mean_delta) |>
  pull(trait_value)

foraging_driver_df <- cluster2_shift_trait_deltas |>
  filter(
    comparison == cluster2_shift_focal_comparison,
    trait_family == "Foraging guild"
  ) |>
  mutate(
    trait_value = factor(trait_value, levels = foraging_driver_order),
    transition = factor(transition_display_labels[.data$transition], levels = transition_display_levels)
  )

# ---- Foraging-guild trait annotations ----
# These side panels describe the species composition of each guild, not the
# size of the modelled shift itself.
guild_driver_annotations <- species_trait_table |>
  filter(!is.na(foraging_guild), foraging_guild != "") |>
  mutate(
    foraging_habitat = if_else(
      foraging_guild %in% aquatic_wetland_guilds,
      "Aquatic/wetland",
      "Terrestrial/aerial"
    )
  ) |>
  group_by(foraging_guild) |>
  summarise(
    trait_value = first(foraging_guild),
    foraging_habitat = if_else(
      mean(foraging_habitat == "Aquatic/wetland", na.rm = TRUE) >= 0.5,
      "Aquatic/wetland",
      "Terrestrial/aerial"
    ),
    n_species = n(),
    .groups = "drop"
  ) |>
  filter(trait_value %in% foraging_driver_order) |>
  mutate(
    trait_value = factor(trait_value, levels = foraging_driver_order),
    foraging_habitat = factor(foraging_habitat, levels = c("Aquatic/wetland", "Terrestrial/aerial"))
  )

guild_migration_composition <- species_trait_table |>
  filter(
    !is.na(foraging_guild), foraging_guild != "",
    !is.na(migratory_strategy), migratory_strategy != "",
    foraging_guild %in% foraging_driver_order
  ) |>
  count(foraging_guild, migratory_strategy, name = "n_species") |>
  group_by(foraging_guild) |>
  mutate(proportion = n_species / sum(n_species, na.rm = TRUE)) |>
  ungroup() |>
  complete(
    foraging_guild = foraging_driver_order,
    migratory_strategy = migration_levels,
    fill = list(n_species = 0, proportion = 0)
  ) |>
  mutate(
    trait_value = factor(foraging_guild, levels = foraging_driver_order),
    migratory_strategy = factor(migratory_strategy, levels = migration_levels)
  )

thermal_display_levels <- rev(thermal_levels)

guild_thermal_composition <- species_trait_table |>
  filter(
    !is.na(foraging_guild), foraging_guild != "",
    !is.na(thermal_group),
    foraging_guild %in% foraging_driver_order
  ) |>
  count(foraging_guild, thermal_group, name = "n_species") |>
  group_by(foraging_guild) |>
  mutate(proportion = n_species / sum(n_species, na.rm = TRUE)) |>
  ungroup() |>
  complete(
    foraging_guild = foraging_driver_order,
    thermal_group = thermal_display_levels,
    fill = list(n_species = 0, proportion = 0)
  ) |>
  mutate(
    trait_value = factor(foraging_guild, levels = foraging_driver_order),
    thermal_group = factor(thermal_group, levels = thermal_display_levels)
  )

overall_mean_sti <- mean(species_trait_table$species_thermal_index, na.rm = TRUE)

# Mean and 95% interval of species thermal index within each foraging guild.
guild_sti_summary <- species_trait_table |>
  filter(
    !is.na(foraging_guild), foraging_guild != "",
    !is.na(species_thermal_index),
    foraging_guild %in% foraging_driver_order
  ) |>
  group_by(foraging_guild) |>
  summarise(
    trait_value = first(foraging_guild),
    mean_sti = mean(species_thermal_index, na.rm = TRUE),
    sti_q025 = quantile(species_thermal_index, 0.025, na.rm = TRUE),
    sti_q975 = quantile(species_thermal_index, 0.975, na.rm = TRUE),
    n_species = n(),
    .groups = "drop"
  ) |>
  mutate(
    trait_value = factor(trait_value, levels = foraging_driver_order)
  )

# ---- Panel A: foraging guild heatmap with trait annotations ----
guild_habitat_strip <- ggplot(
  guild_driver_annotations,
  aes(x = "", y = trait_value, fill = foraging_habitat)
) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_manual(
    values = c("Aquatic/wetland" = "#2b8cbe", "Terrestrial/aerial" = "#7f8f3a"),
    name = "Foraging environment",
    guide = guide_legend(ncol = 1)
  ) +
  labs(x = 'Foraging environment', y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks = element_blank(),
    legend.position = "bottom"
  )

guild_migration_bar <- ggplot(
  guild_migration_composition,
  aes(x = proportion, y = trait_value, fill = migratory_strategy)
) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    breaks = c(0, 0.5, 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(
    values = migration_palette,
    breaks = migration_levels,
    drop = FALSE,
    name = "Migratory strategy",
    guide = guide_legend(ncol = 1, byrow = TRUE)
  ) +
  labs(x = "Migratory ttrategy", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 14, 5.5, 5.5),
    legend.position = "bottom"
  )

# Dashed line marks the mean STI across all species in the trait table.
guild_sti_forest_plot <- ggplot(
  guild_sti_summary,
  aes(y = trait_value)
) +
  geom_vline(xintercept = overall_mean_sti, colour = "grey35", linewidth = 0.35, linetype = "dashed") +
  geom_segment(
    aes(x = sti_q025, xend = sti_q975, yend = trait_value),
    colour = "grey45",
    linewidth = 0.55
  ) +
  geom_point(aes(x = mean_sti, colour = mean_sti), size = 2.4) +
  scale_x_continuous(
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  scale_colour_gradientn(
    colours = thermal_colours,
    name = "Species Thermal Index",
    guide = guide_colourbar(
      direction = "vertical",
      barheight = unit(38, "mm"),
      barwidth = unit(4, "mm")
    )
  ) +
  labs(x = "Species Thermal Index", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 14, 5.5, 5.5),
    legend.position = "bottom"
  )

# Main Panel A response surface: positive values indicate modelled increases in
# expected richness per species within that foraging guild.
cluster2_shift_foraging_heatmap <- ggplot(
  foraging_driver_df,
  aes(x = transition, y = trait_value, fill = mean_expected_richness_delta_per_species)
) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#d95f02",
    mid = "white",
    high = "#1b9e77",
    midpoint = 0,
    name = "Mean expected\nrichness delta\nper species",
    guide = guide_colourbar(
      direction = "vertical",
      barheight = unit(42, "mm"),
      barwidth = unit(4, "mm")
    )
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

cluster2_shift_foraging_plot <- (
  guild_habitat_strip + guild_migration_bar + guild_sti_forest_plot + cluster2_shift_foraging_heatmap +
    plot_layout(widths = c(1.05, 2.1, 2.1, 3.4), guides = "collect")
) &
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.box = "vertical",
    legend.direction = "vertical"
  )

# ---- Panels B and C: compact trait-driver heatmaps ----
# Use a common colour scale so migration and thermal heatmaps are directly
# comparable as mean expected-richness delta per species.
driver_delta_limit <- cluster2_shift_trait_deltas |>
  filter(
    comparison == cluster2_shift_focal_comparison,
    trait_family %in% c("Migratory strategy", "Thermal affinity")
  ) |>
  pull(mean_expected_richness_delta_per_species) |>
  abs() |>
  max(na.rm = TRUE)
driver_delta_limit <- driver_delta_limit * 1.08

driver_panel_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0)
  )

# Panel B keeps the migration-group colours as a left strip, and the response
# heatmap combines the two focal transitions into one axis.
migration_driver_heatmap_df <- thermal_migration_driver_df |>
  filter(trait_family == "Migratory strategy") |>
  mutate(
    trait_value = factor(as.character(.data$trait_value), levels = rev(migration_levels)),
    trait_colour_group = as.character(.data$trait_value)
  )

migration_driver_strip <- migration_driver_heatmap_df |>
  distinct(trait_value, trait_colour_group) |>
  ggplot(aes(x = "Group", y = trait_value, fill = trait_colour_group)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = migration_palette, drop = FALSE, guide = "none") +
  labs(x = NULL, y = NULL) +
  driver_panel_theme +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 9),
    plot.margin = margin(5.5, 2, 5.5, 5.5)
  )

migration_driver_heatmap <- migration_driver_heatmap_df |>
  ggplot(aes(
    x = transition,
    y = trait_value,
    fill = mean_expected_richness_delta_per_species
  )) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#d95f02",
    mid = "white",
    high = "#1b9e77",
    midpoint = 0,
    limits = c(-driver_delta_limit, driver_delta_limit),
    guide = "none"
  ) +
  scale_x_discrete(labels = transition_axis_labels) +
  labs(
    title = "Migratory strategies",
    x = NULL,
    y = NULL
  ) +
  driver_panel_theme +
  theme(
    axis.text.x = element_text(size = 8.2, lineheight = 0.9),
    axis.text.y = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 0)
  )

migration_driver_panel <- migration_driver_strip + migration_driver_heatmap +
  plot_layout(widths = c(0.62, 1.55))

# Panel C mirrors Panel B, but uses STI thermal groups and the Fig. 4 palette.
thermal_driver_heatmap_df <- thermal_migration_driver_df |>
  filter(trait_family == "Thermal affinity") |>
  mutate(
    trait_value = factor(as.character(.data$trait_value), levels = thermal_levels),
    trait_colour_group = as.character(.data$trait_value)
  )

thermal_driver_strip <- thermal_driver_heatmap_df |>
  distinct(trait_value, trait_colour_group) |>
  ggplot(aes(x = "Group", y = trait_value, fill = trait_colour_group)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = thermal_palette, drop = FALSE, guide = "none") +
  labs(x = NULL, y = NULL) +
  driver_panel_theme +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 9),
    plot.margin = margin(5.5, 2, 5.5, 5.5)
  )

thermal_driver_heatmap <- thermal_driver_heatmap_df |>
  ggplot(aes(
    x = transition,
    y = trait_value,
    fill = mean_expected_richness_delta_per_species
  )) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#d95f02",
    mid = "white",
    high = "#1b9e77",
    midpoint = 0,
    limits = c(-driver_delta_limit, driver_delta_limit),
    guide = "none"
  ) +
  scale_x_discrete(labels = transition_axis_labels) +
  labs(
    title = "Thermal affinity groups",
    x = NULL,
    y = NULL
  ) +
  driver_panel_theme +
  theme(
    axis.text.x = element_text(size = 8.2, lineheight = 0.9),
    axis.text.y = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 0)
  )

thermal_driver_panel <- thermal_driver_strip + thermal_driver_heatmap +
  plot_layout(widths = c(0.55, 1.55))

# ---- Assemble manuscript composite ----
# cowplot coordinates keep the mixed patchwork/ggplot objects in a stable layout:
# Panel A spans the top row; Panels B, C, and D sit in the bottom row.
trait_driver_panel_plot <- ggdraw() +
  draw_plot(cluster2_shift_foraging_plot, x = 0.01, y = 0.34, width = 0.98, height = 0.64) +
  draw_plot(migration_driver_panel, x = 0.02, y = 0.025, width = 0.31, height = 0.29) +
  draw_plot(thermal_driver_panel, x = 0.36, y = 0.025, width = 0.29, height = 0.29) +
  draw_plot(atlas1_bioregion_map, x = 0.69, y = 0.00, width = 0.29, height = 0.325) +
  draw_label("A", x = 0.012, y = 0.985, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("B", x = 0.012, y = 0.315, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("C", x = 0.342, y = 0.315, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("D", x = 0.662, y = 0.315, hjust = 0, vjust = 1, fontface = "bold", size = 14)

# ---- Save outputs ----
thermal_png <- file.path(out_dir, paste0(pattern, "-fig7-avian-bioregion-thermal-migration-drivers.png"))
thermal_pdf <- file.path(out_dir, paste0(pattern, "-fig7-avian-bioregion-thermal-migration-drivers.pdf"))
foraging_png <- file.path(out_dir, paste0(pattern, "-fig7-avian-bioregion-foraging-guild-drivers.png"))
foraging_pdf <- file.path(out_dir, paste0(pattern, "-fig7-avian-bioregion-foraging-guild-drivers.pdf"))
trait_driver_png <- file.path(out_dir, paste0(pattern, "-fig7-avian-bioregion-trait-driver-panels.png"))
trait_driver_pdf <- file.path(out_dir, paste0(pattern, "-fig7-avian-bioregion-trait-driver-panels.pdf"))

ggsave(thermal_png, cluster2_shift_thermal_migration_plot, width = 10.5, height = 6.5, units = "in", dpi = 300, bg = "white")
ggsave(thermal_pdf, cluster2_shift_thermal_migration_plot, width = 10.5, height = 6.5, units = "in", bg = "white")
ggsave(foraging_png, cluster2_shift_foraging_plot, width = 16, height = 8.8, units = "in", dpi = 300, bg = "white")
ggsave(foraging_pdf, cluster2_shift_foraging_plot, width = 16, height = 8.8, units = "in", bg = "white")
ggsave(trait_driver_png, trait_driver_panel_plot, width = 13, height = 11.5, units = "in", dpi = 300, bg = "white")
ggsave(trait_driver_pdf, trait_driver_panel_plot, width = 13, height = 11.5, units = "in", bg = "white")

message("Saved: ", thermal_png)
message("Saved: ", thermal_pdf)
message("Saved: ", foraging_png)
message("Saved: ", foraging_pdf)
message("Saved: ", trait_driver_png)
message("Saved: ", trait_driver_pdf)
