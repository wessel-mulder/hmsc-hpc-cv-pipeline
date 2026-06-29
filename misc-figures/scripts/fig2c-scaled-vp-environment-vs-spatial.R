rm(list = ls())

# Companion to Figure 2: scaled environmental versus spatial variance.
#
# The variance-partitioning CSVs store relative VP values for each species, plus
# a TjurR2 row. `load_vp_estimates(..., scaled = TRUE)` multiplies each VP column
# by its species-level TjurR2, so every species' stacked components sum to the
# absolute variance explained by the model for that species.

library(tidyverse)
library(scales)

source(file.path("support_scripts", "figure_data_helpers.R"))

model_dir <- "HmscOutputs"
model_pattern <- "2026-03-13"
output_dir <- file.path("misc-figures", "outputs", "main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure_slug <- paste0(model_pattern, "-fig2c-scaled-vp-environment-vs-spatial")
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))
output_species_csv <- file.path(output_dir, paste0(figure_slug, "-species.csv"))
output_summary_csv <- file.path(output_dir, paste0(figure_slug, "-summary.csv"))

atlas_lookup <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

atlas_order <- unname(atlas_lookup)

component_order <- c(
  "Climate",
  "Land-use",
  "Spatial"
)

component_colors <- c(
  "Climate" = "firebrick3",
  "Land-use" = "goldenrod1",
  "Spatial" = "ivory2"
)

climate_rows <- c("tmean_breeding", "prec_breeding")

model_folders <- figure_model_folders(pattern = model_pattern, base_dir = model_dir)
vp_scaled <- load_vp_estimates(model_folders, base_dir = model_dir, scaled = TRUE)
vp_raw <- load_vp_estimates(model_folders, base_dir = model_dir, scaled = FALSE)

read_scaled_components <- function(scaled_df, raw_df, atlas_id) {
  random_row <- "Random: site"

  if (!random_row %in% rownames(scaled_df)) {
    stop("Missing `Random: site` in scaled VP matrix for atlas ", atlas_id, call. = FALSE)
  }

  if (!"TjurR2" %in% rownames(raw_df)) {
    stop("Missing `TjurR2` in raw VP matrix for atlas ", atlas_id, call. = FALSE)
  }

  measured_rows <- setdiff(rownames(scaled_df), random_row)
  available_climate_rows <- intersect(measured_rows, climate_rows)
  land_use_rows <- setdiff(measured_rows, available_climate_rows)

  climate_scaled <- colSums(scaled_df[available_climate_rows, , drop = FALSE], na.rm = TRUE)
  land_use_scaled <- colSums(scaled_df[land_use_rows, , drop = FALSE], na.rm = TRUE)
  environment_scaled <- climate_scaled + land_use_scaled
  spatial_scaled <- as.numeric(scaled_df[random_row, ])
  names(spatial_scaled) <- colnames(scaled_df)

  tibble(
    atlas = atlas_lookup[[atlas_id]],
    atlas_index = as.integer(atlas_id),
    species = colnames(scaled_df),
    climate = climate_scaled,
    land_use = land_use_scaled,
    environmental_variables = environment_scaled,
    spatial_biotic_dependencies = spatial_scaled,
    total_explained = environment_scaled + spatial_scaled,
    tjur_r2 = as.numeric(raw_df["TjurR2", colnames(scaled_df)]),
    scaling_check_abs_diff = abs(total_explained - tjur_r2)
  )
}

species_scaled_vp <- imap_dfr(vp_scaled, function(scaled_df, atlas_id) {
  read_scaled_components(
    scaled_df = scaled_df,
    raw_df = vp_raw[[atlas_id]],
    atlas_id = atlas_id
  )
}) |>
  mutate(
    atlas = factor(.data$atlas, levels = atlas_order)
  )

if (any(species_scaled_vp$scaling_check_abs_diff > 1e-8, na.rm = TRUE)) {
  stop("Scaled VP components do not sum back to TjurR2 for at least one species.", call. = FALSE)
}

species_long <- species_scaled_vp |>
  select(
    atlas,
    atlas_index,
    species,
    total_explained,
    climate,
    land_use,
    spatial_biotic_dependencies
  ) |>
  group_by(.data$atlas) |>
  arrange(.data$total_explained, .by_group = TRUE) |>
  mutate(species_rank = row_number()) |>
  ungroup() |>
  pivot_longer(
    cols = c(climate, land_use, spatial_biotic_dependencies),
    names_to = "component",
    values_to = "scaled_variance_explained"
  ) |>
  mutate(
    component = recode(
      .data$component,
      climate = "Climate",
      land_use = "Land-use",
      spatial_biotic_dependencies = "Spatial"
    ),
    component = factor(.data$component, levels = component_order)
  )

summary_long <- species_long |>
  group_by(.data$atlas, .data$atlas_index, .data$component) |>
  summarise(
    n_species = n_distinct(.data$species),
    mean_scaled_variance_explained = mean(.data$scaled_variance_explained, na.rm = TRUE),
    median_scaled_variance_explained = median(.data$scaled_variance_explained, na.rm = TRUE),
    q25_scaled_variance_explained = quantile(.data$scaled_variance_explained, 0.25, na.rm = TRUE),
    q75_scaled_variance_explained = quantile(.data$scaled_variance_explained, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

total_summary <- species_scaled_vp |>
  group_by(.data$atlas, .data$atlas_index) |>
  summarise(
    n_species = n_distinct(.data$species),
    mean_total_explained = mean(.data$total_explained, na.rm = TRUE),
    median_total_explained = median(.data$total_explained, na.rm = TRUE),
    q25_total_explained = quantile(.data$total_explained, 0.25, na.rm = TRUE),
    q75_total_explained = quantile(.data$total_explained, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

summary_export <- summary_long |>
  left_join(total_summary, by = c("atlas", "atlas_index", "n_species")) |>
  mutate(
    mean_component_share_of_total = .data$mean_scaled_variance_explained / .data$mean_total_explained
  )

inset_summary <- summary_export |>
  select("atlas", "component", "mean_scaled_variance_explained") |>
  pivot_wider(
    names_from = "component",
    values_from = "mean_scaled_variance_explained"
  ) |>
  left_join(
    total_summary |>
      select("atlas", "mean_total_explained"),
    by = "atlas"
  ) |>
  mutate(
    atlas = factor(.data$atlas, levels = atlas_order),
    inset_left = 3,
    inset_right = 68,
    inset_top = max(species_scaled_vp$total_explained, na.rm = TRUE) * 0.98,
    inset_gap = max(species_scaled_vp$total_explained, na.rm = TRUE) * 0.05,
    box_ymax = .data$inset_top + .data$inset_gap * 0.15,
    box_ymin = .data$inset_top - .data$inset_gap * 4.05,
    mean_x = .data$inset_left + 1.8,
    mean_y = .data$inset_top - .data$inset_gap * 0.45,
    mean_label = paste0(
      "mean R2: ",
      percent(.data$mean_total_explained, accuracy = 0.1)
    )
  )

inset_components <- inset_summary |>
  select(
    "atlas",
    "inset_left",
    "inset_top",
    "inset_gap",
    "Spatial",
    "Land-use",
    "Climate"
  ) |>
  pivot_longer(
    cols = c(
      "Spatial",
      "Land-use",
      "Climate"
    ),
    names_to = "component",
    values_to = "mean_scaled_variance_explained"
  ) |>
  mutate(
    component = factor(
      .data$component,
      levels = c("Spatial", "Land-use", "Climate")
    ),
    component_label = recode(
      .data$component,
      "Spatial" = "spatial",
      "Land-use" = "land-use",
      "Climate" = "climate"
    ),
    square_x = .data$inset_left + 3.2,
    text_x = .data$inset_left + 9,
    y = .data$inset_top - .data$inset_gap * (as.integer(.data$component) + 0.45),
    label = paste0(
      .data$component_label,
      ": ",
      percent(.data$mean_scaled_variance_explained, accuracy = 0.1)
    )
  )

write_csv(species_scaled_vp, output_species_csv)
write_csv(summary_export, output_summary_csv)

base_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9),
    plot.background = element_rect(fill = "grey99", color = NA),
    panel.background = element_rect(fill = "grey99", color = NA),
    strip.text = element_text(face = "bold")
  )

species_panel <- ggplot(
  species_long,
  aes(
    x = .data$species_rank,
    y = .data$scaled_variance_explained,
    fill = .data$component
  )
) +
  geom_col(width = 1, color = "grey85", linewidth = 0.02, position = position_stack(reverse = TRUE)) +
  geom_text(
    data = inset_summary,
    aes(
      x = .data$mean_x,
      y = .data$mean_y,
      label = .data$mean_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 3.1,
    fontface = "bold",
    color = "grey15"
  ) +
  geom_point(
    data = inset_components,
    aes(
      x = .data$square_x,
      y = .data$y,
      fill = .data$component
    ),
    inherit.aes = FALSE,
    shape = 22,
    size = 3.4,
    stroke = 0.22,
    color = "grey35",
    show.legend = FALSE
  ) +
  geom_text(
    data = inset_components,
    aes(
      x = .data$text_x,
      y = .data$y,
      label = .data$label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 3,
    color = "grey15"
  ) +
  facet_wrap(vars(.data$atlas), nrow = 1) +
  scale_fill_manual(values = component_colors) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Species-level scaled variance partitioning",
    subtitle = "Each thin bar is one species; total bar height is Tjur R2.",
    x = "Species ordered by total explained variance within atlas",
    y = "Absolute variance explained"
  ) +
  base_theme +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none"
  )

final_plot <- species_panel

ggsave(
  filename = output_png,
  plot = final_plot,
  width = 9,
  height = 4.4,
  dpi = 320
)

message("Wrote: ", output_png)
message("Wrote: ", output_species_csv)
message("Wrote: ", output_summary_csv)
