rm(list = ls())

# Standalone panels for the amount of variation explained by traits.
# VP$R2T$Y gives one atlas-level summary for occurrence variation overall.
# VP$R2T$Beta gives environmental-variable summaries for the species responses
# to each predictor.

library(tidyverse)
library(patchwork)
library(scales)

source(file.path("support_scripts", "project_paths.R"))

model_dir <- "HmscOutputs"
model_pattern <- "2026-03-13"

output_dir <- file.path("misc-figures", "outputs", "main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure_slug <- paste0(model_pattern, "-trait-occurrence-r2t-panel")
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))
output_csv <- file.path(output_dir, paste0(figure_slug, ".csv"))

env_figure_slug <- paste0(model_pattern, "-trait-environment-r2t-panel")
env_output_png <- file.path(output_dir, paste0(env_figure_slug, ".png"))
env_output_csv <- file.path(output_dir, paste0(env_figure_slug, ".csv"))

sum_figure_slug <- paste0(model_pattern, "-trait-environment-r2t-average-panel")
sum_output_png <- file.path(output_dir, paste0(sum_figure_slug, ".png"))
sum_output_csv <- file.path(output_dir, paste0(sum_figure_slug, ".csv"))

combined_figure_slug <- paste0(model_pattern, "-trait-r2t-summary-panels")
combined_output_png <- file.path(output_dir, paste0(combined_figure_slug, ".png"))

var_rename <- c(
  "tmean_breeding" = "Temperature",
  "prec_breeding" = "Precipitation",
  "hh" = "Land-use heterogeneity",
  "perc_urban" = "Urban",
  "perc_cropland" = "Cropland",
  "perc_pasture" = "Pasture",
  "perc_forest" = "Forest",
  "perc_grass_shrub" = "Grass/shrubland"
)

var_order <- unname(var_rename)

env_colors <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban" = "snow3",
  "Cropland" = "goldenrod1",
  "Pasture" = "darkorange",
  "Forest" = "springgreen4",
  "Grass/shrubland" = "springgreen2"
)

atlas_rename <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

read_r2t_y <- function(model_folder) {
  atlas_num <- str_extract(model_folder, "(?<=Atlas)\\d+")
  file_path <- file.path(
    model_dir,
    model_folder,
    "Results",
    sprintf("%sparameter_estimates_VP_R2T_Y.csv", model_folder)
  )

  if (!file.exists(file_path)) {
    stop("Missing VP R2T Y export: ", file_path, call. = FALSE)
  }

  read.csv(file_path, check.names = TRUE) |>
    as_tibble() |>
    rename(response_component = 1, trait_r2t_y = 2) |>
    mutate(
      atlas = atlas_rename[atlas_num],
      atlas = factor(.data$atlas, levels = unname(atlas_rename))
    ) |>
    select(atlas, response_component, trait_r2t_y)
}

read_r2t_beta <- function(model_folder) {
  atlas_num <- str_extract(model_folder, "(?<=Atlas)\\d+")
  file_path <- file.path(
    model_dir,
    model_folder,
    "Results",
    sprintf("%sparameter_estimates_VP_R2T_Beta.csv", model_folder)
  )

  if (!file.exists(file_path)) {
    stop("Missing VP R2T Beta export: ", file_path, call. = FALSE)
  }

  read.csv(file_path, check.names = TRUE) |>
    as_tibble() |>
    rename(variable_raw = 1, trait_r2t_beta = 2) |>
    filter(.data$variable_raw %in% names(var_rename)) |>
    mutate(
      atlas = atlas_rename[atlas_num],
      atlas = factor(.data$atlas, levels = unname(atlas_rename)),
      variable = var_rename[.data$variable_raw],
      variable = factor(.data$variable, levels = var_order)
    ) |>
    select(atlas, variable, variable_raw, trait_r2t_beta)
}

model_folders <- find_model_folders(base_dir = model_dir, pattern = model_pattern)
if (length(model_folders) == 0) {
  stop("No model folders found for pattern: ", model_pattern, call. = FALSE)
}

r2t_y <- map_dfr(model_folders, read_r2t_y) |>
  arrange(.data$atlas)

r2t_beta <- map_dfr(model_folders, read_r2t_beta) |>
  arrange(.data$variable, .data$atlas)

r2t_beta_sum <- r2t_beta |>
  group_by(.data$atlas) |>
  summarise(
    summed_trait_r2t_beta = sum(.data$trait_r2t_beta, na.rm = TRUE),
    mean_trait_r2t_beta = mean(.data$trait_r2t_beta, na.rm = TRUE),
    n_variables = n(),
    .groups = "drop"
  ) |>
  arrange(.data$atlas)

r2t_beta_average_stack <- r2t_beta |>
  left_join(
    r2t_beta_sum |>
      select(atlas, n_variables),
    by = "atlas"
  ) |>
  mutate(
    average_trait_r2t_beta = .data$trait_r2t_beta / .data$n_variables,
    variable = factor(.data$variable, levels = rev(var_order))
  )

write_csv(r2t_y, output_csv)
write_csv(r2t_beta, env_output_csv)
write_csv(r2t_beta_sum, sum_output_csv)

r2t_y_limit <- max(r2t_y$trait_r2t_y, na.rm = TRUE) * 1.18
r2t_beta_limit <- max(r2t_beta$trait_r2t_beta, na.rm = TRUE) * 1.18
r2t_beta_mean_limit <- max(r2t_beta_sum$mean_trait_r2t_beta, na.rm = TRUE) * 1.18

r2t_y_panel <- ggplot(r2t_y, aes(x = atlas, y = trait_r2t_y, group = 1)) +
  geom_col(width = 0.62, fill = "grey45", colour = "white", linewidth = 0.25) +
  geom_line(colour = "grey20", linewidth = 0.35) +
  geom_point(size = 2.4, colour = "grey20") +
  geom_text(
    aes(label = percent(.data$trait_r2t_y, accuracy = 0.1)),
    vjust = -0.55,
    size = 3.2,
    colour = "grey15"
  ) +
  scale_y_continuous(
    limits = c(0, r2t_y_limit),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Occurrence variation explained by traits",
    x = NULL,
    y = "Trait-explained\noccurrence variation"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

ggsave(output_png, r2t_y_panel, width = 4.6, height = 2.6, units = "in", dpi = 300, bg = "white")

env_r2t_panel <- ggplot(r2t_beta, aes(x = atlas, y = trait_r2t_beta, group = 1)) +
  geom_col(aes(fill = variable), width = 0.62, colour = "white", linewidth = 0.2) +
  geom_line(colour = "grey20", linewidth = 0.32) +
  geom_point(size = 1.8, colour = "grey20") +
  geom_text(
    aes(label = percent(.data$trait_r2t_beta, accuracy = 1)),
    vjust = -0.45,
    size = 2.6,
    colour = "grey15"
  ) +
  facet_wrap(vars(variable), ncol = 4) +
  scale_fill_manual(values = env_colors, guide = "none") +
  scale_y_continuous(
    limits = c(0, r2t_beta_limit),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    title = "Species-environment response variation explained by traits",
    x = NULL,
    y = "Trait-explained\nBeta variation"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8.2, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(size = 8.6, face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = "white"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.spacing = grid::unit(0.75, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

ggsave(env_output_png, env_r2t_panel, width = 8.2, height = 4.8, units = "in", dpi = 300, bg = "white")

summed_r2t_panel <- ggplot(
  r2t_beta_average_stack,
  aes(x = atlas, y = average_trait_r2t_beta, fill = variable)
) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.18) +
  geom_line(
    data = r2t_beta_sum,
    aes(x = atlas, y = mean_trait_r2t_beta, group = 1),
    inherit.aes = FALSE,
    colour = "grey20",
    linewidth = 0.35
  ) +
  geom_point(
    data = r2t_beta_sum,
    aes(x = atlas, y = mean_trait_r2t_beta),
    inherit.aes = FALSE,
    size = 2.4,
    colour = "grey20"
  ) +
  geom_text(
    data = r2t_beta_sum,
    aes(
      x = atlas,
      y = mean_trait_r2t_beta,
      label = percent(.data$mean_trait_r2t_beta, accuracy = 1)
    ),
    inherit.aes = FALSE,
    vjust = -0.55,
    size = 3.2,
    colour = "grey15"
  ) +
  scale_fill_manual(
    values = env_colors,
    breaks = var_order,
    name = NULL,
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  scale_y_continuous(
    limits = c(0, r2t_beta_mean_limit),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Average trait-explained environmental responses",
    x = NULL,
    y = "Average trait-explained\nBeta variation"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "bottom",
    legend.text = element_text(size = 7.4),
    legend.key.width = grid::unit(0.9, "lines"),
    legend.key.height = grid::unit(0.75, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

ggsave(sum_output_png, summed_r2t_panel, width = 5.3, height = 3.25, units = "in", dpi = 300, bg = "white")

combined_r2t_panel <- r2t_y_panel + summed_r2t_panel +
  plot_layout(widths = c(1, 1.08)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.01, 0.99)
  )

ggsave(combined_output_png, combined_r2t_panel, width = 10.4, height = 3.55, units = "in", dpi = 300, bg = "white")

message("Wrote: ", output_png)
message("Wrote: ", output_csv)
message("Wrote: ", env_output_png)
message("Wrote: ", env_output_csv)
message("Wrote: ", sum_output_png)
message("Wrote: ", sum_output_csv)
message("Wrote: ", combined_output_png)
