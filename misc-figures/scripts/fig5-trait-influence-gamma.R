rm(list = ls())

library(tidyverse)
library(readxl)
library(patchwork)

source(file.path("support_scripts", "project_paths.R"))

model_dir <- "HmscOutputs"
model_pattern <- "2026-03-13"
min_species_per_trait_level <- 5
support_level <- 0.95

output_dir <- file.path("misc-figures", "outputs", "main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure_slug <- paste0(model_pattern, "-fig5-trait-influence-gamma")
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))
output_pdf <- file.path(output_dir, paste0(figure_slug, ".pdf"))

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

atlas_rename <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

atlas_sizes <- c(
  "1970s" = 1.8,
  "1990s" = 3.4,
  "2010s" = 5.3
)

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

category_colors <- c(
  "Thermal index" = "#3B6D11",
  "Migration" = "#7F77DD",
  "Foraging guild" = "#D85A30"
)

reference_categories <- paste(
  "Reference categories:",
  "migration = long-distance migrants;",
  "foraging guild = aerial insectivores."
)

migration_order <- c(
  "sedentary",
  "sedentary and short-distance",
  "short-distance",
  "short-and long-distance"
)

foraging_order <- c(
  "Tit-like birds",
  "Low flycatching feeders",
  "Dabbling ducks",
  "Scolopacids",
  "Plovers"
)

trait_order_requested <- c(
  foraging_order,
  rev(migration_order),
  "Thermal index"
)

clean_trait_name <- function(trait) {
  trait |>
    str_remove("^Migration_a3_DOF") |>
    str_remove("^foraging_guild_consensus") |>
    str_replace("^species_thermal_index$", "Thermal index") |>
    str_trim()
}

trait_category <- function(trait) {
  case_when(
    str_starts(trait, "Migration") ~ "Migration",
    str_starts(trait, "foraging_guild") ~ "Foraging guild",
    str_starts(trait, "species_therm") ~ "Thermal index",
    TRUE ~ NA_character_
  )
}

read_gamma_effects <- function(model_folder) {
  atlas_num <- str_extract(model_folder, "(?<=Atlas)\\d+")
  file_path <- file.path(
    model_dir,
    model_folder,
    "Results",
    sprintf("%sparameter_estimates_Gamma_.xlsx", model_folder)
  )

  if (!file.exists(file_path)) {
    stop("Missing Gamma export: ", file_path, call. = FALSE)
  }

  means <- read_excel(file_path, sheet = "Posterior mean")
  support_pos <- read_excel(file_path, sheet = "Pr(x>0)")

  means_long <- means |>
    pivot_longer(
      cols = -Traits,
      names_to = "variable",
      values_to = "effect_size"
    )

  support_long <- support_pos |>
    pivot_longer(
      cols = -Traits,
      names_to = "variable",
      values_to = "support_pos"
    )

  means_long |>
    left_join(support_long, by = c("Traits", "variable")) |>
    mutate(
      atlas = atlas_num,
      significant = support_pos >= support_level | support_pos <= 1 - support_level,
      direction = case_when(
        support_pos >= support_level ~ "Positive",
        support_pos <= 1 - support_level ~ "Negative",
        TRUE ~ "Not supported"
      )
    )
}

model_folders <- find_model_folders(base_dir = model_dir, pattern = model_pattern)
if (length(model_folders) == 0) {
  stop("No model folders found for pattern: ", model_pattern, call. = FALSE)
}

gamma_raw <- map_dfr(model_folders, read_gamma_effects)

preprocessed_file <- "Data/preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData"
if (!file.exists(preprocessed_file)) {
  stop("Missing preprocessed trait data: ", preprocessed_file, call. = FALSE)
}

load(preprocessed_file)

migration_keep <- names(which(table(Tr$Migration_a3_DOF) >= min_species_per_trait_level))
foraging_keep <- names(which(table(Tr$foraging_guild_consensus) >= min_species_per_trait_level))

trait_counts <- bind_rows(
  tibble(
    trait_category = "Migration",
    trait_clean = names(table(Tr$Migration_a3_DOF)),
    n_species = as.integer(table(Tr$Migration_a3_DOF))
  ),
  tibble(
    trait_category = "Foraging guild",
    trait_clean = names(table(Tr$foraging_guild_consensus)),
    n_species = as.integer(table(Tr$foraging_guild_consensus))
  ),
  tibble(
    trait_category = "Thermal index",
    trait_clean = "Thermal index",
    n_species = nrow(Tr)
  )
)

gamma_plot <- gamma_raw |>
  filter( 
    Traits != "(Intercept)",
    variable != "(Intercept)",
    variable %in% names(var_rename)
  ) |>
  mutate(
    trait_category = trait_category(Traits),
    trait_clean = clean_trait_name(Traits),
    variable = var_rename[variable],
    atlas = atlas_rename[atlas]
  ) |>
  filter(!is.na(trait_category), !is.na(variable), !is.na(atlas)) |>
  filter(
    (trait_category == "Migration" & trait_clean %in% migration_keep) |
      (trait_category == "Foraging guild" & trait_clean %in% foraging_keep) |
      trait_category == "Thermal index"
  ) |>
  left_join(trait_counts, by = c("trait_category", "trait_clean")) |>
  mutate(
    trait_category = factor(
      trait_category,
      levels = c("Foraging guild", "Migration", "Thermal index")
    ),
    variable = factor(variable, levels = var_order),
    atlas = factor(atlas, levels = names(atlas_sizes)),
    direction = factor(direction, levels = c("Positive", "Negative", "Not supported"))
  )

traits_with_supported_effects <- gamma_plot |>
  filter(significant) |>
  distinct(trait_clean) |>
  pull(trait_clean)

gamma_plot <- gamma_plot |>
  filter(trait_clean %in% traits_with_supported_effects)

trait_order <- trait_order_requested[trait_order_requested %in% unique(gamma_plot$trait_clean)]

plot_data <- gamma_plot |>
  mutate(
    trait_clean = factor(trait_clean, levels = trait_order),
    x = as.numeric(variable),
    y = as.numeric(trait_clean)
  )

grid_cells <- plot_data |>
  distinct(trait_clean, variable, x, y)

category_bands <- plot_data |>
  distinct(trait_clean, trait_category) |>
  mutate(y = as.numeric(trait_clean)) |>
  group_by(trait_category) |>
  summarise(
    ymin = min(y) - 0.5,
    ymax = max(y) + 0.5,
    ymid = mean(y),
    .groups = "drop"
  ) |>
  mutate(
    category_fill = category_colors[as.character(trait_category)],
    category_text = category_colors[as.character(trait_category)]
  )

plot_points <- plot_data |>
  filter(significant) |>
  arrange(desc(atlas))

legend_envs <- intersect(var_order, as.character(unique(plot_points$variable)))

n_envs <- length(var_order)
divider_lines <- category_bands |>
  filter(trait_category != last(levels(plot_data$trait_category)))

make_direction_panel <- function(direction_label, show_environment_legend = TRUE) {
  points <- plot_points |>
    filter(direction == direction_label)
  environment_guide <- if (show_environment_legend) {
    guide_legend(
      order = 2,
      nrow = 2,
      override.aes = list(shape = 21, fill = NA, size = 3)
    )
  } else {
    "none"
  }

  ggplot() +
    geom_rect(
      data = grid_cells,
      aes(xmin = x - 0.5, xmax = x + 0.5, ymin = y - 0.5, ymax = y + 0.5),
      fill = "grey98",
      color = "grey88",
      linewidth = 0.18
    ) +
    geom_hline(
      data = divider_lines,
      aes(yintercept = ymax),
      color = "grey55",
      linewidth = 0.35,
      linetype = "dashed"
    ) +
    geom_point(
      data = points |> filter(atlas == "2010s"),
      aes(
        x = x,
        y = y,
        size = atlas,
        color = variable
      ),
      shape = 21,
      fill = NA,
      stroke = 1.25,
      alpha = 0.95
    ) +
    geom_point(
      data = points |> filter(atlas == "1990s"),
      aes(
        x = x,
        y = y,
        size = atlas,
        color = variable
      ),
      shape = 21,
      fill = NA,
      stroke = 1.05,
      alpha = 0.95
    ) +
    geom_point(
      data = points |> filter(atlas == "1970s"),
      aes(
        x = x,
        y = y,
        size = atlas,
        color = variable,
        fill = variable
      ),
      shape = 21,
      stroke = 0.35,
      alpha = 0.95
    ) +
    geom_rect(
      data = category_bands,
      aes(
        xmin = n_envs + 0.58,
        xmax = n_envs + 0.83,
        ymin = ymin,
        ymax = ymax
      ),
      fill = category_bands$category_fill,
      color = NA
    ) +
    geom_text(
      data = category_bands,
      aes(
        x = n_envs + 1.0,
        y = ymid,
        label = trait_category
      ),
      color = category_bands$category_text,
      hjust = 0,
      size = 3,
      fontface = "bold"
    ) +
    scale_size_manual(
      values = atlas_sizes,
      breaks = names(atlas_sizes),
      limits = names(atlas_sizes),
      name = "Atlas period",
      guide = guide_legend(override.aes = list(fill = "grey65"))
    ) +
    scale_color_manual(
      values = env_colors,
      breaks = legend_envs,
      limits = legend_envs,
      name = "Environmental variable"
    ) +
    scale_fill_manual(
      values = env_colors,
      breaks = legend_envs,
      limits = legend_envs,
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = seq_along(var_order),
      labels = var_order,
      expand = expansion(add = c(0.5, 3.15))
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(plot_data$trait_clean)),
      labels = levels(plot_data$trait_clean),
      expand = expansion(add = 0.5)
    ) +
    coord_fixed(clip = "off") +
    labs(
      title = paste(direction_label, "trait moderation"),
      subtitle = paste0(
        "Circles mark Gamma coefficients with Pr(>0) >= ",
        support_level,
        " or <= ",
        1 - support_level,
        ". 1990s and 2010s are rings."
      ),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, color = "grey35"),
      plot.margin = margin(8, 78, 8, 8)
    ) +
    guides(
      size = guide_legend(
        order = 1,
        override.aes = list(
          shape = 21,
          color = "grey30",
          fill = c("grey65", NA, NA),
          stroke = c(0.35, 1.05, 1.25)
        )
      ),
      color = environment_guide
    )
}

p_positive <- make_direction_panel("Positive", show_environment_legend = TRUE)
p_negative <- make_direction_panel("Negative", show_environment_legend = FALSE)

p_combined <- p_positive + p_negative +
  plot_layout(guides = "collect") +
  plot_annotation(caption = reference_categories) &
  theme(
    legend.position = "bottom",
    plot.caption = element_text(size = 8.5, color = "grey35", hjust = 0)
  )

ggsave(output_png, p_combined, width = 18, height = 10, dpi = 300)
ggsave(output_pdf, p_combined, width = 18, height = 10)

# Compatibility alias for the old draft filename. The publication filename above
# is the canonical output.
ggsave(file.path(output_dir, "gamma_circles.png"), p_combined, width = 18, height = 10, dpi = 300)
ggsave(file.path(output_dir, "gamma_circles.pdf"), p_combined, width = 18, height = 10)

message("Wrote: ", output_png)
message("Wrote: ", output_pdf)
