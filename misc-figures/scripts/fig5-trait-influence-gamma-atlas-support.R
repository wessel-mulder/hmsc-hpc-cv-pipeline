rm(list = ls())

# Atlas-split Gamma support figure.
#
# This is an alternate view of the trait-moderation Gamma coefficients. Each
# tile asks whether a trait level has atlas-specific posterior support for a
# positive or negative moderation of an environmental response. Environmental
# variables are shown as separate panels; atlas periods are on the x-axis within
# each panel.

library(tidyverse)
library(readxl)
library(patchwork)

source(file.path("support_scripts", "project_paths.R"))

model_dir <- "HmscOutputs"
model_pattern <- "2026-03-13"
min_species_per_trait_level <- 2
support_level <- 0.95

output_dir <- file.path("misc-figures", "outputs", "main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure_slug <- paste0(model_pattern, "-fig5-trait-influence-gamma-atlas-support")
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))
output_csv <- file.path(output_dir, paste0(figure_slug, ".csv"))

var_rename <- c(
  "tmean_breeding" = "Temperature",
  "prec_breeding" = "Precipitation",
  "hh" = "Land-Use Heterogeneity",
  "perc_urban" = "Urban",
  "perc_cropland" = "Cropland",
  "perc_pasture" = "Pasture",
  "perc_forest" = "Forest",
  "perc_grass_shrub" = "Grass/Shrubland"
)

var_order <- unname(var_rename)

atlas_rename <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

migration_order <- c(
  "long-distance",
  "short-and long-distance",
  "short-distance",
  "sedentary and short-distance",
  "sedentary"
)

foraging_order <- c(
  "Passerine seedeaters",
  "Omnivorous corvids",
  "Tit-like birds",
  "Foliage gleaners",
  "Low flycatching feeders",
  "Marsh warblers",
  "Open-land insectivores",
  "Aerial insectivores",
  "Thrushes",
  "Buntings",
  "Flycatchers",
  "Shrikes",
  "Diurnal raptors",
  "Owls",
  "Woodpeckers",
  "Columbids",
  "Gallinaceous birds",
  "Scolopacids",
  "Dabbling ducks",
  "Aquatic pursuers",
  "Plovers",
  "Gulls",
  "Rails",
  "Terns",
  "Diving ducks",
  "Grazing waterfowl",
  "Wading birds"
)

trait_category_labels <- c(
  "Species thermal index" = "Species thermal index",
  "Migration" = "Migratory strategy",
  "Foraging guild" = "Foraging guild"
)

gamma_colors <- c(
  "Negative" = "#d95f02",
  "No supported effect" = "grey98",
  "Positive" = "#1b9e77"
)

figure_caption <- paste0(
  "Tiles show atlas-specific supported Gamma coefficients for each trait level ",
  "and environmental variable. Species counts are shown in parentheses ",
  "(green = positive; red = negative; white = no ",
  "supported effect). Supported Gamma coefficients have Pr(x > 0) >= ",
  support_level, " or Pr(x > 0) <= ", 1 - support_level,
  ". Foraging guilds are shown only when they have at least one supported ",
  "Gamma coefficient. Section titles identify treatment-coded reference ",
  "categories."
)

clean_trait_name <- function(trait) {
  trait |>
    str_remove("^Migration_a3_DOF") |>
    str_remove("^foraging_guild_consensus") |>
    str_replace("^species_thermal_index$", "Species thermal index") |>
    str_trim()
}

format_trait_label <- function(trait) {
  recode(
    trait,
    "long-distance" = "Long-Distance",
    "short-and long-distance" = "Short- and Long-Distance",
    "short-distance" = "Short-Distance",
    "sedentary and short-distance" = "Sedentary and Short-Distance",
    "sedentary" = "Sedentary",
    .default = trait
  )
}

trait_category <- function(trait) {
  case_when(
    str_starts(trait, "Migration") ~ "Migration",
    str_starts(trait, "foraging_guild") ~ "Foraging guild",
    str_starts(trait, "species_therm") ~ "Species thermal index",
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
      significant = .data$support_pos >= support_level |
        .data$support_pos <= 1 - support_level,
      direction = case_when(
        .data$support_pos >= support_level ~ "Positive",
        .data$support_pos <= 1 - support_level ~ "Negative",
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
    trait_category = "Species thermal index",
    trait_clean = "Species thermal index",
    n_species = nrow(Tr)
  )
)

gamma_plot <- gamma_raw |>
  filter(
    .data$Traits != "(Intercept)",
    .data$variable != "(Intercept)",
    .data$variable %in% names(var_rename)
  ) |>
  mutate(
    trait_category = trait_category(.data$Traits),
    trait_clean = clean_trait_name(.data$Traits),
    variable = var_rename[.data$variable],
    atlas = atlas_rename[.data$atlas]
  ) |>
  filter(!is.na(.data$trait_category), !is.na(.data$variable), !is.na(.data$atlas)) |>
  filter(
      (.data$trait_category == "Migration" & .data$trait_clean %in% migration_keep) |
      (.data$trait_category == "Foraging guild" & .data$trait_clean %in% foraging_keep) |
      .data$trait_category == "Species thermal index"
  ) |>
  left_join(trait_counts, by = c("trait_category", "trait_clean")) |>
  mutate(
    trait_category = factor(
      .data$trait_category,
      levels = c("Species thermal index", "Migration", "Foraging guild")
    ),
    variable = factor(.data$variable, levels = var_order),
    atlas = factor(.data$atlas, levels = unname(atlas_rename)),
    direction = factor(.data$direction, levels = c("Positive", "Negative", "Not supported"))
  )

thermal_order <- "Species thermal index"
trait_levels_with_enough_species <- trait_counts |>
  filter(.data$n_species >= min_species_per_trait_level | .data$trait_category == "Species thermal index") |>
  pull(.data$trait_clean)

significant_foraging_levels <- gamma_plot |>
  filter(
    .data$trait_category == "Foraging guild",
    .data$significant,
    .data$direction %in% c("Positive", "Negative")
  ) |>
  distinct(.data$trait_clean) |>
  pull(.data$trait_clean)

migration_reference_levels <- setdiff(
  migration_keep,
  gamma_plot |>
    filter(.data$trait_category == "Migration") |>
    distinct(.data$trait_clean) |>
    pull(.data$trait_clean)
)

foraging_reference_levels <- setdiff(
  foraging_keep,
  gamma_plot |>
    filter(.data$trait_category == "Foraging guild") |>
    distinct(.data$trait_clean) |>
    pull(.data$trait_clean)
)

migration_reference_label <- paste0(
  "Migratory strategy - relative to long-distance migrants (n = ",
  trait_counts$n_species[match(migration_reference_levels[[1]], trait_counts$trait_clean)],
  ")"
)

foraging_reference_label <- paste0(
  "Foraging guild - relative to ",
  format_trait_label(foraging_reference_levels[[1]]),
  " (n = ",
  trait_counts$n_species[match(foraging_reference_levels[[1]], trait_counts$trait_clean)],
  ")"
)

foraging_order_plotted <- foraging_order[foraging_order %in% significant_foraging_levels]

trait_level_order <- c(
  thermal_order,
  setdiff(migration_order, migration_reference_levels),
  foraging_order_plotted
)
trait_level_order <- trait_level_order[trait_level_order %in% trait_levels_with_enough_species]

trait_lookup <- trait_counts |>
  filter(.data$trait_clean %in% trait_level_order) |>
  distinct(.data$trait_clean, .data$trait_category, .data$n_species) |>
  mutate(
    trait_label = paste0(format_trait_label(.data$trait_clean), " (n = ", .data$n_species, ")")
  )

tile_values <- gamma_plot |>
  mutate(
    gamma_support = case_when(
      .data$significant & .data$direction == "Positive" ~ "Positive",
      .data$significant & .data$direction == "Negative" ~ "Negative",
      TRUE ~ "No supported effect"
    )
  ) |>
  select(
    atlas,
    trait_clean,
    variable,
    trait_category,
    n_species,
    effect_size,
    support_pos,
    gamma_support
  )

tile_df <- expand_grid(
  atlas = unname(atlas_rename),
  trait_clean = trait_level_order,
  variable = var_order
) |>
  left_join(tile_values, by = c("atlas", "trait_clean", "variable")) |>
  mutate(
    trait_category = coalesce(.data$trait_category, trait_lookup$trait_category[match(.data$trait_clean, trait_lookup$trait_clean)]),
    n_species = coalesce(.data$n_species, trait_lookup$n_species[match(.data$trait_clean, trait_lookup$trait_clean)]),
    trait_label = trait_lookup$trait_label[match(.data$trait_clean, trait_lookup$trait_clean)],
    gamma_support = replace_na(.data$gamma_support, "No supported effect"),
    gamma_support = factor(
      .data$gamma_support,
      levels = c("Negative", "No supported effect", "Positive")
    ),
    atlas = factor(.data$atlas, levels = unname(atlas_rename)),
    trait_clean = factor(.data$trait_clean, levels = rev(trait_level_order)),
    variable = factor(.data$variable, levels = var_order),
    trait_category = factor(
      .data$trait_category,
      levels = c("Species thermal index", "Migration", "Foraging guild")
    )
  )

write_csv(
  tile_df |>
    mutate(
      atlas = as.character(.data$atlas),
      trait_clean = as.character(.data$trait_clean),
      trait_label = as.character(.data$trait_label),
      variable = as.character(.data$variable),
      trait_category = as.character(.data$trait_category),
      gamma_support = as.character(.data$gamma_support)
    ),
  output_csv
)

panel_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 7.1, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6.8),
    axis.title = element_blank(),
    strip.text = element_text(size = 8.8, face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = "white"),
    legend.position = "bottom",
    legend.title = element_text(size = 8.8, face = "bold"),
    legend.text = element_text(size = 8.2),
    panel.spacing.x = grid::unit(0.42, "lines"),
    panel.spacing.y = grid::unit(0.65, "lines"),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin = margin(5.5, 5.5, 4.5, 5.5)
  )

make_trait_panel <- function(data, trait_group, title, trait_order, show_legend = FALSE) {
  trait_label_order <- data |>
    filter(.data$trait_category == trait_group) |>
    distinct(.data$trait_clean, .data$trait_label) |>
    arrange(match(as.character(.data$trait_clean), trait_order)) |>
    pull(.data$trait_label)

  panel_data <- data |>
    filter(.data$trait_category == trait_group) |>
    mutate(
      trait_label = factor(as.character(.data$trait_label), levels = rev(trait_label_order))
    )

  ggplot(
    panel_data,
    aes(x = .data$atlas, y = .data$trait_label, fill = .data$gamma_support)
  ) +
    geom_tile(colour = "grey88", linewidth = 0.22, width = 0.94, height = 0.94) +
    facet_wrap(vars(.data$variable), ncol = length(var_order)) +
    scale_fill_manual(
      values = gamma_colors,
      drop = FALSE,
      name = "Supported\nGamma Sign"
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = TRUE) +
    labs(title = title) +
    panel_theme +
    theme(legend.position = if (show_legend) "bottom" else "none")
}

thermal_panel <- make_trait_panel(
  tile_df,
  "Species thermal index",
  trait_category_labels[["Species thermal index"]],
  thermal_order
)

migration_panel <- make_trait_panel(
  tile_df,
  "Migration",
  migration_reference_label,
  setdiff(migration_order, migration_reference_levels)
)

foraging_panel <- make_trait_panel(
  tile_df,
  "Foraging guild",
  foraging_reference_label,
  foraging_order_plotted,
  show_legend = TRUE
)

foraging_height <- max(2.2, length(foraging_order_plotted) * 0.42)
figure_height <- 5.4 + foraging_height

gamma_atlas_support_plot <- thermal_panel / migration_panel / foraging_panel +
  plot_layout(heights = c(1, 2.5, foraging_height)) &
  theme(
    plot.caption = element_blank()
  )

ggsave(output_png, gamma_atlas_support_plot, width = 17.2, height = figure_height, units = "in", dpi = 300, bg = "white")

message("Wrote: ", output_png)
message("Wrote: ", output_csv)
