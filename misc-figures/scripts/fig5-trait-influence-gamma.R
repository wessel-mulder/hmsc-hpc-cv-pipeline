rm(list = ls())


# GETTING STARTED ---------------------------------------------------------


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
output_atlas_r2t_png <- file.path(output_dir, paste0(figure_slug, "-atlas-r2t.png"))
output_atlas_r2t_pdf <- file.path(output_dir, paste0(figure_slug, "-atlas-r2t.pdf"))

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

figure_caption <- paste0(
  "Top bars show VP$R2T$Beta, the fraction of species-environment responses ",
  "explained by traits for each environmental variable and atlas period. ",
  "Tiles show atlas-specific supported Gamma coefficients ",
  "(green = positive; red = negative).\nSupported Gamma coefficients have ",
  "Pr(x > 0) >= ", support_level, " or Pr(x > 0) <= ", 1 - support_level, "."
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

  read.csv(file_path, check.names = FALSE) |>
    as_tibble() |>
    rename(variable_raw = 1, r2t_beta = 2) |>
    filter(.data$variable_raw %in% names(var_rename)) |>
    mutate(
      atlas = atlas_rename[atlas_num],
      variable = var_rename[.data$variable_raw]
    ) |>
    select(.data$atlas, .data$variable, .data$r2t_beta)
}


# LOADING DATA  -----------------------------------------------------------
model_folders <- find_model_folders(base_dir = model_dir, pattern = model_pattern)
if (length(model_folders) == 0) {
  stop("No model folders found for pattern: ", model_pattern, call. = FALSE)
}

gamma_raw <- map_dfr(model_folders, read_gamma_effects)
r2t_beta <- map_dfr(model_folders, read_r2t_beta)

preprocessed_file <- "Data/preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData"
if (!file.exists(preprocessed_file)) {
  stop("Missing preprocessed trait data: ", preprocessed_file, call. = FALSE)
}

load(preprocessed_file)

migration_keep <- names(which(table(Tr$Migration_a3_DOF) >= min_species_per_trait_level))
foraging_keep <- names(which(table(Tr$foraging_guild_consensus) >= min_species_per_trait_level))
# CLEAN AND PROCESS  ------------------------------------------------------
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
      levels = c("Thermal index", "Migration", "Foraging guild")
    ),
    variable = factor(variable, levels = var_order),
    atlas = factor(atlas, levels = unname(atlas_rename)),
    direction = factor(direction, levels = c("Positive", "Negative", "Not supported"))
  )

# ---- Temporal sign check and atlas-specific signed tile plot ----
# Keep only coefficients with posterior support for a positive or negative
# Gamma effect. Unsupported coefficients become white zero-count tiles below.
supported_gamma <- gamma_plot |>
  filter(.data$significant, .data$direction %in% c("Positive", "Negative"))

# Atlas-specific panels can show temporal sign changes directly. Still write an
# audit table if a trait-variable pair has supported positive and negative Gamma
# coefficients in different atlas periods, because those are useful to inspect.
conflict_check <- supported_gamma |>
  group_by(.data$trait_category, .data$trait_clean, .data$variable) |>
  summarise(
    n_positive = sum(.data$direction == "Positive", na.rm = TRUE),
    n_negative = sum(.data$direction == "Negative", na.rm = TRUE),
    has_positive = n_positive > 0,
    has_negative = n_negative > 0,
    has_both = has_positive & has_negative,
    .groups = "drop"
  )

conflicting_effects <- conflict_check |>
  filter(.data$has_both)

if (nrow(conflicting_effects) > 0) {
  conflict_path <- file.path(output_dir, paste0(figure_slug, "-conflicting-effects.csv"))
  write_csv(conflicting_effects, conflict_path)
  message("Temporal Gamma sign-change audit written to: ", conflict_path)
}

thermal_order <- "Thermal index"
trait_level_order <- c(
  thermal_order,
  rev(migration_order),
  foraging_order
)
trait_level_order <- trait_level_order[trait_level_order %in% unique(gamma_plot$trait_clean)]

trait_lookup <- gamma_plot |>
  distinct(.data$trait_clean, .data$trait_category)

tile_values <- supported_gamma |>
  mutate(
    signed_gamma = case_when(
      .data$direction == "Positive" ~ 1,
      .data$direction == "Negative" ~ -1,
      TRUE ~ 0
    )
  ) |>
  group_by(.data$atlas, .data$trait_clean, .data$variable) |>
  summarise(signed_gamma = first(.data$signed_gamma), .groups = "drop")

tile_df <- expand_grid(
  atlas = unname(atlas_rename),
  trait_clean = trait_level_order,
  variable = var_order
) |>
  left_join(tile_values, by = c("atlas", "trait_clean", "variable")) |>
  left_join(trait_lookup, by = "trait_clean") |>
  mutate(
    signed_gamma = replace_na(.data$signed_gamma, 0),
    gamma_support = case_when(
      .data$signed_gamma > 0 ~ "Positive",
      .data$signed_gamma < 0 ~ "Negative",
      TRUE ~ "No supported effect"
    ),
    gamma_support = factor(
      .data$gamma_support,
      levels = c("Negative", "No supported effect", "Positive")
    ),
    trait_clean = factor(.data$trait_clean, levels = rev(trait_level_order)),
    variable = factor(.data$variable, levels = var_order),
    atlas = factor(.data$atlas, levels = unname(atlas_rename)),
    trait_category = factor(
      .data$trait_category,
      levels = c("Thermal index", "Migration", "Foraging guild")
    )
  )

r2t_beta <- r2t_beta |>
  mutate(
    atlas = factor(.data$atlas, levels = unname(atlas_rename)),
    variable = factor(.data$variable, levels = var_order)
  )

r2t_limit <- max(r2t_beta$r2t_beta, na.rm = TRUE) * 1.08
gamma_colors <- c(
  "Negative" = "#d95f02",
  "No supported effect" = "grey98",
  "Positive" = "#1b9e77"
)
trait_category_labels <- c(
  "Thermal index" = "Species thermal index",
  "Migration" = "Migratory strategy",
  "Foraging guild" = "Foraging guild"
)

driver_panel_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", colour = NA)
  )

r2t_strip <- ggplot(r2t_beta, aes(x = variable, y = r2t_beta)) +
  geom_col(width = 0.82, fill = "grey40", colour = "white", linewidth = 0.18) +
  facet_grid(. ~ atlas) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    limits = c(0, r2t_limit),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(x = NULL, y = "Trait-explained\nBeta variation") +
  driver_panel_theme +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7.5),
    axis.title.y = element_text(size = 8.2, face = "bold"),
    strip.text.x = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = "white"),
    plot.margin = margin(5.5, 5.5, 1.5, 5.5)
  )

gamma_heatmap <- ggplot(tile_df, aes(x = variable, y = trait_clean, fill = gamma_support)) +
  geom_tile(colour = "grey90", linewidth = 0.25) +
  facet_grid(
    trait_category ~ atlas,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(trait_category = as_labeller(trait_category_labels))
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(
    values = gamma_colors,
    drop = FALSE,
    name = "Supported\nGamma sign"
  ) +
  labs(x = NULL, y = NULL) +
  driver_panel_theme +
  theme(
    axis.text.x = element_text(size = 8.2, angle = 35, hjust = 1),
    axis.text.y = element_text(size = 8.8),
    legend.position = "right",
    legend.title = element_text(size = 8.8, face = "bold"),
    legend.text = element_text(size = 8),
    strip.text.x = element_blank(),
    strip.background.x = element_blank(),
    strip.text.y = element_text(size = 9, face = "bold"),
    strip.background.y = element_rect(fill = "grey92", colour = "white"),
    panel.spacing.x = grid::unit(0.9, "lines"),
    panel.spacing.y = grid::unit(0.65, "lines"),
    plot.margin = margin(1.5, 5.5, 5.5, 5.5)
  )

p_gamma <- r2t_strip / gamma_heatmap +
  plot_layout(heights = c(1.0, 5.4), guides = "collect") +
  plot_annotation(caption = figure_caption) &
  theme(
    legend.position = "right",
    plot.caption = element_text(size = 8.2, colour = "grey25", hjust = 0),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

ggsave(output_png, p_gamma, width = 14.6, height = 8.2, units = "in", dpi = 300, bg = "white")
ggsave(output_pdf, p_gamma, width = 14.6, height = 8.2, units = "in", bg = "white")
ggsave(output_atlas_r2t_png, p_gamma, width = 14.6, height = 8.2, units = "in", dpi = 300, bg = "white")
ggsave(output_atlas_r2t_pdf, p_gamma, width = 14.6, height = 8.2, units = "in", bg = "white")

# Compatibility alias for the old draft filename. The publication filename above
# is the canonical output.
ggsave(file.path(output_dir, "gamma_circles.png"), p_gamma, width = 14.6, height = 8.2, units = "in", dpi = 300, bg = "white")
ggsave(file.path(output_dir, "gamma_circles.pdf"), p_gamma, width = 14.6, height = 8.2, units = "in", bg = "white")

message("Wrote: ", output_png)
message("Wrote: ", output_pdf)
message("Wrote: ", output_atlas_r2t_png)
message("Wrote: ", output_atlas_r2t_pdf)
