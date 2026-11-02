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

source(file.path("scripts", "support_scripts", "project_paths.R"))

model_dir <- file.path("intermediate-outputs", "hmsc-output")
model_pattern <- "2026-07-27"
min_species_per_trait_level <- 3

args <- commandArgs(trailingOnly = TRUE)
# Approved layouts: default is compact with outer labels only.
# Use --repeat-labels for environmental headings and years in every group.
compact_rows <- TRUE
spaced_rows <- FALSE
left_group_headings <- TRUE
contained_headings <- TRUE
repeat_labels <- TRUE
aligned_headings <- TRUE
square_tiles <- TRUE
outer_labels_only <- !("--repeat-labels" %in% args)
tight_groups <- outer_labels_only
row_height_in <- 0.38
figure_width <- 17.2
text_size <- 12
support_level <- 0.95
support_level_arg <- args[grepl("^--support-level=", args)]
if (length(support_level_arg) > 1) {
  stop("Supply at most one --support-level=<number> argument.", call. = FALSE)
}
if (length(support_level_arg) == 1) {
  support_level <- as.numeric(sub("^--support-level=", "", support_level_arg))
}
if (!is.finite(support_level) || support_level <= 0.5 || support_level >= 1) {
  stop("--support-level must be a number greater than 0.5 and less than 1.", call. = FALSE)
}

output_dir <- file.path("manuscript-outputs", "figures", "main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

support_suffix <- if (isTRUE(all.equal(support_level, 0.95))) {
  ""
} else {
  paste0("-p", formatC(support_level * 100, format = "f", digits = 0))
}

figure_slug <- paste0(model_pattern, "-fig2-trait-influence-gamma-atlas-support", support_suffix)
if (!outer_labels_only) figure_slug <- paste0(figure_slug, "-repeated-labels")
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))

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
    "short-and long-distance" = "Short- and long-distance",
    "short-distance" = "Short-distance",
    "sedentary and short-distance" = "Sedentary and short-distance",
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

preprocessed_file <- file.path("data-new", "derived", "hmsc-inputs", "preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData")
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
  str_to_lower(format_trait_label(foraging_reference_levels[[1]])),
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


section_title_lookup <- c(
  "Species thermal index" = trait_category_labels[["Species thermal index"]],
  "Migration" = migration_reference_label,
  "Foraging guild" = foraging_reference_label
)

trait_labels_for_section <- function(data, trait_group, trait_order) {
  data |>
    filter(.data$trait_category == trait_group) |>
    distinct(.data$trait_clean, .data$trait_label) |>
    arrange(match(as.character(.data$trait_clean), trait_order)) |>
    pull(.data$trait_label)
}

panel_trait_rows <- function(data, trait_group) {
  data |>
    filter(.data$trait_category == trait_group) |>
    distinct(.data$trait_clean) |>
    nrow()
}

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
    legend.direction = "horizontal",
    legend.title = element_text(size = 8.8, face = "bold"),
    legend.text = element_text(size = 8.2),
    legend.box = "horizontal",
    panel.spacing.x = grid::unit(0.42, "lines"),
    panel.spacing.y = grid::unit(0.65, "lines"),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin = margin(5.5, 5.5, 4.5, 5.5)
  )

if (!is.null(text_size)) {
  panel_theme <- panel_theme + theme(
    axis.text.x = element_text(size = text_size, angle = 45, hjust = 1),
    axis.text.y = element_text(size = text_size),
    strip.text = element_text(size = text_size, face = "bold"),
    legend.title = element_text(size = text_size, face = "bold"),
    legend.text = element_text(size = text_size),
    plot.title = element_text(size = text_size, face = "bold", hjust = 0)
  )
}

variable_labeller <- function(labels) {
  if (square_tiles) labels[labels == "Grass/Shrubland"] <- "Grass/\nShrubland"
  if (compact_rows) {
    labels[labels == "Land-Use Heterogeneity"] <- "Land use\nheterogeneity"
  } else if (!is.null(text_size)) {
    labels[labels == "Land-Use Heterogeneity"] <- "Land-Use\nHeterogeneity"
  }
  labels
}

make_trait_panel <- function(data, trait_group, title, trait_order, show_legend = FALSE) {
  if (left_group_headings && trait_group == "Species thermal index") title <- "Thermal affinity"
  if (contained_headings) title <- NULL
  trait_label_order <- trait_labels_for_section(data, trait_group, trait_order)

  panel_data <- data |>
    filter(.data$trait_category == trait_group) |>
    mutate(
      trait_label = factor(as.character(.data$trait_label), levels = rev(trait_label_order))
    )

  ggplot(
    panel_data,
    aes(x = .data$atlas, y = .data$trait_label, fill = .data$gamma_support)
  ) +
    geom_tile(colour = "grey88", linewidth = 0.22, width = 0.94,
              height = if (spaced_rows) 0.94 * 0.25 / row_height_in else 0.94) +
    facet_wrap(vars(.data$variable), ncol = length(var_order),
               labeller = as_labeller(variable_labeller)) +
    scale_fill_manual(
      values = gamma_colors,
      drop = FALSE,
      name = if (compact_rows) "Trait - environment association" else "Gamma sign",
      labels = if (compact_rows) c("Negative", "No association", "Positive") else waiver()
    ) +
    scale_x_discrete(drop = FALSE, expand = expansion(add = 0.6)) +
    scale_y_discrete(drop = TRUE,
                     expand = if (compact_rows) expansion(add = 0.5) else waiver()) +
    labs(title = title, x = NULL, y = NULL) +
    panel_theme +
    theme(legend.position = if (show_legend) "bottom" else "none") +
    theme(
      plot.title.position = if (left_group_headings) "plot" else "panel",
      plot.margin = if (left_group_headings && trait_group != "Species thermal index") margin(19, 5.5, 4.5, 5.5) else margin(5.5, 5.5, 4.5, 5.5),
      axis.text.x = if (compact_rows && !repeat_labels && trait_group != "Foraging guild") element_blank() else panel_theme$axis.text.x,
      strip.text = if (compact_rows && !repeat_labels && trait_group != "Species thermal index") element_blank() else panel_theme$strip.text,
      strip.background = if (compact_rows && !repeat_labels && trait_group != "Species thermal index") element_blank() else panel_theme$strip.background
    ) +
    guides(
      fill = guide_legend(
        title.position = "left",
        title.hjust = 1,
        nrow = 1,
        override.aes = list(colour = "grey88", linewidth = 0.22)
      )
    )
}

thermal_panel <- make_trait_panel(
  tile_df,
  "Species thermal index",
  section_title_lookup[["Species thermal index"]],
  thermal_order
)

migration_panel <- make_trait_panel(
  tile_df,
  "Migration",
  section_title_lookup[["Migration"]],
  setdiff(migration_order, migration_reference_levels)
)

foraging_panel <- make_trait_panel(
  tile_df,
  "Foraging guild",
  section_title_lookup[["Foraging guild"]],
  foraging_order_plotted
)

legend_grob <- cowplot::get_legend(
  make_trait_panel(
    tile_df,
    "Foraging guild",
    section_title_lookup[["Foraging guild"]],
    foraging_order_plotted,
    show_legend = TRUE
  )
)

panel_grobs <- lapply(
  list(thermal_panel, migration_panel, foraging_panel),
  ggplotGrob
)

max_panel_widths <- do.call(
  grid::unit.pmax,
  lapply(panel_grobs, function(grob) grob$widths)
)

panel_grobs <- lapply(panel_grobs, function(grob) {
  grob$widths <- max_panel_widths
  grob
})

panel_rows <- c(
  panel_trait_rows(tile_df, "Species thermal index"),
  panel_trait_rows(tile_df, "Migration"),
  panel_trait_rows(tile_df, "Foraging guild")
)
panel_overhead <- 1.7
legend_height <- 0.65
figure_height <- max(8.2, 0.58 * sum(panel_rows) + 3.4)

if (compact_rows) {
  # Fix the physical row pitch, rather than allocating height to whole plots:
  # their titles, strips and axes have different amounts of overhead.
  panel_grobs <- Map(function(grob, n_rows) {
    panel_indices <- unique(grob$layout$t[grepl("^panel", grob$layout$name)])
    grob$heights[panel_indices] <- grid::unit(n_rows * row_height_in, "in")
    grob
  }, panel_grobs, panel_rows)
  if (contained_headings) {
    heading_size <- if (is.null(text_size)) 12 else text_size
    headings <- c("Thermal affinity", "Migratory strategy", "Foraging guild")
    references <- c("", sub("^Migratory strategy - ", "", migration_reference_label),
                    sub("^Foraging guild - ", "", foraging_reference_label))
    reference_width <- max(vapply(references[-1], function(label) {
      grid::convertWidth(grid::grobWidth(grid::textGrob(label,
        gp = grid::gpar(fontsize = heading_size, fontface = "italic"))), "in", valueOnly = TRUE)
    }, numeric(1))) + 0.15
    panel_grobs <- lapply(seq_along(panel_grobs), function(i) {
      grob <- panel_grobs[[i]]
      label_col <- min(grob$layout$l[grepl("^axis-l", grob$layout$name)])
      grob$widths[label_col] <- grid::unit(reference_width, "in")
      panel_row <- min(grob$layout$t[grepl("^panel", grob$layout$name)])
      if (aligned_headings) {
        panel_row <- min(grob$layout$t[grepl("^strip-t", grob$layout$name)])
        grob$heights[panel_row] <- grid::unit(0.48, "in")
      } else {
        grob <- gtable::gtable_add_rows(grob, grid::unit(if (i == 1) 0.25 else 0.48, "in"), pos = panel_row - 1)
      }
      heading <- grid::textGrob(headings[i], x = grid::unit(1, "npc") - grid::unit(0.06, "in"),
        y = if (aligned_headings && i == 1) grid::unit(0.5, "npc") else grid::unit(1, "npc") - grid::unit(0.02, "in"),
        just = c("right", if (aligned_headings && i == 1) "centre" else "top"),
        gp = grid::gpar(fontsize = heading_size, fontface = "bold"))
      if (i > 1) {
        heading <- grid::grobTree(heading, grid::textGrob(references[i],
          x = grid::unit(1, "npc") - grid::unit(0.06, "in"),
          y = grid::unit(1, "npc") - grid::unit(0.25, "in"), just = c("right", "top"),
          gp = grid::gpar(fontsize = heading_size, fontface = "italic")))
      }
      gtable::gtable_add_grob(grob, heading, t = panel_row, l = label_col, clip = "on", name = "contained-heading")
    })
  }
  if (outer_labels_only) {
    # Preserve the approved geometry and spacing; hide only repeated labels.
    panel_grobs <- lapply(seq_along(panel_grobs), function(i) {
      grob <- panel_grobs[[i]]
      hide <- (i > 1 & grepl("^strip-t", grob$layout$name)) |
        (i < length(panel_grobs) & grepl("^axis-b", grob$layout$name))
      grob$grobs[hide] <- lapply(which(hide), function(j) grid::nullGrob())
      if (tight_groups) {
        if (i < length(panel_grobs)) {
          axis_rows <- unique(grob$layout$t[grepl("^axis-b", grob$layout$name)])
          grob$heights[axis_rows] <- grid::unit(0, "in")
        }
        # Keep only a small gap and the two-line heading above each group.
        grob$heights[1] <- grid::unit(if (i > 1) 0.06 else 0.04, "in")
        grob$heights[length(grob$heights)] <- grid::unit(0.02, "in")
      }
      grob
    })
  }
  panel_heights_in <- vapply(panel_grobs, function(grob) {
    grid::convertHeight(sum(grob$heights), "in", valueOnly = TRUE)
  }, numeric(1))
  if (square_tiles) {
    # Three x levels plus 0.6 expansion at each end span 3.2 units.
    # Matching their physical scale to y makes each 0.94-unit tile square.
    panel_grobs <- lapply(panel_grobs, function(grob) {
      panel_cols <- unique(grob$layout$l[grepl("^panel", grob$layout$name)])
      grob$widths[panel_cols] <- grid::unit(3.2 * row_height_in, "in")
      grob
    })
    figure_width <- grid::convertWidth(sum(panel_grobs[[1]]$widths), "in", valueOnly = TRUE)
  }
  legend_height_in <- grid::convertHeight(sum(legend_grob$heights), "in", valueOnly = TRUE) + 0.12
  figure_height <- sum(panel_heights_in) + legend_height_in
  gamma_atlas_support_plot <- gridExtra::arrangeGrob(
    grobs = c(panel_grobs, list(legend_grob)), ncol = 1,
    heights = grid::unit(c(panel_heights_in, legend_height_in), "in")
  )
} else {
gamma_atlas_support_plot <- gridExtra::arrangeGrob(
  grobs = c(panel_grobs, list(legend_grob)),
  ncol = 1,
  heights = c(panel_rows + panel_overhead, legend_height)
)
}

ggsave(output_png, gamma_atlas_support_plot, width = figure_width, height = figure_height, units = "in", dpi = 300, bg = "white")

message("Wrote: ", output_png)
