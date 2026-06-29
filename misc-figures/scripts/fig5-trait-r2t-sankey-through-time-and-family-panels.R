rm(list = ls())

# Alternate combined figure for trait-explained species-environment responses.
#
# Panel A uses a Sankey-bump display for the same full VP$R2T$Beta values shown
# in the original top panel. Ribbon widths are absolute R2T_Beta values divided
# by the number of environmental variables. Because Sankey-bump ribbons are
# stacked, the total stack height is the average full R2T_Beta across the eight
# environmental variables rather than a composition or raw sum.
#
# Panel B keeps the post-hoc drop-one trait-family contribution bars, with black
# ticks marking the full VP$R2T$Beta value for each environmental variable and
# atlas period.

library(tidyverse)
library(patchwork)
library(scales)

if (!requireNamespace("ggbump", quietly = TRUE) || !requireNamespace("ggsankey", quietly = TRUE)) {
  stop(
    "Install the `ggbump` and `ggsankey` packages before running this script. ",
    "They are used for the Sankey-bump R2T panel.",
    call. = FALSE
  )
}

library(ggbump)
library(ggsankey)

model_pattern <- "2026-03-13"

output_dir <- file.path("misc-figures", "outputs", "main")
full_r2t_csv <- file.path(output_dir, paste0(model_pattern, "-trait-environment-r2t-panel.csv"))
drop_one_csv <- file.path(output_dir, paste0(model_pattern, "-trait-family-drop-one-r2t.csv"))

figure_slug <- paste0(model_pattern, "-trait-r2t-sankey-through-time-and-family-panels")
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))

standalone_slug <- paste0(model_pattern, "-trait-environment-r2t-sankey-bump")
standalone_png <- file.path(output_dir, paste0(standalone_slug, ".png"))

var_order <- c(
  "Temperature",
  "Precipitation",
  "Land-use heterogeneity",
  "Urban",
  "Cropland",
  "Pasture",
  "Forest",
  "Grass/shrubland"
)

atlas_order <- c("1970s", "1990s", "2010s")

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

trait_family_order <- c(
  "Species thermal index",
  "Migratory strategy",
  "Foraging guild"
)

# Use colors already present in the Figure 2 palette, now assigned to trait
# groups rather than environmental variables.
trait_family_colors <- c(
  "Species thermal index" = "firebrick3",
  "Migratory strategy" = "dodgerblue3",
  "Foraging guild" = "springgreen4"
)

n_environment_variables <- length(var_order)

if (!file.exists(full_r2t_csv)) {
  stop("Missing full R2T_Beta CSV: ", full_r2t_csv, call. = FALSE)
}

if (!file.exists(drop_one_csv)) {
  stop("Missing drop-one trait-family CSV: ", drop_one_csv, call. = FALSE)
}

full_r2t <- read_csv(
  full_r2t_csv,
  col_types = cols(atlas = col_character())
) |>
  mutate(
    atlas = factor(.data$atlas, levels = atlas_order),
    variable = factor(.data$variable, levels = var_order),
    # Keep ribbons visible if a future R2T estimate is exactly zero, without
    # changing the underlying values used in labels or saved CSVs.
    plot_r2t_beta = pmax(.data$trait_r2t_beta / n_environment_variables, 0.0001)
  )

drop_one_summary <- read_csv(
  drop_one_csv,
  col_types = cols(atlas = col_character())
) |>
  mutate(
    atlas = factor(.data$atlas, levels = atlas_order),
    variable = factor(.data$variable, levels = var_order),
    trait_family = factor(.data$trait_family, levels = trait_family_order)
  )

bar_totals <- drop_one_summary |>
  group_by(.data$atlas, .data$variable) |>
  summarise(
    full_r2t_beta = first(.data$full_r2t_beta),
    stacked_drop_one_contribution = sum(.data$contribution_mean, na.rm = TRUE),
    .groups = "drop"
  )

expected_full_rows <- length(atlas_order) * length(var_order)
expected_drop_one_rows <- expected_full_rows * length(trait_family_order)

if (nrow(full_r2t) != expected_full_rows) {
  stop("Unexpected full R2T row count: ", nrow(full_r2t), call. = FALSE)
}

if (nrow(drop_one_summary) != expected_drop_one_rows) {
  stop("Unexpected drop-one row count: ", nrow(drop_one_summary), call. = FALSE)
}

max_full_value <- max(
  full_r2t$trait_r2t_beta,
  bar_totals$full_r2t_beta,
  bar_totals$stacked_drop_one_contribution,
  drop_one_summary$contribution_mean,
  na.rm = TRUE
)

y_limits <- c(0, max_full_value * 1.18)

base_panel_theme <- theme_minimal(base_size = 10) +
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

sankey_panel <- ggplot(
  full_r2t,
  aes(
    x = .data$atlas,
    node = .data$variable,
    group = .data$variable,
    fill = .data$variable,
    value = .data$plot_r2t_beta
  )
) +
  geom_sankey_bump(
    space = 0,
    type = "alluvial",
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = env_colors,
    breaks = var_order,
    drop = FALSE,
    name = "Environmental variable"
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, by = 0.05),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Average trait-explained species-environment variation",
    x = NULL,
    y = "Average absolute\nR2T_Beta"
  ) +
  base_panel_theme +
  theme(
    legend.position = "right",
    legend.key.size = grid::unit(0.36, "cm"),
    legend.text = element_text(size = 8.3),
    legend.title = element_text(size = 8.7, face = "bold")
  )

stacked_family_panel <- ggplot(
  drop_one_summary,
  aes(x = .data$atlas, y = .data$contribution_mean, fill = .data$trait_family)
) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey25") +
  geom_col(width = 0.66, colour = "white", linewidth = 0.18) +
  geom_point(
    data = bar_totals,
    aes(x = .data$atlas, y = .data$full_r2t_beta),
    inherit.aes = FALSE,
    shape = 95,
    size = 8.5,
    stroke = 1.25,
    colour = "grey8"
  ) +
  geom_text(
    data = bar_totals,
    aes(
      x = .data$atlas,
      y = pmax(.data$full_r2t_beta, .data$stacked_drop_one_contribution),
      label = percent(.data$full_r2t_beta, accuracy = 1)
    ),
    inherit.aes = FALSE,
    vjust = -0.45,
    size = 2.35,
    colour = "grey12"
  ) +
  facet_wrap(vars(.data$variable), ncol = 4) +
  scale_fill_manual(values = trait_family_colors, drop = FALSE, name = NULL) +
  scale_y_continuous(
    limits = y_limits,
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    title = "Post-hoc trait-family contributions",
    x = NULL,
    y = "Drop-one contribution\nto R2T_Beta"
  ) +
  base_panel_theme +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8.6),
    legend.key.width = grid::unit(1.2, "lines")
  )

combined_plot <- sankey_panel / stacked_family_panel +
  plot_layout(heights = c(0.82, 1.35)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 13),
    plot.tag.position = c(0.005, 0.99)
  )

ggsave(standalone_png, sankey_panel, width = 8.8, height = 3.6, units = "in", dpi = 300, bg = "white")
ggsave(output_png, combined_plot, width = 8.8, height = 10.4, units = "in", dpi = 300, bg = "white")

message("Wrote: ", standalone_png)
message("Wrote: ", output_png)
