rm(list = ls())

# Combined figure for trait-explained species-environment responses.
#
# Panel A shows how the full VP$R2T$Beta value changes through atlas periods for
# each environmental variable. Panel B shows the same bars decomposed into the
# three post-hoc drop-one trait-family contributions. The black ticks in panel B
# mark the full VP$R2T$Beta values from panel A.

library(tidyverse)
library(patchwork)
library(scales)

model_pattern <- "2026-03-13"

output_dir <- file.path("misc-figures", "outputs", "main")
full_r2t_csv <- file.path(output_dir, paste0(model_pattern, "-trait-environment-r2t-panel.csv"))
drop_one_csv <- file.path(output_dir, paste0(model_pattern, "-trait-family-drop-one-r2t.csv"))

figure_slug <- paste0(model_pattern, "-trait-r2t-through-time-and-family-panels")
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))

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
    variable = factor(.data$variable, levels = var_order)
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

through_time_panel <- ggplot(
  full_r2t,
  aes(x = atlas, y = trait_r2t_beta, fill = variable)
) +
  geom_col(width = 0.66, colour = "white", linewidth = 0.18) +
  geom_line(aes(group = 1), colour = "grey20", linewidth = 0.3) +
  geom_point(size = 1.7, colour = "grey20") +
  geom_text(
    aes(label = percent(.data$trait_r2t_beta, accuracy = 1)),
    vjust = -0.45,
    size = 2.35,
    colour = "grey12"
  ) +
  facet_wrap(vars(variable), ncol = 4) +
  scale_fill_manual(values = env_colors, guide = "none") +
  scale_y_continuous(
    limits = y_limits,
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    title = "Full trait-explained species-environment variation",
    x = NULL,
    y = "R2T_Beta"
  ) +
  base_panel_theme

stacked_family_panel <- ggplot(
  drop_one_summary,
  aes(x = atlas, y = contribution_mean, fill = trait_family)
) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey25") +
  geom_col(width = 0.66, colour = "white", linewidth = 0.18) +
  geom_point(
    data = bar_totals,
    aes(x = atlas, y = full_r2t_beta),
    inherit.aes = FALSE,
    shape = 95,
    size = 8.5,
    stroke = 1.25,
    colour = "grey8"
  ) +
  geom_text(
    data = bar_totals,
    aes(
      x = atlas,
      y = pmax(full_r2t_beta, stacked_drop_one_contribution),
      label = percent(full_r2t_beta, accuracy = 1)
    ),
    inherit.aes = FALSE,
    vjust = -0.45,
    size = 2.35,
    colour = "grey12"
  ) +
  facet_wrap(vars(variable), ncol = 4) +
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

combined_plot <- through_time_panel / stacked_family_panel +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 13),
    plot.tag.position = c(0.005, 0.99)
  )

ggsave(output_png, combined_plot, width = 8.8, height = 10.3, units = "in", dpi = 300, bg = "white")

message("Wrote: ", output_png)
