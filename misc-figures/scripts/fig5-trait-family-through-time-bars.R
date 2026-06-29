rm(list = ls())

# Quick through-time summary of trait-family influence.
#
# This aggregates the post-hoc drop-one contributions across the eight
# environmental variables, giving one atlas-level influence score per trait
# family. The plotted value is the mean contribution across variables; the CSV
# also keeps the summed contribution for reference.

library(tidyverse)
library(scales)

model_pattern <- "2026-03-13"

output_dir <- file.path("misc-figures", "outputs", "main")
input_csv <- file.path(output_dir, paste0(model_pattern, "-trait-family-drop-one-r2t.csv"))

figure_slug <- paste0(model_pattern, "-trait-family-through-time-bars")
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))
output_csv <- file.path(output_dir, paste0(figure_slug, ".csv"))

atlas_order <- c("1970s", "1990s", "2010s")

trait_family_order <- c(
  "Species thermal index",
  "Migratory strategy",
  "Foraging guild"
)

# Reuse colors from the Figure 2 palette so the trait-family panels remain in
# the same visual language as the rest of the figure set.
trait_family_colors <- c(
  "Species thermal index" = "firebrick3",
  "Migratory strategy" = "dodgerblue3",
  "Foraging guild" = "springgreen4"
)

if (!file.exists(input_csv)) {
  stop(
    "Missing drop-one summary CSV. Run misc-figures/scripts/fig5-trait-family-drop-one-r2t.R first: ",
    input_csv,
    call. = FALSE
  )
}

drop_one_summary <- read_csv(
  input_csv,
  col_types = cols(atlas = col_character())
) |>
  mutate(
    atlas = factor(.data$atlas, levels = atlas_order),
    trait_family = factor(.data$trait_family, levels = trait_family_order)
  )

trait_family_time <- drop_one_summary |>
  group_by(.data$atlas, .data$trait_family) |>
  summarise(
    mean_contribution = mean(.data$contribution_mean, na.rm = TRUE),
    sd_contribution = sd(.data$contribution_mean, na.rm = TRUE),
    summed_contribution = sum(.data$contribution_mean, na.rm = TRUE),
    n_variables = n(),
    .groups = "drop"
  ) |>
  arrange(.data$trait_family, .data$atlas)

expected_rows <- length(atlas_order) * length(trait_family_order)
if (nrow(trait_family_time) != expected_rows) {
  stop(
    "Unexpected trait-family summary row count. Expected ",
    expected_rows,
    " but found ",
    nrow(trait_family_time),
    ".",
    call. = FALSE
  )
}

write_csv(
  trait_family_time |>
    mutate(
      atlas = as.character(.data$atlas),
      trait_family = as.character(.data$trait_family)
    ),
  output_csv
)

y_limit <- max(
  trait_family_time$mean_contribution + trait_family_time$sd_contribution,
  na.rm = TRUE
) * 1.18
dodge <- position_dodge(width = 0.72)

trait_family_time_plot <- ggplot(
  trait_family_time,
  aes(x = atlas, y = mean_contribution, fill = trait_family, group = trait_family)
) +
  geom_col(position = dodge, width = 0.62, colour = "white", linewidth = 0.22) +
  geom_errorbar(
    aes(
      ymin = pmax(0, .data$mean_contribution - .data$sd_contribution),
      ymax = .data$mean_contribution + .data$sd_contribution
    ),
    position = dodge,
    width = 0.18,
    linewidth = 0.45,
    colour = "grey18"
  ) +
  geom_line(
    aes(colour = trait_family),
    position = dodge,
    linewidth = 0.45,
    alpha = 0.85
  ) +
  geom_point(
    aes(colour = trait_family),
    position = dodge,
    size = 1.9,
    alpha = 0.95
  ) +
  geom_text(
    aes(
      y = .data$mean_contribution + .data$sd_contribution,
      label = percent(.data$mean_contribution, accuracy = 0.1)
    ),
    position = dodge,
    vjust = -0.45,
    size = 3,
    colour = "grey12"
  ) +
  scale_fill_manual(values = trait_family_colors, name = NULL) +
  scale_colour_manual(values = trait_family_colors, guide = "none") +
  scale_y_continuous(
    limits = c(0, y_limit),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Trait-family influence through time",
    x = NULL,
    y = "Mean drop-one contribution\nacross environmental variables"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 8.8),
    legend.key.width = grid::unit(1.2, "lines"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

ggsave(output_png, trait_family_time_plot, width = 5.8, height = 3.55, units = "in", dpi = 300, bg = "white")

message("Wrote: ", output_png)
message("Wrote: ", output_csv)
