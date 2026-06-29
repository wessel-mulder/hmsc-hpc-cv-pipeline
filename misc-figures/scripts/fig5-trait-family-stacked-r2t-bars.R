rm(list = ls())

# Stacked trait-family version of the species-environment R2T_Beta barplots.
#
# This script is intentionally lightweight: it reads the post-hoc drop-one
# summary produced by fig5-trait-family-drop-one-r2t.R, then redraws the
# species-environment bars with the three trait families as stacked components.
# The black tick marks the actual full VP$R2T$Beta value, because drop-one
# contributions can overlap and therefore do not need to sum exactly to the
# fitted full trait-explained Beta variation.

library(tidyverse)
library(scales)

model_pattern <- "2026-03-13"

output_dir <- file.path("misc-figures", "outputs", "main")
input_csv <- file.path(output_dir, paste0(model_pattern, "-trait-family-drop-one-r2t.csv"))

figure_slug <- paste0(model_pattern, "-trait-family-stacked-r2t-bars")
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))
output_csv <- file.path(output_dir, paste0(figure_slug, ".csv"))

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

trait_family_order <- c(
  "Species thermal index",
  "Migratory strategy",
  "Foraging guild"
)

# Reuse colors from the figure 2 palette so this panel still belongs to the same
# visual family, while assigning them here to trait-family components.
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
    variable = factor(.data$variable, levels = var_order),
    trait_family = factor(.data$trait_family, levels = trait_family_order)
  )

expected_rows <- length(atlas_order) * length(var_order) * length(trait_family_order)
if (nrow(drop_one_summary) != expected_rows) {
  stop(
    "Unexpected input row count. Expected ",
    expected_rows,
    " rows but found ",
    nrow(drop_one_summary),
    ".",
    call. = FALSE
  )
}

bar_totals <- drop_one_summary |>
  group_by(.data$atlas, .data$variable) |>
  summarise(
    full_r2t_beta = first(.data$full_r2t_beta),
    stacked_drop_one_contribution = sum(.data$contribution_mean, na.rm = TRUE),
    shared_or_overlapping_contribution = .data$full_r2t_beta - .data$stacked_drop_one_contribution,
    .groups = "drop"
  )

write_csv(
  bar_totals |>
    mutate(
      atlas = as.character(.data$atlas),
      variable = as.character(.data$variable)
    ),
  output_csv
)

y_min <- min(0, drop_one_summary$contribution_mean, na.rm = TRUE)
y_max <- max(
  bar_totals$full_r2t_beta,
  bar_totals$stacked_drop_one_contribution,
  drop_one_summary$contribution_mean,
  na.rm = TRUE
)

stacked_r2t_plot <- ggplot(
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
    limits = c(y_min * 1.1, y_max * 1.18),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    title = "Trait-family contributions to species-environment response variation",
    x = NULL,
    y = "Drop-one contribution\nto R2T_Beta"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8.2, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(size = 8.6, face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = "white"),
    legend.position = "bottom",
    legend.text = element_text(size = 8.6),
    legend.key.width = grid::unit(1.2, "lines"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.spacing = grid::unit(0.75, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

ggsave(output_png, stacked_r2t_plot, width = 8.4, height = 5.25, units = "in", dpi = 300, bg = "white")

message("Wrote: ", output_png)
message("Wrote: ", output_csv)
