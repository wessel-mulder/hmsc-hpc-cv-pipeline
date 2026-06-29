rm(list = ls())

# Plot the posterior distribution of the HMSC phylogenetic signal parameter.
#
# `rho` measures how strongly species-level residual/latent structure follows
# the supplied phylogeny. This script keeps the posterior draws so the figure
# shows uncertainty through time, rather than only plotting a point estimate.

library(tidyverse)
library(Hmsc)
library(coda)
library(scales)

source(file.path("support_scripts", "figure_data_helpers.R"))

model_pattern <- "2026-03-13"
model_dir <- "HmscOutputs"

output_dir <- file.path("misc-figures", "outputs", "main")
figure_slug <- paste0(model_pattern, "-phylogeny-rho-through-time")

output_png <- file.path(output_dir, paste0(figure_slug, ".png"))
output_draws_csv <- file.path(output_dir, paste0(figure_slug, "-posterior-draws.csv"))
output_summary_csv <- file.path(output_dir, paste0(figure_slug, ".csv"))
output_contrast_csv <- file.path(output_dir, paste0(figure_slug, "-contrasts.csv"))

atlas_rename <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

atlas_order <- unname(atlas_rename)

# A quiet, single-parameter palette. This avoids competing with the red/green
# effect-size palette used for signed trait moderation figures.
rho_fill <- "#7aa6c2"
rho_line <- "#1f4e6a"

rho_draws_from_model <- function(model, atlas_num) {
  if (is.null(model$phyloTree)) {
    stop("Atlas ", atlas_num, " does not contain a phyloTree.", call. = FALSE)
  }

  rho_mcmc <- convertToCodaObject(model)$Rho

  if (is.null(rho_mcmc) || length(rho_mcmc) == 0) {
    stop("Atlas ", atlas_num, " does not contain posterior Rho draws.", call. = FALSE)
  }

  imap_dfr(rho_mcmc, function(chain_draws, chain_id) {
    tibble(
      atlas = atlas_rename[[as.character(atlas_num)]],
      atlas_num = as.integer(atlas_num),
      chain = as.integer(chain_id),
      iteration = seq_along(as.numeric(chain_draws)),
      rho = as.numeric(chain_draws)
    )
  }) |>
    group_by(.data$atlas, .data$atlas_num) |>
    mutate(draw_id = row_number()) |>
    ungroup()
}

rho_summary <- function(draws) {
  draws |>
    group_by(.data$atlas, .data$atlas_num) |>
    summarise(
      mean_rho = mean(.data$rho, na.rm = TRUE),
      median_rho = median(.data$rho, na.rm = TRUE),
      lower_95 = quantile(.data$rho, 0.025, na.rm = TRUE),
      upper_95 = quantile(.data$rho, 0.975, na.rm = TRUE),
      prob_rho_positive = mean(.data$rho > 0, na.rm = TRUE),
      n_draws = n(),
      .groups = "drop"
    ) |>
    mutate(atlas = factor(.data$atlas, levels = atlas_order)) |>
    arrange(.data$atlas)
}

rho_contrasts <- function(draws) {
  wide_draws <- draws |>
    select(draw_id, atlas, rho) |>
    mutate(atlas = factor(.data$atlas, levels = atlas_order)) |>
    pivot_wider(names_from = atlas, values_from = rho)

  contrast_specs <- tribble(
    ~contrast, ~later, ~earlier,
    "1990s - 1970s", "1990s", "1970s",
    "2010s - 1990s", "2010s", "1990s",
    "2010s - 1970s", "2010s", "1970s"
  )

  pmap_dfr(contrast_specs, function(contrast, later, earlier) {
    difference <- wide_draws[[later]] - wide_draws[[earlier]]

    tibble(
      contrast = contrast,
      mean_difference = mean(difference, na.rm = TRUE),
      median_difference = median(difference, na.rm = TRUE),
      lower_95 = quantile(difference, 0.025, na.rm = TRUE),
      upper_95 = quantile(difference, 0.975, na.rm = TRUE),
      prob_increase = mean(difference > 0, na.rm = TRUE),
      n_draws = sum(!is.na(difference))
    )
  })
}

model_folders <- figure_model_folders(pattern = model_pattern, base_dir = model_dir)
models <- load_hmsc_posteriors(model_folders, base_dir = model_dir)

rho_draws <- imap_dfr(models, rho_draws_from_model) |>
  mutate(atlas = factor(.data$atlas, levels = atlas_order))

if (any(is.na(rho_draws$rho))) {
  stop("Rho posterior draws include NA values.", call. = FALSE)
}

rho_summary_df <- rho_summary(rho_draws)
rho_contrast_df <- rho_contrasts(rho_draws)

write_csv(
  rho_draws |>
    mutate(atlas = as.character(.data$atlas)),
  output_draws_csv
)

write_csv(
  rho_summary_df |>
    mutate(atlas = as.character(.data$atlas)),
  output_summary_csv
)

write_csv(rho_contrast_df, output_contrast_csv)

rho_plot <- ggplot(rho_draws, aes(x = atlas, y = rho)) +
  geom_violin(
    fill = rho_fill,
    colour = "grey20",
    linewidth = 0.35,
    alpha = 0.72,
    width = 0.72,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    colour = "grey18",
    linewidth = 0.35,
    alpha = 0.85
  ) +
  geom_line(
    data = rho_summary_df,
    aes(x = atlas, y = mean_rho, group = 1),
    inherit.aes = FALSE,
    colour = rho_line,
    linewidth = 0.65
  ) +
  geom_point(
    data = rho_summary_df,
    aes(x = atlas, y = mean_rho),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.7,
    stroke = 0.55,
    fill = "white",
    colour = rho_line
  ) +
  geom_text(
    data = rho_summary_df,
    aes(
      x = atlas,
      y = mean_rho,
      label = number(.data$mean_rho, accuracy = 0.01)
    ),
    inherit.aes = FALSE,
    vjust = -1.25,
    size = 3.1,
    colour = "grey12"
  ) +
  scale_y_continuous(
    breaks = seq(0.5, 1, by = 0.1),
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.01, 0.08))
  ) +
  coord_cartesian(ylim = c(0.45, 1)) +
  labs(
    title = "Phylogenetic signal through time",
    x = NULL,
    y = expression(rho)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 9, colour = "grey20"),
    axis.title.y = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.margin = margin(8, 12, 8, 10)
  )

ggsave(
  output_png,
  rho_plot,
  width = 4.8,
  height = 3.6,
  dpi = 320,
  bg = "white"
)

message("Wrote: ", output_summary_csv)
message("Wrote: ", output_draws_csv)
message("Wrote: ", output_contrast_csv)
message("Wrote: ", output_png)
