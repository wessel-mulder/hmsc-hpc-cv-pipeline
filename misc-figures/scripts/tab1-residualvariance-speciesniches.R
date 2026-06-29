rm(list = ls())

# Residual species-niche variance-covariance from the HMSC V matrix.
#
# convertToCodaObject(model)$V gives posterior draws of the covariance matrix for
# species-specific responses to the environmental covariates. The diagonal terms
# are how much species vary in their response to each covariate. The off-diagonal
# terms are covariance in species responses to pairs of covariates; for plotting
# these are converted to correlations so every pair is on a comparable scale.

library(tidyverse)
library(Hmsc)
library(coda)
library(patchwork)
library(scales)

source(file.path("support_scripts", "figure_data_helpers.R"))

model_dir <- "HmscOutputs"
model_pattern <- "2026-03-13"

output_dir <- file.path("misc-figures", "outputs", "main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure_slug <- paste0(model_pattern, "-tab1-residualvariance-speciesniches")
diag_csv <- file.path(output_dir, paste0(figure_slug, "-diagonal-variance.csv"))
corr_csv <- file.path(output_dir, paste0(figure_slug, "-response-correlations.csv"))
diag_png <- file.path(output_dir, paste0(figure_slug, "-diagonal-variance.png"))
corr_png <- file.path(output_dir, paste0(figure_slug, "-response-correlations.png"))
combined_png <- file.path(output_dir, paste0(figure_slug, ".png"))

atlas_lookup <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

atlas_order <- unname(atlas_lookup)

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

v_name_pattern <- "^V\\[(.*) \\(C\\d+\\), (.*) \\(C\\d+\\)\\]$"

combine_v_chains <- function(v_chains, atlas_label) {
  map2_dfr(v_chains, seq_along(v_chains), function(chain, chain_id) {
    as.matrix(chain) |>
      as_tibble(.name_repair = "minimal") |>
      mutate(
        atlas = atlas_label,
        chain = chain_id,
        draw = row_number(),
        .before = 1
      )
  })
}

parse_v_draws <- function(model, atlas_num) {
  atlas_label <- atlas_lookup[[atlas_num]]
  message("Parsing V draws for atlas: ", atlas_label)

  codam <- convertToCodaObject(model)

  combine_v_chains(codam$V, atlas_label = atlas_label) |>
    pivot_longer(
      cols = starts_with("V["),
      names_to = "term",
      values_to = "covariance"
    ) |>
    extract(
      term,
      into = c("row_var_raw", "col_var_raw"),
      regex = v_name_pattern,
      remove = TRUE
    )
}

make_correlation_draw <- function(draw_df) {
  covariates <- unique(c(draw_df$row_var_raw, draw_df$col_var_raw))
  cov_mat <- matrix(
    NA_real_,
    nrow = length(covariates),
    ncol = length(covariates),
    dimnames = list(covariates, covariates)
  )

  cov_mat[cbind(draw_df$row_var_raw, draw_df$col_var_raw)] <- draw_df$covariance

  # The posterior object stores both triangles. Average tiny numerical
  # asymmetries away before converting covariance to correlation.
  cov_mat <- (cov_mat + t(cov_mat)) / 2
  corr_mat <- cov2cor(cov_mat)

  as.data.frame(as.table(corr_mat), stringsAsFactors = FALSE) |>
    as_tibble() |>
    rename(
      row_var_raw = Var1,
      col_var_raw = Var2,
      correlation = Freq
    )
}

model_folders <- figure_model_folders(pattern = model_pattern, base_dir = model_dir)
models <- load_hmsc_posteriors(model_folders, base_dir = model_dir)

v_long <- imap_dfr(models, parse_v_draws)

environment_v_long <- v_long |>
  filter(
    .data$row_var_raw %in% names(var_rename),
    .data$col_var_raw %in% names(var_rename)
  ) |>
  mutate(
    atlas = factor(.data$atlas, levels = atlas_order),
    row_var = var_rename[.data$row_var_raw],
    col_var = var_rename[.data$col_var_raw],
    row_var = factor(.data$row_var, levels = var_order),
    col_var = factor(.data$col_var, levels = var_order)
  )

diagonal_summary <- environment_v_long |>
  filter(.data$row_var_raw == .data$col_var_raw) |>
  group_by(.data$atlas, .data$row_var_raw, .data$row_var) |>
  summarise(
    mean_variance = mean(.data$covariance, na.rm = TRUE),
    median_variance = median(.data$covariance, na.rm = TRUE),
    lower_95 = quantile(.data$covariance, 0.025, na.rm = TRUE),
    upper_95 = quantile(.data$covariance, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(.data$row_var, .data$atlas)

correlation_draws <- environment_v_long |>
  group_by(.data$atlas, .data$chain, .data$draw) |>
  group_modify(~ make_correlation_draw(.x)) |>
  ungroup() |>
  filter(
    .data$row_var_raw %in% names(var_rename),
    .data$col_var_raw %in% names(var_rename)
  ) |>
  mutate(
    row_var = var_rename[.data$row_var_raw],
    col_var = var_rename[.data$col_var_raw],
    row_var = factor(.data$row_var, levels = rev(var_order)),
    col_var = factor(.data$col_var, levels = var_order)
  )

correlation_summary <- correlation_draws |>
  filter(.data$row_var_raw != .data$col_var_raw) |>
  group_by(.data$atlas, .data$row_var_raw, .data$col_var_raw, .data$row_var, .data$col_var) |>
  summarise(
    mean_correlation = mean(.data$correlation, na.rm = TRUE),
    median_correlation = median(.data$correlation, na.rm = TRUE),
    lower_95 = quantile(.data$correlation, 0.025, na.rm = TRUE),
    upper_95 = quantile(.data$correlation, 0.975, na.rm = TRUE),
    prob_positive = mean(.data$correlation > 0, na.rm = TRUE),
    supported = .data$prob_positive >= 0.95 | .data$prob_positive <= 0.05,
    .groups = "drop"
  ) |>
  arrange(.data$atlas, .data$row_var, .data$col_var)

write_csv(
  diagonal_summary |>
    mutate(
      atlas = as.character(.data$atlas),
      row_var = as.character(.data$row_var)
    ),
  diag_csv
)

write_csv(
  correlation_summary |>
    mutate(
      atlas = as.character(.data$atlas),
      row_var = as.character(.data$row_var),
      col_var = as.character(.data$col_var)
    ),
  corr_csv
)

variance_limit <- max(diagonal_summary$upper_95, na.rm = TRUE) * 1.06

diag_plot <- ggplot(
  diagonal_summary,
  aes(x = atlas, y = mean_variance, colour = row_var, group = row_var)
) +
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    width = 0.13,
    linewidth = 0.32,
    colour = "grey35"
  ) +
  geom_line(linewidth = 0.48) +
  geom_point(size = 1.7) +
  facet_wrap(vars(row_var), ncol = 4, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  guides(colour = "none") +
  labs(
    title = "Residual variance in species responses",
    x = NULL,
    y = "Posterior variance"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 7.8, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(size = 8.3, face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = "white"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.spacing = grid::unit(0.7, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

correlation_limit <- max(abs(correlation_summary$mean_correlation), na.rm = TRUE)

corr_plot <- ggplot(
  correlation_summary,
  aes(x = col_var, y = row_var, fill = mean_correlation)
) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = if_else(.data$supported, "*", "")),
    size = 3.2,
    colour = "grey10"
  ) +
  facet_wrap(vars(atlas), nrow = 1) +
  scale_fill_gradient2(
    low = "#d95f02",
    mid = "grey98",
    high = "#1b9e77",
    midpoint = 0,
    limits = c(-correlation_limit, correlation_limit),
    name = "Response\ncorrelation"
  ) +
  coord_equal() +
  labs(
    title = "Covariation among species responses to environmental gradients",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 7.4, angle = 35, hjust = 1),
    axis.text.y = element_text(size = 7.4),
    strip.text = element_text(size = 9.2, face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = "white"),
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

combined_plot <- diag_plot / corr_plot +
  plot_layout(heights = c(1.05, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.01, 0.99)
  )

ggsave(diag_png, diag_plot, width = 8.4, height = 4.8, units = "in", dpi = 300, bg = "white")
ggsave(corr_png, corr_plot, width = 10.2, height = 4.25, units = "in", dpi = 300, bg = "white")
ggsave(combined_png, combined_plot, width = 10.2, height = 9.2, units = "in", dpi = 300, bg = "white")

message("Wrote: ", diag_csv)
message("Wrote: ", corr_csv)
message("Wrote: ", diag_png)
message("Wrote: ", corr_png)
message("Wrote: ", combined_png)
