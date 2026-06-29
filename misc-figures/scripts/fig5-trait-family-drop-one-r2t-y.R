rm(list = ls())

# Post-hoc drop-one trait-family decomposition of VP$R2T$Y.
#
# This script does not refit the HMSC models. It reuses posterior draws from the
# fitted models and asks how much trait-explained occurrence variation is lost
# when each trait family is removed from the fitted trait design matrix.

library(tidyverse)
library(Hmsc)
library(scales)

source(file.path("support_scripts", "figure_data_helpers.R"))

model_dir <- "HmscOutputs"
model_pattern <- "2026-03-13"
posterior_start <- 1
match_tolerance <- 1e-10

output_dir <- file.path("misc-figures", "outputs", "main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure_slug <- paste0(model_pattern, "-trait-family-drop-one-r2t-y")
output_csv <- file.path(output_dir, paste0(figure_slug, ".csv"))
output_draws_csv <- file.path(output_dir, paste0(figure_slug, "-posterior-draws.csv"))
output_check_csv <- file.path(output_dir, paste0(figure_slug, "-checks.csv"))
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))

atlas_rename <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

trait_family_order <- c(
  "Species thermal index",
  "Migratory strategy",
  "Foraging guild"
)

trait_family_patterns <- c(
  "Species thermal index" = "^species_thermal_index$",
  "Migratory strategy" = "^Migration_a3_DOF",
  "Foraging guild" = "^foraging_guild_consensus"
)

# Reuse colors from the Figure 2 palette so the trait-family panels remain in
# the same visual language as the rest of the figure set.
trait_family_colors <- c(
  "Species thermal index" = "firebrick3",
  "Migratory strategy" = "dodgerblue3",
  "Foraging guild" = "springgreen4"
)

trait_family_columns <- function(trait_names) {
  family_columns <- imap(
    trait_family_patterns,
    ~ which(str_detect(trait_names, .x))
  )

  empty_families <- names(family_columns)[map_int(family_columns, length) == 0]
  if (length(empty_families) > 0) {
    stop(
      "No fitted trait-design columns found for: ",
      paste(empty_families, collapse = ", "),
      call. = FALSE
    )
  }

  family_columns
}

read_exported_r2t_y <- function(model_folder) {
  atlas_num <- str_extract(model_folder, "(?<=Atlas)\\d+")
  file_path <- file.path(
    model_dir,
    model_folder,
    "Results",
    sprintf("%sparameter_estimates_VP_R2T_Y.csv", model_folder)
  )

  if (!file.exists(file_path)) {
    stop("Missing VP R2T Y export: ", file_path, call. = FALSE)
  }

  read.csv(file_path, check.names = TRUE) |>
    as_tibble() |>
    rename(exported_full_r2t_y = 2) |>
    mutate(atlas = atlas_rename[atlas_num]) |>
    select(atlas, exported_full_r2t_y)
}

trait_predicted_beta <- function(trait_matrix, gamma_matrix) {
  t(trait_matrix %*% t(gamma_matrix))
}

fixed_prediction <- function(model, beta_matrix) {
  switch(
    class(model$X)[1L],
    matrix = {
      model$X %*% beta_matrix
    },
    list = {
      prediction <- matrix(NA_real_, nrow = model$ny, ncol = model$ns)
      for (species_id in seq_len(model$ns)) {
        prediction[, species_id] <- model$X[[species_id]] %*% beta_matrix[, species_id]
      }
      prediction
    },
    stop("Unsupported model$X class: ", class(model$X)[1L], call. = FALSE)
  )
}

center_species_within_sites <- function(prediction_matrix) {
  prediction_matrix - matrix(
    rep(rowMeans(prediction_matrix, na.rm = TRUE), ncol(prediction_matrix)),
    ncol = ncol(prediction_matrix)
  )
}

r2t_y_from_predictions <- function(trait_prediction, full_prediction, n_species) {
  trait_centered <- center_species_within_sites(trait_prediction)
  full_centered <- center_species_within_sites(full_prediction)

  covariance_by_site <- rowSums(trait_centered * full_centered, na.rm = TRUE) /
    (n_species - 1)
  trait_variance_by_site <- rowSums(trait_centered * trait_centered, na.rm = TRUE) /
    (n_species - 1)
  full_variance_by_site <- rowSums(full_centered * full_centered, na.rm = TRUE) /
    (n_species - 1)

  numerator <- sum(covariance_by_site^2, na.rm = TRUE)
  denominator <- sum(trait_variance_by_site * full_variance_by_site, na.rm = TRUE)

  if (denominator == 0) {
    return(NA_real_)
  }

  numerator / denominator
}

compute_drop_one_draws <- function(model, atlas_label) {
  trait_names <- colnames(model$Tr)
  trait_families <- trait_family_columns(trait_names)
  posterior_draws <- Hmsc:::poolMcmcChains(model$postList, start = posterior_start)

  map_dfr(seq_along(posterior_draws), function(draw_id) {
    draw <- posterior_draws[[draw_id]]

    full_prediction <- fixed_prediction(model, draw$Beta)
    full_trait_beta <- trait_predicted_beta(model$Tr, draw$Gamma)
    full_trait_prediction <- fixed_prediction(model, full_trait_beta)

    full_r2t_y <- r2t_y_from_predictions(
      trait_prediction = full_trait_prediction,
      full_prediction = full_prediction,
      n_species = model$ns
    )

    map_dfr(names(trait_families), function(trait_family) {
      keep_traits <- setdiff(seq_along(trait_names), trait_families[[trait_family]])
      drop_trait_beta <- trait_predicted_beta(
        trait_matrix = model$Tr[, keep_traits, drop = FALSE],
        gamma_matrix = draw$Gamma[, keep_traits, drop = FALSE]
      )
      drop_trait_prediction <- fixed_prediction(model, drop_trait_beta)
      drop_r2t_y <- r2t_y_from_predictions(
        trait_prediction = drop_trait_prediction,
        full_prediction = full_prediction,
        n_species = model$ns
      )

      tibble(
        atlas = atlas_label,
        draw = draw_id,
        trait_family = trait_family,
        full_r2t_y = full_r2t_y,
        drop_r2t_y = drop_r2t_y,
        contribution = full_r2t_y - drop_r2t_y
      )
    })
  })
}

model_folders <- figure_model_folders(pattern = model_pattern, base_dir = model_dir)
models <- load_hmsc_posteriors(model_folders, base_dir = model_dir)

drop_one_draws <- imap_dfr(models, function(model, atlas_num) {
  compute_drop_one_draws(
    model = model,
    atlas_label = atlas_rename[[atlas_num]]
  )
})

exported_r2t_y <- map_dfr(model_folders, read_exported_r2t_y)

drop_one_summary <- drop_one_draws |>
  group_by(.data$atlas, .data$trait_family) |>
  summarise(
    full_r2t_y = mean(.data$full_r2t_y, na.rm = TRUE),
    drop_r2t_y = mean(.data$drop_r2t_y, na.rm = TRUE),
    contribution_mean = mean(.data$contribution, na.rm = TRUE),
    contribution_median = median(.data$contribution, na.rm = TRUE),
    contribution_lower_95 = quantile(.data$contribution, 0.025, na.rm = TRUE),
    contribution_upper_95 = quantile(.data$contribution, 0.975, na.rm = TRUE),
    prob_positive = mean(.data$contribution > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(exported_r2t_y, by = "atlas") |>
  mutate(
    full_export_difference = .data$full_r2t_y - .data$exported_full_r2t_y,
    atlas = factor(.data$atlas, levels = unname(atlas_rename)),
    trait_family = factor(.data$trait_family, levels = trait_family_order)
  ) |>
  arrange(.data$trait_family, .data$atlas)

drop_one_checks <- drop_one_summary |>
  group_by(.data$atlas) |>
  summarise(
    n_trait_families = n(),
    max_abs_full_export_difference = max(abs(.data$full_export_difference), na.rm = TRUE),
    .groups = "drop"
  )

expected_rows <- length(atlas_rename) * length(trait_family_order)
if (nrow(drop_one_summary) != expected_rows) {
  stop(
    "Unexpected drop-one R2T_Y summary row count. Expected ",
    expected_rows,
    " but found ",
    nrow(drop_one_summary),
    ".",
    call. = FALSE
  )
}

if (any(drop_one_checks$n_trait_families != length(trait_family_order))) {
  stop("At least one atlas is missing trait-family rows.", call. = FALSE)
}

if (max(drop_one_checks$max_abs_full_export_difference, na.rm = TRUE) > match_tolerance) {
  write_csv(drop_one_checks, output_check_csv)
  stop(
    "Recomputed full R2T_Y does not match exported VP$R2T$Y within ",
    match_tolerance,
    ". See: ",
    output_check_csv,
    call. = FALSE
  )
}

drop_one_summary_output <- drop_one_summary |>
  transmute(
    atlas = as.character(.data$atlas),
    trait_family = as.character(.data$trait_family),
    full_r2t_y = .data$full_r2t_y,
    drop_r2t_y = .data$drop_r2t_y,
    contribution_mean = .data$contribution_mean,
    contribution_median = .data$contribution_median,
    contribution_lower_95 = .data$contribution_lower_95,
    contribution_upper_95 = .data$contribution_upper_95,
    prob_positive = .data$prob_positive
  )

write_csv(drop_one_summary_output, output_csv)
write_csv(drop_one_draws, output_draws_csv)
write_csv(drop_one_checks, output_check_csv)

drop_one_draws_plot <- drop_one_draws |>
  mutate(
    atlas = factor(.data$atlas, levels = unname(atlas_rename)),
    trait_family = factor(.data$trait_family, levels = trait_family_order)
  )

y_limit <- max(
  abs(drop_one_draws_plot$contribution),
  na.rm = TRUE
) * 1.12

drop_one_plot <- ggplot(
  drop_one_draws_plot,
  aes(x = trait_family, y = contribution, fill = trait_family)
) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.35) +
  geom_violin(
    colour = "grey20",
    linewidth = 0.35,
    width = 0.82,
    alpha = 0.72,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    colour = "grey18",
    fill = "white",
    linewidth = 0.35,
    alpha = 0.8
  ) +
  geom_point(
    data = drop_one_summary,
    aes(x = trait_family, y = contribution_mean),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.2,
    stroke = 0.45,
    fill = "white",
    colour = "grey12"
  ) +
  geom_text(
    data = drop_one_summary,
    aes(
      x = trait_family,
      y = contribution_mean,
      label = percent(.data$contribution_mean, accuracy = 0.1)
    ),
    inherit.aes = FALSE,
    vjust = -1.05,
    size = 2.8,
    colour = "grey12",
    show.legend = FALSE
  ) +
  facet_wrap(vars(atlas), nrow = 1) +
  scale_fill_manual(values = trait_family_colors, name = NULL) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.03, 0.1))
  ) +
  coord_cartesian(ylim = c(-y_limit, y_limit)) +
  labs(
    title = "Posterior trait-family contributions to trait-explained occurrence variation",
    x = NULL,
    y = "Drop-one contribution\nto R2T_Y"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = "white"),
    legend.position = "bottom",
    legend.text = element_text(size = 8.8),
    legend.key.width = grid::unit(1.2, "lines"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

ggsave(output_png, drop_one_plot, width = 7.4, height = 3.9, units = "in", dpi = 300, bg = "white")

message("Wrote: ", output_csv)
message("Wrote: ", output_draws_csv)
message("Wrote: ", output_check_csv)
message("Wrote: ", output_png)
