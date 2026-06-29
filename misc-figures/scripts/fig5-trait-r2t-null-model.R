rm(list = ls())

# Trait-profile shuffle null model for VP$R2T$Beta.
#
# This is a post-hoc null model: the fitted HMSC models and posterior draws are
# held fixed, and only the species-to-trait-profile matching is randomized. Each
# null iteration permutes whole rows of the fitted trait design matrix within an
# atlas, preserving the trait distributions and correlations among trait columns.

library(tidyverse)
library(Hmsc)
library(scales)

source(file.path("support_scripts", "figure_data_helpers.R"))

model_dir <- "HmscOutputs"
model_pattern <- "2026-03-13"
posterior_start <- 1
n_permutations <- 999
null_seed <- 20260521
match_tolerance <- 1e-10

output_dir <- file.path("misc-figures", "outputs", "main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure_slug <- paste0(model_pattern, "-trait-r2t-profile-shuffle-null")
output_csv <- file.path(output_dir, paste0(figure_slug, ".csv"))
output_permutation_csv <- file.path(output_dir, paste0(figure_slug, "-permutations.csv"))
output_check_csv <- file.path(output_dir, paste0(figure_slug, "-checks.csv"))
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))

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

center_rows <- function(x) {
  sweep(x, 1, rowMeans(x, na.rm = TRUE), check.margin = FALSE)
}

cor2_from_centered_rows <- function(beta_centered, pred_centered, denom, perm = NULL) {
  if (!is.null(perm)) {
    pred_centered <- pred_centered[, perm, drop = FALSE]
  }

  corr <- rowSums(beta_centered * pred_centered, na.rm = TRUE) / denom
  corr[denom == 0] <- NA_real_
  corr^2
}

read_exported_r2t_beta <- function(model_folder) {
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

  read.csv(file_path, check.names = TRUE) |>
    as_tibble() |>
    rename(variable_raw = 1, exported_r2t_beta = 2) |>
    filter(.data$variable_raw %in% names(var_rename)) |>
    mutate(
      atlas = atlas_rename[atlas_num],
      variable = var_rename[.data$variable_raw]
    ) |>
    select(atlas, variable_raw, variable, exported_r2t_beta)
}

trait_distributions_preserved <- function(trait_matrix, perm) {
  if (!identical(sort(perm), seq_len(nrow(trait_matrix)))) {
    return(FALSE)
  }

  all(vapply(
    seq_len(ncol(trait_matrix)),
    function(column_id) {
      isTRUE(all.equal(
        sort(trait_matrix[, column_id]),
        sort(trait_matrix[perm, column_id]),
        check.attributes = FALSE
      ))
    },
    logical(1)
  ))
}

prepare_variable_matrices <- function(model, posterior_draws, variable_raw, variable_index) {
  n_species <- nrow(model$Tr)

  map2(variable_raw, variable_index, function(raw_name, beta_row) {
    beta_matrix <- t(vapply(
      posterior_draws,
      function(draw) draw$Beta[beta_row, ],
      numeric(n_species)
    ))

    # For a whole-row trait-profile shuffle, Tr_perm %*% Gamma is equivalent to
    # permuting the species order of this fitted trait prediction. Precomputing
    # this once keeps the 999-permutation null computationally cheap.
    pred_matrix <- t(vapply(
      posterior_draws,
      function(draw) (model$Tr %*% t(draw$Gamma))[, beta_row],
      numeric(n_species)
    ))

    beta_centered <- center_rows(beta_matrix)
    pred_centered <- center_rows(pred_matrix)

    list(
      variable_raw = raw_name,
      variable = var_rename[[raw_name]],
      beta_centered = beta_centered,
      pred_centered = pred_centered,
      denom = sqrt(rowSums(beta_centered^2) * rowSums(pred_centered^2))
    )
  })
}

compute_atlas_null <- function(model, atlas_label) {
  posterior_draws <- Hmsc:::poolMcmcChains(model$postList, start = posterior_start)
  variable_raw <- intersect(names(var_rename), model$covNames)
  variable_index <- match(variable_raw, model$covNames)
  variable_matrices <- prepare_variable_matrices(
    model = model,
    posterior_draws = posterior_draws,
    variable_raw = variable_raw,
    variable_index = variable_index
  )

  observed <- map_dfr(variable_matrices, function(variable_data) {
    observed_draws <- cor2_from_centered_rows(
      beta_centered = variable_data$beta_centered,
      pred_centered = variable_data$pred_centered,
      denom = variable_data$denom
    )

    tibble(
      atlas = atlas_label,
      variable_raw = variable_data$variable_raw,
      variable = variable_data$variable,
      observed_r2t_beta = mean(observed_draws, na.rm = TRUE)
    )
  })

  preservation_checks <- logical(n_permutations)

  null_permutations <- map_dfr(seq_len(n_permutations), function(permutation_id) {
    perm <- sample.int(nrow(model$Tr))
    preservation_checks[[permutation_id]] <<- trait_distributions_preserved(model$Tr, perm)

    map_dfr(variable_matrices, function(variable_data) {
      null_draws <- cor2_from_centered_rows(
        beta_centered = variable_data$beta_centered,
        pred_centered = variable_data$pred_centered,
        denom = variable_data$denom,
        perm = perm
      )

      tibble(
        atlas = atlas_label,
        permutation = permutation_id,
        variable_raw = variable_data$variable_raw,
        variable = variable_data$variable,
        null_r2t_beta = mean(null_draws, na.rm = TRUE)
      )
    })
  })

  list(
    observed = observed,
    null_permutations = null_permutations,
    atlas_check = tibble(
      atlas = atlas_label,
      atlas_n_permutations = n_permutations,
      trait_distributions_preserved = all(preservation_checks)
    )
  )
}

set.seed(null_seed)

model_folders <- figure_model_folders(pattern = model_pattern, base_dir = model_dir)
models <- load_hmsc_posteriors(model_folders, base_dir = model_dir)

atlas_results <- imap(models, function(model, atlas_num) {
  message("Computing trait-shuffle null for atlas: ", atlas_rename[[atlas_num]])
  compute_atlas_null(model = model, atlas_label = atlas_rename[[atlas_num]])
})

observed_r2t <- map_dfr(atlas_results, "observed")
null_permutations <- map_dfr(atlas_results, "null_permutations")
atlas_checks <- map_dfr(atlas_results, "atlas_check")
exported_r2t <- map_dfr(model_folders, read_exported_r2t_beta)

p_values <- null_permutations |>
  left_join(observed_r2t, by = c("atlas", "variable_raw", "variable")) |>
  group_by(.data$atlas, .data$variable_raw, .data$variable) |>
  summarise(
    empirical_p_upper = (
      1 + sum(.data$null_r2t_beta >= .data$observed_r2t_beta, na.rm = TRUE)
    ) / (1 + n_permutations),
    .groups = "drop"
  )

null_summary <- null_permutations |>
  group_by(.data$atlas, .data$variable_raw, .data$variable) |>
  summarise(
    null_mean = mean(.data$null_r2t_beta, na.rm = TRUE),
    null_median = median(.data$null_r2t_beta, na.rm = TRUE),
    null_lower_95 = quantile(.data$null_r2t_beta, 0.025, na.rm = TRUE),
    null_upper_95 = quantile(.data$null_r2t_beta, 0.975, na.rm = TRUE),
    null_sd = sd(.data$null_r2t_beta, na.rm = TRUE),
    n_permutations = n(),
    .groups = "drop"
  ) |>
  left_join(observed_r2t, by = c("atlas", "variable_raw", "variable")) |>
  left_join(p_values, by = c("atlas", "variable_raw", "variable")) |>
  mutate(
    observed_minus_null_mean = .data$observed_r2t_beta - .data$null_mean,
    standardized_effect = .data$observed_minus_null_mean / .data$null_sd,
    atlas = factor(.data$atlas, levels = unname(atlas_rename)),
    variable = factor(.data$variable, levels = var_order)
  ) |>
  arrange(.data$variable, .data$atlas)

checks <- null_summary |>
  left_join(exported_r2t, by = c("atlas", "variable_raw", "variable")) |>
  mutate(
    observed_export_difference = .data$observed_r2t_beta - .data$exported_r2t_beta,
    n_permutations_expected = n_permutations
  ) |>
  select(
    atlas,
    variable_raw,
    variable,
    n_permutations,
    n_permutations_expected,
    observed_r2t_beta,
    exported_r2t_beta,
    observed_export_difference
  ) |>
  left_join(atlas_checks, by = "atlas")

expected_rows <- length(atlas_rename) * length(var_order)
if (nrow(null_summary) != expected_rows) {
  stop(
    "Unexpected null summary row count. Expected ",
    expected_rows,
    " but found ",
    nrow(null_summary),
    ".",
    call. = FALSE
  )
}

if (any(checks$n_permutations != n_permutations)) {
  stop("At least one atlas-variable combination has the wrong number of null permutations.", call. = FALSE)
}

if (!all(checks$trait_distributions_preserved)) {
  stop("At least one null permutation failed the trait-distribution preservation check.", call. = FALSE)
}

if (max(abs(checks$observed_export_difference), na.rm = TRUE) > match_tolerance) {
  stop(
    "Recomputed observed R2T_Beta does not match exported VP$R2T$Beta within ",
    match_tolerance,
    ". See: ",
    output_check_csv,
    call. = FALSE
  )
}

summary_output <- null_summary |>
  transmute(
    atlas = as.character(.data$atlas),
    variable = as.character(.data$variable),
    observed_r2t_beta = .data$observed_r2t_beta,
    null_mean = .data$null_mean,
    null_median = .data$null_median,
    null_lower_95 = .data$null_lower_95,
    null_upper_95 = .data$null_upper_95,
    empirical_p_upper = .data$empirical_p_upper,
    observed_minus_null_mean = .data$observed_minus_null_mean,
    standardized_effect = .data$standardized_effect
  )

write_csv(summary_output, output_csv)
write_csv(null_permutations, output_permutation_csv)
write_csv(checks, output_check_csv)

plot_data <- null_summary |>
  mutate(
    above_null = .data$empirical_p_upper <= 0.05,
    p_label = if_else(
      .data$empirical_p_upper <= 0.001,
      "p <= 0.001",
      paste0("p = ", number(.data$empirical_p_upper, accuracy = 0.001))
    )
  )

y_limit <- max(plot_data$observed_r2t_beta, plot_data$null_upper_95, na.rm = TRUE) * 1.2

null_plot <- ggplot(plot_data, aes(x = atlas)) +
  geom_col(
    aes(y = observed_r2t_beta, fill = variable, alpha = above_null),
    width = 0.62,
    colour = "white",
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(ymin = null_lower_95, ymax = null_upper_95),
    width = 0.18,
    linewidth = 0.52,
    colour = "grey15"
  ) +
  geom_point(
    aes(y = null_mean),
    shape = 21,
    size = 1.8,
    stroke = 0.45,
    fill = "white",
    colour = "grey15"
  ) +
  geom_text(
    aes(
      y = pmax(observed_r2t_beta, null_upper_95),
      label = p_label
    ),
    vjust = -0.45,
    size = 2.25,
    colour = "grey12"
  ) +
  facet_wrap(vars(variable), ncol = 4) +
  scale_fill_manual(values = env_colors, guide = "none") +
  scale_alpha_manual(
    values = c("TRUE" = 1, "FALSE" = 0.38),
    breaks = c(TRUE, FALSE),
    labels = c("Observed above null", "Observed within null"),
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(0, y_limit),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    title = "Trait-profile shuffle null for species-environment response variation",
    x = NULL,
    y = "R2T_Beta"
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
    legend.text = element_text(size = 8.4),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.spacing = grid::unit(0.75, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

ggsave(output_png, null_plot, width = 8.4, height = 5.1, units = "in", dpi = 300, bg = "white")

message("Wrote: ", output_csv)
message("Wrote: ", output_permutation_csv)
message("Wrote: ", output_check_csv)
message("Wrote: ", output_png)
