rm(list = ls())

# Posterior uncertainty for the scaled variance-partitioning summaries.
#
# The companion figure `fig2c-scaled-vp-environment-vs-spatial.R` shows the
# posterior-mean scaled variance partitioning. This script asks whether the
# differences among atlas periods are large relative to the fitted HMSC
# posterior uncertainty. It keeps HMSC's own grouped variance-partitioning
# calculation, but repeats it one posterior draw at a time so the draw-level
# uncertainty is not averaged away.

library(tidyverse)
library(Hmsc)
library(scales)

source(file.path("support_scripts", "figure_data_helpers.R"))

# HMSC's predict.Hmsc() checks the available core count before entering its
# single-core path. In this non-interactive script that check can resolve to NA
# inside the namespace, so we clone the HMSC method and replace only the initial
# core-count assignment. Everything downstream remains HMSC's own prediction
# code.
predict_hmsc_single_core <- getS3method("predict", "Hmsc")
predict_hmsc_single_core_body <- body(predict_hmsc_single_core)
predict_hmsc_single_core_body[[2]] <- quote(nParallel <- 1L)
body(predict_hmsc_single_core) <- predict_hmsc_single_core_body

model_dir <- "HmscOutputs"
model_pattern <- "2026-03-13"
output_dir <- file.path("misc-figures", "outputs", "main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure_slug <- paste0(model_pattern, "-fig2c-scaled-vp-posterior-uncertainty")
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))
output_atlas_csv <- file.path(output_dir, paste0(figure_slug, "-atlas-summary.csv"))
output_contrast_csv <- file.path(output_dir, paste0(figure_slug, "-contrasts.csv"))

current_summary_csv <- file.path(
  output_dir,
  paste0(model_pattern, "-fig2c-scaled-vp-environment-vs-spatial-summary.csv")
)

atlas_lookup <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

atlas_order <- unname(atlas_lookup)

vp_component_order <- c("Climate", "Land-use", "Spatial")
plot_component_order <- c("mean R2", "Spatial", "Land-use", "Climate")

component_colors <- c(
  "mean R2" = "grey25",
  "Climate" = "firebrick3",
  "Land-use" = "goldenrod1",
  "Spatial" = "ivory2"
)

# HMSC groups the fixed-effect design matrix. The intercept has zero variance,
# but assigning it to climate matches HMSC's default handling of the intercept
# and keeps every design-matrix column assigned to a group.
fixed_effect_group_names <- c("Climate", "Land-use")
fixed_effect_group_by_covariate <- c(
  "(Intercept)" = 1,
  "tmean_breeding" = 1,
  "prec_breeding" = 1,
  "hh" = 2,
  "perc_urban" = 2,
  "perc_cropland" = 2,
  "perc_pasture" = 2,
  "perc_forest" = 2,
  "perc_grass_shrub" = 2
)

# Prediction is run in chunks because all posterior predictions at once would
# be much larger than needed. A chunk size below two can trigger a one-draw edge
# case in HMSC's random-effect prediction path, so keep this >= 2.
prediction_chunk_size <- 40
posterior_start <- 1
posterior_thin <- 1

model_folders <- figure_model_folders(pattern = model_pattern, base_dir = model_dir)

compute_tjur_r2 <- function(y, pred_y) {
  vapply(seq_len(ncol(y)), function(species_index) {
    observed <- y[, species_index]
    predicted <- pred_y[, species_index]

    presences <- observed == 1
    absences <- observed == 0

    if (!any(presences, na.rm = TRUE) || !any(absences, na.rm = TRUE)) {
      return(NA_real_)
    }

    mean(predicted[presences], na.rm = TRUE) -
      mean(predicted[absences], na.rm = TRUE)
  }, numeric(1))
}

load_hmsc_model <- function(model_folder) {
  model_path <- file.path(model_dir, model_folder, "Models", "Fitted", hmsc_fitted_file)
  model_env <- new.env(parent = emptyenv())
  load(model_path, envir = model_env)

  if (!exists("fitted_model", envir = model_env, inherits = FALSE)) {
    stop("Expected `fitted_model` in ", model_path, call. = FALSE)
  }

  model_env$fitted_model$posteriors
}

pooled_posterior_draws <- function(model) {
  poolMcmcChains(
    model$postList,
    start = posterior_start,
    thin = posterior_thin
  )
}

prediction_chunks <- function(draw_count) {
  starts <- seq(1, draw_count, by = prediction_chunk_size)
  purrr::map(starts, ~ .x:min(.x + prediction_chunk_size - 1, draw_count))
}

compute_tjur_draw_matrix <- function(model, draws, atlas_label) {
  draw_count <- length(draws)
  tjur_by_draw <- matrix(
    NA_real_,
    nrow = draw_count,
    ncol = model$ns,
    dimnames = list(seq_len(draw_count), model$spNames)
  )

  for (draw_indices in prediction_chunks(draw_count)) {
    message(
      "  ", atlas_label, ": posterior predictions for draws ",
      min(draw_indices), "-", max(draw_indices), " / ", draw_count
    )

    pred_list <- predict_hmsc_single_core(
      model,
      post = draws[draw_indices],
      expected = TRUE,
      nParallel = 1
    )

    for (i in seq_along(draw_indices)) {
      draw_id <- draw_indices[[i]]
      tjur_by_draw[draw_id, ] <- compute_tjur_r2(model$Y, pred_list[[i]])
    }
  }

  tjur_by_draw
}

validate_group_vector <- function(model) {
  missing_covariates <- setdiff(model$covNames, names(fixed_effect_group_by_covariate))

  if (length(missing_covariates) > 0) {
    stop(
      "The posterior VP script does not know how to group these covariates: ",
      paste(missing_covariates, collapse = ", "),
      call. = FALSE
    )
  }

  unname(fixed_effect_group_by_covariate[model$covNames])
}

compute_scaled_vp_draws <- function(model, draws, tjur_by_draw, atlas_label) {
  fixed_effect_group <- validate_group_vector(model)
  draw_count <- length(draws)

  purrr::map_dfr(seq_len(draw_count), function(draw_id) {
    if (draw_id %% 100 == 0) {
      message("  ", atlas_label, ": grouped VP for draw ", draw_id, " / ", draw_count)
    }

    # HMSC's computeVariancePartitioning() averages over hM$postList. Supplying
    # a one-draw postList lets us retain that draw's grouped VP values.
    one_draw_model <- model
    one_draw_model$postList <- list(list(draws[[draw_id]]))

    vp <- computeVariancePartitioning(
      one_draw_model,
      group = fixed_effect_group,
      groupnames = fixed_effect_group_names,
      start = 1
    )

    vp_vals <- vp$vals
    rownames(vp_vals) <- recode(
      rownames(vp_vals),
      "Random: site" = "Spatial"
    )
    vp_vals <- vp_vals[vp_component_order, , drop = FALSE]

    tjur_r2 <- tjur_by_draw[draw_id, colnames(vp_vals)]
    scaled_vp <- sweep(vp_vals, 2, tjur_r2, `*`)
    scaling_check <- max(abs(colSums(scaled_vp, na.rm = TRUE) - tjur_r2), na.rm = TRUE)
    component_means <- rowMeans(scaled_vp, na.rm = TRUE)

    tibble(
      draw_id = draw_id,
      component = c("mean R2", names(component_means)),
      value = c(mean(tjur_r2, na.rm = TRUE), unname(component_means)),
      max_scaling_abs_diff = scaling_check
    )
  })
}

posterior_draw_summary <- purrr::map_dfr(model_folders, function(model_folder) {
  atlas_index <- atlas_numbers(model_folder)
  atlas_label <- atlas_lookup[[as.character(atlas_index)]]

  message("Processing ", atlas_label, " from ", model_folder)
  model <- load_hmsc_model(model_folder)
  draws <- pooled_posterior_draws(model)

  message("  ", atlas_label, ": using ", length(draws), " pooled posterior draws")
  tjur_by_draw <- compute_tjur_draw_matrix(model, draws, atlas_label)
  scaled_vp_draws <- compute_scaled_vp_draws(model, draws, tjur_by_draw, atlas_label)

  rm(model, draws, tjur_by_draw)
  gc()

  scaled_vp_draws |>
    mutate(
      atlas = factor(atlas_label, levels = atlas_order),
      atlas_index = atlas_index,
      .before = 1
    )
})

max_scaling_error <- max(posterior_draw_summary$max_scaling_abs_diff, na.rm = TRUE)
if (max_scaling_error > 1e-8) {
  stop(
    "Scaled VP components do not sum to draw-specific Tjur R2. Max difference: ",
    signif(max_scaling_error, 4),
    call. = FALSE
  )
}

if (file.exists(current_summary_csv)) {
  current_component_summary <- read_csv(
    current_summary_csv,
    col_types = cols(atlas = col_character(), .default = col_guess())
  ) |>
    mutate(
      component = recode(.data$component, "Spatial / biotic dependencies" = "Spatial")
    )

  current_figure_summary <- bind_rows(
    current_component_summary |>
      distinct(.data$atlas, .data$mean_total_explained) |>
      transmute(
        atlas = factor(.data$atlas, levels = atlas_order),
        component = "mean R2",
        current_figure_mean = .data$mean_total_explained
      ),
    current_component_summary |>
      transmute(
        atlas = factor(.data$atlas, levels = atlas_order),
        component = .data$component,
        current_figure_mean = .data$mean_scaled_variance_explained
      )
  )
} else {
  warning("Current summary CSV not found, so current-figure comparison columns will be NA.")
  current_figure_summary <- posterior_draw_summary |>
    distinct(.data$atlas, .data$component) |>
    mutate(current_figure_mean = NA_real_)
}

atlas_summary <- posterior_draw_summary |>
  group_by(.data$atlas, .data$atlas_index, .data$component) |>
  summarise(
    n_draws = n(),
    posterior_mean = mean(.data$value, na.rm = TRUE),
    posterior_median = median(.data$value, na.rm = TRUE),
    q025 = quantile(.data$value, 0.025, na.rm = TRUE),
    q10 = quantile(.data$value, 0.10, na.rm = TRUE),
    q90 = quantile(.data$value, 0.90, na.rm = TRUE),
    q975 = quantile(.data$value, 0.975, na.rm = TRUE),
    pr_gt0 = mean(.data$value > 0, na.rm = TRUE),
    max_scaling_abs_diff = max(.data$max_scaling_abs_diff, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(current_figure_summary, by = c("atlas", "component")) |>
  mutate(
    posterior_mean_pp = 100 * .data$posterior_mean,
    posterior_median_pp = 100 * .data$posterior_median,
    q025_pp = 100 * .data$q025,
    q10_pp = 100 * .data$q10,
    q90_pp = 100 * .data$q90,
    q975_pp = 100 * .data$q975,
    current_figure_mean_pp = 100 * .data$current_figure_mean,
    posterior_minus_current_pp = 100 * (.data$posterior_mean - .data$current_figure_mean)
  )

contrast_specs <- tribble(
  ~contrast, ~later_atlas, ~earlier_atlas,
  "1990s - 1970s", "1990s", "1970s",
  "2010s - 1990s", "2010s", "1990s",
  "2010s - 1970s", "2010s", "1970s"
)

posterior_draw_values <- posterior_draw_summary |>
  select("atlas", "draw_id", "component", "value")

contrast_draws <- purrr::pmap_dfr(contrast_specs, function(contrast, later_atlas, earlier_atlas) {
  later <- posterior_draw_values |>
    filter(.data$atlas == later_atlas) |>
    rename(later_value = value)

  earlier <- posterior_draw_values |>
    filter(.data$atlas == earlier_atlas) |>
    rename(earlier_value = value)

  inner_join(
    later,
    earlier,
    by = c("draw_id", "component"),
    suffix = c("_later", "_earlier")
  ) |>
    transmute(
      contrast = contrast,
      draw_id = .data$draw_id,
      component = .data$component,
      diff = .data$later_value - .data$earlier_value
    )
})

contrast_summary <- contrast_draws |>
  group_by(.data$contrast, .data$component) |>
  summarise(
    n_draws = n(),
    diff_mean = mean(.data$diff, na.rm = TRUE),
    diff_median = median(.data$diff, na.rm = TRUE),
    q025 = quantile(.data$diff, 0.025, na.rm = TRUE),
    q10 = quantile(.data$diff, 0.10, na.rm = TRUE),
    q90 = quantile(.data$diff, 0.90, na.rm = TRUE),
    q975 = quantile(.data$diff, 0.975, na.rm = TRUE),
    pr_diff_gt0 = mean(.data$diff > 0, na.rm = TRUE),
    pr_diff_lt0 = mean(.data$diff < 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    contrast = factor(.data$contrast, levels = contrast_specs$contrast),
    diff_mean_pp = 100 * .data$diff_mean,
    diff_median_pp = 100 * .data$diff_median,
    q025_pp = 100 * .data$q025,
    q10_pp = 100 * .data$q10,
    q90_pp = 100 * .data$q90,
    q975_pp = 100 * .data$q975,
    component = factor(.data$component, levels = rev(plot_component_order)),
    component_y = as.numeric(.data$component),
    pr_label = paste0("Pr(>0)=", number(.data$pr_diff_gt0, accuracy = 0.01))
  ) |>
  arrange(.data$contrast, desc(.data$component_y))

label_offset <- max(abs(c(contrast_summary$q025_pp, contrast_summary$q975_pp)), na.rm = TRUE) * 0.05
label_offset <- max(label_offset, 0.15)

contrast_summary <- contrast_summary |>
  mutate(
    label_x = if_else(
      .data$diff_median_pp >= 0,
      .data$q975_pp + label_offset,
      .data$q025_pp - label_offset
    ),
    label_hjust = if_else(.data$diff_median_pp >= 0, 0, 1)
  )

x_limit <- max(abs(contrast_summary$label_x), na.rm = TRUE) + label_offset

write_csv(atlas_summary, output_atlas_csv)
write_csv(contrast_summary, output_contrast_csv)

posterior_plot <- ggplot(contrast_summary) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey45") +
  geom_segment(
    aes(
      x = .data$q025_pp,
      xend = .data$q975_pp,
      y = .data$component_y,
      yend = .data$component_y
    ),
    linewidth = 1.2,
    color = "grey55",
    lineend = "round"
  ) +
  geom_segment(
    aes(
      x = .data$q025_pp,
      xend = .data$q975_pp,
      y = .data$component_y,
      yend = .data$component_y,
      color = .data$component
    ),
    linewidth = 0.75,
    alpha = 0.75,
    lineend = "round"
  ) +
  geom_segment(
    aes(
      x = .data$q10_pp,
      xend = .data$q90_pp,
      y = .data$component_y,
      yend = .data$component_y
    ),
    linewidth = 4,
    color = "grey45",
    lineend = "round"
  ) +
  geom_segment(
    aes(
      x = .data$q10_pp,
      xend = .data$q90_pp,
      y = .data$component_y,
      yend = .data$component_y,
      color = .data$component
    ),
    linewidth = 2.8,
    lineend = "round"
  ) +
  geom_point(
    aes(
      x = .data$diff_median_pp,
      y = .data$component_y,
      fill = .data$component
    ),
    shape = 21,
    size = 3.2,
    stroke = 0.35,
    color = "grey20"
  ) +
  geom_text(
    aes(
      x = .data$label_x,
      y = .data$component_y,
      label = .data$pr_label,
      hjust = .data$label_hjust
    ),
    size = 2.8,
    color = "grey20"
  ) +
  facet_wrap(vars(.data$contrast), nrow = 1) +
  scale_y_continuous(
    breaks = seq_along(rev(plot_component_order)),
    labels = rev(plot_component_order),
    expand = expansion(add = 0.55)
  ) +
  scale_x_continuous(
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  scale_color_manual(values = component_colors) +
  scale_fill_manual(values = component_colors) +
  coord_cartesian(xlim = c(-x_limit, x_limit), clip = "off") +
  labs(
    title = "Posterior uncertainty in scaled variance-partitioning differences",
    subtitle = "Points are posterior medians; thick and thin intervals show 80% and 95% posterior intervals.",
    x = "Difference in mean scaled variance explained (percentage points)",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "grey99", color = NA),
    panel.background = element_rect(fill = "grey99", color = NA),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 10),
    plot.margin = margin(8, 26, 8, 8)
  )

ggsave(
  filename = output_png,
  plot = posterior_plot,
  width = 9,
  height = 4.2,
  dpi = 320
)

message("Wrote: ", output_png)
message("Wrote: ", output_atlas_csv)
message("Wrote: ", output_contrast_csv)
