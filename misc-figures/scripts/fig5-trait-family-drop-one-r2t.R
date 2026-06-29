rm(list = ls())

# Post-hoc drop-one trait-family decomposition of VP$R2T$Beta.
#
# This script does not refit the HMSC models. It reuses posterior draws from the
# fitted models and asks how much of the trait-explained species-environment
# response variation is lost when each trait family is removed from the fitted
# trait design matrix.

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

figure_slug <- paste0(model_pattern, "-trait-family-drop-one-r2t")
output_csv <- file.path(output_dir, paste0(figure_slug, ".csv"))
output_draws_csv <- file.path(output_dir, paste0(figure_slug, "-posterior-draws.csv"))
output_check_csv <- file.path(output_dir, paste0(figure_slug, "-checks.csv"))
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))
output_pdf <- file.path(output_dir, paste0(figure_slug, ".pdf"))

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

safe_cor2 <- function(x, y) {
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
    return(NA_real_)
  }

  cor(x, y, use = "pairwise.complete.obs")^2
}

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
    rename(variable_raw = 1, exported_full_r2t_beta = 2) |>
    filter(.data$variable_raw %in% names(var_rename)) |>
    mutate(
      atlas = atlas_rename[atlas_num],
      variable = var_rename[.data$variable_raw]
    ) |>
    select(atlas, variable_raw, variable, exported_full_r2t_beta)
}

compute_drop_one_draws <- function(model, atlas_label) {
  trait_names <- colnames(model$Tr)
  trait_families <- trait_family_columns(trait_names)
  posterior_draws <- Hmsc:::poolMcmcChains(model$postList, start = posterior_start)
  variable_raw <- intersect(names(var_rename), model$covNames)
  variable_index <- match(variable_raw, model$covNames)

  map_dfr(seq_along(posterior_draws), function(draw_id) {
    draw <- posterior_draws[[draw_id]]
    full_trait_prediction <- t(model$Tr %*% t(draw$Gamma))

    map_dfr(seq_along(variable_raw), function(variable_id) {
      raw_name <- variable_raw[[variable_id]]
      beta_row <- variable_index[[variable_id]]
      beta_values <- draw$Beta[beta_row, ]
      full_r2t <- safe_cor2(beta_values, full_trait_prediction[beta_row, ])

      map_dfr(names(trait_families), function(trait_family) {
        keep_traits <- setdiff(seq_along(trait_names), trait_families[[trait_family]])
        drop_trait_prediction <- t(
          model$Tr[, keep_traits, drop = FALSE] %*%
            t(draw$Gamma[, keep_traits, drop = FALSE])
        )
        drop_r2t <- safe_cor2(beta_values, drop_trait_prediction[beta_row, ])

        tibble(
          atlas = atlas_label,
          draw = draw_id,
          variable_raw = raw_name,
          variable = var_rename[[raw_name]],
          trait_family = trait_family,
          full_r2t_beta = full_r2t,
          drop_r2t_beta = drop_r2t,
          contribution = full_r2t - drop_r2t
        )
      })
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

exported_r2t_beta <- map_dfr(model_folders, read_exported_r2t_beta)

drop_one_summary <- drop_one_draws |>
  group_by(.data$atlas, .data$variable_raw, .data$variable, .data$trait_family) |>
  summarise(
    full_r2t_beta = mean(.data$full_r2t_beta, na.rm = TRUE),
    drop_r2t_beta = mean(.data$drop_r2t_beta, na.rm = TRUE),
    contribution_mean = mean(.data$contribution, na.rm = TRUE),
    contribution_median = median(.data$contribution, na.rm = TRUE),
    contribution_lower_95 = quantile(.data$contribution, 0.025, na.rm = TRUE),
    contribution_upper_95 = quantile(.data$contribution, 0.975, na.rm = TRUE),
    prob_positive = mean(.data$contribution > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    exported_r2t_beta,
    by = c("atlas", "variable_raw", "variable")
  ) |>
  mutate(
    full_export_difference = .data$full_r2t_beta - .data$exported_full_r2t_beta,
    atlas = factor(.data$atlas, levels = unname(atlas_rename)),
    variable = factor(.data$variable, levels = var_order),
    trait_family = factor(.data$trait_family, levels = rev(trait_family_order))
  ) |>
  arrange(.data$atlas, .data$trait_family, .data$variable)

drop_one_checks <- drop_one_summary |>
  group_by(.data$atlas, .data$variable_raw, .data$variable) |>
  summarise(
    n_trait_families = n(),
    max_abs_full_export_difference = max(abs(.data$full_export_difference), na.rm = TRUE),
    .groups = "drop"
  )

expected_rows <- length(atlas_rename) * length(var_order) * length(trait_family_order)
if (nrow(drop_one_summary) != expected_rows) {
  stop(
    "Unexpected drop-one summary row count. Expected ",
    expected_rows,
    " but found ",
    nrow(drop_one_summary),
    ".",
    call. = FALSE
  )
}

if (any(drop_one_checks$n_trait_families != length(trait_family_order))) {
  stop("At least one atlas-variable combination is missing trait-family rows.", call. = FALSE)
}

if (max(drop_one_checks$max_abs_full_export_difference, na.rm = TRUE) > match_tolerance) {
  stop(
    "Recomputed full R2T_Beta does not match exported VP$R2T$Beta within ",
    match_tolerance,
    ". See: ",
    output_check_csv,
    call. = FALSE
  )
}

drop_one_summary_output <- drop_one_summary |>
  transmute(
    atlas = as.character(.data$atlas),
    variable = as.character(.data$variable),
    trait_family = as.character(.data$trait_family),
    full_r2t_beta = .data$full_r2t_beta,
    drop_r2t_beta = .data$drop_r2t_beta,
    contribution_mean = .data$contribution_mean,
    contribution_median = .data$contribution_median,
    contribution_lower_95 = .data$contribution_lower_95,
    contribution_upper_95 = .data$contribution_upper_95,
    prob_positive = .data$prob_positive
  )

write_csv(drop_one_summary_output, output_csv)
write_csv(drop_one_draws, output_draws_csv)
write_csv(drop_one_checks, output_check_csv)

contribution_limit <- max(abs(drop_one_summary$contribution_mean), na.rm = TRUE)

drop_one_plot <- ggplot(
  drop_one_summary,
  aes(x = variable, y = trait_family, fill = contribution_mean)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(
    aes(label = percent(.data$contribution_mean, accuracy = 1)),
    size = 3,
    colour = "grey12"
  ) +
  facet_wrap(vars(atlas), nrow = 1) +
  scale_fill_gradient2(
    low = "#d95f02",
    mid = "grey98",
    high = "#1b9e77",
    midpoint = 0,
    limits = c(-contribution_limit, contribution_limit),
    labels = percent_format(accuracy = 1),
    name = "Drop-one\ncontribution"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(
    title = "Post-hoc trait-family contributions to species-environment responses",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8.2, angle = 35, hjust = 1),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = "white"),
    legend.title = element_text(size = 8.8, face = "bold"),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

ggsave(output_png, drop_one_plot, width = 11.2, height = 4.8, units = "in", dpi = 300, bg = "white")
ggsave(output_pdf, drop_one_plot, width = 11.2, height = 4.8, units = "in", bg = "white")

message("Wrote: ", output_csv)
message("Wrote: ", output_draws_csv)
message("Wrote: ", output_check_csv)
message("Wrote: ", output_png)
message("Wrote: ", output_pdf)
