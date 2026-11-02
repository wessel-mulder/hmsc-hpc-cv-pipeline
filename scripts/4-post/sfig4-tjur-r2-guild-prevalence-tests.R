#!/usr/bin/env Rscript

# Focused Tjur R2 temporal-change tests for diurnal raptors and four
# increasing-prevalence comparison guilds.
#
# Inputs:
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig4-guild-cv-accuracy-temporal-change-species-values.csv
#   manuscript-outputs/tables/methods-audit/modelled_species_traits_prevalence.csv
#
# Outputs:
#   manuscript-outputs/figures/supplement/<model-pattern>-sfig4-tjur-r2-guild-prevalence-change.png
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig4-tjur-r2-guild-prevalence-change-values.csv
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig4-tjur-r2-guild-increase-tests.csv
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig4-tjur-r2-raptor-vs-guild-change-tests.csv
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig4-tjur-r2-prevalence-adjusted-models.csv
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig4-tjur-r2-guild-slope-estimates.csv
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig4-tjur-r2-guild-slope-pairwise-tests.csv
#   manuscript-outputs/tables/supplement/<model-pattern>-sfig4-tjur-r2-guild-slope-interaction-tests.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

model_pattern <- "2026-07-27"
n_permutations <- 99999
random_seed <- 42

input_change_path <- file.path(
  "manuscript-outputs",
  "tables",
  "supplement",
  paste0(model_pattern, "-sfig4-guild-cv-accuracy-temporal-change-species-values.csv")
)
trait_table_path <- file.path(
  "manuscript-outputs",
  "tables",
  "methods-audit",
  "modelled_species_traits_prevalence.csv"
)

output_figure_dir <- file.path("manuscript-outputs", "figures", "supplement")
output_table_dir <- file.path("manuscript-outputs", "tables", "supplement")

output_figure_path <- file.path(
  output_figure_dir,
  paste0(model_pattern, "-sfig4-tjur-r2-guild-prevalence-change.png")
)
analysis_values_path <- file.path(
  output_table_dir,
  paste0(model_pattern, "-sfig4-tjur-r2-guild-prevalence-change-values.csv")
)
guild_increase_tests_path <- file.path(
  output_table_dir,
  paste0(model_pattern, "-sfig4-tjur-r2-guild-increase-tests.csv")
)
raptor_vs_guild_tests_path <- file.path(
  output_table_dir,
  paste0(model_pattern, "-sfig4-tjur-r2-raptor-vs-guild-change-tests.csv")
)
prevalence_model_path <- file.path(
  output_table_dir,
  paste0(model_pattern, "-sfig4-tjur-r2-prevalence-adjusted-models.csv")
)
guild_slope_estimates_path <- file.path(
  output_table_dir,
  paste0(model_pattern, "-sfig4-tjur-r2-guild-slope-estimates.csv")
)
guild_slope_pairwise_tests_path <- file.path(
  output_table_dir,
  paste0(model_pattern, "-sfig4-tjur-r2-guild-slope-pairwise-tests.csv")
)
guild_slope_interaction_tests_path <- file.path(
  output_table_dir,
  paste0(model_pattern, "-sfig4-tjur-r2-guild-slope-interaction-tests.csv")
)

selected_guilds <- c(
  "Diurnal raptors",
  "Tit-like birds",
  "Passerine seedeaters",
  "Low flycatching feeders",
  "Foliage gleaners"
)

contrast_levels <- c("1990s - 1970s", "2010s - 1970s")

guild_cols <- c(
  "Diurnal raptors" = "#b23a48",
  "Tit-like birds" = "#26847d",
  "Passerine seedeaters" = "#b9832d",
  "Low flycatching feeders" = "#4b6fb3",
  "Foliage gleaners" = "#4f8f3a"
)

if (!file.exists(input_change_path)) {
  stop(
    "Missing species-level accuracy-change table: ",
    input_change_path,
    "\nRun scripts/4-data-visualisation/sfig4-guild-cv-accuracy-temporal-change.R first.",
    call. = FALSE
  )
}
if (!file.exists(trait_table_path)) {
  stop("Missing species trait/prevalence table: ", trait_table_path, call. = FALSE)
}

accuracy_changes <- read.csv(input_change_path, check.names = FALSE)
species_traits <- read.csv(trait_table_path, check.names = FALSE)

required_trait_columns <- c(
  "species",
  "atlas1_prevalence",
  "atlas2_prevalence",
  "atlas3_prevalence"
)
missing_trait_columns <- setdiff(required_trait_columns, colnames(species_traits))
if (length(missing_trait_columns) > 0) {
  stop("Missing prevalence columns: ", paste(missing_trait_columns, collapse = ", "), call. = FALSE)
}

if (anyDuplicated(species_traits$species) > 0) {
  stop("Species names are not unique in the trait/prevalence table.", call. = FALSE)
}

tjur_analysis <- accuracy_changes |>
  filter(.data$metric == "TjurR2", .data$contrast %in% contrast_levels) |>
  left_join(
    species_traits |>
      select(
        "species",
        "atlas1_prevalence",
        "atlas2_prevalence",
        "atlas3_prevalence"
      ),
    by = "species"
  ) |>
  mutate(
    contrast = factor(.data$contrast, levels = contrast_levels),
    guild = factor(.data$guild, levels = selected_guilds),
    tjur_start = .data$atlas1,
    tjur_end = if_else(.data$contrast == "1990s - 1970s", .data$atlas2, .data$atlas3),
    tjur_change = .data$change,
    prevalence_start = .data$atlas1_prevalence,
    prevalence_end = if_else(
      .data$contrast == "1990s - 1970s",
      .data$atlas2_prevalence,
      .data$atlas3_prevalence
    ),
    prevalence_change = .data$prevalence_end - .data$prevalence_start,
    is_raptor = .data$guild == "Diurnal raptors"
  ) |>
  select(
    "species",
    "guild",
    "contrast",
    "tjur_start",
    "tjur_end",
    "tjur_change",
    "prevalence_start",
    "prevalence_end",
    "prevalence_change",
    "is_raptor"
  )

exact_sign_flip_median_test <- function(values) {
  values <- values[!is.na(values)]
  n_values <- length(values)
  observed_median <- median(values)

  if (n_values == 0) {
    stop("Cannot test an empty vector.", call. = FALSE)
  }

  signs <- expand.grid(rep(list(c(-1, 1)), n_values))
  null_medians <- apply(as.matrix(signs), 1, \(x) median(values * x))
  tolerance <- sqrt(.Machine$double.eps)

  p_greater <- mean(null_medians >= observed_median - tolerance)
  p_less <- mean(null_medians <= observed_median + tolerance)

  tibble(
    median_change = observed_median,
    sign_flip_null_size = length(null_medians),
    p_sign_flip_median_one_sided_increase = p_greater,
    p_sign_flip_median_two_sided = min(1, 2 * min(p_greater, p_less))
  )
}

within_guild_test <- function(data) {
  values <- data$tjur_change
  nonzero_values <- values[!is.na(values) & values != 0]
  n_increased <- sum(nonzero_values > 0)
  n_nonzero <- length(nonzero_values)
  sign_test_p <- if (n_nonzero > 0) {
    binom.test(n_increased, n_nonzero, p = 0.5, alternative = "greater")$p.value
  } else {
    NA_real_
  }

  exact_sign_flip_median_test(values) |>
    mutate(
      n_species = sum(!is.na(values)),
      mean_change = mean(values, na.rm = TRUE),
      n_species_increased = n_increased,
      prop_species_increased = if_else(n_nonzero > 0, n_increased / n_nonzero, NA_real_),
      p_sign_test_one_sided_increase = sign_test_p,
      .before = "median_change"
    )
}

exact_label_permutation_median_test <- function(raptor_values, comparison_values) {
  raptor_values <- raptor_values[!is.na(raptor_values)]
  comparison_values <- comparison_values[!is.na(comparison_values)]
  all_values <- c(raptor_values, comparison_values)
  n_raptors <- length(raptor_values)
  observed_difference <- median(raptor_values) - median(comparison_values)
  combinations <- combn(seq_along(all_values), n_raptors)
  null_differences <- apply(combinations, 2, \(idx) {
    median(all_values[idx]) - median(all_values[-idx])
  })
  tolerance <- sqrt(.Machine$double.eps)

  tibble(
    median_change_raptors = median(raptor_values),
    median_change_comparison = median(comparison_values),
    median_difference_raptors_minus_comparison = observed_difference,
    null_size = length(null_differences),
    method = "exact_label_permutation",
    p_permutation_two_sided = mean(abs(null_differences) >= abs(observed_difference) - tolerance)
  )
}

monte_carlo_label_permutation_median_test <- function(raptor_values, comparison_values, n_perm = n_permutations) {
  raptor_values <- raptor_values[!is.na(raptor_values)]
  comparison_values <- comparison_values[!is.na(comparison_values)]
  all_values <- c(raptor_values, comparison_values)
  n_raptors <- length(raptor_values)
  observed_difference <- median(raptor_values) - median(comparison_values)
  n_total <- length(all_values)
  tolerance <- sqrt(.Machine$double.eps)

  null_differences <- replicate(n_perm, {
    pseudo_raptors <- sample.int(n_total, size = n_raptors, replace = FALSE)
    median(all_values[pseudo_raptors]) - median(all_values[-pseudo_raptors])
  })

  tibble(
    median_change_raptors = median(raptor_values),
    median_change_comparison = median(comparison_values),
    median_difference_raptors_minus_comparison = observed_difference,
    null_size = n_perm,
    method = "monte_carlo_label_permutation",
    p_permutation_two_sided = (sum(abs(null_differences) >= abs(observed_difference) - tolerance) + 1) /
      (n_perm + 1)
  )
}

between_guild_tests_for_contrast <- function(data) {
  comparison_guilds <- setdiff(selected_guilds, "Diurnal raptors")
  raptor_values <- data |>
    filter(.data$guild == "Diurnal raptors") |>
    pull(.data$tjur_change)

  pairwise_tests <- purrr::map_dfr(comparison_guilds, \(comparison_guild) {
    comparison_values <- data |>
      filter(.data$guild == comparison_guild) |>
      pull(.data$tjur_change)

    exact_label_permutation_median_test(raptor_values, comparison_values) |>
      mutate(
        comparison = paste("Diurnal raptors vs", comparison_guild),
        comparison_guild = comparison_guild,
        n_raptors = length(raptor_values),
        n_comparison = length(comparison_values),
        .before = "median_change_raptors"
      )
  })

  pooled_values <- data |>
    filter(.data$guild != "Diurnal raptors") |>
    pull(.data$tjur_change)

  pooled_test <- monte_carlo_label_permutation_median_test(raptor_values, pooled_values) |>
    mutate(
      comparison = "Diurnal raptors vs selected non-raptor guilds pooled",
      comparison_guild = "Selected non-raptor guilds pooled",
      n_raptors = length(raptor_values),
      n_comparison = length(pooled_values),
      .before = "median_change_raptors"
    )

  bind_rows(pairwise_tests, pooled_test)
}

tidy_lm <- function(model, contrast_label, model_id, model_formula) {
  model_summary <- summary(model)
  coefficients <- as.data.frame(model_summary$coefficients)

  tibble(
    contrast = contrast_label,
    model_id = model_id,
    model_formula = model_formula,
    term = rownames(coefficients),
    estimate = coefficients[, "Estimate"],
    std_error = coefficients[, "Std. Error"],
    statistic = coefficients[, "t value"],
    p_value = coefficients[, "Pr(>|t|)"],
    r_squared = model_summary$r.squared,
    adjusted_r_squared = model_summary$adj.r.squared,
    n_species = nobs(model)
  )
}

fit_prevalence_models <- function(data, contrast_label) {
  model_specs <- list(
    prevalence_only = tjur_change ~ prevalence_change,
    raptor_only = tjur_change ~ is_raptor,
    raptor_plus_prevalence_change = tjur_change ~ is_raptor + prevalence_change,
    raptor_plus_starting_prevalence_and_change = tjur_change ~ is_raptor + prevalence_start + prevalence_change
  )

  purrr::imap_dfr(model_specs, \(formula, model_id) {
    tidy_lm(
      lm(formula, data = data),
      contrast_label = contrast_label,
      model_id = model_id,
      model_formula = paste(deparse(formula), collapse = " ")
    )
  })
}

slope_basis <- function(model, guild_name) {
  coefficient_names <- names(coef(model))
  basis <- setNames(rep(0, length(coefficient_names)), coefficient_names)
  basis[["prevalence_change"]] <- 1

  if (guild_name != selected_guilds[[1]]) {
    interaction_names <- c(
      paste0("prevalence_change:guild", guild_name),
      paste0("guild", guild_name, ":prevalence_change")
    )
    interaction_name <- interaction_names[interaction_names %in% coefficient_names][[1]]
    basis[[interaction_name]] <- 1
  }

  basis
}

estimate_from_basis <- function(model, basis) {
  coefficients <- coef(model)
  covariance <- vcov(model)
  estimate <- sum(basis * coefficients)
  std_error <- sqrt(as.numeric(t(basis) %*% covariance %*% basis))
  df <- df.residual(model)
  statistic <- estimate / std_error

  tibble(
    estimate = estimate,
    std_error = std_error,
    df = df,
    statistic = statistic,
    p_value = 2 * pt(abs(statistic), df = df, lower.tail = FALSE),
    conf_low = estimate - qt(0.975, df = df) * std_error,
    conf_high = estimate + qt(0.975, df = df) * std_error
  )
}

fit_guild_slope_models <- function(data, contrast_label) {
  slope_model <- lm(tjur_change ~ prevalence_change * guild, data = data)
  additive_model <- lm(tjur_change ~ prevalence_change + guild, data = data)
  interaction_test <- anova(additive_model, slope_model)
  model_summary <- summary(slope_model)

  slope_estimates <- purrr::map_dfr(selected_guilds, \(guild_name) {
    basis <- slope_basis(slope_model, guild_name)
    estimate_from_basis(slope_model, basis) |>
      mutate(
        contrast = contrast_label,
        guild = guild_name,
        n_species = sum(data$guild == guild_name),
        r_squared = model_summary$r.squared,
        adjusted_r_squared = model_summary$adj.r.squared,
        .before = "estimate"
      )
  })

  pairwise_tests <- purrr::map_dfr(combn(selected_guilds, 2, simplify = FALSE), \(guild_pair) {
    basis_difference <- slope_basis(slope_model, guild_pair[[1]]) -
      slope_basis(slope_model, guild_pair[[2]])

    estimate_from_basis(slope_model, basis_difference) |>
      mutate(
        contrast = contrast_label,
        guild_1 = guild_pair[[1]],
        guild_2 = guild_pair[[2]],
        slope_difference_guild_1_minus_guild_2 = .data$estimate,
        .before = "estimate"
      )
  })

  omnibus_test <- tibble(
    contrast = contrast_label,
    model_comparison = "tjur_change ~ prevalence_change + guild vs tjur_change ~ prevalence_change * guild",
    df_additive = interaction_test$Res.Df[[1]],
    df_interaction = interaction_test$Res.Df[[2]],
    df_difference = interaction_test$Df[[2]],
    f_statistic = interaction_test$F[[2]],
    p_value = interaction_test$`Pr(>F)`[[2]],
    r_squared_interaction_model = model_summary$r.squared,
    adjusted_r_squared_interaction_model = model_summary$adj.r.squared,
    n_species = nobs(slope_model)
  )

  list(
    slope_estimates = slope_estimates,
    pairwise_tests = pairwise_tests,
    omnibus_test = omnibus_test
  )
}

guild_increase_tests <- tjur_analysis |>
  group_by(.data$guild, .data$contrast) |>
  group_modify(~ within_guild_test(.x)) |>
  ungroup() |>
  mutate(
    p_fdr_bh_sign_flip_increase = p.adjust(.data$p_sign_flip_median_one_sided_increase, method = "BH"),
    p_fdr_bh_sign_test_increase = p.adjust(.data$p_sign_test_one_sided_increase, method = "BH")
  ) |>
  arrange(.data$contrast, .data$guild)

set.seed(random_seed)

raptor_vs_guild_tests <- tjur_analysis |>
  group_by(.data$contrast) |>
  group_modify(~ between_guild_tests_for_contrast(.x)) |>
  ungroup() |>
  mutate(
    p_fdr_bh = p.adjust(.data$p_permutation_two_sided, method = "BH"),
    random_seed = if_else(.data$method == "monte_carlo_label_permutation", random_seed, NA_integer_)
  ) |>
  arrange(.data$contrast, .data$p_permutation_two_sided)

prevalence_models <- split(tjur_analysis, tjur_analysis$contrast) |>
  purrr::imap_dfr(\(data, contrast_label) {
    fit_prevalence_models(data, contrast_label = contrast_label)
  })

guild_slope_model_results <- split(tjur_analysis, tjur_analysis$contrast) |>
  purrr::imap(\(data, contrast_label) {
    fit_guild_slope_models(data, contrast_label = contrast_label)
  })

guild_slope_estimates <- purrr::map_dfr(guild_slope_model_results, "slope_estimates") |>
  arrange(.data$contrast, .data$guild)

guild_slope_pairwise_tests <- purrr::map_dfr(guild_slope_model_results, "pairwise_tests") |>
  group_by(.data$contrast) |>
  mutate(p_fdr_bh = p.adjust(.data$p_value, method = "BH")) |>
  ungroup() |>
  arrange(.data$contrast, .data$p_value)

guild_slope_interaction_tests <- purrr::map_dfr(guild_slope_model_results, "omnibus_test") |>
  mutate(p_fdr_bh = p.adjust(.data$p_value, method = "BH"))

slope_interaction_annotations <- guild_slope_interaction_tests |>
  transmute(
    contrast = factor(.data$contrast, levels = contrast_levels),
    x = Inf,
    y = Inf,
    label = paste0(
      "Guild slope interaction: R2 = ",
      sprintf("%.2f", .data$r_squared_interaction_model),
      ", p = ",
      format.pval(.data$p_value, digits = 2, eps = 0.001)
    )
  )

tjur_prevalence_plot <- ggplot(
  tjur_analysis,
  aes(x = .data$prevalence_change, y = .data$tjur_change)
) +
  geom_hline(yintercept = 0, linewidth = 0.28, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.28, colour = "grey55") +
  geom_smooth(
    aes(colour = .data$guild),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.62
  ) +
  geom_point(
    aes(fill = .data$guild),
    shape = 21,
    size = 2.2,
    stroke = 0.20,
    colour = "grey18",
    alpha = 0.82
  ) +
  geom_text(
    data = slope_interaction_annotations,
    aes(x = .data$x, y = .data$y, label = .data$label),
    inherit.aes = FALSE,
    hjust = 1.02,
    vjust = 1.30,
    size = 2.85,
    colour = "grey20"
  ) +
  facet_wrap(vars(.data$contrast), nrow = 1) +
  scale_colour_manual(values = guild_cols, name = NULL) +
  scale_fill_manual(values = guild_cols, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.10))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.14))) +
  labs(
    x = "Change in species prevalence",
    y = "Change in 10-fold cross-validated Tjur R2"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.70, "lines"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text = element_text(face = "bold", colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    axis.title = element_text(size = 9),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.width = unit(0.42, "lines"),
    plot.background = element_rect(fill = "white", colour = NA)
  ) +
  guides(colour = "none")

dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  output_figure_path,
  tjur_prevalence_plot,
  width = 7.7,
  height = 3.9,
  units = "in",
  dpi = 300,
  bg = "white"
)

write.csv(tjur_analysis, analysis_values_path, row.names = FALSE)
write.csv(guild_increase_tests, guild_increase_tests_path, row.names = FALSE)
write.csv(raptor_vs_guild_tests, raptor_vs_guild_tests_path, row.names = FALSE)
write.csv(prevalence_models, prevalence_model_path, row.names = FALSE)
write.csv(guild_slope_estimates, guild_slope_estimates_path, row.names = FALSE)
write.csv(guild_slope_pairwise_tests, guild_slope_pairwise_tests_path, row.names = FALSE)
write.csv(guild_slope_interaction_tests, guild_slope_interaction_tests_path, row.names = FALSE)

print(guild_increase_tests, n = Inf)
print(raptor_vs_guild_tests, n = Inf)
print(prevalence_models, n = Inf)
print(guild_slope_estimates, n = Inf)
print(guild_slope_pairwise_tests, n = Inf)
print(guild_slope_interaction_tests, n = Inf)
message("Wrote: ", output_figure_path)
message("Wrote: ", analysis_values_path)
message("Wrote: ", guild_increase_tests_path)
message("Wrote: ", raptor_vs_guild_tests_path)
message("Wrote: ", prevalence_model_path)
message("Wrote: ", guild_slope_estimates_path)
message("Wrote: ", guild_slope_pairwise_tests_path)
message("Wrote: ", guild_slope_interaction_tests_path)
