# Test whether 2010s foraging-guild richness has quadratic habitat relationships.
#
# Purpose:
#   At the shared Atlas 2 Good / Atlas 3 available grid cells used in the
#   guild-richness maps, test whether Atlas 3 guild richness is associated with:
#     1. habitat heterogeneity,
#     2. percentage forest cover, and
#     3. percentage agricultural cover (cropland + pasture).
#
# Analysis:
#   - Each habitat predictor is represented by a linear and squared term.
#   - All three quadratic relationships are estimated simultaneously, so each
#     curve is adjusted for the other two habitat variables.
#   - Predictors are mean-centred before squaring. This improves numerical
#     stability without changing the fitted curves.
#   - A spatial-error version tests whether the curves persist after accounting
#     for spatial autocorrelation among neighbouring 5 km atlas cells.
#   - Turning points are calculated from the fitted coefficients. Their
#     curvature determines whether each relationship is hump-shaped or U-shaped.
#   - These are cross-sectional associations for the 2010s, not causal effects.
#   - Plot output is PNG only.

rm(list = ls())

suppressPackageStartupMessages({
  library(broom)
  library(car)
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(spatialreg)
  library(spdep)
  library(tidyr)
})

source(file.path("support_scripts", "data_helpers.R"))

#### PATHS ####

guild_richness_path <- file.path(
  "notebooks",
  "exploratory",
  "outputs",
  "atlas-guild-richness",
  "atlas2-atlas3-good-shared-guild-richness-by-site.csv"
)
environment_path <- file.path(
  "Data",
  "data",
  "1_preprocessing",
  "X_environmental",
  "X_Environmental.csv"
)
grid_path <- file.path(
  "Data",
  "data",
  "1_preprocessing",
  "atlas-grids",
  "DOF_Shapefiles_",
  "DK5km_ED50grid_approx_kvadrkod_DOF.shp"
)

out_dir <- file.path(
  "notebooks",
  "exploratory",
  "outputs",
  "atlas-guild-richness",
  "atlas3-habitat-test"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_data_path <- file.path(
  out_dir,
  "atlas3-guild-richness-habitat-data.csv"
)
coefficient_path <- file.path(
  out_dir,
  "atlas3-guild-richness-quadratic-coefficients.csv"
)
model_summary_path <- file.path(
  out_dir,
  "atlas3-guild-richness-quadratic-model-summary.csv"
)
variable_tests_path <- file.path(
  out_dir,
  "atlas3-guild-richness-quadratic-variable-tests.csv"
)
turning_points_path <- file.path(
  out_dir,
  "atlas3-guild-richness-quadratic-turning-points.csv"
)
predictions_path <- file.path(
  out_dir,
  "atlas3-guild-richness-quadratic-predictions.csv"
)
diagnostics_path <- file.path(
  out_dir,
  "atlas3-guild-richness-quadratic-diagnostics.csv"
)
model_comparison_path <- file.path(
  out_dir,
  "atlas3-guild-richness-linear-vs-quadratic-comparison.csv"
)
relationships_png_path <- file.path(
  out_dir,
  "atlas3-guild-richness-quadratic-habitat-relationships.png"
)

#### CONSTANTS ####

expected_atlas3_sites <- 1465
expected_complete_sites <- 1462

predictor_terms <- c(
  "hh",
  "forest_percent",
  "agriculture_percent"
)

predictor_labels <- c(
  "hh" = "Habitat heterogeneity",
  "forest_percent" = "Forest cover (%)",
  "agriculture_percent" = "Agricultural cover (%)"
)

predictor_colours <- c(
  "Habitat heterogeneity" = "#8c510a",
  "Forest cover (%)" = "#1b7837",
  "Agricultural cover (%)" = "#2166ac"
)

centred_terms <- c(
  "hh" = "hh_c",
  "forest_percent" = "forest_c",
  "agriculture_percent" = "agriculture_c"
)

quadratic_formula <- guild_richness ~
  hh_c + I(hh_c^2) +
  forest_c + I(forest_c^2) +
  agriculture_c + I(agriculture_c^2)

linear_formula <- guild_richness ~
  hh_c +
  forest_c +
  agriculture_c

#### HELPERS ####

# Convert spatialreg coefficient output to a tidy table matching broom::tidy().
tidy_spatial_error <- function(model, model_name) {
  coefficient_matrix <- summary(model)$Coef

  data.frame(
    model = model_name,
    term = rownames(coefficient_matrix),
    estimate = coefficient_matrix[, "Estimate"],
    std.error = coefficient_matrix[, "Std. Error"],
    statistic = coefficient_matrix[, "z value"],
    p.value = coefficient_matrix[, "Pr(>|z|)"],
    conf.low = coefficient_matrix[, "Estimate"] -
      qnorm(0.975) * coefficient_matrix[, "Std. Error"],
    conf.high = coefficient_matrix[, "Estimate"] +
      qnorm(0.975) * coefficient_matrix[, "Std. Error"],
    row.names = NULL,
    check.names = FALSE
  )
}

# Build adjusted quadratic curves. The focal variable spans the central 98% of
# its observed values; the other two variables are held at their observed means.
make_polynomial_predictions <- function(
    model,
    data,
    variable_centres,
    model_name,
    spatial = FALSE) {
  bind_rows(lapply(predictor_terms, function(term) {
    prediction_data <- data.frame(
      hh_c = 0,
      forest_c = 0,
      agriculture_c = 0
    )
    prediction_data <- prediction_data[rep(1, 150), , drop = FALSE]

    x <- seq(
      quantile(data[[term]], 0.01, na.rm = TRUE),
      quantile(data[[term]], 0.99, na.rm = TRUE),
      length.out = 150
    )
    prediction_data[[centred_terms[[term]]]] <- x - variable_centres[[term]]

    if (spatial) {
      design_matrix <- model.matrix(
        delete.response(terms(quadratic_formula)),
        prediction_data
      )
      beta <- coef(model)
      beta <- beta[names(beta) != "lambda"]
      beta_vcov <- vcov(model)[names(beta), names(beta), drop = FALSE]
      fit <- as.numeric(design_matrix[, names(beta), drop = FALSE] %*% beta)
      se_fit <- sqrt(rowSums(
        (design_matrix[, names(beta), drop = FALSE] %*% beta_vcov) *
          design_matrix[, names(beta), drop = FALSE]
      ))
    } else {
      prediction <- predict(
        model,
        newdata = prediction_data,
        se.fit = TRUE
      )
      fit <- as.numeric(prediction$fit)
      se_fit <- as.numeric(prediction$se.fit)
    }

    data.frame(
      model = model_name,
      term = term,
      variable = predictor_labels[[term]],
      x = x,
      fit = fit,
      conf.low = fit - qnorm(0.975) * se_fit,
      conf.high = fit + qnorm(0.975) * se_fit
    )
  }))
}

# Calculate each quadratic turning point, returning it to the original scale.
calculate_turning_points <- function(model, model_name, centres, data) {
  beta <- coef(model)

  bind_rows(lapply(predictor_terms, function(term) {
    linear_term <- centred_terms[[term]]
    squared_term <- paste0("I(", linear_term, "^2)")
    quadratic_coefficient <- unname(beta[[squared_term]])
    turning_point <- centres[[term]] -
      unname(beta[[linear_term]]) / (2 * quadratic_coefficient)

    data.frame(
      model = model_name,
      term = term,
      variable = predictor_labels[[term]],
      quadratic_coefficient = quadratic_coefficient,
      curvature = ifelse(
        quadratic_coefficient < 0,
        "concave: local maximum",
        "convex: local minimum"
      ),
      turning_point = turning_point,
      observed_min = min(data[[term]], na.rm = TRUE),
      observed_max = max(data[[term]], na.rm = TRUE),
      turning_point_within_observed_range =
        turning_point >= min(data[[term]], na.rm = TRUE) &&
        turning_point <= max(data[[term]], na.rm = TRUE)
    )
  }))
}

# Compare the full quadratic model with models omitting both polynomial terms
# for one habitat variable. This tests each variable as a complete curve.
calculate_variable_tests <- function(full_model, data) {
  full_r2 <- summary(full_model)$r.squared
  full_sse <- deviance(full_model)

  bind_rows(lapply(predictor_terms, function(term) {
    keep_terms <- setdiff(predictor_terms, term)
    reduced_formula <- as.formula(paste(
      "guild_richness ~",
      paste(
        unlist(lapply(keep_terms, function(kept_term) {
          centred <- centred_terms[[kept_term]]
          c(centred, paste0("I(", centred, "^2)"))
        })),
        collapse = " + "
      )
    ))
    reduced_model <- lm(reduced_formula, data = data)
    comparison <- anova(reduced_model, full_model)
    reduced_sse <- deviance(reduced_model)

    data.frame(
      term = term,
      variable = predictor_labels[[term]],
      degrees_of_freedom = comparison$Df[2],
      F_statistic = comparison$F[2],
      p_value = comparison$`Pr(>F)`[2],
      unique_r_squared = full_r2 - summary(reduced_model)$r.squared,
      partial_r_squared = (reduced_sse - full_sse) / reduced_sse
    )
  }))
}

#### LOAD 2010s GUILD RICHNESS ####

guild_richness <- read.csv(
  guild_richness_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  filter(.data$period == "2010s")

if (nrow(guild_richness) != expected_atlas3_sites) {
  stop(
    "Expected ",
    expected_atlas3_sites,
    " Atlas 3 rows, but found ",
    nrow(guild_richness),
    ".",
    call. = FALSE
  )
}

if (anyDuplicated(guild_richness$site)) {
  stop("Atlas 3 guild-richness data contain duplicated sites.", call. = FALSE)
}

#### LOAD AND PREPARE 2010s HABITAT DATA ####

environment_raw <- read.csv(
  environment_path,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
environment_raw$survey <- rownames(environment_raw)

environment <- environment_raw |>
  clean_lulc_columns() |>
  mutate(
    perc_agriculture = .data$perc_cropland + .data$perc_pasture,
    forest_percent = 100 * .data$perc_forest,
    agriculture_percent = 100 * .data$perc_agriculture
  ) |>
  select(
    "survey",
    "hh",
    "forest_percent",
    "agriculture_percent"
  )

analysis_data_all <- guild_richness |>
  select(
    "survey",
    "site",
    "guild_richness",
    "species_richness"
  ) |>
  left_join(environment, by = "survey")

missing_habitat_sites <- analysis_data_all |>
  filter(if_any(all_of(predictor_terms), is.na)) |>
  pull(.data$site)

analysis_data <- analysis_data_all |>
  drop_na(all_of(predictor_terms))

if (nrow(analysis_data) != expected_complete_sites) {
  stop(
    "Expected ",
    expected_complete_sites,
    " complete Atlas 3 sites, but found ",
    nrow(analysis_data),
    ".",
    call. = FALSE
  )
}

variable_centres <- vapply(
  analysis_data[predictor_terms],
  mean,
  numeric(1),
  na.rm = TRUE
)

analysis_data <- analysis_data |>
  mutate(
    hh_c = .data$hh - variable_centres[["hh"]],
    forest_c = .data$forest_percent - variable_centres[["forest_percent"]],
    agriculture_c =
      .data$agriculture_percent - variable_centres[["agriculture_percent"]]
  )

#### FIT LINEAR AND QUADRATIC ORDINARY MODELS ####

linear_model <- lm(linear_formula, data = analysis_data)
quadratic_model <- lm(quadratic_formula, data = analysis_data)

linear_vs_quadratic <- anova(linear_model, quadratic_model)

model_comparison <- data.frame(
  comparison = "Linear habitat model vs quadratic habitat model",
  linear_AIC = AIC(linear_model),
  quadratic_AIC = AIC(quadratic_model),
  AIC_improvement = AIC(linear_model) - AIC(quadratic_model),
  additional_degrees_of_freedom = linear_vs_quadratic$Df[2],
  F_statistic = linear_vs_quadratic$F[2],
  p_value = linear_vs_quadratic$`Pr(>F)`[2]
)

#### FIT SPATIAL-ERROR VERSIONS ####

grid_sf <- st_read(grid_path, quiet = TRUE) |>
  filter(.data$kvadratkod %in% analysis_data$site) |>
  left_join(analysis_data, by = c("kvadratkod" = "site"))

if (nrow(grid_sf) != nrow(analysis_data) || anyNA(grid_sf$guild_richness)) {
  stop("Not all analysis sites joined to the atlas grid.", call. = FALSE)
}

neighbours <- suppressWarnings(poly2nb(grid_sf, queen = TRUE))
weights <- nb2listw(neighbours, style = "W", zero.policy = TRUE)
grid_data <- st_drop_geometry(grid_sf)

linear_spatial_model <- errorsarlm(
  linear_formula,
  data = grid_data,
  listw = weights,
  zero.policy = TRUE,
  method = "Matrix"
)

quadratic_spatial_model <- errorsarlm(
  quadratic_formula,
  data = grid_data,
  listw = weights,
  zero.policy = TRUE,
  method = "Matrix"
)

linear_residual_moran <- moran.test(
  residuals(lm(linear_formula, data = grid_data)),
  weights,
  zero.policy = TRUE
)
quadratic_residual_moran <- moran.test(
  residuals(lm(quadratic_formula, data = grid_data)),
  weights,
  zero.policy = TRUE
)
quadratic_spatial_residual_moran <- moran.test(
  residuals(quadratic_spatial_model),
  weights,
  zero.policy = TRUE
)

quadratic_spatial_lambda_test <- LR1.Sarlm(quadratic_spatial_model)

#### SUMMARISE MODELS ####

coefficient_output <- bind_rows(
  broom::tidy(quadratic_model, conf.int = TRUE) |>
    mutate(model = "Quadratic linear model"),
  tidy_spatial_error(
    quadratic_spatial_model,
    "Quadratic spatial-error model"
  )
) |>
  select(
    "model",
    "term",
    "estimate",
    "std.error",
    "statistic",
    "p.value",
    "conf.low",
    "conf.high"
  )

model_summary <- data.frame(
  model = c(
    "Linear model",
    "Quadratic model",
    "Linear spatial-error model",
    "Quadratic spatial-error model"
  ),
  n_sites = nrow(analysis_data),
  adjusted_r_squared = c(
    summary(linear_model)$adj.r.squared,
    summary(quadratic_model)$adj.r.squared,
    NA_real_,
    NA_real_
  ),
  residual_sd = c(
    summary(linear_model)$sigma,
    summary(quadratic_model)$sigma,
    sqrt(linear_spatial_model$s2),
    sqrt(quadratic_spatial_model$s2)
  ),
  AIC = c(
    AIC(linear_model),
    AIC(quadratic_model),
    AIC(linear_spatial_model),
    AIC(quadratic_spatial_model)
  ),
  spatial_lambda = c(
    NA_real_,
    NA_real_,
    linear_spatial_model$lambda,
    quadratic_spatial_model$lambda
  )
)

variable_tests <- calculate_variable_tests(quadratic_model, analysis_data)

turning_points <- bind_rows(
  calculate_turning_points(
    quadratic_model,
    "Quadratic linear model",
    variable_centres,
    analysis_data
  ),
  calculate_turning_points(
    quadratic_spatial_model,
    "Quadratic spatial-error model",
    variable_centres,
    analysis_data
  )
)

vif_output <- data.frame(
  term = names(car::vif(quadratic_model)),
  vif = as.numeric(car::vif(quadratic_model))
)

diagnostics <- bind_rows(
  vif_output |>
    transmute(
      diagnostic = "variance_inflation_factor",
      item = .data$term,
      value = .data$vif
    ),
  data.frame(
    diagnostic = c(
      "linear_model_residual_Moran_I",
      "linear_model_residual_Moran_p_value",
      "quadratic_model_residual_Moran_I",
      "quadratic_model_residual_Moran_p_value",
      "quadratic_spatial_residual_Moran_I",
      "quadratic_spatial_residual_Moran_p_value",
      "quadratic_spatial_lambda",
      "quadratic_spatial_lambda_p_value",
      "neighbour_graph_components",
      "isolated_grid_cells"
    ),
    item = NA_character_,
    value = c(
      unname(linear_residual_moran$estimate[["Moran I statistic"]]),
      linear_residual_moran$p.value,
      unname(quadratic_residual_moran$estimate[["Moran I statistic"]]),
      quadratic_residual_moran$p.value,
      unname(
        quadratic_spatial_residual_moran$estimate[["Moran I statistic"]]
      ),
      quadratic_spatial_residual_moran$p.value,
      quadratic_spatial_model$lambda,
      quadratic_spatial_lambda_test$p.value,
      n.comp.nb(neighbours)$nc,
      sum(card(neighbours) == 0)
    )
  )
)

#### CALCULATE ADJUSTED QUADRATIC CURVES ####

polynomial_predictions <- bind_rows(
  make_polynomial_predictions(
    quadratic_model,
    analysis_data,
    variable_centres,
    "Quadratic linear model",
    spatial = FALSE
  ),
  make_polynomial_predictions(
    quadratic_spatial_model,
    analysis_data,
    variable_centres,
    "Quadratic spatial-error model",
    spatial = TRUE
  )
)

#### WRITE TABLE OUTPUTS ####

write.csv(analysis_data, analysis_data_path, row.names = FALSE, na = "")
write.csv(coefficient_output, coefficient_path, row.names = FALSE, na = "")
write.csv(model_summary, model_summary_path, row.names = FALSE, na = "")
write.csv(variable_tests, variable_tests_path, row.names = FALSE, na = "")
write.csv(turning_points, turning_points_path, row.names = FALSE, na = "")
write.csv(polynomial_predictions, predictions_path, row.names = FALSE, na = "")
write.csv(diagnostics, diagnostics_path, row.names = FALSE, na = "")
write.csv(model_comparison, model_comparison_path, row.names = FALSE, na = "")

#### PNG: ADJUSTED QUADRATIC HABITAT RELATIONSHIPS ####

raw_relationship_data <- analysis_data |>
  select("guild_richness", all_of(predictor_terms)) |>
  pivot_longer(
    cols = all_of(predictor_terms),
    names_to = "term",
    values_to = "x"
  ) |>
  mutate(variable = unname(predictor_labels[.data$term]))

spatial_predictions <- polynomial_predictions |>
  filter(.data$model == "Quadratic spatial-error model")

relationship_plot <- ggplot(
  raw_relationship_data,
  aes(x = .data$x, y = .data$guild_richness)
) +
  geom_point(
    color = "grey35",
    alpha = 0.14,
    size = 0.75
  ) +
  geom_ribbon(
    data = spatial_predictions,
    aes(
      x = .data$x,
      ymin = .data$conf.low,
      ymax = .data$conf.high,
      fill = .data$variable
    ),
    inherit.aes = FALSE,
    alpha = 0.18,
    color = NA
  ) +
  geom_line(
    data = polynomial_predictions,
    aes(
      x = .data$x,
      y = .data$fit,
      color = .data$variable,
      linetype = .data$model
    ),
    inherit.aes = FALSE,
    linewidth = 1.05
  ) +
  facet_wrap(~variable, scales = "free_x") +
  scale_color_manual(values = predictor_colours, guide = "none") +
  scale_fill_manual(values = predictor_colours, guide = "none") +
  scale_linetype_manual(
    values = c(
      "Quadratic linear model" = "dashed",
      "Quadratic spatial-error model" = "solid"
    ),
    name = NULL
  ) +
  labs(
    title = "Quadratic habitat associations with 2010s guild richness",
    subtitle = paste0(
      "Adjusted second-order polynomial curves at ",
      nrow(analysis_data),
      " atlas cells; shaded intervals show the spatial-error model"
    ),
    x = NULL,
    y = "Number of foraging guilds represented"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = relationships_png_path,
  plot = relationship_plot,
  width = 11,
  height = 4.3,
  units = "in",
  dpi = 300,
  bg = "white"
)

#### CONSOLE SUMMARY ####

message("Complete Atlas 3 sites: ", nrow(analysis_data))
message(
  "Sites omitted because habitat data were missing: ",
  paste(missing_habitat_sites, collapse = ", ")
)
message(
  "Quadratic-model adjusted R-squared: ",
  format(round(summary(quadratic_model)$adj.r.squared, 4), nsmall = 4)
)
message(
  "Quadratic AIC improvement over linear model: ",
  format(round(model_comparison$AIC_improvement, 1), nsmall = 1),
  " (p = ",
  format.pval(model_comparison$p_value, digits = 3),
  ")"
)
message(
  "Quadratic spatial residual Moran's I: ",
  format(
    round(
      quadratic_spatial_residual_moran$estimate[["Moran I statistic"]],
      4
    ),
    nsmall = 4
  ),
  " (p = ",
  format.pval(quadratic_spatial_residual_moran$p.value, digits = 3),
  ")"
)
message("Wrote quadratic coefficients: ", coefficient_path)
message("Wrote model summary: ", model_summary_path)
message("Wrote variable tests: ", variable_tests_path)
message("Wrote turning points: ", turning_points_path)
message("Wrote polynomial predictions: ", predictions_path)
message("Wrote diagnostics: ", diagnostics_path)
message("Wrote model comparison: ", model_comparison_path)
message("Wrote quadratic relationship PNG: ", relationships_png_path)
message("Estimated turning points:")
print(turning_points)
