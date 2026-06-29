rm(list = ls())

.libPaths(c("~/Rlibs", .libPaths()))

library(tidyverse)
library(Hmsc)
library(scales)
library(patchwork)

source(file.path("support_scripts", "figure_data_helpers.R"))

base_dir <- "HmscOutputs"
pattern <- "2026-03-13"
min_species_per_group <- 5
atlas_from <- "1"
atlas_to <- "3"
period_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")

predictor_labels <- c(
  "tmean_breeding" = "Temperature",
  "prec_breeding" = "Precipitation",
  "hh" = "Land-use heterogeneity",
  "perc_urban" = "Urban",
  "perc_cropland" = "Cropland",
  "perc_pasture" = "Pasture",
  "perc_forest" = "Forest",
  "perc_grass_shrub" = "Grass/Shrubland"
)

predictor_colours <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban" = "#4d4d4d",
  "Cropland" = "goldenrod1",
  "Pasture" = "darkorange",
  "Forest" = "springgreen4",
  "Grass/Shrubland" = "springgreen2"
)

out_dir <- file.path("notebooks", "exploratory", "outputs", "trait-probability-change")
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading model predictions and traits...")
model_folders <- figure_model_folders(pattern = pattern, base_dir = base_dir)
atlas_ids <- as.character(atlas_numbers(model_folders))
mods <- load_hmsc_posteriors(model_folders, base_dir = base_dir)
designs <- load_hmsc_study_designs(mods)
preds_y <- load_or_compute_site_predictions(mods, model_folders, base_dir = base_dir)
names(mods) <- names(designs) <- names(preds_y) <- atlas_ids

focal_atlases <- c(atlas_from, atlas_to)
mods <- mods[focal_atlases]
designs <- designs[focal_atlases]
preds_y <- preds_y[focal_atlases]

load("Data/preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData")
traits <- Tr |>
  as.data.frame() |>
  rownames_to_column("species") |>
  select(species, Migration_a3_DOF, foraging_guild_consensus)

message("Loading variance partitioning estimates for trait-group composition...")
vp_scaled <- load_vp_estimates(model_folders, base_dir = base_dir, scaled = TRUE)
vp_long <- imap_dfr(vp_scaled[focal_atlases], function(df, atlas) {
  df |>
    rownames_to_column("variable") |>
    pivot_longer(-variable, names_to = "species", values_to = "VP") |>
    mutate(atlas = as.character(atlas))
}) |>
  filter(variable %in% names(predictor_labels)) |>
  mutate(
    variable_clean = recode(variable, !!!predictor_labels),
    variable_clean = factor(variable_clean, levels = names(predictor_colours))
  )

get_base_site <- function(x) sub("_[123]$", "", as.character(x))

group_species_table <- function(traits, trait_col, trait_family, min_species = 1) {
  traits |>
    filter(!is.na(.data[[trait_col]]), .data[[trait_col]] != "") |>
    group_by(trait_family = trait_family, trait_value = .data[[trait_col]]) |>
    summarise(n_species = n_distinct(species), species = list(species), .groups = "drop") |>
    filter(n_species >= min_species)
}

trait_groups_all <- bind_rows(
  group_species_table(traits, "Migration_a3_DOF", "Migratory strategy", min_species = 1),
  group_species_table(traits, "foraging_guild_consensus", "Foraging guild", min_species = 1)
)

trait_groups_plot <- trait_groups_all |>
  filter(n_species >= min_species_per_group)

write_csv(
  trait_groups_all |> select(trait_family, trait_value, n_species),
  file.path(out_dir, paste0(pattern, "-trait-probability-change-groups-all.csv"))
)
write_csv(
  trait_groups_plot |> select(trait_family, trait_value, n_species),
  file.path(out_dir, paste0(pattern, "-trait-probability-change-groups-n-ge-", min_species_per_group, ".csv"))
)

trait_probability_surface <- function(pred_y, design, trait_groups, traits) {
  model_species <- colnames(pred_y)
  denominator_species <- intersect(model_species, traits$species)
  expected_total_richness <- rowSums(pred_y[, denominator_species, drop = FALSE], na.rm = TRUE)

  pmap_dfr(trait_groups, function(trait_family, trait_value, species, n_species) {
    numerator_species <- intersect(model_species, species)
    expected_group_richness <- rowSums(pred_y[, numerator_species, drop = FALSE], na.rm = TRUE)

    tibble(
      survey = rownames(pred_y),
      trait_family = trait_family,
      trait_value = trait_value,
      n_species = n_species,
      trait_probability = ifelse(
        expected_total_richness > 0,
        expected_group_richness / expected_total_richness,
        NA_real_
      ),
      expected_group_richness = expected_group_richness,
      expected_total_richness = expected_total_richness
    )
  }) |>
    left_join(design, by = "survey")
}

message("Computing site-level trait probabilities...")
trait_surfaces_all <- imap_dfr(preds_y, function(pred_y, atlas) {
  trait_probability_surface(pred_y, designs[[atlas]], trait_groups_all, traits) |>
    mutate(
      atlas = as.character(.env$atlas),
      period = period_lookup[[as.character(.env$atlas)]],
      base_site = get_base_site(site),
      .before = 1
    )
})

probability_sum_checks <- trait_surfaces_all |>
  group_by(atlas, period, trait_family, survey, site) |>
  summarise(probability_sum = sum(trait_probability, na.rm = TRUE), .groups = "drop") |>
  group_by(atlas, period, trait_family) |>
  summarise(
    n_sites = n(),
    min_sum = min(probability_sum, na.rm = TRUE),
    mean_sum = mean(probability_sum, na.rm = TRUE),
    max_sum = max(probability_sum, na.rm = TRUE),
    max_abs_deviation_from_one = max(abs(probability_sum - 1), na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  probability_sum_checks,
  file.path(out_dir, paste0(pattern, "-trait-probability-change-sum-checks.csv"))
)

message("Computing ", period_lookup[[atlas_to]], " - ", period_lookup[[atlas_from]], " changes...")
trait_deltas_all <- trait_surfaces_all |>
  filter(atlas %in% focal_atlases) |>
  select(
    atlas, period, trait_family, trait_value, n_species, base_site, site, survey,
    trait_probability, expected_group_richness, expected_total_richness
  ) |>
  pivot_wider(
    id_cols = c(trait_family, trait_value, n_species, base_site),
    names_from = atlas,
    values_from = c(site, survey, period, trait_probability, expected_group_richness, expected_total_richness),
    names_sep = "_"
  ) |>
  filter(!is.na(.data[[paste0("trait_probability_", atlas_from)]]),
         !is.na(.data[[paste0("trait_probability_", atlas_to)]])) |>
  transmute(
    trait_family,
    trait_value,
    n_species,
    base_site,
    site_from = .data[[paste0("site_", atlas_from)]],
    site_to = .data[[paste0("site_", atlas_to)]],
    period_from = .data[[paste0("period_", atlas_from)]],
    period_to = .data[[paste0("period_", atlas_to)]],
    probability_from = .data[[paste0("trait_probability_", atlas_from)]],
    probability_to = .data[[paste0("trait_probability_", atlas_to)]],
    probability_delta = probability_to - probability_from,
    expected_group_richness_from = .data[[paste0("expected_group_richness_", atlas_from)]],
    expected_group_richness_to = .data[[paste0("expected_group_richness_", atlas_to)]],
    expected_group_richness_delta = expected_group_richness_to - expected_group_richness_from,
    expected_total_richness_from = .data[[paste0("expected_total_richness_", atlas_from)]],
    expected_total_richness_to = .data[[paste0("expected_total_richness_", atlas_to)]]
  )

trait_change_summary_all <- trait_deltas_all |>
  group_by(trait_family, trait_value, n_species) |>
  summarise(
    n_sites = n(),
    mean_probability_from = mean(probability_from, na.rm = TRUE),
    mean_probability_to = mean(probability_to, na.rm = TRUE),
    mean_probability_delta = mean(probability_delta, na.rm = TRUE),
    median_probability_delta = median(probability_delta, na.rm = TRUE),
    q10_probability_delta = quantile(probability_delta, 0.10, na.rm = TRUE),
    q25_probability_delta = quantile(probability_delta, 0.25, na.rm = TRUE),
    q75_probability_delta = quantile(probability_delta, 0.75, na.rm = TRUE),
    q90_probability_delta = quantile(probability_delta, 0.90, na.rm = TRUE),
    sd_probability_delta = sd(probability_delta, na.rm = TRUE),
    prop_sites_increase = mean(probability_delta > 0, na.rm = TRUE),
    mean_group_richness_delta = mean(expected_group_richness_delta, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(trait_family) |>
  mutate(
    rank_by_change = dense_rank(desc(mean_probability_delta)),
    direction = ifelse(mean_probability_delta >= 0, "Increase", "Decrease"),
    eligible_for_main_plot = n_species >= min_species_per_group
  ) |>
  ungroup()

write_csv(
  trait_surfaces_all,
  file.path(out_dir, paste0(pattern, "-trait-probability-surfaces-", period_lookup[[atlas_from]], "-", period_lookup[[atlas_to]], ".csv"))
)
write_csv(
  trait_deltas_all,
  file.path(out_dir, paste0(pattern, "-trait-probability-deltas-", period_lookup[[atlas_to]], "-minus-", period_lookup[[atlas_from]], ".csv"))
)
write_csv(
  trait_change_summary_all,
  file.path(out_dir, paste0(pattern, "-trait-probability-change-summary-", period_lookup[[atlas_to]], "-minus-", period_lookup[[atlas_from]], ".csv"))
)

trait_vp_composition <- bind_rows(
  traits |>
    select(species, trait_value = foraging_guild_consensus) |>
    mutate(trait_family = "Foraging guild"),
  traits |>
    select(species, trait_value = Migration_a3_DOF) |>
    mutate(trait_family = "Migratory strategy")
) |>
  filter(!is.na(trait_value), trait_value != "") |>
  inner_join(vp_long, by = "species", relationship = "many-to-many") |>
  group_by(trait_family, trait_value, variable_clean) |>
  summarise(
    mean_vp = mean(VP, na.rm = TRUE),
    median_vp = median(VP, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(trait_family, trait_value) |>
  mutate(
    total_mean_vp = sum(mean_vp, na.rm = TRUE),
    relative_vp = ifelse(total_mean_vp > 0, mean_vp / total_mean_vp, NA_real_)
  ) |>
  ungroup() |>
  left_join(
    trait_change_summary_all |> select(trait_family, trait_value, n_species, eligible_for_main_plot),
    by = c("trait_family", "trait_value")
  )

write_csv(
  trait_vp_composition,
  file.path(out_dir, paste0(pattern, "-trait-probability-change-vp-composition-", period_lookup[[atlas_to]], "-minus-", period_lookup[[atlas_from]], ".csv"))
)

plot_ranked_change <- function(summary_df, min_species = min_species_per_group) {
  plot_df <- summary_df |>
    filter(n_species >= min_species) |>
    group_by(trait_family) |>
    arrange(mean_probability_delta, .by_group = TRUE) |>
    mutate(
      trait_label = paste0(trait_value, " (n=", n_species, ")"),
      trait_label = factor(trait_label, levels = unique(trait_label))
    ) |>
    ungroup()

  ggplot(plot_df, aes(mean_probability_delta, trait_label, fill = direction)) +
    geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.35) +
    geom_errorbar(
      aes(xmin = q10_probability_delta, xmax = q90_probability_delta),
      orientation = "y",
      width = 0,
      colour = "grey45",
      linewidth = 0.45
    ) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
    facet_grid(trait_family ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
    scale_fill_manual(values = c("Increase" = "#1b9e77", "Decrease" = "#d95f02"), name = NULL) +
    labs(
      title = paste0("Trait-probability change: ", period_lookup[[atlas_to]], " minus ", period_lookup[[atlas_from]]),
      subtitle = paste0(
        "Bars show mean change across shared grid cells; grey lines are 10th-90th percentile spatial intervals.\nCategories with n >= ",
        min_species
      ),
      x = "Change in predicted trait probability",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.y = element_text(face = "bold", angle = 0),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey40", size = 10)
    )
}

p_ranked <- plot_ranked_change(trait_change_summary_all, min_species_per_group)
ggsave(
  file.path(plot_dir, paste0(pattern, "-trait-probability-change-ranked-", period_lookup[[atlas_to]], "-minus-", period_lookup[[atlas_from]], ".png")),
  p_ranked,
  width = 9,
  height = 9.5,
  dpi = 300
)

plot_change_with_vp_family <- function(summary_df, vp_df, family, min_species = min_species_per_group) {
  order_df <- summary_df |>
    filter(trait_family == family, n_species >= min_species) |>
    arrange(mean_probability_delta) |>
    mutate(
      trait_label = paste0(trait_value, " (n=", n_species, ")"),
      trait_label = factor(trait_label, levels = trait_label)
    )

  if (nrow(order_df) == 0) {
    stop("No eligible trait groups for ", family, call. = FALSE)
  }

  change_df <- order_df |>
    mutate(direction = factor(direction, levels = c("Decrease", "Increase")))

  vp_plot_df <- vp_df |>
    filter(trait_family == family, n_species >= min_species) |>
    mutate(trait_label = paste0(trait_value, " (n=", n_species, ")")) |>
    semi_join(order_df |> select(trait_label), by = "trait_label") |>
    mutate(
      trait_label = factor(trait_label, levels = levels(order_df$trait_label)),
      variable_clean = factor(variable_clean, levels = names(predictor_colours))
    )

  change_plot <- ggplot(change_df, aes(mean_probability_delta, trait_label, fill = direction)) +
    geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.35) +
    geom_errorbar(
      aes(xmin = q10_probability_delta, xmax = q90_probability_delta),
      orientation = "y",
      width = 0,
      colour = "grey45",
      linewidth = 0.45
    ) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
    scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
    scale_fill_manual(values = c("Increase" = "#1b9e77", "Decrease" = "#d95f02"), name = NULL) +
    labs(
      title = family,
      subtitle = "Probability change",
      x = paste0(period_lookup[[atlas_to]], " - ", period_lookup[[atlas_from]]),
      y = NULL
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(colour = "grey40", size = 9)
    )

  vp_plot <- ggplot(vp_plot_df, aes(relative_vp, trait_label, fill = variable_clean)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
    scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.01))) +
    scale_fill_manual(values = predictor_colours, name = NULL, drop = FALSE) +
    labs(
      subtitle = "Environmental VP composition",
      x = "Relative VP share",
      y = NULL
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.subtitle = element_text(colour = "grey40", size = 9)
    )

  change_plot + vp_plot + plot_layout(widths = c(1.25, 1), guides = "collect") &
    theme(legend.position = "bottom")
}

p_change_vp <- plot_change_with_vp_family(trait_change_summary_all, trait_vp_composition, "Foraging guild") /
  plot_change_with_vp_family(trait_change_summary_all, trait_vp_composition, "Migratory strategy") +
  plot_annotation(
    title = paste0("Trait-probability shifts and environmental VP composition: ", period_lookup[[atlas_to]], " minus ", period_lookup[[atlas_from]]),
    subtitle = paste0(
      "Left: mean probability change across shared grid cells with 10th-90th percentile spatial intervals. ",
      "Right: relative contribution of environmental variables to species-level VP within each group."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(colour = "grey40", size = 10)
    )
  )

ggsave(
  file.path(plot_dir, paste0(pattern, "-trait-probability-change-with-vp-composition-", period_lookup[[atlas_to]], "-minus-", period_lookup[[atlas_from]], ".png")),
  p_change_vp,
  width = 12,
  height = 12,
  dpi = 300
)

top_bottom_summary <- trait_change_summary_all |>
  filter(n_species >= min_species_per_group) |>
  group_by(trait_family) |>
  arrange(mean_probability_delta, .by_group = TRUE) |>
  summarise(
    strongest_decreases = paste0(head(trait_value, 5), " (", percent(head(mean_probability_delta, 5), accuracy = 0.1), ")", collapse = "; "),
    strongest_increases = paste0(rev(tail(trait_value, 5)), " (", percent(rev(tail(mean_probability_delta, 5)), accuracy = 0.1), ")", collapse = "; "),
    .groups = "drop"
  )
write_csv(
  top_bottom_summary,
  file.path(out_dir, paste0(pattern, "-trait-probability-change-top-bottom-", period_lookup[[atlas_to]], "-minus-", period_lookup[[atlas_from]], ".csv"))
)

print(probability_sum_checks)
print(top_bottom_summary)

message("Finished trait-probability change bars.")
