rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  GWmodel, sf, sp, tidyverse, Hmsc, vegan,
  terra, cowplot, RColorBrewer, scales
)

source(file.path("support_scripts", "figure_data_helpers.R"))

base_dir <- "./HmscOutputs"
pattern <- "2026-03-13"
bandwidth_m <- 100000
turnover_neighbour_m <- 200000
min_species_per_group <- 5
focal_atlases <- c("1", "3")
year_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")

out_dir <- file.path("notebooks", "exploratory", "outputs", "trait-change-turnover-gwr")
model_dir <- file.path(out_dir, "models")
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

grid_shape_path <- path.expand("~/box/PhD/logistics/data/distributions/DK5km_ED50grid_approx_kvadrkod_DOF/DK5km_ED50grid_approx_kvadrkod_DOF.shp")
shape_sf <- if (file.exists(grid_shape_path)) st_as_sf(vect(grid_shape_path)) else NULL

matching_folders <- figure_model_folders(pattern = pattern, base_dir = base_dir)
models_nums <- as.character(atlas_numbers(matching_folders))
mods <- load_hmsc_posteriors(matching_folders, base_dir = base_dir)
designs <- load_hmsc_study_designs(mods)
predsY <- load_or_compute_site_predictions(mods, matching_folders, base_dir = base_dir)
names(mods) <- names(designs) <- names(predsY) <- models_nums

mods <- mods[focal_atlases]
designs <- designs[focal_atlases]
predsY <- predsY[focal_atlases]

load("Data/preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData")
traits <- Tr |>
  rownames_to_column("species") |>
  select(species, Migration_a3_DOF, foraging_guild_consensus, species_thermal_index)

get_base_site <- function(x) sub("_[123]$", "", as.character(x))

thermal_groups <- function(sti) {
  breaks <- quantile(sti, c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE)
  cut(
    sti,
    breaks = breaks,
    include.lowest = TRUE,
    labels = c("very cold", "cold", "medium", "warm", "very warm")
  )
}

group_species_table <- function(traits, trait_col, family, min_species = 1) {
  traits |>
    filter(!is.na(.data[[trait_col]]), .data[[trait_col]] != "") |>
    group_by(trait_family = family, trait_value = .data[[trait_col]]) |>
    summarise(species = list(species), n_species = n(), .groups = "drop") |>
    filter(n_species >= min_species)
}

thermal_species <- traits |>
  filter(!is.na(species_thermal_index)) |>
  mutate(thermal_group = thermal_groups(species_thermal_index)) |>
  group_by(trait_family = "Thermal", trait_value = paste0("Thermal group: ", thermal_group)) |>
  summarise(species = list(species), n_species = n(), .groups = "drop")

migration_species_all <- group_species_table(traits, "Migration_a3_DOF", "Migration", min_species = 1)
guild_species_all <- group_species_table(traits, "foraging_guild_consensus", "Foraging guild", min_species = 1)

migration_species <- migration_species_all |> filter(n_species >= min_species_per_group)
guild_species <- guild_species_all |> filter(n_species >= min_species_per_group)

trait_groups_all <- bind_rows(thermal_species, migration_species_all, guild_species_all)
trait_groups_modelled <- bind_rows(
  tibble(trait_family = "Thermal", trait_value = "CWM STI", species = list(traits$species), n_species = nrow(traits)),
  thermal_species,
  migration_species,
  guild_species
)

write_csv(
  trait_groups_modelled |> select(trait_family, trait_value, n_species),
  file.path(out_dir, paste0(pattern, "-trait-groups-modelled.csv"))
)

trait_probability_surface <- function(pred_y, design, trait_groups, traits) {
  model_species <- colnames(pred_y)
  denominator_species <- intersect(model_species, traits$species)
  total_richness <- rowSums(pred_y[, denominator_species, drop = FALSE], na.rm = TRUE)

  pmap_dfr(trait_groups, function(trait_family, trait_value, species, n_species) {
    numerator_species <- intersect(model_species, species)
    expected_group_richness <- rowSums(pred_y[, numerator_species, drop = FALSE], na.rm = TRUE)

    tibble(
      survey = rownames(pred_y),
      trait_family = trait_family,
      trait_value = trait_value,
      n_species = n_species,
      trait_probability = ifelse(total_richness > 0, expected_group_richness / total_richness, NA_real_),
      expected_group_richness = expected_group_richness,
      expected_total_richness = total_richness
    )
  }) |>
    left_join(design, by = "survey")
}

cwm_sti_surface <- function(pred_y, design, traits) {
  sti <- traits$species_thermal_index
  names(sti) <- traits$species
  species <- intersect(colnames(pred_y), names(sti)[!is.na(sti)])
  richness <- rowSums(pred_y[, species, drop = FALSE], na.rm = TRUE)
  cwm <- as.numeric(pred_y[, species, drop = FALSE] %*% sti[species]) / richness

  tibble(
    survey = rownames(pred_y),
    trait_family = "Thermal",
    trait_value = "CWM STI",
    n_species = length(species),
    trait_probability = cwm,
    expected_group_richness = NA_real_,
    expected_total_richness = richness
  ) |>
    left_join(design, by = "survey")
}

trait_surfaces <- imap_dfr(predsY, function(pred_y, atlas) {
  bind_rows(
    cwm_sti_surface(pred_y, designs[[atlas]], traits),
    trait_probability_surface(pred_y, designs[[atlas]], bind_rows(thermal_species, migration_species, guild_species), traits)
  ) |>
    mutate(
      atlas = .env$atlas,
      period = year_lookup[[as.character(.env$atlas)]],
      base_site = get_base_site(site),
      .before = 1
    )
})

trait_surfaces_full_family <- imap_dfr(predsY, function(pred_y, atlas) {
  trait_probability_surface(pred_y, designs[[atlas]], trait_groups_all, traits) |>
    mutate(
      atlas = .env$atlas,
      period = year_lookup[[as.character(.env$atlas)]],
      base_site = get_base_site(site),
      .before = 1
    )
})

family_probability_checks <- trait_surfaces_full_family |>
  filter(trait_family %in% c("Migration", "Foraging guild")) |>
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

write_csv(family_probability_checks, file.path(out_dir, paste0(pattern, "-trait-family-probability-sum-checks.csv")))

make_trait_deltas <- function(surfaces) {
  start <- surfaces |> filter(atlas == "1")
  end <- surfaces |> filter(atlas == "3")

  joined <- inner_join(
    start |>
      select(trait_family, trait_value, base_site, survey_start = survey, site_start = site,
             X_start = X, Y_start = Y, value_start = trait_probability),
    end |>
      select(trait_family, trait_value, base_site, survey_end = survey, site_end = site,
             X = X, Y = Y, value_end = trait_probability),
    by = c("trait_family", "trait_value", "base_site")
  ) |>
    mutate(
      delta = value_end - value_start,
      direction = ifelse(delta >= 0, "Increase", "Decrease")
    )

  scale_lookup <- joined |>
    group_by(trait_family, trait_value) |>
    summarise(delta_sd = sd(delta, na.rm = TRUE), .groups = "drop") |>
    mutate(delta_sd = ifelse(is.na(delta_sd) | delta_sd <= 0, NA_real_, delta_sd))

  zero_variance <- scale_lookup |> filter(is.na(delta_sd))
  if (nrow(zero_variance) > 0) {
    warning("Dropping zero-variance trait changes: ", paste(zero_variance$trait_value, collapse = ", "))
  }

  joined |>
    left_join(scale_lookup, by = c("trait_family", "trait_value")) |>
    filter(!is.na(delta_sd)) |>
    mutate(
      standardized_change = delta / delta_sd,
      abs_standardized_change = abs(standardized_change)
    )
}

trait_deltas <- make_trait_deltas(trait_surfaces)

biggest_change_winners <- trait_deltas |>
  group_by(base_site) |>
  slice_max(abs_standardized_change, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    winner_label = paste(trait_family, trait_value, sep = ": "),
    signed_winner_change = ifelse(direction == "Increase", abs_standardized_change, -abs_standardized_change)
  )

write_csv(trait_surfaces, file.path(out_dir, paste0(pattern, "-trait-surfaces-1970s-2010s.csv")))
write_csv(trait_deltas, file.path(out_dir, paste0(pattern, "-trait-deltas-2010s-minus-1970s.csv")))
write_csv(biggest_change_winners, file.path(out_dir, paste0(pattern, "-biggest-trait-change-winners.csv")))

predicted_richness <- predicted_richness_frames(predsY, designs)

site_beta_turnover <- function(pred_y, design, distance_m = 200000) {
  design <- design[match(rownames(pred_y), design$survey), ]
  coords <- as.matrix(design[, c("X", "Y")])
  community_dissimilarity <- as.matrix(vegan::vegdist(pred_y, method = "bray"))
  geographic_distance <- as.matrix(dist(coords))
  diag(geographic_distance) <- NA_real_

  neighbour_mask <- geographic_distance <= distance_m
  neighbour_mask[is.na(neighbour_mask)] <- FALSE

  tibble(
    survey = rownames(pred_y),
    beta_turnover = map_dbl(seq_len(nrow(pred_y)), function(i) {
      neighbours <- neighbour_mask[i, ]
      if (!any(neighbours)) return(NA_real_)
      mean(community_dissimilarity[i, neighbours], na.rm = TRUE)
    })
  ) |>
    left_join(design, by = "survey")
}

beta_turnover <- map2(predsY, designs, site_beta_turnover, distance_m = turnover_neighbour_m)

responses <- bind_rows(
  imap_dfr(predicted_richness, ~ .x |> mutate(atlas = .y, period = year_lookup[[as.character(.y)]], response = "richness", value = richness)),
  imap_dfr(beta_turnover, ~ .x |> mutate(atlas = .y, period = year_lookup[[as.character(.y)]], response = "beta_turnover", value = beta_turnover))
) |>
  mutate(base_site = get_base_site(site)) |>
  select(atlas, period, response, survey, site, base_site, X, Y, value)

fit_single_trait_gwr <- function(df, response_col = "value", predictor_col = "trait_probability") {
  working <- df |>
    select(survey, site, X, Y, response = all_of(response_col), predictor = all_of(predictor_col)) |>
    drop_na()

  if (nrow(working) < 40 || sd(working$predictor, na.rm = TRUE) <= 0) {
    return(NULL)
  }

  working <- working |>
    mutate(trait_probability = as.numeric(scale(predictor)))
  metadata <- working |>
    select(survey, site)
  coordinates(working) <- c("X", "Y")

  formula <- response ~ trait_probability
  model <- gwr.basic(
    formula,
    data = working,
    bw = bandwidth_m,
    kernel = "bisquare",
    adaptive = FALSE
  )

  list(mod = model, bandwidth = bandwidth_m, n_rows = nrow(working), metadata = metadata)
}

safe_fit_single_trait_gwr <- function(df) {
  tryCatch(
    fit_single_trait_gwr(df),
    error = function(e) {
      warning(
        "GWR failed for ",
        paste(unique(df$response), unique(df$atlas), unique(df$trait_family), unique(df$trait_value), sep = " / "),
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
}

gwr_input <- responses |>
  inner_join(
    trait_surfaces |> select(atlas, survey, trait_family, trait_value, trait_probability),
    by = c("atlas", "survey"),
    relationship = "many-to-many"
  )

gwr_splits <- gwr_input |>
  group_by(response, atlas, period, trait_family, trait_value) |>
  group_split()

gwr_key <- function(df) {
  paste(unique(df$response), unique(df$atlas), unique(df$trait_family), unique(df$trait_value), sep = "__")
}
names(gwr_splits) <- map_chr(gwr_splits, gwr_key)

gwr_models <- map(gwr_splits, safe_fit_single_trait_gwr)

summarise_single_gwr <- function(fit, df) {
  if (is.null(fit)) {
    return(tibble(
      response = unique(df$response),
      atlas = unique(df$atlas),
      period = unique(df$period),
      trait_family = unique(df$trait_family),
      trait_value = unique(df$trait_value),
      n_rows = nrow(df),
      bandwidth_m = bandwidth_m,
      aicc = NA_real_,
      rss = NA_real_,
      local_r2_mean = NA_real_,
      local_r2_median = NA_real_,
      local_t_abs_max = NA_real_,
      status = "failed"
    ))
  }

  sdf <- as.data.frame(fit$mod$SDF)
  diagnostics <- fit$mod$GW.diagnostic
  t_col <- "trait_probability_TV"

  tibble(
    response = unique(df$response),
    atlas = unique(df$atlas),
    period = unique(df$period),
    trait_family = unique(df$trait_family),
    trait_value = unique(df$trait_value),
    n_rows = fit$n_rows,
    bandwidth_m = fit$bandwidth,
    aicc = diagnostics$AICc,
    rss = diagnostics$RSS.gw,
    local_r2_mean = mean(sdf$Local_R2, na.rm = TRUE),
    local_r2_median = median(sdf$Local_R2, na.rm = TRUE),
    local_t_abs_max = max(abs(sdf[[t_col]]), na.rm = TRUE),
    status = "ok"
  )
}

gwr_summary <- map2_dfr(gwr_models, gwr_splits, summarise_single_gwr)

extract_single_gwr_local <- function(fit, df) {
  if (is.null(fit)) return(tibble())

  sdf <- as.data.frame(fit$mod$SDF)
  tibble(
    survey = fit$metadata$survey,
    site = fit$metadata$site,
    X = coordinates(fit$mod$SDF)[, 1],
    Y = coordinates(fit$mod$SDF)[, 2],
    coefficient = sdf$trait_probability,
    local_t = sdf$trait_probability_TV,
    abs_local_t = abs(local_t),
    local_r2 = sdf$Local_R2,
    response = unique(df$response),
    atlas = unique(df$atlas),
    period = unique(df$period),
    trait_family = unique(df$trait_family),
    trait_value = unique(df$trait_value),
    direction = ifelse(local_t >= 0, "Positive", "Negative")
  )
}

gwr_local <- map2_dfr(gwr_models, gwr_splits, extract_single_gwr_local)

gwr_dominant <- gwr_local |>
  group_by(response, atlas, period, trait_family, survey) |>
  slice_max(abs_local_t, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    supported = abs_local_t >= 1.96,
    dominant_trait_value = ifelse(supported, trait_value, "No supported trait predictor"),
    dominant_direction = ifelse(supported, direction, "Not supported")
  )

write_csv(responses, file.path(out_dir, paste0(pattern, "-richness-beta-turnover-responses.csv")))
write_csv(gwr_summary, file.path(out_dir, paste0(pattern, "-trait-probability-gwr-summary-100km.csv")))
write_csv(gwr_local, file.path(out_dir, paste0(pattern, "-trait-probability-gwr-local-results-100km.csv")))
write_csv(gwr_dominant, file.path(out_dir, paste0(pattern, "-trait-probability-gwr-dominant-local-predictors-100km.csv")))
saveRDS(gwr_models, file.path(model_dir, paste0(pattern, "-trait-probability-single-predictor-gwr-models-100km.rds")))

join_shape <- function(df, by_site = "site") {
  if (is.null(shape_sf)) return(NULL)
  shape_sf |> left_join(df, by = c("kvadratkod" = by_site))
}

map_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "white", colour = NA)
  )

plot_discrete_map <- function(df, fill_col, title, legend_title, palette = NULL) {
  if (is.null(palette)) {
    vals <- sort(unique(na.omit(df[[fill_col]])))
    palette <- setNames(hue_pal()(length(vals)), vals)
  }

  if (!is.null(shape_sf)) {
    plot_data <- join_shape(df, by_site = "site_end")
    ggplot(plot_data) +
      geom_sf(aes(fill = .data[[fill_col]]), colour = "grey35", linewidth = 0.08) +
      scale_fill_manual(values = palette, na.value = "transparent", name = legend_title) +
      labs(title = title) +
      map_theme
  } else {
    ggplot(df, aes(X, Y, colour = .data[[fill_col]])) +
      geom_point(size = 0.9) +
      coord_equal() +
      scale_colour_manual(values = palette, na.value = "transparent", name = legend_title) +
      labs(title = title) +
      map_theme
  }
}

plot_continuous_map <- function(df, value_col, title, legend_title) {
  extent <- max(abs(df[[value_col]]), na.rm = TRUE)
  if (!is.finite(extent) || extent == 0) extent <- 1

  if (!is.null(shape_sf)) {
    plot_data <- join_shape(df, by_site = "site_end")
    ggplot(plot_data) +
      geom_sf(aes(fill = .data[[value_col]]), colour = "grey35", linewidth = 0.08) +
      scale_fill_gradient2(
        low = "#2166ac",
        mid = "white",
        high = "#b2182b",
        midpoint = 0,
        limits = c(-extent, extent),
        name = legend_title,
        na.value = "transparent"
      ) +
      labs(title = title) +
      map_theme
  } else {
    ggplot(df, aes(X, Y, colour = .data[[value_col]])) +
      geom_point(size = 0.9) +
      coord_equal() +
      scale_colour_gradient2(
        low = "#2166ac",
        mid = "white",
        high = "#b2182b",
        midpoint = 0,
        limits = c(-extent, extent),
        name = legend_title,
        na.value = "transparent"
      ) +
      labs(title = title) +
      map_theme
  }
}

family_colours <- c("Thermal" = "#d73027", "Migration" = "#7F77DD", "Foraging guild" = "#1b9e77")

family_map <- plot_discrete_map(
  biggest_change_winners,
  fill_col = "trait_family",
  title = "Largest standardized trait change, 2010s minus 1970s",
  legend_title = "Winning trait family",
  palette = family_colours
)

trait_value_map <- plot_discrete_map(
  biggest_change_winners,
  fill_col = "trait_value",
  title = "Trait value with largest standardized change",
  legend_title = "Winning trait value"
)

signed_change_map <- plot_continuous_map(
  biggest_change_winners,
  value_col = "signed_winner_change",
  title = "Signed magnitude of largest standardized trait change",
  legend_title = "Signed |z change|"
)

ggsave(file.path(plot_dir, paste0(pattern, "-biggest-trait-change-family-map.png")), family_map, width = 7.5, height = 6.5, dpi = 300)
ggsave(file.path(plot_dir, paste0(pattern, "-biggest-trait-change-value-map.png")), trait_value_map, width = 9.5, height = 7.5, dpi = 300)
ggsave(file.path(plot_dir, paste0(pattern, "-biggest-trait-change-signed-magnitude-map.png")), signed_change_map, width = 7.5, height = 6.5, dpi = 300)

plot_gwr_dominant_panel <- function(response_name) {
  df <- gwr_dominant |>
    filter(response == response_name) |>
    mutate(
      dominant_trait_value = factor(dominant_trait_value),
      direction_shape = factor(dominant_direction, levels = c("Positive", "Negative", "Not supported"))
    )

  if (!is.null(shape_sf)) {
    plot_data <- shape_sf |>
      left_join(df, by = c("kvadratkod" = "site")) |>
      mutate(direction_shape = factor(dominant_direction, levels = c("Positive", "Negative", "Not supported")))

    fill_levels <- sort(unique(as.character(df$dominant_trait_value)))
    fill_values <- setNames(hue_pal()(length(fill_levels)), fill_levels)
    fill_values["No supported trait predictor"] <- "#d9d9d9"
    direction_points <- plot_data |>
      filter(direction_shape %in% c("Positive", "Negative")) |>
      st_point_on_surface()

    ggplot(plot_data) +
      geom_sf(aes(fill = dominant_trait_value), colour = "grey35", linewidth = 0.05) +
      geom_sf(data = direction_points, aes(shape = direction_shape), size = 0.45, colour = "black", stroke = 0.25) +
      facet_grid(trait_family ~ period) +
      scale_fill_manual(values = fill_values, name = "Dominant trait", na.value = "transparent") +
      scale_shape_manual(values = c("Positive" = 16, "Negative" = 4), name = "Direction") +
      labs(title = paste("Dominant trait-probability predictor of", recode(response_name, richness = "richness", beta_turnover = "local beta turnover"))) +
      map_theme +
      theme(legend.text = element_text(size = 7), strip.text = element_text(face = "bold"))
  } else {
    ggplot(df, aes(X, Y, colour = dominant_trait_value, shape = direction_shape)) +
      geom_point(size = 0.8) +
      coord_equal() +
      facet_grid(trait_family ~ period) +
      labs(title = paste("Dominant trait-probability predictor of", response_name), colour = "Dominant trait", shape = "Direction") +
      map_theme
  }
}

richness_gwr_map <- plot_gwr_dominant_panel("richness")
turnover_gwr_map <- plot_gwr_dominant_panel("beta_turnover")

ggsave(file.path(plot_dir, paste0(pattern, "-dominant-trait-probability-gwr-richness-map.png")), richness_gwr_map, width = 11, height = 8, dpi = 300)
ggsave(file.path(plot_dir, paste0(pattern, "-dominant-trait-probability-gwr-turnover-map.png")), turnover_gwr_map, width = 11, height = 8, dpi = 300)

message("Finished trait-change and trait-probability GWR analysis.")
message("Outputs written to: ", out_dir)
