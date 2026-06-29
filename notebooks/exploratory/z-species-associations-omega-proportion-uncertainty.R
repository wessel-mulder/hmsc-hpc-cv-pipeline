#### Species-pair association-class proportion uncertainty ####
#
# Purpose:
#   Lightweight add-on to the Omega association notebook. This script propagates
#   posterior uncertainty into the aggregate proportions of species pairs
#   classified as supported positive, supported negative, or unsupported. It
#   also asks whether foraging-guild and migratory-strategy pairs show more or
#   less supported association than expected under a group-label permutation
#   null, and whether those group-pair signals persist through time.
#
# Important:
#   This script does not rerender the long R Markdown notebook. It reads the
#   fitted HMSC posterior objects directly, bootstraps posterior Omega draws,
#   writes compact CSV summaries, and saves PNG-only figures.

required_packages <- c("tidyverse", "Hmsc", "coda", "ggplot2", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install required packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(Hmsc)
  library(coda)
  library(ggplot2)
  library(scales)
})

source(file.path("support_scripts", "figure_data_helpers.R"))

base_dir <- "HmscOutputs"
pattern <- "2026-03-13"
fitted_file <- hmsc_fitted_file
support_level <- 0.95
negative_support_cutoff <- 0.05
n_bootstrap <- 1000
bootstrap_draws_per_replicate <- 1000
n_group_permutations <- 999
min_group_pair_count_for_plots <- 10
min_guild_species_for_plots <- 3
set.seed(20260521)

period_order <- c("1970s", "1990s", "2010s")
year_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")
support_levels <- c("Supported positive", "Supported negative", "Unsupported")

out_dir <- file.path(
  "notebooks", "exploratory", "outputs", "species-associations-omega"
)
plot_dir <- file.path(out_dir, "plots")
cache_dir <- file.path(out_dir, "posterior-draw-cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

existing_pair_path <- file.path(out_dir, paste0(pattern, "-omega-species-pairs.csv"))
if (!file.exists(existing_pair_path)) {
  stop(
    "Missing existing pair table: ", existing_pair_path,
    "\nRun the Omega association notebook once before this add-on script.",
    call. = FALSE
  )
}

existing_pairs <- read_csv(
  existing_pair_path,
  col_types = cols(
    atlas = col_character(),
    period = col_character(),
    support_class = col_character(),
    .default = col_guess()
  )
) |>
  mutate(
    atlas = as.character(.data$atlas),
    period = factor(.data$period, levels = period_order),
    support_class = factor(.data$support_class, levels = support_levels)
  )

model_metadata <- tibble(
  model_folder = figure_model_folders(pattern = pattern, base_dir = base_dir)
) |>
  mutate(
    atlas = as.character(atlas_numbers(.data$model_folder)),
    period = unname(year_lookup[.data$atlas])
  ) |>
  arrange(as.numeric(.data$atlas))

save_png <- function(plot, filename, width, height, dpi = 300) {
  # Project convention: plot files from this workflow must be PNG unless the
  # user explicitly asks otherwise.
  if (tolower(tools::file_ext(filename)) != "png") {
    stop("Plot filename must end in .png: ", filename, call. = FALSE)
  }

  path <- file.path(plot_dir, filename)
  ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  path
}

classify_support <- function(pr_positive) {
  case_when(
    pr_positive > support_level ~ "Supported positive",
    pr_positive < negative_support_cutoff ~ "Supported negative",
    TRUE ~ "Unsupported"
  )
}

load_or_cache_omega_pair_signs <- function(model_folder, atlas, period) {
  cache_path <- file.path(
    cache_dir,
    paste0(pattern, "-atlas", atlas, "-site-omega-pair-positive-draws.rds")
  )
  if (file.exists(cache_path)) {
    return(readRDS(cache_path))
  }

  model_path <- file.path(base_dir, model_folder, "Models", "Fitted", fitted_file)
  env <- new.env(parent = emptyenv())
  load(model_path, envir = env)
  if (!exists("fitted_model", envir = env, inherits = FALSE)) {
    stop("Expected `fitted_model` in ", model_path, call. = FALSE)
  }

  model <- env$fitted_model$posteriors
  if (!identical(names(model$ranLevels), "site")) {
    stop("Expected only the site random level in ", model_folder, call. = FALSE)
  }

  # convertToCodaObject returns one column for every cell of the Omega matrix.
  # Its columns are flattened in matrix column-major order. The lower triangle
  # therefore maps to pair IDs species[col]__species[row], matching the existing
  # notebook's upper-triangle species_a__species_b convention.
  omega_draw_matrix <- as.matrix(convertToCodaObject(model)$Omega[[1]])
  species_names <- model$spNames
  lower_index <- which(lower.tri(matrix(NA, nrow = model$ns, ncol = model$ns)), arr.ind = TRUE)
  lower_columns <- lower_index[, "row"] + (lower_index[, "col"] - 1) * model$ns
  pair_ids <- paste(
    species_names[lower_index[, "col"]],
    species_names[lower_index[, "row"]],
    sep = "__"
  )

  out <- list(
    atlas = atlas,
    period = period,
    pair_id = pair_ids,
    positive_draws = omega_draw_matrix[, lower_columns, drop = FALSE] > 0
  )
  saveRDS(out, cache_path)
  out
}

omega_draws_by_atlas <- pmap(
  model_metadata,
  function(model_folder, atlas, period) {
    load_or_cache_omega_pair_signs(
      model_folder = model_folder,
      atlas = atlas,
      period = period
    )
  }
)
names(omega_draws_by_atlas) <- model_metadata$atlas

validate_against_existing_summary <- function(draw_object) {
  existing_atlas <- existing_pairs |>
    filter(.data$atlas == draw_object$atlas) |>
    arrange(.data$pair_id)

  draw_order <- order(draw_object$pair_id)
  draw_pair_ids <- draw_object$pair_id[draw_order]
  if (!identical(existing_atlas$pair_id, draw_pair_ids)) {
    stop("Posterior draw pair IDs do not match existing pair table for atlas ", draw_object$atlas, call. = FALSE)
  }

  draw_pr_positive <- round(colMeans(draw_object$positive_draws[, draw_order, drop = FALSE]), 10)
  max_abs_difference <- max(abs(draw_pr_positive - existing_atlas$pr_positive), na.rm = TRUE)
  if (max_abs_difference > 1e-8) {
    stop(
      "Posterior draws do not reproduce existing Pr(Omega > 0) for atlas ",
      draw_object$atlas, ". Max absolute difference: ", max_abs_difference,
      call. = FALSE
    )
  }

  tibble(
    atlas = draw_object$atlas,
    period = draw_object$period,
    n_pairs = length(draw_pair_ids),
    n_draws = nrow(draw_object$positive_draws),
    max_abs_pr_positive_difference = max_abs_difference
  )
}

validation_summary <- map_dfr(omega_draws_by_atlas, validate_against_existing_summary)
write_csv(
  validation_summary,
  file.path(out_dir, paste0(pattern, "-omega-class-proportion-uncertainty-validation.csv"))
)

bootstrap_class_proportions <- function(draw_object) {
  positive_draws <- draw_object$positive_draws
  n_draws <- nrow(positive_draws)
  n_pair_total <- ncol(positive_draws)

  # Draw-row bootstrap weights let BLAS do the heavy matrix multiply in one
  # shot, which is much quicker than looping over 1,000 colMeans calls.
  bootstrap_weights <- t(rmultinom(
    n = n_bootstrap,
    size = bootstrap_draws_per_replicate,
    prob = rep(1 / n_draws, n_draws)
  ))
  pr_positive_boot <- (bootstrap_weights %*% (positive_draws + 0)) /
    bootstrap_draws_per_replicate

  n_supported_positive <- rowSums(pr_positive_boot > support_level)
  n_supported_negative <- rowSums(pr_positive_boot < negative_support_cutoff)
  n_unsupported <- n_pair_total - n_supported_positive - n_supported_negative

  tibble(
    bootstrap = rep(seq_len(n_bootstrap), each = length(support_levels)),
    atlas = draw_object$atlas,
    period = draw_object$period,
    support_class = factor(rep(support_levels, times = n_bootstrap), levels = support_levels),
    n_pairs = as.integer(c(rbind(
      n_supported_positive,
      n_supported_negative,
      n_unsupported
    ))),
    prop_pairs = .data$n_pairs / n_pair_total
  )
}

bootstrap_proportions <- map_dfr(omega_draws_by_atlas, bootstrap_class_proportions) |>
  mutate(
    atlas = as.character(.data$atlas),
    period = as.character(.data$period),
    support_class = as.character(.data$support_class)
  )

observed_proportions <- existing_pairs |>
  count(.data$atlas, .data$period, .data$support_class, name = "n_pairs") |>
  complete(
    nesting(atlas, period),
    support_class = factor(support_levels, levels = support_levels),
    fill = list(n_pairs = 0)
  ) |>
  group_by(.data$atlas, .data$period) |>
  mutate(prop_pairs = .data$n_pairs / sum(.data$n_pairs)) |>
  ungroup() |>
  mutate(
    atlas = as.character(.data$atlas),
    period = as.character(.data$period),
    support_class = as.character(.data$support_class)
  )

proportion_uncertainty <- bootstrap_proportions |>
  group_by(.data$atlas, .data$period, .data$support_class) |>
  summarise(
    mean = mean(.data$prop_pairs, na.rm = TRUE),
    median = median(.data$prop_pairs, na.rm = TRUE),
    lower_80 = quantile(.data$prop_pairs, 0.10, na.rm = TRUE),
    upper_80 = quantile(.data$prop_pairs, 0.90, na.rm = TRUE),
    lower_95 = quantile(.data$prop_pairs, 0.025, na.rm = TRUE),
    upper_95 = quantile(.data$prop_pairs, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    observed_proportions |>
      select(atlas, period, support_class, observed_n_pairs = n_pairs, observed_prop_pairs = prop_pairs),
    by = c("atlas", "period", "support_class")
  ) |>
  mutate(
    period = factor(.data$period, levels = period_order),
    support_class = factor(.data$support_class, levels = support_levels)
  ) |>
  arrange(as.numeric(.data$atlas), .data$support_class)

comparison_lookup <- tribble(
  ~comparison, ~atlas_to, ~atlas_from,
  "Atlas 2 - Atlas 1", "2", "1",
  "Atlas 3 - Atlas 1", "3", "1",
  "Atlas 3 - Atlas 2", "3", "2"
)

bootstrap_change_draws <- pmap_dfr(
  comparison_lookup,
  function(comparison, atlas_to, atlas_from) {
    to_df <- bootstrap_proportions |>
      filter(.data$atlas == atlas_to) |>
      select(bootstrap, support_class, prop_to = prop_pairs)
    from_df <- bootstrap_proportions |>
      filter(.data$atlas == atlas_from) |>
      select(bootstrap, support_class, prop_from = prop_pairs)

    inner_join(to_df, from_df, by = c("bootstrap", "support_class")) |>
      mutate(
        comparison = comparison,
        atlas_to = atlas_to,
        atlas_from = atlas_from,
        prop_change = .data$prop_to - .data$prop_from
      )
  }
)

observed_change <- pmap_dfr(
  comparison_lookup,
  function(comparison, atlas_to, atlas_from) {
    to_df <- observed_proportions |>
      filter(.data$atlas == atlas_to) |>
      select(support_class, observed_prop_to = prop_pairs)
    from_df <- observed_proportions |>
      filter(.data$atlas == atlas_from) |>
      select(support_class, observed_prop_from = prop_pairs)

    inner_join(to_df, from_df, by = "support_class") |>
      mutate(
        comparison = comparison,
        atlas_to = atlas_to,
        atlas_from = atlas_from,
        observed_prop_change = .data$observed_prop_to - .data$observed_prop_from
      )
  }
)

change_uncertainty <- bootstrap_change_draws |>
  group_by(.data$comparison, .data$atlas_to, .data$atlas_from, .data$support_class) |>
  summarise(
    mean = mean(.data$prop_change, na.rm = TRUE),
    median = median(.data$prop_change, na.rm = TRUE),
    lower_80 = quantile(.data$prop_change, 0.10, na.rm = TRUE),
    upper_80 = quantile(.data$prop_change, 0.90, na.rm = TRUE),
    lower_95 = quantile(.data$prop_change, 0.025, na.rm = TRUE),
    upper_95 = quantile(.data$prop_change, 0.975, na.rm = TRUE),
    pr_change_gt_0 = mean(.data$prop_change > 0, na.rm = TRUE),
    pr_change_lt_0 = mean(.data$prop_change < 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ungroup() |>
  left_join(
    observed_change |>
      select(comparison, atlas_to, atlas_from, support_class, observed_prop_change),
    by = c("comparison", "atlas_to", "atlas_from", "support_class")
  ) |>
  arrange(.data$comparison, .data$support_class)

write_csv(
  proportion_uncertainty,
  file.path(out_dir, paste0(pattern, "-omega-class-proportion-uncertainty.csv"))
)
write_csv(
  change_uncertainty,
  file.path(out_dir, paste0(pattern, "-omega-class-proportion-change-uncertainty.csv"))
)
write_csv(
  bootstrap_proportions,
  file.path(out_dir, paste0(pattern, "-omega-class-proportion-bootstrap-draws.csv"))
)

support_colours <- c(
  "Supported positive" = "#b2182b",
  "Supported negative" = "#2166ac",
  "Unsupported" = "grey55"
)

p_proportions <- ggplot(
  proportion_uncertainty,
  aes(x = .data$period, y = .data$observed_prop_pairs, colour = .data$support_class)
) +
  geom_errorbar(
    aes(
      ymin = .data$lower_95,
      ymax = .data$upper_95,
      group = .data$support_class
    ),
    width = 0.05,
    linewidth = 0.7,
    alpha = 0.35
  ) +
  geom_errorbar(
    aes(
      ymin = .data$lower_80,
      ymax = .data$upper_80,
      group = .data$support_class
    ),
    width = 0.12,
    linewidth = 1.1
  ) +
  geom_line(aes(group = .data$support_class), linewidth = 0.6) +
  geom_point(size = 2) +
  scale_colour_manual(values = support_colours, name = NULL) +
  scale_y_continuous(labels = percent_format(), limits = c(0, NA)) +
  labs(x = NULL, y = "Species-pair proportion") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

p_changes <- ggplot(
  change_uncertainty,
  aes(x = .data$observed_prop_change, y = .data$support_class, colour = .data$support_class)
) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = .data$lower_95, xmax = .data$upper_95),
    orientation = "y",
    width = 0,
    linewidth = 1,
    alpha = 0.45
  ) +
  geom_errorbar(
    aes(xmin = .data$lower_80, xmax = .data$upper_80),
    orientation = "y",
    width = 0,
    linewidth = 2
  ) +
  geom_point(size = 2.6) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_colour_manual(values = support_colours, guide = "none") +
  scale_x_continuous(labels = percent_format()) +
  labs(x = "Change in species-pair proportion", y = NULL) +
  theme_minimal(base_size = 11)

saved_plot_paths <- c(
  save_png(
    p_proportions,
    paste0(pattern, "-omega-class-proportion-uncertainty.png"),
    width = 7,
    height = 4.5
  ),
  save_png(
    p_changes,
    paste0(pattern, "-omega-class-proportion-change-uncertainty.png"),
    width = 10,
    height = 4
  )
)

plot_check <- tibble(path = saved_plot_paths) |>
  mutate(
    exists = file.exists(.data$path),
    extension = tools::file_ext(.data$path),
    is_png = tolower(.data$extension) == "png"
  )
stopifnot(all(plot_check$exists), all(plot_check$is_png))

#### Group-pair association patterns ####

load(file.path("Data", "preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData"))
if (!exists("Tr")) {
  stop("Expected object `Tr` in preprocessed trait data.", call. = FALSE)
}

species_groups <- Tr |>
  as.data.frame() |>
  rownames_to_column(var = "species") |>
  transmute(
    species = .data$species,
    foraging_guild = .data$foraging_guild_consensus,
    migratory_strategy = .data$Migration_a3_DOF
  )

group_count_check <- species_groups |>
  summarise(
    n_foraging_guilds = n_distinct(.data$foraging_guild),
    n_migratory_strategies = n_distinct(.data$migratory_strategy)
  )
stopifnot(
  group_count_check$n_foraging_guilds == 34,
  group_count_check$n_migratory_strategies == 5
)

group_specs <- tibble(
  grouping_type = c("foraging_guild", "migratory_strategy"),
  group_column = c("foraging_guild", "migratory_strategy"),
  display_name = c("Foraging guild", "Migratory strategy")
)

pair_species <- existing_pairs |>
  filter(.data$atlas == "1") |>
  distinct(.data$pair_id, .data$species_a, .data$species_b) |>
  arrange(.data$pair_id)

make_group_pair_map <- function(grouping_type, group_column) {
  labels <- species_groups |>
    select(species, group = all_of(group_column))

  out <- pair_species |>
    left_join(labels |> rename(species_a = species, group_a = group), by = "species_a") |>
    left_join(labels |> rename(species_b = species, group_b = group), by = "species_b") |>
    mutate(
      grouping_type = grouping_type,
      group_1 = pmin(.data$group_a, .data$group_b),
      group_2 = pmax(.data$group_a, .data$group_b),
      group_pair = if_else(
        .data$group_1 == .data$group_2,
        .data$group_1,
        paste(.data$group_1, .data$group_2, sep = " | ")
      ),
      group_pair_type = if_else(
        .data$group_1 == .data$group_2,
        "Within group",
        "Between group"
      )
    ) |>
    select(
      .data$grouping_type, .data$pair_id, .data$species_a, .data$species_b,
      .data$group_a, .data$group_b, .data$group_1, .data$group_2,
      .data$group_pair, .data$group_pair_type
    )

  if (any(is.na(out$group_a)) || any(is.na(out$group_b))) {
    stop("Missing group labels for grouping type: ", grouping_type, call. = FALSE)
  }

  out
}

group_pair_maps <- pmap_dfr(
  group_specs,
  function(grouping_type, group_column, display_name) {
    make_group_pair_map(grouping_type, group_column)
  }
)

group_sizes <- group_specs |>
  pmap_dfr(function(grouping_type, group_column, display_name) {
    species_groups |>
      count(group = .data[[group_column]], name = "n_species") |>
      mutate(
        grouping_type = grouping_type,
        display_name = display_name
      ) |>
      select(.data$grouping_type, .data$display_name, .data$group, .data$n_species)
  })

make_observed_group_proportions <- function(grouping_type) {
  group_map <- group_pair_maps |>
    filter(.data$grouping_type == .env$grouping_type)

  pair_groups <- existing_pairs |>
    select(.data$atlas, .data$period, .data$pair_id, .data$support_class) |>
    left_join(
      group_map |>
        select(
          .data$pair_id, .data$group_1, .data$group_2,
          .data$group_pair, .data$group_pair_type
        ),
      by = "pair_id"
    ) |>
    mutate(
      grouping_type = grouping_type,
      period = as.character(.data$period),
      support_class = as.character(.data$support_class)
    )

  group_keys <- pair_groups |>
    distinct(
      .data$grouping_type, .data$atlas, .data$period,
      .data$group_pair_type, .data$group_1, .data$group_2,
      .data$group_pair
    )

  counts <- pair_groups |>
    count(
      .data$grouping_type, .data$atlas, .data$period,
      .data$group_pair_type, .data$group_1, .data$group_2,
      .data$group_pair, .data$support_class,
      name = "n_class_pairs"
    )

  expand_grid(group_keys, support_class = support_levels) |>
    left_join(
      counts,
      by = c(
        "grouping_type", "atlas", "period", "group_pair_type",
        "group_1", "group_2", "group_pair", "support_class"
      )
    ) |>
    mutate(n_class_pairs = replace_na(.data$n_class_pairs, 0L)) |>
    group_by(
      .data$grouping_type, .data$atlas, .data$period,
      .data$group_pair_type, .data$group_1, .data$group_2,
      .data$group_pair
    ) |>
    mutate(
      n_group_pairs = sum(.data$n_class_pairs),
      observed_prop = .data$n_class_pairs / .data$n_group_pairs
    ) |>
    ungroup() |>
    left_join(
      observed_proportions |>
        select(
          .data$atlas, .data$period, .data$support_class,
          atlas_prop = .data$prop_pairs
        ),
      by = c("atlas", "period", "support_class")
    ) |>
    mutate(
      observed_minus_atlas_prop = .data$observed_prop - .data$atlas_prop,
      unstable_small_pair_count =
        .data$n_group_pairs < min_group_pair_count_for_plots
    )
}

observed_group_proportions <- map_dfr(
  group_specs$grouping_type,
  make_observed_group_proportions
)

group_pair_count_validation <- observed_group_proportions |>
  distinct(
    .data$grouping_type, .data$atlas, .data$period,
    .data$group_pair, .data$n_group_pairs
  ) |>
  group_by(.data$grouping_type, .data$atlas, .data$period) |>
  summarise(total_pairs = sum(.data$n_group_pairs), .groups = "drop")
stopifnot(all(group_pair_count_validation$total_pairs == 12246))

bootstrap_group_class_proportions <- function(draw_object, grouping_type) {
  group_map <- group_pair_maps |>
    filter(.data$grouping_type == .env$grouping_type) |>
    arrange(.data$pair_id)

  draw_order <- order(draw_object$pair_id)
  if (!identical(group_map$pair_id, draw_object$pair_id[draw_order])) {
    stop("Group map pair IDs do not match posterior draw pair IDs.", call. = FALSE)
  }

  positive_draws <- draw_object$positive_draws[, draw_order, drop = FALSE]
  n_draws <- nrow(positive_draws)
  bootstrap_weights <- t(rmultinom(
    n = n_bootstrap,
    size = bootstrap_draws_per_replicate,
    prob = rep(1 / n_draws, n_draws)
  ))
  pr_positive_boot <- (bootstrap_weights %*% (positive_draws + 0)) /
    bootstrap_draws_per_replicate

  class_index <- matrix(3L, nrow = n_bootstrap, ncol = ncol(pr_positive_boot))
  class_index[pr_positive_boot > support_level] <- 1L
  class_index[pr_positive_boot < negative_support_cutoff] <- 2L

  group_keys <- group_map |>
    distinct(
      .data$grouping_type, .data$group_pair_type,
      .data$group_1, .data$group_2, .data$group_pair
    ) |>
    arrange(.data$group_1, .data$group_2)

  map_dfr(seq_len(nrow(group_keys)), function(i) {
    group_info <- group_keys[i, ]
    pair_index <- which(group_map$group_pair == group_info$group_pair)
    n_group_pairs <- length(pair_index)

    class_props <- cbind(
      rowSums(class_index[, pair_index, drop = FALSE] == 1L),
      rowSums(class_index[, pair_index, drop = FALSE] == 2L),
      rowSums(class_index[, pair_index, drop = FALSE] == 3L)
    ) / n_group_pairs
    colnames(class_props) <- support_levels

    tibble(
      grouping_type = grouping_type,
      atlas = draw_object$atlas,
      period = draw_object$period,
      group_pair_type = group_info$group_pair_type,
      group_1 = group_info$group_1,
      group_2 = group_info$group_2,
      group_pair = group_info$group_pair,
      n_group_pairs = n_group_pairs,
      support_class = support_levels,
      bootstrap_mean = colMeans(class_props, na.rm = TRUE),
      bootstrap_median = apply(class_props, 2, median, na.rm = TRUE),
      bootstrap_lower_80 = apply(class_props, 2, quantile, probs = 0.10, na.rm = TRUE),
      bootstrap_upper_80 = apply(class_props, 2, quantile, probs = 0.90, na.rm = TRUE),
      bootstrap_lower_95 = apply(class_props, 2, quantile, probs = 0.025, na.rm = TRUE),
      bootstrap_upper_95 = apply(class_props, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
  })
}

group_pair_bootstrap_uncertainty <- pmap_dfr(
  expand_grid(
    draw_name = names(omega_draws_by_atlas),
    grouping_type = group_specs$grouping_type
  ),
  function(draw_name, grouping_type) {
    bootstrap_group_class_proportions(
      draw_object = omega_draws_by_atlas[[draw_name]],
      grouping_type = grouping_type
    )
  }
)

group_pair_class_proportions <- observed_group_proportions |>
  left_join(
    group_pair_bootstrap_uncertainty,
    by = c(
      "grouping_type", "atlas", "period", "group_pair_type",
      "group_1", "group_2", "group_pair", "n_group_pairs",
      "support_class"
    )
  ) |>
  arrange(
    .data$grouping_type, .data$atlas, .data$support_class,
    .data$group_1, .data$group_2
  )

make_group_permutation_null <- function(grouping_type, group_column) {
  group_map <- group_pair_maps |>
    filter(.data$grouping_type == .env$grouping_type) |>
    arrange(.data$pair_id)

  group_levels <- group_map |>
    distinct(
      .data$group_pair_type, .data$group_1,
      .data$group_2, .data$group_pair
    ) |>
    arrange(.data$group_1, .data$group_2)
  group_pair_levels <- group_levels$group_pair
  n_group_pairs <- as.integer(table(factor(group_map$group_pair, levels = group_pair_levels)))

  species_labels <- species_groups[[group_column]]
  names(species_labels) <- species_groups$species

  support_by_atlas <- map(
    model_metadata$atlas,
    function(atlas_id) {
      existing_pairs |>
        filter(.data$atlas == atlas_id) |>
        arrange(.data$pair_id) |>
        pull(.data$support_class) |>
        as.character()
    }
  )
  names(support_by_atlas) <- model_metadata$atlas

  null_keys <- expand_grid(
    tibble(atlas = model_metadata$atlas, period = model_metadata$period),
    group_levels,
    support_class = support_levels
  ) |>
    mutate(
      grouping_type = grouping_type,
      n_group_pairs = rep(
        rep(n_group_pairs, each = length(support_levels)),
        times = nrow(model_metadata)
      )
    ) |>
    select(
      .data$grouping_type, .data$atlas, .data$period,
      .data$group_pair_type, .data$group_1, .data$group_2,
      .data$group_pair, .data$n_group_pairs, .data$support_class
    )

  null_values <- matrix(
    NA_real_,
    nrow = n_group_permutations,
    ncol = nrow(null_keys)
  )

  for (iteration in seq_len(n_group_permutations)) {
    permuted_labels <- sample(species_labels)
    names(permuted_labels) <- names(species_labels)
    perm_group_a <- unname(permuted_labels[group_map$species_a])
    perm_group_b <- unname(permuted_labels[group_map$species_b])
    perm_group_pair <- if_else(
      pmin(perm_group_a, perm_group_b) == pmax(perm_group_a, perm_group_b),
      pmin(perm_group_a, perm_group_b),
      paste(pmin(perm_group_a, perm_group_b), pmax(perm_group_a, perm_group_b), sep = " | ")
    )

    cursor <- 1L
    for (atlas_id in model_metadata$atlas) {
      class_table <- table(
        factor(perm_group_pair, levels = group_pair_levels),
        factor(support_by_atlas[[atlas_id]], levels = support_levels)
      )
      class_props <- sweep(class_table, 1, rowSums(class_table), "/")
      block_values <- as.vector(t(class_props[, support_levels, drop = FALSE]))
      null_values[iteration, cursor:(cursor + length(block_values) - 1L)] <- block_values
      cursor <- cursor + length(block_values)
    }
  }

  observed_vector <- group_pair_class_proportions |>
    filter(.data$grouping_type == grouping_type) |>
    select(
      .data$grouping_type, .data$atlas, .data$period,
      .data$group_pair_type, .data$group_1, .data$group_2,
      .data$group_pair, .data$n_group_pairs, .data$support_class,
      .data$observed_prop
    )

  null_summary <- null_keys |>
    left_join(
      observed_vector,
      by = c(
        "grouping_type", "atlas", "period", "group_pair_type",
        "group_1", "group_2", "group_pair", "n_group_pairs",
        "support_class"
      )
    ) |>
    mutate(
      null_mean = colMeans(null_values, na.rm = TRUE),
      null_sd = apply(null_values, 2, sd, na.rm = TRUE),
      null_lower_95 = apply(null_values, 2, quantile, probs = 0.025, na.rm = TRUE),
      null_upper_95 = apply(null_values, 2, quantile, probs = 0.975, na.rm = TRUE),
      observed_minus_null = .data$observed_prop - .data$null_mean,
      null_z = if_else(
        .data$null_sd > 0,
        .data$observed_minus_null / .data$null_sd,
        NA_real_
      ),
      empirical_p_upper =
        (colSums(sweep(null_values, 2, .data$observed_prop, `>=`), na.rm = TRUE) + 1) /
        (n_group_permutations + 1),
      empirical_p_lower =
        (colSums(sweep(null_values, 2, .data$observed_prop, `<=`), na.rm = TRUE) + 1) /
        (n_group_permutations + 1),
      empirical_p_two_sided = pmin(
        1,
        2 * pmin(.data$empirical_p_upper, .data$empirical_p_lower)
      ),
      unstable_small_pair_count =
        .data$n_group_pairs < min_group_pair_count_for_plots
    ) |>
    group_by(.data$grouping_type, .data$support_class) |>
    mutate(
      fdr_p = p.adjust(.data$empirical_p_two_sided, method = "fdr"),
      fdr_0_05 = .data$fdr_p < 0.05
    ) |>
    ungroup()

  null_summary
}

group_pair_null_results <- pmap_dfr(
  group_specs,
  function(grouping_type, group_column, display_name) {
    make_group_permutation_null(grouping_type, group_column)
  }
) |>
  left_join(
    group_specs |> select(.data$grouping_type, .data$display_name),
    by = "grouping_type"
  ) |>
  arrange(
    .data$grouping_type, .data$support_class,
    .data$atlas, desc(abs(.data$observed_minus_null))
  )

group_pair_temporal_consistency <- group_pair_null_results |>
  filter(.data$support_class %in% c("Supported positive", "Supported negative")) |>
  group_by(
    .data$grouping_type, .data$display_name, .data$group_pair_type,
    .data$group_1, .data$group_2, .data$group_pair,
    .data$n_group_pairs, .data$support_class
  ) |>
  summarise(
    n_periods = n_distinct(.data$period),
    mean_observed_minus_null = mean(.data$observed_minus_null, na.rm = TRUE),
    min_observed_minus_null = min(.data$observed_minus_null, na.rm = TRUE),
    max_observed_minus_null = max(.data$observed_minus_null, na.rm = TRUE),
    min_abs_null_z = min(abs(.data$null_z), na.rm = TRUE),
    max_abs_observed_minus_null = max(abs(.data$observed_minus_null), na.rm = TRUE),
    n_fdr_0_05 = sum(.data$fdr_0_05, na.rm = TRUE),
    all_above_null = all(.data$observed_minus_null > 0, na.rm = TRUE),
    all_below_null = all(.data$observed_minus_null < 0, na.rm = TRUE),
    excess_1970s = .data$observed_minus_null[.data$period == "1970s"][1],
    excess_1990s = .data$observed_minus_null[.data$period == "1990s"][1],
    excess_2010s = .data$observed_minus_null[.data$period == "2010s"][1],
    .groups = "drop"
  ) |>
  mutate(
    delta_2010s_minus_1970s = .data$excess_2010s - .data$excess_1970s,
    sign_flip_1970s_to_2010s =
      sign(.data$excess_1970s) != sign(.data$excess_2010s),
    temporal_pattern = case_when(
      .data$all_above_null ~ "Consistently above null",
      .data$all_below_null ~ "Consistently below null",
      .data$sign_flip_1970s_to_2010s ~ "Direction flip",
      TRUE ~ "Mixed direction"
    ),
    unstable_small_pair_count =
      .data$n_group_pairs < min_group_pair_count_for_plots
  ) |>
  arrange(
    .data$grouping_type, .data$support_class,
    desc(.data$n_fdr_0_05), desc(.data$min_abs_null_z)
  )

top_consistent_group_pairs <- group_pair_temporal_consistency |>
  filter(
    !.data$unstable_small_pair_count,
    .data$temporal_pattern %in% c("Consistently above null", "Consistently below null")
  ) |>
  group_by(.data$grouping_type, .data$support_class) |>
  slice_max(
    order_by = .data$min_abs_null_z,
    n = 30,
    with_ties = FALSE
  ) |>
  ungroup()

top_changing_group_pairs <- group_pair_temporal_consistency |>
  filter(!.data$unstable_small_pair_count) |>
  group_by(.data$grouping_type, .data$support_class) |>
  slice_max(
    order_by = abs(.data$delta_2010s_minus_1970s),
    n = 30,
    with_ties = FALSE
  ) |>
  ungroup()

write_csv(
  group_pair_class_proportions,
  file.path(out_dir, paste0(pattern, "-omega-group-pair-class-proportions-uncertainty.csv"))
)
write_csv(
  group_pair_null_results,
  file.path(out_dir, paste0(pattern, "-omega-group-pair-permutation-null-results.csv"))
)
write_csv(
  group_pair_temporal_consistency,
  file.path(out_dir, paste0(pattern, "-omega-group-pair-temporal-consistency.csv"))
)
write_csv(
  top_consistent_group_pairs,
  file.path(out_dir, paste0(pattern, "-omega-group-pair-top-consistent.csv"))
)
write_csv(
  top_changing_group_pairs,
  file.path(out_dir, paste0(pattern, "-omega-group-pair-top-changing.csv"))
)
write_csv(
  group_pair_count_validation,
  file.path(out_dir, paste0(pattern, "-omega-group-pair-count-validation.csv"))
)

make_symmetric_group_tiles <- function(df, value_col, group_levels, value_name = "value") {
  value_col <- rlang::ensym(value_col)

  bind_rows(
    df |>
      transmute(
        period = .data$period,
        support_class = .data$support_class,
        group_row = .data$group_1,
        group_col = .data$group_2,
        value = !!value_col
      ),
    df |>
      filter(.data$group_1 != .data$group_2) |>
      transmute(
        period = .data$period,
        support_class = .data$support_class,
        group_row = .data$group_2,
        group_col = .data$group_1,
        value = !!value_col
      )
  ) |>
    mutate(
      group_row = factor(.data$group_row, levels = rev(group_levels)),
      group_col = factor(.data$group_col, levels = group_levels),
      period = factor(.data$period, levels = period_order),
      support_class = factor(.data$support_class, levels = support_levels)
    ) |>
    rename(!!value_name := value)
}

make_group_heatmap <- function(grouping_type, filename, width, height) {
  if (grouping_type == "foraging_guild") {
    plot_group_levels <- group_sizes |>
      filter(
        .data$grouping_type == .env$grouping_type,
        .data$n_species >= min_guild_species_for_plots
      ) |>
      arrange(desc(.data$n_species), .data$group) |>
      pull(.data$group)
  } else {
    plot_group_levels <- group_sizes |>
      filter(.data$grouping_type == .env$grouping_type) |>
      arrange(desc(.data$n_species), .data$group) |>
      pull(.data$group)
  }

  plot_data <- group_pair_null_results |>
    filter(
      .data$grouping_type == .env$grouping_type,
      .data$support_class %in% c("Supported positive", "Supported negative"),
      .data$n_group_pairs >= min_group_pair_count_for_plots,
      .data$group_1 %in% plot_group_levels,
      .data$group_2 %in% plot_group_levels
    ) |>
    make_symmetric_group_tiles(
      value_col = observed_minus_null,
      group_levels = plot_group_levels,
      value_name = "observed_minus_null"
    )

  fill_limit <- quantile(abs(plot_data$observed_minus_null), 0.98, na.rm = TRUE)

  p <- ggplot(
    plot_data,
    aes(x = .data$group_col, y = .data$group_row, fill = .data$observed_minus_null)
  ) +
    geom_tile(colour = "white", linewidth = 0.15) +
    facet_grid(.data$support_class ~ .data$period) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-fill_limit, fill_limit),
      oob = squish,
      labels = percent_format(),
      name = "Observed - null"
    ) +
    coord_equal(expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank(),
      strip.text = element_text(size = 9)
    )

  save_png(p, filename, width = width, height = height)
}

make_group_temporal_plot <- function(grouping_type, filename, width, height) {
  selected_pairs <- group_pair_temporal_consistency |>
    filter(
      .data$grouping_type == .env$grouping_type,
      !.data$unstable_small_pair_count,
      .data$support_class %in% c("Supported positive", "Supported negative")
    ) |>
    group_by(.data$support_class) |>
    slice_max(
      order_by = .data$max_abs_observed_minus_null,
      n = 8,
      with_ties = FALSE
    ) |>
    ungroup() |>
    mutate(pair_class_id = paste(.data$group_pair, .data$support_class, sep = " :: "))

  plot_data <- group_pair_null_results |>
    semi_join(
      selected_pairs |>
        select(.data$grouping_type, .data$group_pair, .data$support_class),
      by = c("grouping_type", "group_pair", "support_class")
    ) |>
    mutate(
      period = factor(.data$period, levels = period_order),
      pair_label = fct_reorder(
        paste(.data$group_pair, .data$support_class, sep = " :: "),
        abs(.data$observed_minus_null),
        .fun = max,
        .desc = TRUE
      )
    )

  p <- ggplot(
    plot_data,
    aes(x = .data$period, y = .data$observed_minus_null, group = .data$pair_label)
  ) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.35) +
    geom_line(aes(colour = .data$support_class), linewidth = 0.55) +
    geom_point(aes(colour = .data$support_class), size = 1.8) +
    facet_wrap(~ pair_label, ncol = 4) +
    scale_colour_manual(values = support_colours, name = NULL) +
    scale_y_continuous(labels = percent_format()) +
    labs(x = NULL, y = "Observed - permutation null") +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 7)
    )

  save_png(p, filename, width = width, height = height)
}

group_plot_paths <- c(
  make_group_heatmap(
    grouping_type = "foraging_guild",
    filename = paste0(pattern, "-omega-foraging-guild-group-pair-excess-heatmaps.png"),
    width = 14,
    height = 10
  ),
  make_group_heatmap(
    grouping_type = "migratory_strategy",
    filename = paste0(pattern, "-omega-migratory-strategy-group-pair-excess-heatmaps.png"),
    width = 9,
    height = 5
  ),
  make_group_temporal_plot(
    grouping_type = "foraging_guild",
    filename = paste0(pattern, "-omega-foraging-guild-group-pair-temporal-signals.png"),
    width = 12,
    height = 8
  ),
  make_group_temporal_plot(
    grouping_type = "migratory_strategy",
    filename = paste0(pattern, "-omega-migratory-strategy-group-pair-temporal-signals.png"),
    width = 10,
    height = 5
  )
)

group_plot_check <- tibble(path = group_plot_paths) |>
  mutate(
    exists = file.exists(.data$path),
    extension = tools::file_ext(.data$path),
    is_png = tolower(.data$extension) == "png"
  )
stopifnot(all(group_plot_check$exists), all(group_plot_check$is_png))

message("Wrote uncertainty outputs to: ", out_dir)
message("Wrote PNG plots to: ", plot_dir)
