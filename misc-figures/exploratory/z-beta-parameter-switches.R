# Beta-parameter switching through time
#
# This exploratory figure asks which species changed their environmental
# response profiles most strongly between adjacent atlas periods, and which
# species stayed comparatively stable. The response profile for each species is
# its vector of HMSC beta coefficients across the environmental predictors.
#
# Output policy for this project: this script writes PNG figures only.

remove(list = ls())
.libPaths(c("~/Rlibs", .libPaths()))

required_packages <- c("tidyverse", "readxl", "patchwork", "scales", "Hmsc")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(tidyverse)
library(readxl)
library(patchwork)
library(scales)
library(Hmsc)

source(file.path("support_scripts", "figure_data_helpers.R"))

# ---- Configuration ----
pattern <- "2026-03-13"
base_dir <- "HmscOutputs"
out_dir <- file.path("misc-figures", "outputs", "exploratory", "beta-parameter-switches")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_trait_group_n <- 5
n_species_each_tail <- 8
n_trait_groups_each_tail <- 5

period_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")
transition_pairs <- tibble(
  atlas_from = c("1", "2"),
  atlas_to = c("2", "3"),
  transition = c("1970s to 1990s", "1990s to 2010s")
)

var_labels <- c(
  "tmean_breeding" = "Temperature",
  "prec_breeding" = "Precipitation",
  "hh" = "Land-use heterogeneity",
  "perc_urban" = "Urban",
  "perc_cropland" = "Cropland",
  "perc_pasture" = "Pasture",
  "perc_forest" = "Forest",
  "perc_grass_shrub" = "Grass/Shrubland"
)

driver_colours <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban" = "#4d4d4d",
  "Cropland" = "goldenrod1",
  "Pasture" = "darkorange",
  "Forest" = "springgreen4",
  "Grass/Shrubland" = "springgreen2"
)

thermal_levels <- c("very cold", "cold", "medium", "warm", "very warm")

# ---- Helpers ----
predictor_names <- function(df) {
  c(
    "tmean_breeding",
    "prec_breeding",
    "hh",
    grep("^perc_", colnames(df), value = TRUE)
  )
}

thermal_groups <- function(sti) {
  cut(
    sti,
    breaks = quantile(sti, seq(0, 1, 0.2), na.rm = TRUE),
    include.lowest = TRUE,
    labels = thermal_levels
  )
}

support_state <- function(prob_positive) {
  case_when(
    prob_positive >= 0.95 ~ "Supported positive",
    prob_positive <= 0.05 ~ "Supported negative",
    TRUE ~ "Unsupported"
  )
}

read_beta_parameter_workbook <- function(model_folder, base_dir = "HmscOutputs") {
  beta_path <- file.path(
    base_dir,
    model_folder,
    "Results",
    paste0(model_folder, "parameter_estimates_Beta_.xlsx")
  )

  posterior_mean <- read_excel(beta_path, sheet = "Posterior mean") |>
    pivot_longer(
      cols = -Species,
      names_to = "variable",
      values_to = "beta_mean"
    )

  posterior_positive <- read_excel(beta_path, sheet = "Pr(x>0)") |>
    pivot_longer(
      cols = -Species,
      names_to = "variable",
      values_to = "prob_positive"
    )

  posterior_mean |>
    left_join(posterior_positive, by = c("Species", "variable")) |>
    mutate(
      atlas = as.character(sub(".*Atlas([0-9]+).*", "\\1", model_folder)),
      model_folder = model_folder,
      period = period_lookup[atlas],
      certainty_weight = 2 * abs(prob_positive - 0.5)
    )
}

predictor_sds_by_atlas <- function(models, variables) {
  imap_dfr(models, function(model, atlas) {
    tibble(
      atlas = as.character(atlas),
      variable = variables,
      predictor_sd = map_dbl(variables, ~ sd(model$XData[[.x]], na.rm = TRUE))
    )
  }) |>
    mutate(period = period_lookup[atlas])
}

make_transition_beta <- function(beta_long, atlas_from, atlas_to, transition) {
  beta_from <- beta_long |>
    filter(atlas == atlas_from) |>
    transmute(
      species = Species,
      variable,
      variable_label,
      period_from = period,
      weighted_beta_from = weighted_standardized_beta,
      standardized_beta_from = standardized_beta,
      prob_positive_from = prob_positive,
      support_state_from = beta_support_state
    )

  beta_to <- beta_long |>
    filter(atlas == atlas_to) |>
    transmute(
      species = Species,
      variable,
      variable_label,
      period_to = period,
      weighted_beta_to = weighted_standardized_beta,
      standardized_beta_to = standardized_beta,
      prob_positive_to = prob_positive,
      support_state_to = beta_support_state
    )

  beta_from |>
    inner_join(beta_to, by = c("species", "variable", "variable_label")) |>
    mutate(
      atlas_from = atlas_from,
      atlas_to = atlas_to,
      transition = transition,
      delta_weighted_beta = weighted_beta_to - weighted_beta_from,
      abs_delta_weighted_beta = abs(delta_weighted_beta),
      delta_standardized_beta = standardized_beta_to - standardized_beta_from,
      abs_delta_standardized_beta = abs(delta_standardized_beta),
      support_state_switch = support_state_from != support_state_to,
      supported_sign_flip = (
        support_state_from == "Supported positive" & support_state_to == "Supported negative"
      ) | (
        support_state_from == "Supported negative" & support_state_to == "Supported positive"
      )
    )
}

summarise_trait_group <- function(species_jump_df, trait_col, trait_family) {
  species_jump_df |>
    filter(!is.na(.data[[trait_col]]), .data[[trait_col]] != "") |>
    group_by(transition, trait_value = .data[[trait_col]]) |>
    summarise(
      n_species = n_distinct(species),
      median_total_jump = median(total_beta_jump, na.rm = TRUE),
      mean_total_jump = mean(total_beta_jump, na.rm = TRUE),
      iqr_total_jump = IQR(total_beta_jump, na.rm = TRUE),
      stable_species_share = mean(total_beta_jump <= quantile(total_beta_jump, 0.25, na.rm = TRUE), na.rm = TRUE),
      extreme_species_share = mean(total_beta_jump >= quantile(total_beta_jump, 0.75, na.rm = TRUE), na.rm = TRUE),
      share_with_any_state_switch = mean(n_support_state_switches > 0, na.rm = TRUE),
      share_with_supported_sign_flip = mean(n_supported_sign_flips > 0, na.rm = TRUE),
      median_predictor_jumps = median(n_predictors_with_large_jump, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      trait_family = trait_family,
      reliable_group = n_species >= min_trait_group_n
    )
}

select_extreme_rows <- function(df, value_col, n_each_tail) {
  value_col <- enquo(value_col)

  arranged_df <- df |>
    arrange(!!value_col) |>
    mutate(
      stable_rank = row_number(),
      jump_rank = row_number(desc(!!value_col))
    )

  # When there are only a few groups, using the full tail size would put the
  # same row in both the "stable" and "jump" tails. Shrinking the tail size keeps
  # the two labels mutually exclusive and leaves the middle rows out.
  tail_n <- max(1, min(n_each_tail, floor(nrow(arranged_df) / 2)))

  arranged_df |>
    mutate(
      focus_class = case_when(
        stable_rank <= tail_n ~ "Most stable",
        jump_rank <= tail_n ~ "Largest jumps",
        TRUE ~ "Middle"
      )
    ) |>
    filter(focus_class != "Middle") |>
    arrange(focus_class, !!value_col)
}

# ---- Load model outputs and beta coefficients ----
message("Loading model folders and fitted-model metadata...")
model_folders <- figure_model_folders(pattern = pattern, base_dir = base_dir)
model_atlases <- atlas_numbers(model_folders)

mods <- load_hmsc_posteriors(model_folders, base_dir = base_dir)
names(mods) <- as.character(model_atlases)

beta_predictors <- intersect(
  predictor_names(mods[[1]]$XData),
  names(var_labels)
)

message("Reading beta workbooks...")
beta_raw <- map_dfr(model_folders, read_beta_parameter_workbook, base_dir = base_dir) |>
  filter(variable != "(Intercept)", variable %in% beta_predictors)

beta_sds <- predictor_sds_by_atlas(mods, beta_predictors)

beta_long <- beta_raw |>
  left_join(beta_sds, by = c("atlas", "period", "variable")) |>
  mutate(
    standardized_beta = beta_mean * predictor_sd,
    weighted_standardized_beta = standardized_beta * certainty_weight,
    beta_support_state = support_state(prob_positive),
    variable_label = recode(variable, !!!var_labels),
    variable_label = factor(variable_label, levels = unname(var_labels))
  )

# ---- Load species trait metadata ----
message("Loading species trait metadata...")
load(file.path("Data", "preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData"))

species_traits <- Tr |>
  as.data.frame() |>
  rownames_to_column("species") |>
  transmute(
    species,
    migratory_strategy = Migration_a3_DOF,
    foraging_guild = foraging_guild_consensus,
    species_thermal_index,
    thermal_group = thermal_groups(species_thermal_index),
    thermal_group = factor(thermal_group, levels = thermal_levels)
  )

# ---- Measure beta switches between adjacent atlas periods ----
message("Calculating adjacent-period beta changes...")
transition_beta <- pmap_dfr(
  transition_pairs,
  ~ make_transition_beta(
    beta_long = beta_long,
    atlas_from = ..1,
    atlas_to = ..2,
    transition = ..3
  )
) |>
  mutate(
    transition = factor(transition, levels = transition_pairs$transition)
  )

# A predictor is counted as a large contributor to a species' jump if it is in
# the upper quartile of all absolute predictor-level changes for that transition.
predictor_jump_cutoffs <- transition_beta |>
  group_by(transition) |>
  summarise(
    large_jump_cutoff = quantile(abs_delta_weighted_beta, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

transition_beta <- transition_beta |>
  left_join(predictor_jump_cutoffs, by = "transition") |>
  mutate(predictor_large_jump = abs_delta_weighted_beta >= large_jump_cutoff)

species_jump <- transition_beta |>
  group_by(species, transition, atlas_from, atlas_to, period_from, period_to) |>
  summarise(
    total_beta_jump = sqrt(sum(delta_weighted_beta^2, na.rm = TRUE)),
    mean_abs_predictor_jump = mean(abs_delta_weighted_beta, na.rm = TRUE),
    max_abs_predictor_jump = max(abs_delta_weighted_beta, na.rm = TRUE),
    dominant_driver = as.character(variable_label[which.max(abs_delta_weighted_beta)][1]),
    n_support_state_switches = sum(support_state_switch, na.rm = TRUE),
    n_supported_sign_flips = sum(supported_sign_flip, na.rm = TRUE),
    n_predictors_with_large_jump = sum(predictor_large_jump, na.rm = TRUE),
    n_predictors = n(),
    .groups = "drop"
  ) |>
  left_join(species_traits, by = "species") |>
  group_by(transition) |>
  mutate(
    jump_percentile = percent_rank(total_beta_jump),
    species_stability_rank = row_number(total_beta_jump),
    species_jump_rank = row_number(desc(total_beta_jump))
  ) |>
  ungroup()

trait_group_summary <- bind_rows(
  summarise_trait_group(species_jump, "foraging_guild", "Foraging guild"),
  summarise_trait_group(species_jump, "migratory_strategy", "Migratory strategy"),
  summarise_trait_group(species_jump, "thermal_group", "Thermal affinity")
) |>
  group_by(trait_family, transition) |>
  mutate(
    group_stability_rank = row_number(median_total_jump),
    group_jump_rank = row_number(desc(median_total_jump))
  ) |>
  ungroup()

# ---- Save reusable data tables ----
write_csv(
  beta_long,
  file.path(out_dir, paste0(pattern, "-beta-standardized-effects.csv"))
)
write_csv(
  transition_beta,
  file.path(out_dir, paste0(pattern, "-beta-transition-predictor-jumps.csv"))
)
write_csv(
  species_jump,
  file.path(out_dir, paste0(pattern, "-beta-transition-species-jumps.csv"))
)
write_csv(
  trait_group_summary,
  file.path(out_dir, paste0(pattern, "-beta-transition-trait-group-summary.csv"))
)

# ---- Plot species-level largest jumps and strongest stability ----
species_focus <- species_jump |>
  group_by(transition) |>
  group_modify(~ select_extreme_rows(.x, total_beta_jump, n_species_each_tail)) |>
  ungroup() |>
  mutate(
    focus_class = factor(focus_class, levels = c("Largest jumps", "Most stable")),
    species_label = str_replace_all(species, "_", " "),
    facet_label = paste(species_label, transition, focus_class, sep = "___")
  ) |>
  arrange(transition, focus_class, total_beta_jump) |>
  mutate(facet_label = factor(facet_label, levels = unique(facet_label)))

species_switch_plot <- ggplot(
  species_focus,
  aes(x = total_beta_jump, y = facet_label, fill = dominant_driver)
) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  geom_point(
    aes(shape = n_supported_sign_flips > 0),
    size = 2.2,
    colour = "grey15",
    fill = "white",
    stroke = 0.7
  ) +
  facet_wrap(focus_class ~ transition, ncol = 2, scales = "free_y") +
  scale_y_discrete(labels = ~ str_remove(.x, "___.*$")) +
  scale_fill_manual(values = driver_colours, name = "Largest changing\npredictor") +
  scale_shape_manual(
    values = c(`FALSE` = 21, `TRUE` = 24),
    labels = c(`FALSE` = "No", `TRUE` = "Yes"),
    name = "Supported\nsign flip"
  ) +
  labs(
    title = "Species with the largest and smallest beta-parameter jumps",
    subtitle = "Jump size is Euclidean distance between adjacent certainty-weighted standardized beta vectors.",
    x = "Total beta-vector jump",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "grey35", size = 9)
  )

# ---- Plot trait-level summaries of jumpy versus stable groups ----
trait_focus <- trait_group_summary |>
  filter(reliable_group) |>
  group_by(trait_family, transition) |>
  group_modify(~ select_extreme_rows(.x, median_total_jump, n_trait_groups_each_tail)) |>
  ungroup() |>
  mutate(
    focus_class = factor(focus_class, levels = c("Largest jumps", "Most stable")),
    trait_value = as.character(trait_value),
    trait_label = paste0(trait_value, " (n=", n_species, ")"),
    facet_label = paste(trait_label, trait_family, transition, focus_class, sep = "___")
  ) |>
  arrange(trait_family, transition, focus_class, median_total_jump) |>
  mutate(facet_label = factor(facet_label, levels = unique(facet_label)))

trait_switch_plot <- ggplot(
  trait_focus,
  aes(x = median_total_jump, y = facet_label, fill = focus_class)
) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  geom_point(
    aes(size = share_with_any_state_switch),
    shape = 21,
    colour = "grey15",
    fill = "white",
    stroke = 0.45
  ) +
  facet_wrap(trait_family ~ transition, ncol = 2, scales = "free_y") +
  scale_y_discrete(labels = ~ str_remove(.x, "___.*$")) +
  scale_fill_manual(
    values = c("Largest jumps" = "#b2182b", "Most stable" = "#2166ac"),
    name = NULL
  ) +
  scale_size_continuous(
    range = c(1.8, 5.5),
    labels = percent_format(accuracy = 1),
    name = "Species with any\nsupport-state switch"
  ) +
  labs(
    title = "Trait groups separating unstable and stable beta responses",
    subtitle = paste0("Groups with fewer than ", min_trait_group_n, " species are excluded from this panel."),
    x = "Median species beta-vector jump",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "grey35", size = 9)
  )

final_plot <- species_switch_plot / trait_switch_plot +
  plot_layout(heights = c(1.05, 1.2), guides = "collect") +
  plot_annotation(
    title = "Beta-parameter switching across atlas periods",
    subtitle = "Environmental beta changes are standardized by within-atlas predictor SD and weighted by posterior certainty.",
    tag_levels = "A"
  ) &
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(colour = "grey30", size = 10),
    plot.tag = element_text(face = "bold", size = 13)
  )

ggsave(
  filename = file.path(out_dir, paste0(pattern, "-fig2b-beta-parameter-switches.png")),
  plot = final_plot,
  width = 13,
  height = 13.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("Finished beta-parameter switching figure.")
