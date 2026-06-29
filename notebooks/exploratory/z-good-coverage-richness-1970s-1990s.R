# Shared-site richness and species-turnover summaries for the 1970s and 1990s.
#
# Purpose:
#   Compare observed species richness between atlas 1 (1970s) and atlas 2
#   (1990s) after using a site and species pool that is comparable across atlas
#   periods.
#
# Filtering order:
#   1. Keep grid cells with pct_lnd >= 25, matching S0_Data_Definitions.R.
#   2. Keep only "Good" effort rows for atlas 1 and atlas 2.
#   3. Ignore effort for atlas 3 because that coverage table is not available.
#   4. Keep only sites present after those filters in atlas 1, 2, and 3.
#   5. Keep only species with at least 5 occurrences in each of atlas 1, 2, and
#      3 after shared-site filtering.
#   6. Use atlas 1 and atlas 2 only for richness maps, delta maps, Moran's I,
#      and species turnover summaries.
#
# Plot output is PNG only.

rm(list = ls())

library(dplyr)
library(ggplot2)
library(sf)
library(spdep)
library(stringr)

source(file.path("support_scripts", "data_helpers.R"))

#### PATHS ####

grid_threshold_path <- file.path(
  "Data", "data", "1_preprocessing", "atlas-grids",
  "grids-ocean-thresholds", "grids_ocean_thresholds.shp"
)

occurrence_path <- file.path(
  "Data", "data", "1_preprocessing", "Y_occurrences", "Y_occurrences.csv"
)

traits_path <- file.path(
  "Data", "data", "1_preprocessing", "Tr_aits", "traits-guild_migration.csv"
)

out_dir <- file.path(
  "notebooks", "exploratory", "outputs",
  "good-coverage-richness-1970s-1990s"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

filter_summary_path <- file.path(out_dir, "shared-site-filter-summary.csv")
richness_csv_path <- file.path(out_dir, "shared-site-richness-by-grid.csv")
delta_csv_path <- file.path(out_dir, "shared-site-richness-delta.csv")
moran_csv_path <- file.path(out_dir, "moran-richness-delta-summary.csv")
species_change_csv_path <- file.path(out_dir, "species-occupancy-change-by-guild.csv")
species_top_change_csv_path <- file.path(out_dir, "species-top-gains-losses.csv")
species_class_csv_path <- file.path(out_dir, "species-change-class-proportions.csv")
guild_class_csv_path <- file.path(out_dir, "guild-normalized-change-class-proportions.csv")

richness_plot_path <- file.path(out_dir, "shared-site-richness-1970s-1990s.png")
delta_plot_path <- file.path(out_dir, "shared-site-delta-richness-1990s-minus-1970s.png")
species_class_plot_path <- file.path(out_dir, "species-change-class-proportions.png")
guild_class_plot_path <- file.path(out_dir, "guild-normalized-change-class-proportions.png")

#### CONSTANTS ####

land_threshold_pct <- 25
min_species_occurrences <- 5
species_change_threshold <- 0.10

all_atlas_lookup <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

plot_atlas_lookup <- all_atlas_lookup[c("1", "2")]
plot_period_levels <- unname(plot_atlas_lookup)
species_change_levels <- c("Decrease", "Stable", "Increase")

#### SMALL HELPERS ####

survey_ids <- function(sites, atlas) {
  paste(sites, atlas, sep = "_")
}

percent_label <- function(x, accuracy = 1) {
  paste0(round(100 * x, accuracy), "%")
}

run_global_moran <- function(values, response_name, nb, listw) {
  # zero.policy = TRUE keeps isolated shared sites in the analysis object while
  # assigning them zero spatial lag contribution.
  moran_result <- spdep::moran.test(
    x = values,
    listw = listw,
    zero.policy = TRUE,
    na.action = na.fail
  )

  data.frame(
    response = response_name,
    moran_i = unname(moran_result$estimate[["Moran I statistic"]]),
    expected_i = unname(moran_result$estimate[["Expectation"]]),
    variance_i = unname(moran_result$estimate[["Variance"]]),
    p_value = moran_result$p.value,
    n_sites = length(values),
    n_no_neighbour_sites = sum(spdep::card(nb) == 0),
    stringsAsFactors = FALSE
  )
}

#### LOAD INPUTS ####

grid_thresholds <- st_read(grid_threshold_path, quiet = TRUE)

raw_y <- read.csv(
  occurrence_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

traits_raw <- read.csv(
  traits_path,
  row.names = 2,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

traits <- data.frame(
  species = rownames(traits_raw),
  foraging_guild_consensus = traits_raw$foraging_guild_consensus,
  stringsAsFactors = FALSE
)

effort <- read_effort_coverage(file.path("Data"))

#### FILTER SITES BEFORE FILTERING SPECIES ####

land_sites <- grid_thresholds |>
  st_drop_geometry() |>
  filter(.data$pct_lnd >= land_threshold_pct) |>
  pull(.data$kvdrtkd)

effort_12 <- effort |>
  mutate(
    site = str_remove(.data$survey, "_[12]$"),
    atlas = str_extract(.data$survey, "[12]$")
  ) |>
  filter(
    .data$atlas %in% names(plot_atlas_lookup),
    .data$coverage == "Good",
    .data$site %in% land_sites
  )

missing_effort_surveys <- setdiff(effort_12$survey, rownames(raw_y))

if (length(missing_effort_surveys) > 0) {
  stop(
    "Some Good-coverage atlas 1/2 surveys are absent from the occurrence matrix: ",
    paste(head(missing_effort_surveys, 10), collapse = ", "),
    call. = FALSE
  )
}

atlas3_land_surveys <- data.frame(
  site = land_sites,
  atlas = "3",
  survey = survey_ids(land_sites, "3"),
  coverage = "Not assessed",
  stringsAsFactors = FALSE
) |>
  filter(.data$survey %in% rownames(raw_y))

sites_by_atlas <- list(
  "1" = effort_12$site[effort_12$atlas == "1"],
  "2" = effort_12$site[effort_12$atlas == "2"],
  "3" = atlas3_land_surveys$site
)

# Keep the shapefile/site order so all later atlas-specific matrices are lined
# up before calculating richness and deltas.
shared_sites <- land_sites[land_sites %in% Reduce(intersect, sites_by_atlas)]

shared_survey_check <- unlist(
  lapply(names(all_atlas_lookup), function(atlas) survey_ids(shared_sites, atlas)),
  use.names = FALSE
)

missing_shared_surveys <- setdiff(shared_survey_check, rownames(raw_y))

if (length(missing_shared_surveys) > 0) {
  stop(
    "Some shared-site surveys are absent from the occurrence matrix: ",
    paste(head(missing_shared_surveys, 10), collapse = ", "),
    call. = FALSE
  )
}

#### FILTER SPECIES AFTER SHARED-SITE FILTERING ####

species_counts_by_atlas <- sapply(
  names(all_atlas_lookup),
  function(atlas) {
    colSums(
      raw_y[survey_ids(shared_sites, atlas), , drop = FALSE],
      na.rm = TRUE
    )
  }
)

retained_species <- rownames(species_counts_by_atlas)[
  apply(species_counts_by_atlas >= min_species_occurrences, 1, all)
]

if (length(retained_species) == 0) {
  stop("No species passed the shared-site >=5 occurrence filter.", call. = FALSE)
}

if (!all(retained_species %in% traits$species)) {
  missing_traits <- setdiff(retained_species, traits$species)
  warning(
    "Some retained species are missing guild metadata and will be labelled Unknown: ",
    paste(missing_traits, collapse = ", ")
  )
}

#### BUILD RICHNESS TABLES FOR ATLAS 1 AND 2 ####

richness_df <- bind_rows(
  lapply(names(plot_atlas_lookup), function(atlas) {
    y_atlas <- raw_y[
      survey_ids(shared_sites, atlas),
      retained_species,
      drop = FALSE
    ]

    data.frame(
      site = shared_sites,
      survey = survey_ids(shared_sites, atlas),
      atlas = atlas,
      period = unname(plot_atlas_lookup[[atlas]]),
      species_richness = rowSums(y_atlas, na.rm = TRUE),
      retained_species_pool_size = length(retained_species),
      stringsAsFactors = FALSE
    )
  })
) |>
  mutate(period = factor(.data$period, levels = plot_period_levels)) |>
  left_join(
    grid_thresholds |>
      st_drop_geometry() |>
      select("kvdrtkd", "pct_lnd"),
    by = c("site" = "kvdrtkd")
  ) |>
  arrange(.data$atlas, .data$site)

richness_1970s <- richness_df |>
  filter(.data$atlas == "1") |>
  arrange(.data$site)

richness_1990s <- richness_df |>
  filter(.data$atlas == "2") |>
  arrange(.data$site)

if (!identical(richness_1970s$site, richness_1990s$site)) {
  stop("1970s and 1990s richness rows do not have identical site order.", call. = FALSE)
}

delta_df <- richness_1970s |>
  transmute(
    site = .data$site,
    richness_1970s = .data$species_richness
  ) |>
  left_join(
    richness_1990s |>
      transmute(
        site = .data$site,
        richness_1990s = .data$species_richness
      ),
    by = "site"
  ) |>
  mutate(delta_richness_1990s_minus_1970s = .data$richness_1990s - .data$richness_1970s) |>
  left_join(
    grid_thresholds |>
      st_drop_geometry() |>
      select("kvdrtkd", "pct_lnd"),
    by = c("site" = "kvdrtkd")
  )

#### MORAN'S I FOR RICHNESS AND DELTA RICHNESS ####

shared_grid_sf <- grid_thresholds |>
  select("kvdrtkd") |>
  filter(.data$kvdrtkd %in% shared_sites) |>
  mutate(kvdrtkd = factor(.data$kvdrtkd, levels = shared_sites)) |>
  arrange(.data$kvdrtkd)

if (!identical(as.character(shared_grid_sf$kvdrtkd), shared_sites)) {
  stop("Shared grid polygons are not ordered like the shared site vector.", call. = FALSE)
}

shared_nb <- spdep::poly2nb(shared_grid_sf, queen = TRUE)
shared_listw <- spdep::nb2listw(shared_nb, style = "W", zero.policy = TRUE)

moran_summary <- bind_rows(
  run_global_moran(richness_1970s$species_richness, "1970s richness", shared_nb, shared_listw),
  run_global_moran(richness_1990s$species_richness, "1990s richness", shared_nb, shared_listw),
  run_global_moran(
    delta_df$delta_richness_1990s_minus_1970s,
    "1990s minus 1970s richness",
    shared_nb,
    shared_listw
  )
)

#### SPECIES OCCUPANCY CHANGE AND GUILD SUMMARIES ####

occupancy_counts <- sapply(
  names(plot_atlas_lookup),
  function(atlas) {
    y_atlas <- raw_y[
      survey_ids(shared_sites, atlas),
      retained_species,
      drop = FALSE
    ]
    colSums(y_atlas > 0, na.rm = TRUE)
  }
)

species_record_counts <- data.frame(
  species = retained_species,
  records_1970s = as.integer(species_counts_by_atlas[retained_species, "1"]),
  records_1990s = as.integer(species_counts_by_atlas[retained_species, "2"]),
  records_2010s = as.integer(species_counts_by_atlas[retained_species, "3"]),
  stringsAsFactors = FALSE
)

species_occupancy <- data.frame(
  species = retained_species,
  occupied_grid_cells_1970s = as.integer(occupancy_counts[retained_species, "1"]),
  occupied_grid_cells_1990s = as.integer(occupancy_counts[retained_species, "2"]),
  n_shared_grid_cells = length(shared_sites),
  stringsAsFactors = FALSE
) |>
  mutate(
    prop_occupied_1970s = .data$occupied_grid_cells_1970s / .data$n_shared_grid_cells,
    prop_occupied_1990s = .data$occupied_grid_cells_1990s / .data$n_shared_grid_cells,
    delta_occupied_grid_cells = .data$occupied_grid_cells_1990s - .data$occupied_grid_cells_1970s,
    delta_occupancy_proportion = .data$prop_occupied_1990s - .data$prop_occupied_1970s,
    delta_occupancy_percentage_points = 100 * .data$delta_occupancy_proportion,
    change_class = case_when(
      .data$delta_occupancy_proportion >= species_change_threshold ~ "Increase",
      .data$delta_occupancy_proportion <= -species_change_threshold ~ "Decrease",
      TRUE ~ "Stable"
    ),
    change_class = factor(.data$change_class, levels = species_change_levels)
  ) |>
  left_join(species_record_counts, by = "species") |>
  left_join(traits, by = "species") |>
  mutate(
    foraging_guild_consensus = ifelse(
      is.na(.data$foraging_guild_consensus) | .data$foraging_guild_consensus == "",
      "Unknown",
      .data$foraging_guild_consensus
    )
  ) |>
  arrange(.data$change_class, desc(abs(.data$delta_occupancy_percentage_points)), .data$species)

species_top_change <- species_occupancy |>
  filter(.data$change_class != "Stable") |>
  arrange(desc(.data$delta_occupancy_percentage_points))

species_class_proportions <- data.frame(
  change_class = species_change_levels,
  stringsAsFactors = FALSE
) |>
  left_join(
    species_occupancy |>
      count(.data$change_class, name = "n_species") |>
      mutate(change_class = as.character(.data$change_class)),
    by = "change_class"
  ) |>
  mutate(
    n_species = ifelse(is.na(.data$n_species), 0L, .data$n_species),
    total_species = sum(.data$n_species),
    proportion_species = .data$n_species / .data$total_species,
    percentage_species = 100 * .data$proportion_species,
    change_class = factor(.data$change_class, levels = species_change_levels)
  )

guild_class_grid <- expand.grid(
  foraging_guild_consensus = sort(unique(species_occupancy$foraging_guild_consensus)),
  change_class = species_change_levels,
  stringsAsFactors = FALSE
)

guild_class_proportions <- guild_class_grid |>
  left_join(
    species_occupancy |>
      count(.data$foraging_guild_consensus, .data$change_class, name = "n_species") |>
      mutate(change_class = as.character(.data$change_class)),
    by = c("foraging_guild_consensus", "change_class")
  ) |>
  mutate(n_species = ifelse(is.na(.data$n_species), 0L, .data$n_species)) |>
  group_by(.data$foraging_guild_consensus) |>
  mutate(
    guild_species_count = sum(.data$n_species),
    proportion_within_guild = .data$n_species / .data$guild_species_count,
    percentage_within_guild = 100 * .data$proportion_within_guild
  ) |>
  ungroup() |>
  mutate(
    change_class = factor(.data$change_class, levels = species_change_levels)
  )

guild_order <- guild_class_proportions |>
  group_by(.data$foraging_guild_consensus) |>
  summarise(
    guild_species_count = max(.data$guild_species_count),
    net_change_signal = sum(
      ifelse(.data$change_class == "Increase", .data$proportion_within_guild, 0) -
        ifelse(.data$change_class == "Decrease", .data$proportion_within_guild, 0)
    ),
    .groups = "drop"
  ) |>
  arrange(.data$net_change_signal, .data$foraging_guild_consensus) |>
  pull(.data$foraging_guild_consensus)

guild_class_proportions <- guild_class_proportions |>
  mutate(
    foraging_guild_consensus = factor(.data$foraging_guild_consensus, levels = guild_order)
  )

#### INTERNAL CHECKS ####

if (nrow(richness_1970s) != length(shared_sites) || nrow(richness_1990s) != length(shared_sites)) {
  stop("Richness tables do not contain exactly one row per shared site and plot atlas.", call. = FALSE)
}

if (!all(species_counts_by_atlas[retained_species, ] >= min_species_occurrences)) {
  stop("At least one retained species failed the >=5 occurrence filter.", call. = FALSE)
}

guild_sum_check <- guild_class_proportions |>
  group_by(.data$foraging_guild_consensus) |>
  summarise(total = sum(.data$proportion_within_guild), .groups = "drop")

if (any(abs(guild_sum_check$total - 1) > 1e-8)) {
  stop("At least one guild-normalized class proportion does not sum to 1.", call. = FALSE)
}

#### WRITE TABLE OUTPUTS ####

filter_summary <- data.frame(
  metric = c(
    "land_threshold_pct",
    "land_threshold_sites",
    "good_effort_atlas1_sites",
    "good_effort_atlas2_sites",
    "atlas3_land_threshold_sites_effort_ignored",
    "shared_sites_after_site_filters",
    "raw_species_pool_size",
    "retained_species_pool_size",
    "min_occurrences_required_per_atlas",
    "species_change_threshold_percentage_points"
  ),
  value = as.character(c(
    land_threshold_pct,
    length(land_sites),
    length(sites_by_atlas[["1"]]),
    length(sites_by_atlas[["2"]]),
    length(sites_by_atlas[["3"]]),
    length(shared_sites),
    ncol(raw_y),
    length(retained_species),
    min_species_occurrences,
    100 * species_change_threshold
  )),
  stringsAsFactors = FALSE
)

write.csv(filter_summary, filter_summary_path, row.names = FALSE, na = "")
write.csv(richness_df, richness_csv_path, row.names = FALSE, na = "")
write.csv(delta_df, delta_csv_path, row.names = FALSE, na = "")
write.csv(moran_summary, moran_csv_path, row.names = FALSE, na = "")
write.csv(species_occupancy, species_change_csv_path, row.names = FALSE, na = "")
write.csv(species_top_change, species_top_change_csv_path, row.names = FALSE, na = "")
write.csv(species_class_proportions, species_class_csv_path, row.names = FALSE, na = "")
write.csv(guild_class_proportions, guild_class_csv_path, row.names = FALSE, na = "")

#### PNG FIGURES ####

richness_plot_sf <- shared_grid_sf |>
  inner_join(richness_df, by = c("kvdrtkd" = "site"))

richness_limit <- max(richness_plot_sf$species_richness, na.rm = TRUE)
richness_breaks <- pretty(c(0, richness_limit), n = 6)
richness_palette <- rev(grDevices::rainbow(14, start = 0, end = 0.78))

richness_map <- ggplot(richness_plot_sf) +
  geom_sf(
    aes(fill = .data$species_richness),
    colour = "white",
    linewidth = 0.04
  ) +
  facet_wrap(vars(.data$period), nrow = 1) +
  coord_sf(
    xlim = c(400000, 910000),
    ylim = c(6000000, 6460000),
    expand = FALSE
  ) +
  scale_fill_gradientn(
    colours = richness_palette,
    limits = c(0, richness_limit),
    breaks = richness_breaks,
    name = "Observed\nrichness"
  ) +
  labs(
    title = "Shared-site species richness in Good-coverage grid cells",
    subtitle = paste0(
      length(shared_sites), " shared sites; ",
      length(retained_species), " species retained with >= ",
      min_species_occurrences, " records in each atlas period."
    )
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.7, "cm"),
    plot.title = element_text(face = "bold", colour = "grey15"),
    plot.subtitle = element_text(colour = "grey30", margin = margin(b = 8)),
    strip.text = element_text(face = "bold", colour = "grey15", size = 12),
    plot.margin = margin(10, 14, 10, 14)
  )

delta_plot_sf <- shared_grid_sf |>
  inner_join(delta_df, by = c("kvdrtkd" = "site"))

delta_limit <- max(abs(delta_plot_sf$delta_richness_1990s_minus_1970s), na.rm = TRUE)

delta_map <- ggplot(delta_plot_sf) +
  geom_sf(
    aes(fill = .data$delta_richness_1990s_minus_1970s),
    colour = "white",
    linewidth = 0.04
  ) +
  coord_sf(
    xlim = c(400000, 910000),
    ylim = c(6000000, 6460000),
    expand = FALSE
  ) +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-delta_limit, delta_limit),
    name = "Delta\nrichness"
  ) +
  labs(
    title = "Delta richness: 1990s minus 1970s",
    subtitle = "Positive values indicate grid cells with higher retained-species richness in the 1990s."
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.7, "cm"),
    plot.title = element_text(face = "bold", colour = "grey15"),
    plot.subtitle = element_text(colour = "grey30", margin = margin(b = 8)),
    plot.margin = margin(10, 14, 10, 14)
  )

class_colours <- c(
  "Decrease" = "#B2182B",
  "Stable" = "grey72",
  "Increase" = "#1A9850"
)

species_class_plot <- ggplot(
  species_class_proportions,
  aes(x = .data$change_class, y = .data$proportion_species, fill = .data$change_class)
) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.25) +
  geom_text(
    aes(label = paste0(.data$n_species, "\n", percent_label(.data$proportion_species))),
    vjust = -0.25,
    size = 3.5,
    lineheight = 0.9
  ) +
  scale_fill_manual(values = class_colours, guide = "none") +
  scale_y_continuous(
    labels = percent_label,
    limits = c(0, max(species_class_proportions$proportion_species) * 1.18),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Proportion of retained species",
    title = "Species occupancy-change classes",
    subtitle = "Classes use a >=10 percentage-point change in occupied shared grid cells."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(colour = "grey20")
  )

guild_class_plot <- ggplot(
  guild_class_proportions,
  aes(
    x = .data$foraging_guild_consensus,
    y = .data$proportion_within_guild,
    fill = .data$change_class
  )
) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  coord_flip() +
  scale_fill_manual(values = class_colours, name = NULL) +
  scale_y_continuous(labels = percent_label, expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = NULL,
    y = "Within-guild proportion of species",
    title = "Occupancy-change classes normalized within foraging guild",
    subtitle = "Each guild sums to 100%, so guilds with more species do not dominate the comparison."
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.y = element_text(colour = "grey20", size = 8.5),
    plot.margin = margin(10, 28, 10, 10)
  )

ggsave(
  filename = richness_plot_path,
  plot = richness_map,
  width = 9.2,
  height = 5.4,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = delta_plot_path,
  plot = delta_map,
  width = 6.5,
  height = 6.2,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = species_class_plot_path,
  plot = species_class_plot,
  width = 7.2,
  height = 4.6,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = guild_class_plot_path,
  plot = guild_class_plot,
  width = 8.2,
  height = max(6.2, 0.26 * length(unique(species_occupancy$foraging_guild_consensus)) + 2.1),
  dpi = 300,
  bg = "white"
)

non_png_plot_outputs <- list.files(
  out_dir,
  pattern = "\\.(pdf|svg|jpg|jpeg|tif|tiff)$",
  ignore.case = TRUE,
  full.names = TRUE
)

if (length(non_png_plot_outputs) > 0) {
  stop(
    "Unexpected non-PNG plot output(s): ",
    paste(non_png_plot_outputs, collapse = ", "),
    call. = FALSE
  )
}

#### CONSOLE SUMMARY FOR QUICK CHECKING ####

print(filter_summary)
print(moran_summary)
print(species_class_proportions)

message("Wrote filter summary to: ", filter_summary_path)
message("Wrote shared-site richness data to: ", richness_csv_path)
message("Wrote delta richness data to: ", delta_csv_path)
message("Wrote Moran's I summary to: ", moran_csv_path)
message("Wrote species occupancy changes to: ", species_change_csv_path)
message("Wrote top gains/losses to: ", species_top_change_csv_path)
message("Wrote species class proportions to: ", species_class_csv_path)
message("Wrote guild-normalized proportions to: ", guild_class_csv_path)
message("Wrote PNG richness map to: ", richness_plot_path)
message("Wrote PNG delta richness map to: ", delta_plot_path)
message("Wrote PNG species class plot to: ", species_class_plot_path)
message("Wrote PNG guild class plot to: ", guild_class_plot_path)
