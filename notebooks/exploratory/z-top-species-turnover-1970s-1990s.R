# Observed species turnover maps between the 1970s and 1990s.
#
# Purpose:
#   Identify the 20 species with the largest observed distribution turnover
#   between Atlas 1 (1970s) and Atlas 2 (1990s), then map where each species
#   was gained, lost, retained, or absent across matched grid cells.
#
# Analysis choices:
#   - Uses the preprocessed good-and-average coverage model-input data, so the
#     species list and site filtering match the main HMSC workflow.
#   - Restricts comparisons to grid cells surveyed in both Atlas 1 and Atlas 2.
#   - Ranks species by changed grid cells: gains + losses.
#   - Writes PNG figures only, following the project plotting instruction.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(sf)
  library(stringr)
  library(tibble)
  library(tidyr)
})

#### PATHS ####

preprocessed_path <- file.path(
  "Data",
  "preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData"
)

grid_path <- file.path(
  "Data",
  "data",
  "1_preprocessing",
  "atlas-grids",
  "grids-ocean-thresholds",
  "grids_ocean_thresholds.shp"
)

out_dir <- file.path(
  "notebooks",
  "exploratory",
  "outputs",
  "top-species-turnover-1970s-1990s"
)
plot_dir <- file.path(out_dir, "plots")
individual_map_dir <- file.path(plot_dir, "individual-species-maps")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_map_dir, recursive = TRUE, showWarnings = FALSE)

ranked_csv_path <- file.path(out_dir, "species-turnover-ranked-1970s-1990s.csv")
top20_csv_path <- file.path(out_dir, "top20-species-turnover-1970s-1990s.csv")
map_data_csv_path <- file.path(out_dir, "top20-species-turnover-map-data-1970s-1990s.csv")
guild_summary_csv_path <- file.path(out_dir, "top20-turnover-guild-summary-1970s-1990s.csv")

faceted_map_png_path <- file.path(plot_dir, "top20-species-turnover-maps-1970s-1990s.png")
guild_summary_png_path <- file.path(plot_dir, "top20-species-turnover-guild-summary-1970s-1990s.png")

#### CONSTANTS ####

comparison_label <- "1970s to 1990s"
atlas_start <- "1"
atlas_end <- "2"

turnover_state_levels <- c("Gained", "Lost", "Retained", "Absent in both")
turnover_state_colours <- c(
  "Gained" = "#2C7FB8",
  "Lost" = "#D95F02",
  "Retained" = "#238443",
  "Absent in both" = "#F0F0F0"
)

direction_colours <- c(
  "More occupied cells in 1990s" = "#2C7FB8",
  "Fewer occupied cells in 1990s" = "#D95F02",
  "Balanced gains and losses" = "#666666"
)

#### HELPERS ####

clean_species_label <- function(species_name) {
  # Species names are stored as Genus_species. The figure labels are easier to
  # scan with spaces, while the CSV keeps the original machine-readable name.
  str_replace_all(species_name, "_", " ")
}

safe_filename <- function(label) {
  # Individual map files use a conservative filename so they are portable across
  # systems and easy to sort.
  label |>
    str_replace_all("[^A-Za-z0-9]+", "-") |>
    str_replace_all("(^-|-$)", "") |>
    str_to_lower()
}

#### LOAD DATA ####

load(preprocessed_path)

required_objects <- c("Y", "Design", "Tr")
missing_objects <- setdiff(required_objects, ls())
if (length(missing_objects) > 0) {
  stop(
    "The preprocessed file is missing object(s): ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

if (!all(colnames(Y) %in% rownames(Tr))) {
  stop("Some species in Y are missing from Tr.", call. = FALSE)
}

if (!all(rownames(Y) %in% rownames(Design))) {
  stop("Some survey rows in Y are missing from Design.", call. = FALSE)
}

#### BUILD PAIRED ATLAS MATRICES ####

design_with_ids <- Design |>
  rownames_to_column(var = "survey") |>
  mutate(atlas = as.character(.data$atlas))

sites_start <- design_with_ids |>
  filter(.data$atlas == atlas_start) |>
  pull(.data$site)

sites_end <- design_with_ids |>
  filter(.data$atlas == atlas_end) |>
  pull(.data$site)

common_sites <- sort(intersect(sites_start, sites_end))

if (length(common_sites) == 0) {
  stop("No matched grid cells were found for Atlas 1 and Atlas 2.", call. = FALSE)
}

design_start <- design_with_ids |>
  filter(.data$atlas == atlas_start, .data$site %in% common_sites) |>
  arrange(.data$site)

design_end <- design_with_ids |>
  filter(.data$atlas == atlas_end, .data$site %in% common_sites) |>
  arrange(.data$site)

if (!identical(design_start$site, design_end$site)) {
  stop("Atlas 1 and Atlas 2 site ordering is not aligned.", call. = FALSE)
}

y_start <- as.matrix(Y[design_start$survey, , drop = FALSE])
y_end <- as.matrix(Y[design_end$survey, , drop = FALSE])

if (!identical(colnames(y_start), colnames(y_end))) {
  stop("Atlas 1 and Atlas 2 species columns are not aligned.", call. = FALSE)
}

#### SUMMARISE SPECIES TURNOVER ####

gained_matrix <- y_start == 0 & y_end == 1
lost_matrix <- y_start == 1 & y_end == 0
retained_matrix <- y_start == 1 & y_end == 1
absent_matrix <- y_start == 0 & y_end == 0

species_turnover <- tibble(
  species = colnames(y_start),
  comparison = comparison_label,
  n_common_sites = nrow(y_start),
  occupied_cells_1970s = colSums(y_start, na.rm = TRUE),
  occupied_cells_1990s = colSums(y_end, na.rm = TRUE),
  gained_cells = colSums(gained_matrix, na.rm = TRUE),
  lost_cells = colSums(lost_matrix, na.rm = TRUE),
  retained_cells = colSums(retained_matrix, na.rm = TRUE),
  absent_in_both_cells = colSums(absent_matrix, na.rm = TRUE)
) |>
  mutate(
    changed_cells = .data$gained_cells + .data$lost_cells,
    turnover_fraction_of_common_sites = .data$changed_cells / .data$n_common_sites,
    net_cell_change = .data$gained_cells - .data$lost_cells,
    absolute_net_cell_change = abs(.data$net_cell_change),
    occupancy_fraction_1970s = .data$occupied_cells_1970s / .data$n_common_sites,
    occupancy_fraction_1990s = .data$occupied_cells_1990s / .data$n_common_sites,
    # Sørensen dissimilarity is included as a relative turnover diagnostic, but
    # the ranking below uses changed grid-cell count to emphasize map footprint.
    sorensen_dissimilarity = if_else(
      .data$occupied_cells_1970s + .data$occupied_cells_1990s > 0,
      1 - (2 * .data$retained_cells) /
        (.data$occupied_cells_1970s + .data$occupied_cells_1990s),
      NA_real_
    ),
    direction = case_when(
      .data$net_cell_change > 0 ~ "More occupied cells in 1990s",
      .data$net_cell_change < 0 ~ "Fewer occupied cells in 1990s",
      TRUE ~ "Balanced gains and losses"
    )
  ) |>
  left_join(
    Tr |>
      as.data.frame() |>
      rownames_to_column(var = "species") |>
      transmute(
        .data$species,
        foraging_guild = as.character(.data$foraging_guild_consensus),
        migration_strategy = as.character(.data$Migration_a3_DOF),
        species_thermal_index = .data$species_thermal_index
      ),
    by = "species"
  ) |>
  mutate(
    species_label = clean_species_label(.data$species),
    foraging_guild = if_else(
      is.na(.data$foraging_guild) | .data$foraging_guild == "",
      "Unknown guild",
      .data$foraging_guild
    )
  ) |>
  arrange(
    desc(.data$changed_cells),
    desc(.data$absolute_net_cell_change),
    .data$species
  )

top20_species <- species_turnover |>
  slice_head(n = 20) |>
  mutate(rank = row_number())

#### BUILD MAP DATA ####

map_state_for_species <- function(species_name) {
  start_values <- y_start[, species_name]
  end_values <- y_end[, species_name]

  tibble(
    site = design_start$site,
    species = species_name,
    observed_1970s = as.integer(start_values),
    observed_1990s = as.integer(end_values),
    turnover_state = case_when(
      start_values == 0 & end_values == 1 ~ "Gained",
      start_values == 1 & end_values == 0 ~ "Lost",
      start_values == 1 & end_values == 1 ~ "Retained",
      TRUE ~ "Absent in both"
    )
  )
}

top20_map_data <- lapply(top20_species$species, map_state_for_species) |>
  bind_rows() |>
  left_join(
    top20_species |>
      select(
        species,
        rank,
        species_label,
        foraging_guild,
        changed_cells,
        gained_cells,
        lost_cells,
        direction
      ),
    by = "species"
  ) |>
  mutate(
    turnover_state = factor(.data$turnover_state, levels = turnover_state_levels),
    facet_label = sprintf(
      "%02d. %s\n%s; %s changed",
      .data$rank,
      .data$species_label,
      .data$foraging_guild,
      .data$changed_cells
    )
  )

grid_sf <- st_read(grid_path, quiet = TRUE) |>
  select(site = kvdrtkd, pct_lnd, geometry)

map_sf <- grid_sf |>
  inner_join(top20_map_data, by = "site")

if (nrow(map_sf) == 0) {
  stop("No grid geometries joined to top-20 turnover data.", call. = FALSE)
}

#### WRITE TABLES ####

guild_summary <- top20_species |>
  group_by(.data$foraging_guild) |>
  summarise(
    n_top20_species = n(),
    total_changed_cells = sum(.data$changed_cells),
    total_gained_cells = sum(.data$gained_cells),
    total_lost_cells = sum(.data$lost_cells),
    mean_turnover_fraction = mean(.data$turnover_fraction_of_common_sites),
    species = paste(.data$species, collapse = "; "),
    .groups = "drop"
  ) |>
  arrange(desc(.data$n_top20_species), desc(.data$total_changed_cells), .data$foraging_guild)

write_csv(species_turnover, ranked_csv_path)
write_csv(top20_species, top20_csv_path)
write_csv(st_drop_geometry(map_sf), map_data_csv_path)
write_csv(guild_summary, guild_summary_csv_path)

#### FACETED MAP PANEL ####

faceted_map <- ggplot(map_sf) +
  geom_sf(
    aes(fill = .data$turnover_state),
    colour = "grey92",
    linewidth = 0.02
  ) +
  facet_wrap(vars(.data$facet_label), ncol = 5) +
  scale_fill_manual(
    values = turnover_state_colours,
    drop = FALSE,
    name = NULL
  ) +
  coord_sf(expand = FALSE) +
  labs(
    title = "Top 20 observed species turnover maps, 1970s to 1990s",
    subtitle = "Species ranked by gained plus lost matched grid cells in the good-and-average coverage model-input data.",
    caption = paste0(
      "Matched grid cells: ", length(common_sites),
      ". Grey cells indicate absence in both atlas periods."
    )
  ) +
  theme_void(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.8, "cm"),
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10.5, colour = "grey25", margin = margin(b = 8)),
    plot.caption = element_text(size = 8.5, colour = "grey35", hjust = 0),
    strip.text = element_text(face = "bold", size = 7.5, lineheight = 0.95),
    panel.spacing = unit(0.6, "lines")
  )

ggsave(
  faceted_map_png_path,
  faceted_map,
  width = 14,
  height = 10.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

#### INDIVIDUAL SPECIES MAPS ####

individual_species_map <- function(species_name) {
  species_info <- top20_species |>
    filter(.data$species == species_name) |>
    slice(1)

  species_sf <- map_sf |>
    filter(.data$species == species_name)

  ggplot(species_sf) +
    geom_sf(
      aes(fill = .data$turnover_state),
      colour = "grey90",
      linewidth = 0.03
    ) +
    scale_fill_manual(
      values = turnover_state_colours,
      drop = FALSE,
      name = NULL
    ) +
    coord_sf(expand = FALSE) +
    labs(
      title = sprintf(
        "%02d. %s",
        species_info$rank,
        species_info$species_label
      ),
      subtitle = sprintf(
        "%s | %s changed cells: %s gained, %s lost",
        species_info$foraging_guild,
        species_info$changed_cells,
        species_info$gained_cells,
        species_info$lost_cells
      ),
      caption = comparison_label
    ) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 15, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 10, colour = "grey25", margin = margin(b = 6)),
      plot.caption = element_text(size = 8.5, colour = "grey35", hjust = 0)
    )
}

invisible(lapply(top20_species$species, function(species_name) {
  species_info <- top20_species |>
    filter(.data$species == species_name) |>
    slice(1)

  output_path <- file.path(
    individual_map_dir,
    sprintf(
      "%02d-%s-turnover-map-1970s-1990s.png",
      species_info$rank,
      safe_filename(species_info$species)
    )
  )

  ggsave(
    output_path,
    individual_species_map(species_name),
    width = 6,
    height = 6.8,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}))

#### GUILD SUMMARY FIGURE ####

guild_levels <- top20_species |>
  count(.data$foraging_guild, wt = .data$changed_cells, name = "changed_cells") |>
  arrange(desc(.data$changed_cells)) |>
  pull(.data$foraging_guild)

top20_plot_data <- top20_species |>
  mutate(
    species_label = fct_reorder(.data$species_label, .data$changed_cells),
    foraging_guild = factor(.data$foraging_guild, levels = guild_levels),
    direction = factor(.data$direction, levels = names(direction_colours)),
    direction_x = .data$changed_cells + max(.data$changed_cells) * 0.035
  )

guild_count_data <- top20_species |>
  mutate(foraging_guild = factor(.data$foraging_guild, levels = guild_levels)) |>
  count(.data$foraging_guild, name = "n_species") |>
  arrange(.data$foraging_guild)

guild_palette <- setNames(
  scales::hue_pal(l = 55, c = 90)(length(guild_levels)),
  guild_levels
)

ranked_species_plot <- ggplot(
  top20_plot_data,
  aes(x = .data$changed_cells, y = .data$species_label)
) +
  geom_col(aes(fill = .data$foraging_guild), width = 0.72) +
  geom_point(
    aes(x = .data$direction_x, colour = .data$direction),
    size = 2.4
  ) +
  scale_fill_manual(values = guild_palette, guide = "none") +
  scale_colour_manual(values = direction_colours, name = "Net direction") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = "Changed grid cells (gains + losses)",
    y = NULL,
    title = "Species with the largest 1970s to 1990s turnover"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 8.8)
  )

guild_count_plot <- ggplot(
  guild_count_data,
  aes(x = .data$n_species, y = fct_rev(.data$foraging_guild), fill = .data$foraging_guild)
) +
  geom_col(width = 0.72) +
  geom_text(aes(label = .data$n_species), hjust = -0.25, size = 3.2) +
  scale_fill_manual(values = guild_palette, guide = "none") +
  scale_x_continuous(
    breaks = seq(0, max(guild_count_data$n_species), by = 1),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    x = "Species in top 20",
    y = NULL,
    title = "Guild composition"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 8.8)
  )

guild_summary_plot <- ranked_species_plot + guild_count_plot +
  plot_layout(widths = c(2.35, 1), guides = "collect") +
  plot_annotation(
    subtitle = paste0(
      "Ranking uses observed gains plus losses across ",
      length(common_sites),
      " matched grid cells; colours show each species' foraging guild."
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.subtitle = element_text(size = 10.5, colour = "grey25")
  )

ggsave(
  guild_summary_png_path,
  guild_summary_plot,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

#### OUTPUT CHECKS ####

non_png_plot_outputs <- list.files(
  plot_dir,
  pattern = "\\.(pdf|svg|jpg|jpeg|tif|tiff)$",
  ignore.case = TRUE,
  full.names = TRUE,
  recursive = TRUE
)

if (length(non_png_plot_outputs) > 0) {
  stop(
    "Unexpected non-PNG plot output(s): ",
    paste(non_png_plot_outputs, collapse = ", "),
    call. = FALSE
  )
}

message("Wrote ranked turnover table: ", ranked_csv_path)
message("Wrote top-20 turnover table: ", top20_csv_path)
message("Wrote top-20 map data: ", map_data_csv_path)
message("Wrote guild summary table: ", guild_summary_csv_path)
message("Wrote faceted map PNG: ", faceted_map_png_path)
message("Wrote guild summary PNG: ", guild_summary_png_path)
message("Wrote individual species map PNGs to: ", individual_map_dir)

print(
  top20_species |>
    select(
      rank,
      species,
      foraging_guild,
      changed_cells,
      gained_cells,
      lost_cells,
      net_cell_change
    )
)
