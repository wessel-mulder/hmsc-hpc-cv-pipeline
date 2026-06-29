#### 2010s species with many supported negative Omega associations ####
#
# Purpose:
#   Rank species by the number of supported negative Omega associations they
#   have in the 2010s Atlas 3 model, then map whether the highest-ranked
#   species occupy the same model-included atlas sites.
#
# Important analysis choices:
#   - "Negative cooccurrence" here means a supported negative residual Omega
#     association in the existing species-pair table.
#   - This script only uses the 2010s period for now.
#   - Occurrence maps use the same filtered model input file as the Omega
#     workflow, so sites and species match the HMSC model inputs.
#   - Plot output is PNG only, following the project convention.

required_packages <- c("tidyverse", "sf", "ggplot2")
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
  library(sf)
  library(ggplot2)
})

#### PATHS ####

pattern <- "2026-03-13"
out_dir <- file.path(
  "notebooks", "exploratory", "outputs", "species-associations-omega"
)
plot_dir <- file.path(out_dir, "plots")

pair_table_path <- file.path(out_dir, paste0(pattern, "-omega-species-pairs.csv"))
model_input_path <- file.path(
  "Data", "preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData"
)
atlas_grid_path <- file.path(
  "Data", "data", "1_preprocessing", "atlas-grids",
  "grids-ocean-thresholds", "grids_ocean_thresholds.shp"
)

species_rank_csv <- file.path(
  out_dir, paste0(pattern, "-omega-2010s-top-negative-degree-species.csv")
)
partner_csv <- file.path(
  out_dir, paste0(pattern, "-omega-2010s-top-negative-degree-partners.csv")
)
site_summary_csv <- file.path(
  out_dir, paste0(pattern, "-omega-2010s-top-negative-species-site-summary.csv")
)
site_species_csv <- file.path(
  out_dir, paste0(pattern, "-omega-2010s-top-negative-species-site-presence.csv")
)

count_map_png <- file.path(
  plot_dir, paste0(pattern, "-omega-2010s-top-negative-species-richness-map.png")
)
facet_map_png <- file.path(
  plot_dir, paste0(pattern, "-omega-2010s-top-negative-species-presence-facets.png")
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

#### CONSTANTS ####

target_period <- "2010s"
target_atlas <- 3
support_class_negative <- "Supported negative"
n_top_species <- 10

#### HELPERS ####

save_png <- function(plot, path, width, height, dpi = 300) {
  if (tolower(tools::file_ext(path)) != "png") {
    stop("Plot path must end in .png: ", path, call. = FALSE)
  }

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

clean_species_label <- function(x) {
  str_replace_all(x, "_", " ")
}

#### LOAD AND RANK SPECIES BY NEGATIVE OMEGA DEGREE ####

if (!file.exists(pair_table_path)) {
  stop("Missing Omega species-pair table: ", pair_table_path, call. = FALSE)
}

omega_pairs <- read_csv(
  pair_table_path,
  col_types = cols(
    period = col_character(),
    support_class = col_character(),
    .default = col_guess()
  )
)

omega_2010s_negative <- omega_pairs |>
  filter(
    .data$period == target_period,
    .data$support_class == support_class_negative
  )

if (nrow(omega_2010s_negative) == 0) {
  stop("No supported negative pairs found for ", target_period, call. = FALSE)
}

# Each Omega pair is symmetric. Bind both directions so species degree is the
# number of distinct supported negative partners for each focal species.
negative_pair_degree_long <- bind_rows(
  omega_2010s_negative |>
    transmute(
      species = .data$species_a,
      partner_species = .data$species_b,
      omega_mean = .data$omega_mean,
      pr_positive = .data$pr_positive,
      pr_negative = .data$pr_negative,
      pair_id = .data$pair_id
    ),
  omega_2010s_negative |>
    transmute(
      species = .data$species_b,
      partner_species = .data$species_a,
      omega_mean = .data$omega_mean,
      pr_positive = .data$pr_positive,
      pr_negative = .data$pr_negative,
      pair_id = .data$pair_id
    )
)

species_degree <- negative_pair_degree_long |>
  group_by(.data$species) |>
  summarise(
    n_supported_negative_partners = n_distinct(.data$partner_species),
    median_negative_omega = median(.data$omega_mean, na.rm = TRUE),
    mean_negative_omega = mean(.data$omega_mean, na.rm = TRUE),
    strongest_negative_omega = min(.data$omega_mean, na.rm = TRUE),
    median_pr_negative = median(.data$pr_negative, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(
    desc(.data$n_supported_negative_partners),
    .data$strongest_negative_omega,
    .data$species
  ) |>
  mutate(rank_negative_degree = row_number())

top_species <- species_degree |>
  slice_head(n = n_top_species)

top_species_names <- top_species$species

partner_table <- negative_pair_degree_long |>
  filter(.data$species %in% top_species_names) |>
  left_join(
    top_species |> select("species", "rank_negative_degree"),
    by = "species"
  ) |>
  arrange(.data$rank_negative_degree, .data$omega_mean, .data$partner_species)

#### LOAD MODEL OCCURRENCES AND SITE METADATA ####

if (!file.exists(model_input_path)) {
  stop("Missing model input file: ", model_input_path, call. = FALSE)
}

load(model_input_path)

if (!all(c("Design", "Y", "Tr") %in% ls())) {
  stop("Expected Design, Y, and Tr in ", model_input_path, call. = FALSE)
}

if (!all(rownames(Design) %in% rownames(Y))) {
  stop("Some Design row names are missing from Y.", call. = FALSE)
}

missing_top_species <- setdiff(top_species_names, colnames(Y))
if (length(missing_top_species) > 0) {
  stop(
    "Top negative-degree species missing from Y: ",
    paste(missing_top_species, collapse = ", "),
    call. = FALSE
  )
}

trait_table <- Tr |>
  rownames_to_column("species")

atlas3_surveys <- rownames(Design)[Design$atlas == target_atlas]
Y_atlas3 <- Y[atlas3_surveys, top_species_names, drop = FALSE]
Design_atlas3 <- Design[atlas3_surveys, , drop = FALSE]

site_metadata <- tibble(
  survey = rownames(Design_atlas3),
  site = as.character(Design_atlas3$site),
  atlas = Design_atlas3$atlas,
  year = Design_atlas3$year,
  easting = Design_atlas3$lon,
  northing = Design_atlas3$lat
)

site_presence <- as_tibble(Y_atlas3, rownames = "survey") |>
  pivot_longer(
    cols = all_of(top_species_names),
    names_to = "species",
    values_to = "present"
  ) |>
  mutate(present = as.integer(.data$present > 0)) |>
  left_join(site_metadata, by = "survey") |>
  left_join(
    top_species |>
      select(
        "species",
        "rank_negative_degree",
        "n_supported_negative_partners"
      ),
    by = "species"
  ) |>
  left_join(trait_table, by = "species") |>
  mutate(
    species_label = paste0(
      .data$rank_negative_degree,
      ". ",
      clean_species_label(.data$species),
      " (",
      .data$n_supported_negative_partners,
      ")"
    ),
    species_label = factor(
      .data$species_label,
      levels = unique(.data$species_label[order(.data$rank_negative_degree)])
    )
  )

site_summary <- site_presence |>
  group_by(.data$site, .data$atlas, .data$year, .data$easting, .data$northing) |>
  summarise(
    n_top_negative_species_present = sum(.data$present, na.rm = TRUE),
    top_negative_species_present = paste(
      clean_species_label(.data$species[.data$present == 1]),
      collapse = "; "
    ),
    .groups = "drop"
  ) |>
  mutate(
    any_top_negative_species_present = .data$n_top_negative_species_present > 0
  )

species_occ_summary <- site_presence |>
  group_by(
    .data$species,
    .data$species_label,
    .data$rank_negative_degree,
    .data$n_supported_negative_partners
  ) |>
  summarise(
    atlas3_sites_present = sum(.data$present, na.rm = TRUE),
    atlas3_sites_total = n_distinct(.data$site),
    atlas3_site_occupancy_percent =
      100 * .data$atlas3_sites_present / .data$atlas3_sites_total,
    .groups = "drop"
  )

site_count_summary <- site_summary |>
  count(.data$n_top_negative_species_present, name = "n_sites") |>
  mutate(percent_sites = 100 * .data$n_sites / sum(.data$n_sites)) |>
  arrange(.data$n_top_negative_species_present)

species_rank_out <- top_species |>
  left_join(trait_table, by = "species") |>
  left_join(
    species_occ_summary |>
      select(
        "species",
        "atlas3_sites_present",
        "atlas3_sites_total",
        "atlas3_site_occupancy_percent"
      ),
    by = "species"
  ) |>
  arrange(.data$rank_negative_degree)

#### WRITE TABLE OUTPUTS ####

write_csv(species_rank_out, species_rank_csv, na = "")
write_csv(partner_table, partner_csv, na = "")
write_csv(site_summary, site_summary_csv, na = "")
write_csv(site_presence, site_species_csv, na = "")

#### MAP OUTPUTS ####

if (!file.exists(atlas_grid_path)) {
  stop("Missing atlas grid shapefile: ", atlas_grid_path, call. = FALSE)
}

atlas_grid <- st_read(atlas_grid_path, quiet = TRUE)

site_map <- atlas_grid |>
  inner_join(site_summary, by = c("kvdrtkd" = "site"))

site_species_map <- atlas_grid |>
  inner_join(site_presence, by = c("kvdrtkd" = "site"))

if (nrow(site_map) != nrow(site_summary)) {
  warning(
    "Site map rows do not exactly match site summary rows: map=",
    nrow(site_map),
    " summary=",
    nrow(site_summary)
  )
}

count_breaks <- seq(
  0,
  max(site_summary$n_top_negative_species_present, na.rm = TRUE),
  by = 1
)

count_colours <- c(
  "#f7f7f7", "#dbe9f6", "#a6bddb", "#74a9cf",
  "#2b8cbe", "#045a8d", "#023858"
)

count_map <- ggplot() +
  geom_sf(data = atlas_grid, fill = "grey95", colour = "white", linewidth = 0.05) +
  geom_sf(
    data = site_map,
    aes(fill = .data$n_top_negative_species_present),
    colour = "grey85",
    linewidth = 0.03
  ) +
  scale_fill_gradientn(
    colours = count_colours,
    limits = c(0, max(count_breaks)),
    breaks = count_breaks,
    name = "Top-ten species\npresent"
  ) +
  coord_sf(datum = NA, expand = FALSE) +
  labs(
    title = "Top negative-degree species, 2010s",
    subtitle = "Colour = number of top-ten species recorded at each model-included atlas site"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot",
    plot.margin = margin(8, 8, 8, 8),
    legend.position = "right"
  )

facet_map <- ggplot() +
  geom_sf(data = atlas_grid, fill = "grey96", colour = "white", linewidth = 0.03) +
  geom_sf(
    data = site_species_map,
    aes(
      fill = factor(
        .data$present,
        levels = c(0, 1),
        labels = c("Absent", "Present")
      )
    ),
    colour = NA
  ) +
  facet_wrap(~ species_label, ncol = 5) +
  scale_fill_manual(
    values = c("Absent" = "grey88", "Present" = "#2166ac"),
    name = NULL
  ) +
  coord_sf(datum = NA, expand = FALSE) +
  labs(
    title = "Occurrence footprints of top negative-degree species, 2010s",
    subtitle = "Facet-title numbers are supported negative Omega partner counts"
  ) +
  theme_void(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.title.position = "plot",
    plot.margin = margin(8, 8, 8, 8),
    strip.text = element_text(face = "bold", size = 7.5),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines")
  )

save_png(count_map, count_map_png, width = 8.8, height = 7.2)
save_png(facet_map, facet_map_png, width = 12.5, height = 6.6)

#### CONSOLE SUMMARY ####

cat(species_rank_csv, "\n")
cat(partner_csv, "\n")
cat(site_summary_csv, "\n")
cat(site_species_csv, "\n")
cat(count_map_png, "\n")
cat(facet_map_png, "\n")

cat("\nTop species:\n")
print(
  species_rank_out |>
    select(
      "rank_negative_degree",
      "species",
      "n_supported_negative_partners",
      "atlas3_sites_present",
      "atlas3_site_occupancy_percent",
      "foraging_guild_consensus",
      "Migration_a3_DOF"
    )
)

cat("\nTop-ten species count per 2010s site:\n")
print(site_count_summary)
