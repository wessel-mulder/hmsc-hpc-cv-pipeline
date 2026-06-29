rm(list = ls())

# Guild-level niche explainability from the HMSC variance-partitioning matrix.
#
# VP$vals is a species-level matrix: rows are variance components and columns
# are species. Because columns sum to 1, the measured-environment niche fraction
# is the sum of the environmental rows, while "Random: site" is the residual
# spatial/unexplained-by-measured-environment fraction.

library(tidyverse)
library(scales)

source(file.path("support_scripts", "figure_data_helpers.R"))

model_dir <- "HmscOutputs"
model_pattern <- "2026-03-13"
min_species_per_guild_plot <- 5

output_dir <- file.path("misc-figures", "outputs", "main")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

figure_slug <- paste0(model_pattern, "-guild-niche-explainability-vp")
output_png <- file.path(output_dir, paste0(figure_slug, ".png"))
output_species_csv <- file.path(output_dir, paste0(figure_slug, "-species.csv"))
output_summary_csv <- file.path(output_dir, paste0(figure_slug, "-summary.csv"))
output_trend_csv <- file.path(output_dir, paste0(figure_slug, "-trend-flags.csv"))

atlas_lookup <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

atlas_order <- unname(atlas_lookup)

atlas_index_lookup <- c(
  "1970s" = 1,
  "1990s" = 2,
  "2010s" = 3
)

read_vp_species <- function(model_folder, model) {
  atlas_num <- str_extract(model_folder, "(?<=Atlas)\\d+")
  atlas_label <- atlas_lookup[[atlas_num]]
  vp_path <- file.path(
    model_dir,
    model_folder,
    "Results",
    paste0(model_folder, "VP.rds")
  )

  if (!file.exists(vp_path)) {
    stop("Missing VP.rds: ", vp_path, call. = FALSE)
  }

  vp <- readRDS(vp_path)
  vp_matrix <- vp$vals
  random_row <- "Random: site"

  if (!random_row %in% rownames(vp_matrix)) {
    stop("VP matrix has no `Random: site` row: ", vp_path, call. = FALSE)
  }

  species_traits <- model$TrData |>
    as_tibble(rownames = "species") |>
    select(
      species,
      foraging_guild = "foraging_guild_consensus"
    )

  environmental_rows <- setdiff(rownames(vp_matrix), random_row)

  tibble(
    species = colnames(vp_matrix),
    atlas = atlas_label,
    measured_environment_explained = colSums(vp_matrix[environmental_rows, , drop = FALSE]),
    unexplained_or_residual_spatial = vp_matrix[random_row, ],
    vp_column_sum = colSums(vp_matrix)
  ) |>
    mutate(
      atlas = atlas_label
    ) |>
    select(
      atlas,
      species,
      measured_environment_explained,
      unexplained_or_residual_spatial,
      vp_column_sum
    ) |>
    left_join(species_traits, by = "species")
}

model_folders <- figure_model_folders(pattern = model_pattern, base_dir = model_dir)
models <- load_hmsc_posteriors(model_folders, base_dir = model_dir)

species_vp <- map2_dfr(model_folders, models, read_vp_species) |>
  mutate(
    atlas = factor(.data$atlas, levels = atlas_order),
    atlas_index = atlas_index_lookup[as.character(.data$atlas)]
  )

if (any(abs(species_vp$vp_column_sum - 1) > 1e-8, na.rm = TRUE)) {
  stop("At least one species VP column does not sum to 1.", call. = FALSE)
}

guild_summary <- species_vp |>
  group_by(.data$atlas, .data$atlas_index, .data$foraging_guild) |>
  summarise(
    n_species = n(),
    mean_environment_explained = mean(.data$measured_environment_explained, na.rm = TRUE),
    sd_environment_explained = sd(.data$measured_environment_explained, na.rm = TRUE),
    mean_unexplained = mean(.data$unexplained_or_residual_spatial, na.rm = TRUE),
    sd_unexplained = sd(.data$unexplained_or_residual_spatial, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(.data$foraging_guild, .data$atlas)

guild_trends <- guild_summary |>
  group_by(.data$foraging_guild) |>
  summarise(
    n_species = first(.data$n_species),
    mean_environment_explained_1970s =
      .data$mean_environment_explained[.data$atlas == "1970s"][[1]],
    mean_environment_explained_2010s =
      .data$mean_environment_explained[.data$atlas == "2010s"][[1]],
    delta_environment_explained =
      .data$mean_environment_explained_2010s - .data$mean_environment_explained_1970s,
    increasingly_explained = .data$delta_environment_explained > 0,
    mean_unexplained_1970s = .data$mean_unexplained[.data$atlas == "1970s"][[1]],
    mean_unexplained_2010s = .data$mean_unexplained[.data$atlas == "2010s"][[1]],
    delta_unexplained = .data$mean_unexplained_2010s - .data$mean_unexplained_1970s,
    .groups = "drop"
  ) |>
  arrange(desc(.data$delta_environment_explained))

write_csv(
  species_vp |>
    mutate(atlas = as.character(.data$atlas)),
  output_species_csv
)
write_csv(
  guild_summary |>
    mutate(atlas = as.character(.data$atlas)),
  output_summary_csv
)
write_csv(guild_trends, output_trend_csv)

plot_guilds <- guild_trends |>
  filter(.data$n_species >= min_species_per_guild_plot) |>
  arrange(desc(.data$delta_environment_explained)) |>
  pull(.data$foraging_guild)

plot_summary <- guild_summary |>
  filter(.data$foraging_guild %in% plot_guilds) |>
  mutate(
    foraging_guild = factor(.data$foraging_guild, levels = rev(plot_guilds)),
    increasingly_explained = .data$foraging_guild %in%
      guild_trends$foraging_guild[guild_trends$increasingly_explained]
  )

y_limit <- max(
  plot_summary$mean_environment_explained + plot_summary$sd_environment_explained,
  na.rm = TRUE
) * 1.08

guild_plot <- ggplot(
  plot_summary,
  aes(
    x = atlas,
    y = mean_environment_explained,
    group = foraging_guild,
    colour = increasingly_explained
  )
) +
  geom_errorbar(
    aes(
      ymin = pmax(0, .data$mean_environment_explained - .data$sd_environment_explained),
      ymax = pmin(1, .data$mean_environment_explained + .data$sd_environment_explained)
    ),
    width = 0.12,
    linewidth = 0.35,
    colour = "grey35"
  ) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.8) +
  facet_wrap(vars(foraging_guild), ncol = 4) +
  scale_colour_manual(
    values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02"),
    breaks = c(TRUE, FALSE),
    labels = c("Yes", "No"),
    name = "Higher in\n2010s?"
  ) +
  scale_y_continuous(
    limits = c(0, y_limit),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    title = "Guild niche variation explained by measured environments",
    x = NULL,
    y = "Mean species-level\nexplained niche variation"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 7.6, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(size = 7.9, face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = "white"),
    legend.position = "bottom",
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 8.3),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.spacing = grid::unit(0.65, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

ggsave(output_png, guild_plot, width = 8.8, height = 6.6, units = "in", dpi = 300, bg = "white")

message("Wrote: ", output_png)
message("Wrote: ", output_species_csv)
message("Wrote: ", output_summary_csv)
message("Wrote: ", output_trend_csv)
