# Quick atlas species-turnover analysis from the raw occurrence matrix.
#
# Purpose:
#   Build national species lists for each Danish atlas period from
#   Data/Y_occurrences/Y_occurrences.csv, count species gained and lost between
#   adjacent atlas periods, and compare the per-decade rates with the published
#   Jarvinen and Ulfstrand 1980 baseline.
#
# Important analysis choices:
#   - Species are treated as present in an atlas if they occur in at least one
#     raw survey row for that atlas.
#   - This is a national species-pool comparison, not a common-site comparison.
#   - Atlas 1 -> 2 and Atlas 2 -> 3 are each treated as rounded two-decade
#     intervals, so rate per decade = species count / 2.
#   - Plot output is PNG only, following the project plotting instruction.

rm(list = ls())

library(ggplot2)

#### PATHS ####

occurrence_path <- file.path("Data", "Y_occurrences", "Y_occurrences.csv")
out_dir <- file.path("notebooks", "exploratory", "outputs", "atlas-species-turnover")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species_lists_path <- file.path(out_dir, "atlas-species-lists.csv")
turnover_summary_path <- file.path(out_dir, "atlas-species-turnover-summary.csv")
turnover_species_path <- file.path(out_dir, "atlas-species-turnover-species.csv")
timeline_png_path <- file.path(out_dir, "atlas-species-turnover-timeline.png")

#### CONSTANTS ####

atlas_labels <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

atlas_years <- c(
  "1" = 1970,
  "2" = 1990,
  "3" = 2010
)

atlas_comparisons <- data.frame(
  comparison = c("1970s to 1990s", "1990s to 2010s"),
  start_atlas = c("1", "2"),
  end_atlas = c("2", "3"),
  decades_used_for_rate = c(2, 2),
  stringsAsFactors = FALSE
)

historical_baseline <- data.frame(
  source = "Jarvinen and Ulfstrand 1980",
  comparison = "1850 to 1970",
  start_label = "1850",
  end_label = "1970",
  start_year = 1850,
  end_year = 1970,
  decades_used_for_rate = NA_real_,
  gained_species = NA_integer_,
  lost_species = NA_integer_,
  net_species_change = NA_integer_,
  total_turnover = NA_integer_,
  gains_per_decade = 2.8,
  losses_per_decade = 0.6,
  rate_source = "published per-decade rate",
  stringsAsFactors = FALSE
)

#### HELPERS ####

parse_atlas_from_survey <- function(survey_id) {
  # Survey IDs are of the form SITE_ATLAS, for example AD90_1.
  # Splitting from the right is safer than relying on a fixed site-code length.
  vapply(
    strsplit(survey_id, "_", fixed = TRUE),
    function(parts) tail(parts, 1),
    character(1)
  )
}

build_species_list <- function(y_atlas, species_columns, atlas_code) {
  # A species is in the national list for an atlas if it is recorded in at
  # least one survey row in that atlas period.
  observed_counts <- colSums(y_atlas[, species_columns, drop = FALSE], na.rm = TRUE)
  present_species <- names(observed_counts)[observed_counts > 0]
  n_surveys <- nrow(y_atlas)

  data.frame(
    atlas = atlas_code,
    period = unname(atlas_labels[atlas_code]),
    species = present_species,
    occupied_surveys = as.integer(observed_counts[present_species]),
    total_surveys = n_surveys,
    occupancy_proportion = as.numeric(observed_counts[present_species] / n_surveys),
    stringsAsFactors = FALSE
  )
}

turnover_species_frame <- function(start_species, end_species, comparison_label) {
  gained <- sort(setdiff(end_species, start_species))
  lost <- sort(setdiff(start_species, end_species))

  data.frame(
    comparison = comparison_label,
    turnover_type = c(rep("gained", length(gained)), rep("lost", length(lost))),
    species = c(gained, lost),
    stringsAsFactors = FALSE
  )
}

#### LOAD AND CHECK RAW OCCURRENCE MATRIX ####

y_raw <- read.csv(occurrence_path, check.names = FALSE, stringsAsFactors = FALSE)

if (!"survey" %in% names(y_raw)) {
  stop("Expected a 'survey' column in ", occurrence_path, call. = FALSE)
}

species_columns <- setdiff(names(y_raw), "survey")

if (length(species_columns) != 201) {
  stop(
    "Expected 201 species columns, but found ",
    length(species_columns),
    ". Check that the raw Y_occurrences file has not changed.",
    call. = FALSE
  )
}

y_raw$atlas <- parse_atlas_from_survey(y_raw$survey)

unexpected_atlas_codes <- setdiff(unique(y_raw$atlas), names(atlas_labels))
if (length(unexpected_atlas_codes) > 0) {
  stop(
    "Unexpected atlas code(s) parsed from survey IDs: ",
    paste(unexpected_atlas_codes, collapse = ", "),
    call. = FALSE
  )
}

#### BUILD NATIONAL ATLAS SPECIES LISTS ####

species_lists <- do.call(
  rbind,
  lapply(names(atlas_labels), function(atlas_code) {
    build_species_list(
      y_atlas = y_raw[y_raw$atlas == atlas_code, , drop = FALSE],
      species_columns = species_columns,
      atlas_code = atlas_code
    )
  })
)

species_by_atlas <- split(species_lists$species, species_lists$atlas)

#### COUNT GAINS AND LOSSES BETWEEN ADJACENT ATLASES ####

turnover_species <- do.call(
  rbind,
  lapply(seq_len(nrow(atlas_comparisons)), function(i) {
    contrast <- atlas_comparisons[i, ]
    turnover_species_frame(
      start_species = species_by_atlas[[contrast$start_atlas]],
      end_species = species_by_atlas[[contrast$end_atlas]],
      comparison_label = contrast$comparison
    )
  })
)

atlas_summary <- do.call(
  rbind,
  lapply(seq_len(nrow(atlas_comparisons)), function(i) {
    contrast <- atlas_comparisons[i, ]
    start_species <- species_by_atlas[[contrast$start_atlas]]
    end_species <- species_by_atlas[[contrast$end_atlas]]

    gained_species <- length(setdiff(end_species, start_species))
    lost_species <- length(setdiff(start_species, end_species))

    data.frame(
      source = "Raw Y_occurrences national list",
      comparison = contrast$comparison,
      start_label = unname(atlas_labels[contrast$start_atlas]),
      end_label = unname(atlas_labels[contrast$end_atlas]),
      start_year = unname(atlas_years[contrast$start_atlas]),
      end_year = unname(atlas_years[contrast$end_atlas]),
      decades_used_for_rate = contrast$decades_used_for_rate,
      gained_species = gained_species,
      lost_species = lost_species,
      net_species_change = gained_species - lost_species,
      total_turnover = gained_species + lost_species,
      gains_per_decade = gained_species / contrast$decades_used_for_rate,
      losses_per_decade = lost_species / contrast$decades_used_for_rate,
      rate_source = "count divided by rounded two-decade interval",
      stringsAsFactors = FALSE
    )
  })
)

turnover_summary <- rbind(historical_baseline, atlas_summary)

#### INTERNAL TESTS ####

if (!all(c("1", "2", "3") %in% names(species_by_atlas))) {
  stop("Expected species lists for atlas codes 1, 2, and 3.", call. = FALSE)
}

for (i in seq_len(nrow(atlas_comparisons))) {
  contrast <- atlas_comparisons[i, ]
  start_species <- species_by_atlas[[contrast$start_atlas]]
  end_species <- species_by_atlas[[contrast$end_atlas]]
  gained <- setdiff(end_species, start_species)
  lost <- setdiff(start_species, end_species)

  if (length(intersect(gained, lost)) > 0) {
    stop("Set-operation check failed for ", contrast$comparison, call. = FALSE)
  }

  summary_row <- atlas_summary[atlas_summary$comparison == contrast$comparison, ]
  if (!isTRUE(all.equal(summary_row$gained_species / 2, summary_row$gains_per_decade))) {
    stop("Gain-rate check failed for ", contrast$comparison, call. = FALSE)
  }
  if (!isTRUE(all.equal(summary_row$lost_species / 2, summary_row$losses_per_decade))) {
    stop("Loss-rate check failed for ", contrast$comparison, call. = FALSE)
  }
}

#### WRITE TABLE OUTPUTS ####

write.csv(species_lists, species_lists_path, row.names = FALSE, na = "")
write.csv(turnover_summary, turnover_summary_path, row.names = FALSE, na = "")
write.csv(turnover_species, turnover_species_path, row.names = FALSE, na = "")

#### PNG TIMELINE FIGURE ####

# The historical 1850-1970 segment is intentionally compressed to keep the two
# modern atlas intervals readable. The plotted x positions therefore represent a
# manuscript sketch timeline rather than a strictly linear year axis.
timeline_positions <- data.frame(
  comparison = c("1850 to 1970", "1970s to 1990s", "1990s to 2010s"),
  x_start = c(0, 3, 4.5),
  x_end = c(2, 4.5, 6),
  interval_group = c("Published baseline", "Danish atlas estimates", "Danish atlas estimates"),
  stringsAsFactors = FALSE
)

timeline_plot_data <- merge(turnover_summary, timeline_positions, by = "comparison", all.x = TRUE)
timeline_plot_data$interval_mid <- (timeline_plot_data$x_start + timeline_plot_data$x_end) / 2

rate_plot_data <- rbind(
  transform(
    timeline_plot_data,
    turnover_type = "Gains",
    rate_per_decade = gains_per_decade
  ),
  transform(
    timeline_plot_data,
    turnover_type = "Losses",
    rate_per_decade = losses_per_decade
  )
)

rate_plot_data$turnover_type <- factor(rate_plot_data$turnover_type, levels = c("Gains", "Losses"))

max_rate <- max(rate_plot_data$rate_per_decade, na.rm = TRUE)

timeline_plot <- ggplot(
  rate_plot_data,
  aes(
    x = x_start,
    xend = x_end,
    y = rate_per_decade,
    yend = rate_per_decade,
    colour = turnover_type
  )
) +
  geom_vline(
    xintercept = c(0, 2, 3, 4.5, 6),
    colour = "grey88",
    linewidth = 0.35
  ) +
  geom_segment(linewidth = 2.2, lineend = "round") +
  geom_point(aes(x = interval_mid), size = 3.2) +
  geom_text(
    aes(
      x = interval_mid,
      y = rate_per_decade + 0.35,
      label = sprintf("%.1f / decade", rate_per_decade)
    ),
    colour = "grey15",
    size = 3.5,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = 1,
    y = max_rate + 1.2,
    label = "Published baseline",
    colour = "grey25",
    fontface = "bold",
    size = 3.7
  ) +
  annotate(
    "text",
    x = 4.5,
    y = max_rate + 1.2,
    label = "Raw atlas estimates",
    colour = "grey25",
    fontface = "bold",
    size = 3.7
  ) +
  annotate(
    "text",
    x = 2.5,
    y = -0.28,
    label = "//",
    colour = "grey35",
    size = 8
  ) +
  scale_colour_manual(
    values = c("Gains" = "#287D60", "Losses" = "#B13B3D"),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(0, 2, 3, 4.5, 6),
    labels = c("1850", "1970", "1970s", "1990s", "2010s"),
    limits = c(-0.15, 6.15),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(-0.5, max_rate + 1.6),
    breaks = seq(0, ceiling(max_rate + 1), by = 1),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Species turnover rate (species / decade)",
    title = "National breeding-bird species gains and losses through time",
    subtitle = "The 1850-1970 baseline is compressed; atlas intervals use rounded two-decade rates."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(colour = "grey20"),
    axis.text.y = element_text(colour = "grey20")
  )

ggsave(
  filename = timeline_png_path,
  plot = timeline_plot,
  width = 8.5,
  height = 4.8,
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

print(turnover_summary)

message("Wrote species lists to: ", species_lists_path)
message("Wrote turnover summary to: ", turnover_summary_path)
message("Wrote turnover species table to: ", turnover_species_path)
message("Wrote PNG timeline to: ", timeline_png_path)
