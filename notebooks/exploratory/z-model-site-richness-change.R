# Species-richness change across atlas periods at model-included sites.
#
# Purpose:
#   Show how observed species richness changes across the Danish atlas periods
#   after restricting the raw occurrence matrix to the survey rows included in
#   the HMSC model input.
#
# Important analysis choices:
#   - Richness is counted from the raw 201-species Y_occurrences matrix.
#   - Sites/survey rows are filtered using Design from
#     Data/preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData.
#     That file is produced by S0_Data_Definitions.R after applying the land
#     threshold, Good/Average effort filtering, environmental NA filtering, and
#     model-preparation steps.
#   - The main "All model-included surveys" summary uses every survey row in the
#     model Design object. A "Common sites in all atlases" subset is also saved
#     and plotted so the temporal comparison can be read without changing site
#     membership across periods.
#   - Plot output is PNG only.

rm(list = ls())

library(ggplot2)

#### PATHS ####

raw_occurrence_path <- file.path("Data", "Y_occurrences", "Y_occurrences.csv")
model_input_path <- file.path("Data", "preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData")
out_dir <- file.path("notebooks", "exploratory", "outputs", "model-site-richness-change")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

survey_richness_path <- file.path(out_dir, "model-site-raw-richness-by-survey.csv")
summary_path <- file.path(out_dir, "model-site-raw-richness-summary.csv")
paired_changes_path <- file.path(out_dir, "model-site-raw-richness-paired-changes.csv")
plot_path <- file.path(out_dir, "model-site-raw-richness-by-atlas.png")

#### CONSTANTS ####

period_lookup <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

period_levels <- unname(period_lookup)

#### HELPERS ####

summarise_richness <- function(df) {
  # Keep this summary in base R so the script is easy to run on a fresh machine.
  split_df <- split(df, list(df$analysis_set, df$atlas), drop = TRUE)

  do.call(
    rbind,
    lapply(split_df, function(x) {
      data.frame(
        analysis_set = unique(x$analysis_set),
        atlas = unique(x$atlas),
        period = unique(as.character(x$period)),
        n_surveys = nrow(x),
        n_sites = length(unique(x$site)),
        mean_richness = mean(x$species_richness, na.rm = TRUE),
        median_richness = median(x$species_richness, na.rm = TRUE),
        sd_richness = sd(x$species_richness, na.rm = TRUE),
        se_richness = sd(x$species_richness, na.rm = TRUE) / sqrt(nrow(x)),
        q25_richness = unname(quantile(x$species_richness, 0.25, na.rm = TRUE)),
        q75_richness = unname(quantile(x$species_richness, 0.75, na.rm = TRUE)),
        min_richness = min(x$species_richness, na.rm = TRUE),
        max_richness = max(x$species_richness, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
}

#### LOAD DATA ####

raw_y <- read.csv(
  raw_occurrence_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

load(model_input_path)

if (!all(c("Y", "Design") %in% ls())) {
  stop("Expected objects Y and Design in ", model_input_path, call. = FALSE)
}

model_surveys <- rownames(Design)
missing_model_surveys <- setdiff(model_surveys, rownames(raw_y))

if (length(missing_model_surveys) > 0) {
  stop(
    "Some model Design surveys are absent from raw Y_occurrences: ",
    paste(head(missing_model_surveys, 10), collapse = ", "),
    call. = FALSE
  )
}

if (ncol(raw_y) != 201) {
  stop(
    "Expected 201 species columns in raw Y_occurrences, but found ",
    ncol(raw_y),
    ".",
    call. = FALSE
  )
}

#### CALCULATE RAW SPECIES RICHNESS AT MODEL-INCLUDED SURVEYS ####

raw_y_model_rows <- raw_y[model_surveys, , drop = FALSE]

survey_richness <- data.frame(
  survey = model_surveys,
  site = as.character(Design$site),
  atlas = as.character(Design$atlas),
  period = factor(unname(period_lookup[as.character(Design$atlas)]), levels = period_levels),
  model_year = Design$year,
  species_richness = rowSums(raw_y_model_rows, na.rm = TRUE),
  raw_species_pool_size = ncol(raw_y_model_rows),
  stringsAsFactors = FALSE
)

survey_richness$analysis_set <- "All model-included surveys"

sites_by_atlas <- split(survey_richness$site, survey_richness$atlas)
common_sites <- Reduce(intersect, sites_by_atlas[names(period_lookup)])

common_site_richness <- survey_richness[survey_richness$site %in% common_sites, ]
common_site_richness$analysis_set <- "Common sites in all atlases"

plot_richness <- rbind(survey_richness, common_site_richness)
plot_richness$analysis_set <- factor(
  plot_richness$analysis_set,
  levels = c("All model-included surveys", "Common sites in all atlases")
)

richness_summary <- summarise_richness(plot_richness)
richness_summary$analysis_set <- factor(
  richness_summary$analysis_set,
  levels = levels(plot_richness$analysis_set)
)
richness_summary$period <- factor(richness_summary$period, levels = period_levels)
richness_summary <- richness_summary[order(richness_summary$analysis_set, richness_summary$atlas), ]

#### PAIRED COMMON-SITE CHANGES ####

common_wide <- reshape(
  common_site_richness[, c("site", "period", "species_richness")],
  idvar = "site",
  timevar = "period",
  direction = "wide"
)

names(common_wide) <- sub("^species_richness\\.", "richness_", names(common_wide))

paired_changes <- data.frame(
  site = common_wide$site,
  richness_1970s = common_wide$richness_1970s,
  richness_1990s = common_wide$richness_1990s,
  richness_2010s = common_wide$richness_2010s,
  delta_1990s_minus_1970s = common_wide$richness_1990s - common_wide$richness_1970s,
  delta_2010s_minus_1990s = common_wide$richness_2010s - common_wide$richness_1990s,
  delta_2010s_minus_1970s = common_wide$richness_2010s - common_wide$richness_1970s,
  stringsAsFactors = FALSE
)

#### INTERNAL CHECKS ####

if (nrow(survey_richness) != nrow(Design)) {
  stop("Survey-richness row count does not match model Design row count.", call. = FALSE)
}

if (!identical(sort(unique(as.character(survey_richness$period))), sort(period_levels))) {
  stop("Expected richness rows for all three atlas periods.", call. = FALSE)
}

if (anyNA(survey_richness$species_richness)) {
  stop("Missing species-richness values found.", call. = FALSE)
}

if (length(common_sites) == 0) {
  stop("No sites are shared across all three atlas periods.", call. = FALSE)
}

#### WRITE TABLE OUTPUTS ####

write.csv(survey_richness, survey_richness_path, row.names = FALSE, na = "")
write.csv(richness_summary, summary_path, row.names = FALSE, na = "")
write.csv(paired_changes, paired_changes_path, row.names = FALSE, na = "")

#### PNG FIGURE ####

mean_line_data <- richness_summary
mean_line_data$mean_plus_se <- mean_line_data$mean_richness + mean_line_data$se_richness
mean_line_data$mean_minus_se <- mean_line_data$mean_richness - mean_line_data$se_richness
mean_line_data$label <- sprintf(
  "mean %.1f\nn=%s",
  mean_line_data$mean_richness,
  format(mean_line_data$n_surveys, big.mark = ",")
)
mean_line_data$label_y <- max(plot_richness$species_richness, na.rm = TRUE) + 9

richness_plot <- ggplot(plot_richness, aes(x = period, y = species_richness)) +
  geom_boxplot(
    aes(fill = period),
    width = 0.58,
    outlier.alpha = 0.08,
    outlier.size = 0.45,
    linewidth = 0.35,
    colour = "grey25"
  ) +
  geom_line(
    data = mean_line_data,
    aes(x = period, y = mean_richness, group = 1),
    inherit.aes = FALSE,
    linewidth = 0.7,
    colour = "grey15"
  ) +
  geom_errorbar(
    data = mean_line_data,
    aes(x = period, ymin = mean_minus_se, ymax = mean_plus_se),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 0.55,
    colour = "grey15"
  ) +
  geom_point(
    data = mean_line_data,
    aes(x = period, y = mean_richness),
    inherit.aes = FALSE,
    size = 2.6,
    colour = "grey15"
  ) +
  geom_text(
    data = mean_line_data,
    aes(x = period, y = label_y, label = label),
    inherit.aes = FALSE,
    size = 3.2,
    lineheight = 0.9,
    colour = "grey12"
  ) +
  facet_wrap(~ analysis_set, ncol = 1) +
  scale_fill_manual(
    values = c("1970s" = "#4477AA", "1990s" = "#CCBB44", "2010s" = "#228833"),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, max(plot_richness$species_richness, na.rm = TRUE) + 18),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Observed species richness per survey row",
    title = "Species richness at model-included sites",
    subtitle = "Raw 201-species occurrence matrix; rows filtered by land threshold and Good/Average effort."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", hjust = 0),
    strip.background = element_rect(fill = "grey92", colour = NA),
    axis.text.x = element_text(colour = "grey20"),
    axis.text.y = element_text(colour = "grey20"),
    plot.margin = margin(12, 18, 12, 18)
  )

ggsave(
  filename = plot_path,
  plot = richness_plot,
  width = 7.4,
  height = 7.2,
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

print(richness_summary)

message("Common sites represented in all three atlas periods: ", length(common_sites))
message("Wrote survey-level richness to: ", survey_richness_path)
message("Wrote richness summary to: ", summary_path)
message("Wrote paired common-site changes to: ", paired_changes_path)
message("Wrote PNG plot to: ", plot_path)
