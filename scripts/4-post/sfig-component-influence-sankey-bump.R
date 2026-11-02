rm(list = ls())

required_packages <- c("tidyverse", "ggplot2", "ggsankey")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install required package(s) before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggsankey)
})

pattern <- "2026-07-27"
scaled <- TRUE

model_dir <- file.path("intermediate-outputs", "hmsc-output")
figure_out_dir <- file.path("manuscript-outputs", "figures", "supplement")
table_out_dir <- file.path("manuscript-outputs", "tables", "supplement")

component_rows <- list(
  "Spatial component" = "Random: site",
  "Land-use component" = c(
    "hh",
    "perc_urban",
    "perc_cropland",
    "perc_pasture",
    "perc_forest",
    "perc_grass_shrub"
  ),
  "Climate component" = c("tmean_breeding", "prec_breeding")
)

component_colours <- c(
  "Spatial component" = "grey68",
  "Land-use component" = "springgreen4",
  "Climate component" = "dodgerblue3"
)

atlas_labels <- c(
  "1" = "1970s",
  "2" = "1990s",
  "3" = "2010s"
)

percent_decimal_labels <- function(x) {
  vals <- format(round(x * 100, 1), trim = TRUE, scientific = FALSE)
  paste0(sub("\\.0$", "", vals), "%")
}

read_vp_csv <- function(path) {
  raw <- read.csv(path, check.names = FALSE)
  rownames(raw) <- raw[[1]]
  raw[-1]
}

model_folders <- list.dirs(model_dir, recursive = FALSE, full.names = FALSE)
model_folders <- model_folders[grepl(pattern, basename(model_folders))]
model_folders <- model_folders[
  order(as.numeric(str_extract(model_folders, "(?<=Atlas)[0-9]+")))
]

if (length(model_folders) == 0) {
  stop(
    "No model folders matching pattern '",
    pattern,
    "' were found under ",
    model_dir,
    call. = FALSE
  )
}

component_by_species <- map_dfr(model_folders, function(model) {
  atlas_num <- str_extract(model, "(?<=Atlas)[0-9]+")
  vp_path <- file.path(
    model_dir,
    model,
    "Results",
    sprintf("%sparameter_estimates_VP_.csv", model)
  )

  vp <- read_vp_csv(vp_path)

  required_rows <- c(unlist(component_rows, use.names = FALSE), "TjurR2")
  missing_rows <- setdiff(required_rows, rownames(vp))
  if (length(missing_rows) > 0) {
    stop(
      "Missing expected row(s) in ",
      vp_path,
      ": ",
      paste(missing_rows, collapse = ", "),
      call. = FALSE
    )
  }

  vp_values <- vp[setdiff(rownames(vp), "TjurR2"), , drop = FALSE]

  if (isTRUE(scaled)) {
    tjur_r2 <- as.numeric(vp["TjurR2", ])
    names(tjur_r2) <- colnames(vp)
    vp_values <- sweep(vp_values, 2, tjur_r2, "*")
  }

  imap_dfr(component_rows, function(rows, component) {
    species_values <- colSums(vp_values[rows, , drop = FALSE], na.rm = TRUE)

    tibble(
      atlas = atlas_labels[[atlas_num]],
      component = component,
      species = names(species_values),
      value = as.numeric(species_values)
    )
  })
})

component_summary <- component_by_species %>%
  group_by(atlas, component) %>%
  summarise(
    n_species = n(),
    mean_vp = mean(value, na.rm = TRUE),
    median_vp = median(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    atlas = factor(atlas, levels = unname(atlas_labels)),
    component = factor(component, levels = names(component_rows))
  ) %>%
  arrange(atlas, component)

component_totals <- component_summary %>%
  group_by(atlas) %>%
  summarise(total_mean_vp = sum(mean_vp), .groups = "drop")

y_max <- ceiling(max(component_totals$total_mean_vp) * 100 / 5) * 5 / 100
y_breaks <- seq(0, y_max, by = 0.025)

y_label <- if (isTRUE(scaled)) {
  "Absolute variance explained (%)"
} else {
  "Relative variance explained (%)"
}

component_sankey <- ggplot(
  component_summary,
  aes(
    x = atlas,
    node = component,
    fill = component,
    value = mean_vp
  )
) +
  geom_sankey_bump(space = 0, type = "alluvial", alpha = 0.8) +
  scale_y_continuous(
    limits = c(0, y_max),
    breaks = y_breaks,
    minor_breaks = seq(0, y_max, by = 0.0125),
    labels = percent_decimal_labels
  ) +
  scale_fill_manual(values = component_colours) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    x = NULL,
    y = y_label,
    fill = "Component"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.25),
    plot.background = element_rect(fill = "grey99", color = NA),
    axis.text.x = element_text(face = "bold", size = 12)
  )

dir.create(figure_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_out_dir, recursive = TRUE, showWarnings = FALSE)

figure_path <- file.path(
  figure_out_dir,
  sprintf("%s-sfig-component-influence-sankey-bump.png", pattern)
)
data_path <- file.path(
  table_out_dir,
  sprintf("%s-sfig-component-influence-sankey-bump-data.csv", pattern)
)

write_csv(
  component_summary %>%
    mutate(
      mean_vp_percent = mean_vp * 100,
      median_vp_percent = median_vp * 100
    ),
  data_path
)

ggsave(
  figure_path,
  component_sankey,
  width = 10,
  height = 6
)

message("Saved figure: ", figure_path)
message("Saved plotted data: ", data_path)
