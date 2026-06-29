rm(list = ls())

#### GETTING STARTED ####
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, Hmsc, RColorBrewer, ggplot2,
  rnaturalearth, rnaturalearthdata,
  gridExtra, patchwork, sf, cowplot,
  terra
)
source(file.path("support_scripts", "figure_data_helpers.R"))

#### LOAD MODELS ####
dir <- "./HmscOutputs"
pattern <- "2026-03-13"

matching_folders <- figure_model_folders(pattern = pattern, base_dir = dir)
models_nums <- atlas_numbers(matching_folders)
mods <- load_hmsc_posteriors(matching_folders, base_dir = dir)
designs <- load_hmsc_study_designs(mods)

#### GET PREDICTIONS & EMPIRICAL DATA ####
predsY <- load_or_compute_site_predictions(mods, matching_folders, base_dir = dir)
emps <- load_empirical_responses(mods)

predicted_richness <- predicted_richness_frames(predsY, designs)

empirical_richness <- map2(emps, designs, function(emp_y, design) {
  emp_y |>
    as.data.frame() |>
    mutate(richness = rowSums(emp_y)) |>
    rownames_to_column(var = "survey") |>
    select(survey, richness) |>
    left_join(design, by = "survey")
})

#### PROBABILISTIC SORENSEN INDEX ####
probsorensen <- function(preds.vector, emp.vector) {
  sp.present <- names(emp.vector)[emp.vector == 1]

  if (length(sp.present) == 0) {
    return(NA_real_)
  }

  sp.present.preds <- preds.vector[names(preds.vector) %in% sp.present]
  sp.present.min.preds <- min(sp.present.preds, na.rm = TRUE)
  preds.greater.than.min.preds <- preds.vector[preds.vector >= sp.present.min.preds]

  numerator <- 2 * sum(sp.present.preds, na.rm = TRUE)
  denominator <- numerator + sum(preds.greater.than.min.preds, na.rm = TRUE)

  if (denominator == 0) {
    return(NA_real_)
  }

  numerator / denominator
}

sorensens <- map2(predsY, emps, function(pred_y, emp_y) {
  vals <- map_dbl(rownames(pred_y), function(row_name) {
    probsorensen(pred_y[row_name, ], emp_y[row_name, ])
  })

  tibble(survey = rownames(pred_y), sorensen = vals)
})

sorensen_frames <- map2(sorensens, designs, function(sorensen, design) {
  left_join(sorensen, design, by = "survey")
})

#### RICHNESS DIFFERENCE ####
richness_difference <- map2(predicted_richness, empirical_richness, function(pred_df, emp_df) {
  pred_df |>
    select(survey, site, X, Y, predicted_richness = richness) |>
    inner_join(
      emp_df |> select(survey, empirical_richness = richness),
      by = "survey"
    ) |>
    mutate(richness_difference = predicted_richness - empirical_richness)
})

#### DIAGNOSTICS ####
diagnostics <- imap_dfr(predicted_richness, function(pred_df, atlas) {
  emp_df <- empirical_richness[[atlas]]
  sor_df <- sorensen_frames[[atlas]]

  richness_check <- pred_df |>
    select(survey, predicted_richness = richness) |>
    inner_join(
      emp_df |> select(survey, empirical_richness = richness),
      by = "survey"
    ) |>
    mutate(richness_difference = predicted_richness - empirical_richness)

  tibble(
    atlas = atlas,
    n_sites = nrow(richness_check),
    predicted_observed_identical = isTRUE(all.equal(
      richness_check$predicted_richness,
      richness_check$empirical_richness,
      tolerance = 0
    )),
    mean_abs_richness_difference = mean(abs(richness_check$richness_difference), na.rm = TRUE),
    max_abs_richness_difference = max(abs(richness_check$richness_difference), na.rm = TRUE),
    richness_correlation = cor(
      richness_check$predicted_richness,
      richness_check$empirical_richness,
      use = "complete.obs"
    ),
    sorensen_min = min(sor_df$sorensen, na.rm = TRUE),
    sorensen_mean = mean(sor_df$sorensen, na.rm = TRUE),
    sorensen_max = max(sor_df$sorensen, na.rm = TRUE),
    sorensen_all_one = all(dplyr::near(sor_df$sorensen, 1))
  )
})

print(diagnostics)

if (any(diagnostics$predicted_observed_identical)) {
  stop("Predicted and observed richness are identical for at least one atlas.", call. = FALSE)
}

if (any(diagnostics$sorensen_all_one)) {
  stop("Probabilistic Sørensen is 1 for every site in at least one atlas.", call. = FALSE)
}

#### SHARED AESTHETICS ####
names(predicted_richness) <- names(empirical_richness) <- names(richness_difference) <-
  names(sorensen_frames) <- c("1", "2", "3")

year_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")

richness_limits <- bind_rows(c(
  map(predicted_richness, select, richness),
  map(empirical_richness, select, richness)
)) |>
  summarise(
    min = min(richness, na.rm = TRUE),
    max = max(richness, na.rm = TRUE)
  ) |>
  unlist(use.names = FALSE)
richness_breaks <- seq(richness_limits[1], richness_limits[2], length.out = 5)

richness_pal <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))
sorensen_pal <- colorRampPalette(c("#d73027", "#fee08b", "#1a9850"))
diff_extent <- max(abs(map_dbl(richness_difference, ~ min(.x$richness_difference, na.rm = TRUE))),
                   abs(map_dbl(richness_difference, ~ max(.x$richness_difference, na.rm = TRUE))))
diff_limits <- c(-diff_extent, diff_extent)
diff_breaks <- seq(diff_limits[1], diff_limits[2], length.out = 5)
diff_pal <- rev(brewer.pal(11, "RdYlBu"))
sorensen_limits <- range(map_dbl(sorensen_frames, ~ min(.x$sorensen, na.rm = TRUE)),
                         map_dbl(sorensen_frames, ~ max(.x$sorensen, na.rm = TRUE)))
sorensen_breaks <- seq(sorensen_limits[1], sorensen_limits[2], length.out = 5)

#### READ ATLAS POLYGON ####
shape <- vect("~/box/PhD/logistics/data/distributions/DK5km_ED50grid_approx_kvadrkod_DOF/DK5km_ED50grid_approx_kvadrkod_DOF.shp")
shape_sf <- st_as_sf(shape)

mainland_bbox <- st_bbox(c(
  xmin = 400000,
  xmax = 750000,
  ymin = 6000000,
  ymax = 6450000
), crs = st_crs(25832))

bornholm_bbox <- st_bbox(c(
  xmin = 855000,
  xmax = 905000,
  ymin = 6100000,
  ymax = 6160000
), crs = st_crs(25832))

mainland_width <- as.numeric(mainland_bbox["xmax"] - mainland_bbox["xmin"])
mainland_height <- as.numeric(mainland_bbox["ymax"] - mainland_bbox["ymin"])
bornholm_width <- as.numeric(bornholm_bbox["xmax"] - bornholm_bbox["xmin"])
bornholm_height <- as.numeric(bornholm_bbox["ymax"] - bornholm_bbox["ymin"])

inset_w <- bornholm_width / mainland_width
inset_h <- bornholm_height / mainland_height

#### MAP HELPERS ####
make_map <- function(df, name, value_col, fill_scale, legend_title,
                     show_title = TRUE, inset_offset = 0.02) {
  year_label <- year_lookup[[name]]

  plot_data <- shape_sf |>
    left_join(df, by = c("kvadratkod" = "site"))

  make_geom <- function(data) {
    list(
      geom_sf(data = data, aes(fill = .data[[value_col]]), color = "grey30", size = 0.1),
      fill_scale
    )
  }

  p_main <- ggplot() +
    make_geom(plot_data) +
    coord_sf(
      xlim = c(mainland_bbox["xmin"], mainland_bbox["xmax"]),
      ylim = c(mainland_bbox["ymin"], mainland_bbox["ymax"]),
      expand = FALSE
    ) +
    labs(x = NULL, y = NULL, title = if (show_title) year_label else NULL, fill = legend_title) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )

  p_inset <- ggplot() +
    make_geom(plot_data) +
    coord_sf(
      xlim = c(bornholm_bbox["xmin"], bornholm_bbox["xmax"]),
      ylim = c(bornholm_bbox["ymin"], bornholm_bbox["ymax"]),
      expand = FALSE
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.5)
    )

  p_final <- ggdraw(p_main) +
    draw_plot(
      p_inset,
      x = 1 - inset_w - inset_offset,
      y = 1 - inset_h - inset_offset,
      width = inset_w,
      height = inset_h
    )

  list(
    grob = as_grob(p_final),
    ggplot = p_main
  )
}

row_label <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, angle = 90, fontface = "bold", size = 4) +
    theme_void()
}

richness_scale <- scale_fill_gradientn(
  colors = richness_pal(10),
  limits = richness_limits,
  breaks = richness_breaks,
  labels = round(richness_breaks, 0),
  name = "Richness",
  na.value = "transparent"
)

diff_scale <- scale_fill_gradient2(
  low = diff_pal[1],
  mid = "white",
  high = diff_pal[11],
  midpoint = 0,
  limits = diff_limits,
  breaks = diff_breaks,
  labels = round(diff_breaks, 0),
  name = "Pred. - obs. richness",
  na.value = "transparent"
)

sorensen_scale <- scale_fill_gradientn(
  colors = sorensen_pal(10),
  limits = sorensen_limits,
  breaks = sorensen_breaks,
  labels = round(sorensen_breaks, 2),
  name = "Prob. Sørensen",
  na.value = "transparent"
)

predicted_maps <- imap(predicted_richness, ~ make_map(
  .x, .y, "richness", richness_scale, "Richness", show_title = TRUE
))

empirical_maps <- imap(empirical_richness, ~ make_map(
  .x, .y, "richness", richness_scale, "Richness", show_title = FALSE
))

difference_maps <- imap(richness_difference, ~ make_map(
  .x, .y, "richness_difference", diff_scale, "Predicted - observed richness", show_title = FALSE
))

sorensen_maps <- imap(sorensen_frames, ~ make_map(
  .x, .y, "sorensen", sorensen_scale, "Probabilistic Sørensen", show_title = FALSE
))

richness_legend <- as_grob(get_legend(
  predicted_maps[[1]]$ggplot + theme(legend.position = "right")
))
diff_legend <- as_grob(get_legend(
  difference_maps[[1]]$ggplot + theme(legend.position = "right")
))
sorensen_legend <- as_grob(get_legend(
  sorensen_maps[[1]]$ggplot + theme(legend.position = "right")
))

#### ASSEMBLE ####
layout <- "
RABCL
SEFGL
THIJM
UKNOP
"

final_layout <- wrap_plots(
  R = row_label("Predicted richness"),
  A = predicted_maps[[1]]$grob,
  B = predicted_maps[[2]]$grob,
  C = predicted_maps[[3]]$grob,
  L = richness_legend,
  S = row_label("Observed richness"),
  E = empirical_maps[[1]]$grob,
  F = empirical_maps[[2]]$grob,
  G = empirical_maps[[3]]$grob,
  T = row_label("Richness difference"),
  H = difference_maps[[1]]$grob,
  I = difference_maps[[2]]$grob,
  J = difference_maps[[3]]$grob,
  M = diff_legend,
  U = row_label("Community similarity"),
  K = sorensen_maps[[1]]$grob,
  N = sorensen_maps[[2]]$grob,
  O = sorensen_maps[[3]]$grob,
  P = sorensen_legend,
  design = layout
) +
  plot_layout(widths = c(0.12, 1, 1, 1, 0.6), heights = c(1, 1, 1, 1))

ggsave(
  sprintf("misc-figures/outputs/supplementary/%s-sfig4-richness-sorenson-similarity.png", pattern),
  final_layout,
  width = 11,
  height = 11
)
