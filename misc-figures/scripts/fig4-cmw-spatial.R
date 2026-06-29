rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, Hmsc, cowplot, ggbeeswarm, scales)

source(file.path("support_scripts", "figure_data_helpers.R"))

base_dir <- "./HmscOutputs"
pattern <- "2026-03-13"
out_dir <- file.path("misc-figures", "outputs", "main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

matching_folders <- figure_model_folders(pattern = pattern, base_dir = base_dir)
mods <- load_hmsc_posteriors(matching_folders, base_dir = base_dir)
designs <- load_hmsc_study_designs(mods)
preds_y <- load_or_compute_site_predictions(mods, matching_folders, base_dir = base_dir)

period_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")
period_levels <- unname(period_lookup)
thermal_levels <- c("very cold", "cold", "medium", "warm", "very warm")

thermal_palette_options <- list(
  "Blue-teal-gold" = c("#244C7A", "#2F86A6", "#5BB98C", "#F0D96A", "#E9913C"),
  "Navy-mint-apricot" = c("#263C6A", "#3E8CBF", "#8FD2B3", "#F6E0A6", "#D9895B"),
  "Cividis thermal" = c("#31446B", "#496C89", "#829E82", "#D2C56F", "#F28E2B"),
  "Blue-cream-copper" = c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#A6611A"),
  "Blue-purple-rose" = c("#313695", "#74ADD1", "#FEE090", "#F46D43", "#9E0142")
)
selected_thermal_palette_name <- "Blue-purple-rose"
use_quantile_colour_stops <- TRUE
cwm_colour_reference_period <- "1970s"
thermal_colours <- thermal_palette_options[[selected_thermal_palette_name]]
thermal_palette <- setNames(
  colorRampPalette(thermal_colours)(length(thermal_levels)),
  thermal_levels
)
delta_colours <- c("blue4", "white", "red2")

base_site <- function(x) sub("_[123]$", "", x)

sti_sp <- mods[[1]]$Tr[, "species_thermal_index"]
species <- Reduce(intersect, c(list(names(sti_sp)), map(preds_y, colnames)))
sti_sp <- sti_sp[species]
preds_y <- map(preds_y, ~ .x[, species, drop = FALSE])

shared_sites <- Reduce(intersect, map(preds_y, ~ base_site(rownames(.x))))

preds_shared <- imap(preds_y, function(pred, atlas_id) {
  base <- base_site(rownames(pred))
  keep <- base %in% shared_sites
  out <- pred[keep, species, drop = FALSE]
  rownames(out) <- base[keep]
  out[order(rownames(out)), , drop = FALSE]
})

design_shared <- imap_dfr(designs, function(design, atlas_id) {
  design |>
    mutate(
      site = base_site(.data$survey),
      atlas_id = atlas_id
    ) |>
    filter(.data$site %in% shared_sites) |>
    select(site, atlas_id, X, Y) |>
    arrange(.data$site)
})

thermal_breaks <- quantile(sti_sp, seq(0, 1, 0.2), na.rm = TRUE)
sp_group <- cut(
  sti_sp,
  breaks = thermal_breaks,
  include.lowest = TRUE,
  labels = thermal_levels
)
names(sp_group) <- names(sti_sp)

rel_abund <- map(preds_shared, function(pred) {
  richness <- rowSums(pred)
  sweep(pred, 1, richness, "/")
})

group_probability <- imap_dfr(rel_abund, function(rel, atlas_id) {
  map_dfc(thermal_levels, function(group) {
    group_species <- names(sp_group)[sp_group == group]
    tibble(!!group := rowSums(rel[, group_species, drop = FALSE]))
  }) |>
    mutate(
      site = rownames(rel),
      atlas_id = atlas_id
    )
})

dominant_group <- group_probability |>
  pivot_longer(
    cols = all_of(thermal_levels),
    names_to = "group",
    values_to = "probability"
  ) |>
  group_by(.data$site, .data$atlas_id) |>
  slice_max(.data$probability, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    period = factor(unname(period_lookup[.data$atlas_id]), levels = period_levels),
    group = factor(.data$group, levels = thermal_levels)
  )

cwm_sti <- imap_dfr(rel_abund, function(rel, atlas_id) {
  tibble(
    site = rownames(rel),
    atlas_id = atlas_id,
    sti = as.numeric(rel %*% sti_sp)
  )
}) |>
  left_join(design_shared, by = c("site", "atlas_id")) |>
  mutate(period = factor(unname(period_lookup[.data$atlas_id]), levels = period_levels))

plot_df <- cwm_sti |>
  left_join(
    dominant_group |> select(site, atlas_id, group),
    by = c("site", "atlas_id")
  )

pct_df <- plot_df |>
  count(.data$period, .data$group, name = "n") |>
  complete(
    period = factor(period_levels, levels = period_levels),
    group = factor(thermal_levels, levels = thermal_levels),
    fill = list(n = 0)
  ) |>
  group_by(.data$period) |>
  mutate(
    pct = 100 * .data$n / sum(.data$n),
    pct_label = case_when(
      .data$n == 0 ~ "0%",
      .data$pct < 1 ~ "<1%",
      TRUE ~ paste0(round(.data$pct), "%")
    )
  ) |>
  ungroup()

delta_df <- cwm_sti |>
  select(site, period, sti, X, Y) |>
  pivot_wider(names_from = period, values_from = sti) |>
  mutate(
    `1990s minus 1970s` = .data$`1990s` - .data$`1970s`,
    `2010s minus 1970s` = .data$`2010s` - .data$`1970s`
  ) |>
  select(site, X, Y, `1990s minus 1970s`, `2010s minus 1970s`) |>
  pivot_longer(
    cols = c(`1990s minus 1970s`, `2010s minus 1970s`),
    names_to = "contrast",
    values_to = "delta_sti"
  ) |>
  mutate(contrast = factor(.data$contrast, levels = c("1990s minus 1970s", "2010s minus 1970s")))

period_summary <- cwm_sti |>
  group_by(.data$period) |>
  summarise(
    median_sti = median(.data$sti, na.rm = TRUE),
    mean_sti = mean(.data$sti, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    median_delta = .data$median_sti - .data$median_sti[.data$period == "1970s"],
    mean_delta = .data$mean_sti - .data$mean_sti[.data$period == "1970s"]
  )

mainland_bbox <- c(xmin = 400000, xmax = 750000, ymin = 6000000, ymax = 6450000)
bornholm_bbox <- c(xmin = 855000, xmax = 905000, ymin = 6100000, ymax = 6160000)
mainland_width <- mainland_bbox[["xmax"]] - mainland_bbox[["xmin"]]
mainland_height <- mainland_bbox[["ymax"]] - mainland_bbox[["ymin"]]
bornholm_width <- bornholm_bbox[["xmax"]] - bornholm_bbox[["xmin"]]
bornholm_height <- bornholm_bbox[["ymax"]] - bornholm_bbox[["ymin"]]
bornholm_inset_width <- bornholm_width / mainland_width
bornholm_inset_height <- bornholm_height / mainland_height

cwm_limits <- range(cwm_sti$sti, na.rm = TRUE)
cwm_breaks <- pretty(cwm_limits, n = 5)
delta_limit <- max(abs(delta_df$delta_sti), na.rm = TRUE)
delta_limits <- c(-delta_limit, delta_limit)
delta_breaks <- pretty(delta_limits, n = 5)

cwm_colour_values <- function(n) {
  if (!use_quantile_colour_stops) {
    return(seq(0, 1, length.out = n))
  }

  reference_sti <- cwm_sti |>
    filter(.data$period == cwm_colour_reference_period) |>
    pull(.data$sti)

  values <- quantile(
    reference_sti,
    probs = seq(0, 1, length.out = n),
    na.rm = TRUE,
    names = FALSE
  )
  values <- c(cwm_limits[[1]], values[-c(1, n)], cwm_limits[[2]])
  values <- rescale(values, from = cwm_limits)

  if (anyDuplicated(round(values, 8)) > 0) {
    seq(0, 1, length.out = n)
  } else {
    values
  }
}

map_theme <- function(legend_position = "none") {
  theme_minimal(base_size = 10) +
    theme(
      legend.position = legend_position,
      legend.title = element_text(face = "bold"),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(0, 0, 0, 0)
    )
}

plot_cwm_base <- function(df, title = NULL, bbox = mainland_bbox, show_legend = FALSE,
                          border = FALSE, cwm_colours = thermal_colours) {
  ggplot(df) +
    geom_point(aes(x = .data$X, y = .data$Y, colour = .data$sti), size = 1.25, alpha = 0.95) +
    scale_colour_gradientn(
      colours = cwm_colours,
      values = cwm_colour_values(length(cwm_colours)),
      limits = cwm_limits,
      breaks = cwm_breaks,
      labels = label_number(accuracy = 0.1),
      name = expression(CWM[STI]),
      guide = guide_colourbar(barwidth = 8, barheight = 0.7)
    ) +
    coord_fixed(
      xlim = c(bbox[["xmin"]], bbox[["xmax"]]),
      ylim = c(bbox[["ymin"]], bbox[["ymax"]]),
      expand = FALSE
    ) +
    labs(title = title) +
    map_theme(if (show_legend) "bottom" else "none") +
    theme(
      panel.border = if (border) {
        element_rect(colour = "grey35", fill = NA, linewidth = 0.45)
      } else {
        element_blank()
      }
    )
}

plot_delta_base <- function(df, title = NULL, bbox = mainland_bbox, show_legend = FALSE, border = FALSE) {
  ggplot(df) +
    geom_point(aes(x = .data$X, y = .data$Y, colour = .data$delta_sti), size = 1.25, alpha = 0.95) +
    scale_colour_gradient2(
      low = delta_colours[[1]],
      mid = delta_colours[[2]],
      high = delta_colours[[3]],
      midpoint = 0,
      limits = delta_limits,
      breaks = delta_breaks,
      labels = label_number(accuracy = 0.01),
      name = expression(Delta~CWM[STI]),
      guide = guide_colourbar(barwidth = 8, barheight = 0.7)
    ) +
    coord_fixed(
      xlim = c(bbox[["xmin"]], bbox[["xmax"]]),
      ylim = c(bbox[["ymin"]], bbox[["ymax"]]),
      expand = FALSE
    ) +
    labs(title = title) +
    map_theme(if (show_legend) "bottom" else "none") +
    theme(
      panel.border = if (border) {
        element_rect(colour = "grey35", fill = NA, linewidth = 0.45)
      } else {
        element_blank()
      }
    )
}

plot_map_with_inset <- function(df, title, base_fun) {
  p_main <- base_fun(df, title = title, bbox = mainland_bbox)
  p_inset <- base_fun(df, bbox = bornholm_bbox, border = TRUE) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.7)
    )

  ggdraw(p_main) +
    draw_plot(
      p_inset,
      x = 1 - bornholm_inset_width - 0.2,
      y = 1 - bornholm_inset_height - 0.2,
      width = bornholm_inset_width,
      height = bornholm_inset_height
    )
}

cwm_period_maps <- map(period_levels, function(period_name) {
  plot_map_with_inset(
    cwm_sti |> filter(.data$period == period_name),
    title = period_name,
    base_fun = plot_cwm_base
  )
})

delta_maps <- levels(delta_df$contrast) |>
  map(function(contrast_name) {
    plot_map_with_inset(
      delta_df |> filter(.data$contrast == contrast_name),
      title = contrast_name,
      base_fun = plot_delta_base
    )
  })

dominance_tiles <- function(period_name) {
  tile_df <- pct_df |>
    filter(.data$period == period_name) |>
    mutate(pct_text_colour = if_else(.data$group == "medium", "grey20", "white"))

  ggplot(tile_df, aes(x = .data$group, y = 1, fill = .data$group)) +
    geom_tile(width = 0.92, height = 0.82, colour = "white", linewidth = 0.35) +
    geom_text(
      aes(label = .data$pct_label, colour = .data$pct_text_colour),
      size = 2.45,
      fontface = "bold"
    ) +
    scale_fill_manual(values = thermal_palette, drop = FALSE, guide = "none") +
    scale_colour_identity() +
    scale_x_discrete(expand = expansion(add = 0.08)) +
    scale_y_continuous(limits = c(0.45, 1.55), expand = expansion(mult = c(0, 0))) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 4, 0, 4))
}

bee_plot <- function(period_name, show_x_axis = TRUE) {
  bee_df <- plot_df |> filter(.data$period == period_name)
  summary_row <- period_summary |> filter(.data$period == period_name)
  title <- if (period_name == "1970s") {
    period_name
  } else {
    paste0(period_name, " (median ", label_number(accuracy = 0.01, style_positive = "plus")(summary_row$median_delta), ")")
  }
  ref_median <- period_summary$median_sti[period_summary$period == "1970s"]

  ggplot(bee_df, aes(x = .data$sti)) +
    geom_vline(xintercept = ref_median, colour = "grey35", linewidth = 0.35, linetype = "dashed") +
    geom_boxplot(
      aes(y = 0.78),
      width = 0.08,
      outlier.shape = NA,
      fill = NA,
      linewidth = 0.3,
      colour = "grey20"
    ) +
    geom_beeswarm(
      aes(y = 1, colour = .data$group),
      size = 0.58,
      alpha = 0.62,
      cex = 0.45,
      orientation = "y"
    ) +
    scale_colour_manual(values = thermal_palette, drop = FALSE, guide = "none") +
    scale_x_continuous(
      limits = cwm_limits,
      breaks = cwm_breaks,
      labels = label_number(accuracy = 0.1),
      expand = expansion(mult = c(0.015, 0.015))
    ) +
    scale_y_continuous(limits = c(0.66, 1.24), expand = expansion(mult = c(0, 0))) +
    labs(title = title, x = if (show_x_axis) expression(CWM[STI]) else NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = if (show_x_axis) element_text(size = 7) else element_blank(),
      axis.ticks.x = if (show_x_axis) element_line(linewidth = 0.25) else element_blank(),
      axis.title.x = if (show_x_axis) element_text(size = 8, margin = margin(t = 2)) else element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
      plot.margin = if (show_x_axis) margin(1, 4, 2, 4) else margin(1, 4, 0, 4)
    )
}

bee_row <- function(period_name) {
  plot_grid(
    dominance_tiles(period_name),
    bee_plot(period_name, show_x_axis = period_name == tail(period_levels, 1)),
    ncol = 1,
    rel_heights = c(0.2, 1),
    align = "v"
  )
}

map_panel_a <- plot_grid(
  plotlist = cwm_period_maps,
  nrow = 1,
  align = "h",
  rel_widths = c(1, 1, 1)
)

change_panel_b <- plot_grid(
  plotlist = delta_maps,
  nrow = 1,
  align = "h",
  rel_widths = c(1, 1)
)

beeswarm_panel_c_body <- plot_grid(
  plotlist = map(period_levels, bee_row),
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1, 1)
)

beeswarm_panel_c <- plot_grid(
  NULL,
  beeswarm_panel_c_body,
  nrow = 1,
  rel_widths = c(0.045, 1)
)

cwm_legend <- get_legend(plot_cwm_base(cwm_sti, show_legend = TRUE))
delta_legend <- get_legend(plot_delta_base(delta_df, show_legend = TRUE))
group_legend <- get_legend(
  ggplot(pct_df, aes(x = .data$group, y = .data$pct, fill = .data$group)) +
    geom_col() +
    scale_fill_manual(
      values = thermal_palette,
      breaks = thermal_levels,
      name = "Dominant thermal group",
      guide = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.margin = margin(0, 0, 0, 0)
    )
)

legend_row <- plot_grid(
  cwm_legend,
  delta_legend,
  group_legend,
  nrow = 1,
  rel_widths = c(0.9, 0.9, 1.25)
)

top_row <- plot_grid(
  map_panel_a,
  labels = "A",
  label_size = 15,
  label_fontface = "bold"
)

bottom_row <- plot_grid(
  change_panel_b,
  beeswarm_panel_c,
  labels = c("B", "C"),
  label_size = 15,
  label_fontface = "bold",
  nrow = 1,
  rel_widths = c(1.18, 1)
)

fig4_cwm_sti <- plot_grid(
  top_row,
  bottom_row,
  legend_row,
  ncol = 1,
  rel_heights = c(1.08, 1, 0.13)
) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

png_path <- file.path(out_dir, paste0(pattern, "-fig4-cwm-maps-dominance.png"))
pdf_path <- file.path(out_dir, paste0(pattern, "-fig4-cwm-maps-dominance.pdf"))

ggsave(png_path, fig4_cwm_sti, width = 13.5, height = 10.2, units = "in", dpi = 300, bg = "white")
ggsave(pdf_path, fig4_cwm_sti, width = 13.5, height = 10.2, units = "in", bg = "white")

message("Saved: ", png_path)
message("Saved: ", pdf_path)

palette_swatch <- function(palette_name, colours) {
  swatch_df <- tibble(
    group = factor(thermal_levels, levels = thermal_levels),
    x = seq_along(thermal_levels),
    label = c("VC", "C", "M", "W", "VW")
  )

  ggplot(swatch_df, aes(x = .data$x, y = 1, fill = .data$group)) +
    geom_tile(width = 0.96, height = 0.5, colour = "white", linewidth = 0.35) +
    geom_text(aes(label = .data$label), colour = "white", size = 2.5, fontface = "bold") +
    scale_fill_manual(values = setNames(colours, thermal_levels), guide = "none") +
    scale_x_continuous(limits = c(0.45, length(thermal_levels) + 0.55), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0.72, 1.28), expand = expansion(mult = c(0, 0))) +
    labs(title = palette_name) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0),
      plot.margin = margin(0, 4, 0, 2)
    )
}

palette_comparison_row <- function(colours, palette_name) {
  palette_base_fun <- function(df, title = NULL, bbox = mainland_bbox, show_legend = FALSE, border = FALSE) {
    plot_cwm_base(
      df = df,
      title = title,
      bbox = bbox,
      show_legend = show_legend,
      border = border,
      cwm_colours = colours
    )
  }

  palette_maps <- map(period_levels, function(period_name) {
    plot_map_with_inset(
      cwm_sti |> filter(.data$period == period_name),
      title = period_name,
      base_fun = palette_base_fun
    )
  })

  plot_grid(
    plotlist = c(list(palette_swatch(palette_name, colours)), palette_maps),
    nrow = 1,
    rel_widths = c(0.58, 1, 1, 1),
    align = "h"
  )
}

palette_comparison_title <- ggdraw() +
  draw_label(
    expression("CWM"[STI] * " palette comparison: 1970s-reference quantile stops over the observed STI range"),
    x = 0,
    hjust = 0,
    fontface = "bold",
    size = 13
  )

palette_comparison <- plot_grid(
  plotlist = c(
    list(palette_comparison_title),
    imap(thermal_palette_options, palette_comparison_row)
  ),
  ncol = 1,
  rel_heights = c(0.12, rep(1, length(thermal_palette_options)))
) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

palette_png_path <- file.path(out_dir, paste0(pattern, "-fig4-cwm-palette-comparison.png"))
palette_pdf_path <- file.path(out_dir, paste0(pattern, "-fig4-cwm-palette-comparison.pdf"))

ggsave(palette_png_path, palette_comparison, width = 13.5, height = 10.8, units = "in", dpi = 300, bg = "white")
ggsave(palette_pdf_path, palette_comparison, width = 13.5, height = 10.8, units = "in", bg = "white")

message("Saved: ", palette_png_path)
message("Saved: ", palette_pdf_path)
