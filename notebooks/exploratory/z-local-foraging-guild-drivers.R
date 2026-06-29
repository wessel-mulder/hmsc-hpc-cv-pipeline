rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  GWmodel, sf, sp, tidyverse, Hmsc,
  spdep, RColorBrewer, cowplot, terra, patchwork
)

source(file.path("support_scripts", "figure_data_helpers.R"))

base_dir <- "./HmscOutputs"
pattern <- "2026-03-13"
bandwidth_m <- 100000
min_species_per_group <- 5
min_expected_total_richness <- 1

grid_shape_path <- path.expand("~/box/PhD/logistics/data/distributions/DK5km_ED50grid_approx_kvadrkod_DOF/DK5km_ED50grid_approx_kvadrkod_DOF.shp")
mainland_bbox <- st_bbox(c(xmin = 400000, xmax = 750000, ymin = 6000000, ymax = 6450000), crs = st_crs(25832))
bornholm_bbox <- st_bbox(c(xmin = 855000, xmax = 905000, ymin = 6100000, ymax = 6160000), crs = st_crs(25832))

year_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")

matching_folders <- figure_model_folders(pattern = pattern, base_dir = base_dir)
models_nums <- atlas_numbers(matching_folders)
mods <- load_hmsc_posteriors(matching_folders, base_dir = base_dir)
designs <- load_hmsc_study_designs(mods)
predsY <- load_or_compute_site_predictions(mods, matching_folders, base_dir = base_dir)
names(mods) <- names(designs) <- names(predsY) <- models_nums

load("Data/preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData")
traits <- Tr |>
  tibble::rownames_to_column("species") |>
  select(species, Migration_a3_DOF, foraging_guild_consensus)

analysis_specs <- list(
  list(
    trait_col = "foraging_guild_consensus",
    trait_label = "Foraging guild",
    group_col = "foraging_guild",
    response_label = "Foraging guild probability",
    output_stub = "foraging-guild",
    out_dir = file.path("notebooks", "exploratory", "outputs", "local-foraging-guild-drivers")
  ),
  list(
    trait_col = "Migration_a3_DOF",
    trait_label = "Migratory strategy",
    group_col = "migratory_strategy",
    response_label = "Migratory strategy probability",
    output_stub = "migratory-strategy",
    out_dir = file.path("notebooks", "exploratory", "outputs", "local-migratory-strategy-drivers")
  )
)

predictor_names <- function(df) {
  c(
    "tmean_breeding",
    "prec_breeding",
    "hh",
    grep("^perc_", colnames(df), value = TRUE)
  )
}

predictor_labels <- c(
  "tmean_breeding" = "Temperature",
  "prec_breeding" = "Precipitation",
  "hh" = "Land-use heterogeneity",
  "perc_urban" = "Urban (% coverage)",
  "perc_cropland" = "Cropland (% coverage)",
  "perc_pasture" = "Pasture (% coverage)",
  "perc_forest" = "Forest (% coverage)",
  "perc_grass_shrub" = "Grass/Shrubland (% coverage)"
)

driver_colours <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban (% coverage)" = "#4d4d4d",
  "Cropland (% coverage)" = "goldenrod1",
  "Pasture (% coverage)" = "darkorange",
  "Forest (% coverage)" = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2",
  "No significant driver" = "#d9d9d9"
)

probability_palette <- colorRampPalette(c(
  "#f7fbff", "#deebf7", "#9ecae1", "#4292c6", "#08519c"
))

shape_sf <- NULL
if (file.exists(grid_shape_path)) {
  shape_sf <- st_as_sf(vect(grid_shape_path))
}

mainland_width <- as.numeric(mainland_bbox["xmax"] - mainland_bbox["xmin"])
mainland_height <- as.numeric(mainland_bbox["ymax"] - mainland_bbox["ymin"])
bornholm_width <- as.numeric(bornholm_bbox["xmax"] - bornholm_bbox["xmin"])
bornholm_height <- as.numeric(bornholm_bbox["ymax"] - bornholm_bbox["ymin"])
inset_w <- bornholm_width / mainland_width
inset_h <- bornholm_height / mainland_height

group_probability <- function(pred_y, traits, trait_col, group_name, species_subset) {
  model_species <- colnames(pred_y)
  denominator_species <- intersect(model_species, traits$species)
  numerator_species <- intersect(model_species, species_subset)

  expected_total_richness <- rowSums(pred_y[, denominator_species, drop = FALSE], na.rm = TRUE)
  expected_group_richness <- rowSums(pred_y[, numerator_species, drop = FALSE], na.rm = TRUE)

  tibble(
    survey = rownames(pred_y),
    group_probability = ifelse(
      expected_total_richness > 0,
      expected_group_richness / expected_total_richness,
      NA_real_
    ),
    expected_group_richness = expected_group_richness,
    expected_total_richness = expected_total_richness,
    n_species = length(numerator_species),
    trait_group = group_name
  )
}

add_environment <- function(response_df, model) {
  response_df |>
    left_join(model$XData |> rownames_to_column(var = "survey"), by = "survey")
}

make_group_species <- function(traits, trait_col) {
  traits |>
    filter(!is.na(.data[[trait_col]]), .data[[trait_col]] != "") |>
    group_by(trait_group = .data[[trait_col]]) |>
    summarise(species = list(species), n_species = n(), .groups = "drop") |>
    filter(n_species >= min_species_per_group) |>
    arrange(desc(n_species), trait_group)
}

make_period_frames <- function(preds, designs, models, traits, group_species, spec) {
  imap_dfr(preds, function(pred_y, atlas) {
    pmap_dfr(group_species, function(trait_group, species, n_species) {
      group_probability(pred_y, traits, spec$trait_col, trait_group, species) |>
        left_join(designs[[atlas]], by = "survey") |>
        add_environment(models[[atlas]]) |>
        mutate(
          response = "group_probability",
          trait = spec$trait_label,
          atlas = .env$atlas,
          period = year_lookup[[as.character(.env$atlas)]],
          .before = 1
        )
    })
  })
}

prepare_gwr_data <- function(df, response, vars) {
  out <- df |>
    select(survey, site, X, Y, all_of(response), all_of(vars)) |>
    drop_na()

  out[vars] <- scale(out[vars])
  coordinates(out) <- c("X", "Y")
  out
}

fit_gwr_model <- function(df, response, vars, bandwidth_m = 100000) {
  if (nrow(df) < 40) {
    stop("Too few rows for local regression after filtering: ", nrow(df), call. = FALSE)
  }

  gwr_data <- prepare_gwr_data(df, response, vars)
  formula <- reformulate(vars, response = response)

  gwr_model <- gwr.basic(
    formula,
    data = gwr_data,
    bw = bandwidth_m,
    kernel = "bisquare",
    adaptive = FALSE
  )

  collin_diag <- tryCatch(
    gwr.collin.diagno(
      formula,
      data = gwr_data,
      bw = bandwidth_m,
      kernel = "bisquare",
      adaptive = FALSE
    ),
    error = function(e) {
      warning("Local collinearity diagnostic failed: ", conditionMessage(e))
      NULL
    }
  )

  list(
    mod = gwr_model,
    collin = collin_diag,
    bandwidth = bandwidth_m,
    adaptive = FALSE,
    vars = vars,
    response = response
  )
}

safe_fit_gwr_model <- possibly(fit_gwr_model, otherwise = NULL)

summarise_gwr_model <- function(fit, response, trait_label, group_name, period_or_comparison, n_rows, label, predictor_set = NA_character_) {
  if (is.null(fit)) {
    return(tibble(
      response = response,
      trait = trait_label,
      trait_group = group_name,
      period_or_comparison = period_or_comparison,
      model_set = label,
      predictor_set = predictor_set,
      n_rows = n_rows,
      bandwidth_value = bandwidth_m,
      aicc = NA_real_,
      rss = NA_real_,
      local_r2_mean = NA_real_,
      local_r2_median = NA_real_,
      local_r2_min = NA_real_,
      local_r2_max = NA_real_,
      local_cn_mean = NA_real_,
      local_cn_max = NA_real_,
      local_cn_over_30 = NA_real_,
      status = "failed"
    ))
  }

  sdf <- as.data.frame(fit$mod$SDF)
  diagnostics <- fit$mod$GW.diagnostic
  local_cn <- if (is.null(fit$collin)) NA_real_ else fit$collin$SDF$local_CN
  safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

  tibble(
    response = response,
    trait = trait_label,
    trait_group = group_name,
    period_or_comparison = period_or_comparison,
    model_set = label,
    predictor_set = predictor_set,
    n_rows = n_rows,
    bandwidth_value = fit$bandwidth,
    aicc = diagnostics$AICc,
    rss = diagnostics$RSS.gw,
    local_r2_mean = mean(sdf$Local_R2, na.rm = TRUE),
    local_r2_median = median(sdf$Local_R2, na.rm = TRUE),
    local_r2_min = min(sdf$Local_R2, na.rm = TRUE),
    local_r2_max = max(sdf$Local_R2, na.rm = TRUE),
    local_cn_mean = safe_mean(local_cn),
    local_cn_max = safe_max(local_cn),
    local_cn_over_30 = if (all(is.na(local_cn))) NA_real_ else mean(local_cn > 30, na.rm = TRUE),
    status = "ok"
  )
}

dominant_driver_table <- function(fit, metadata, labels) {
  if (is.null(fit)) return(tibble())

  sdf <- as.data.frame(fit$mod$SDF)
  tv_cols <- paste0(fit$vars, "_TV")
  tv_cols <- tv_cols[tv_cols %in% names(sdf)]
  tv_matrix <- as.matrix(sdf[, tv_cols, drop = FALSE])
  dominant_idx <- max.col(abs(tv_matrix), ties.method = "first")
  dominant_col <- tv_cols[dominant_idx]
  raw_tv_value <- tv_matrix[cbind(seq_len(nrow(tv_matrix)), dominant_idx)]
  raw_var <- gsub("_TV$", "", dominant_col)

  tibble(
    survey = rownames(sdf),
    X = coordinates(fit$mod$SDF)[, 1],
    Y = coordinates(fit$mod$SDF)[, 2],
    max_tv = abs(raw_tv_value),
    raw_tv_value = raw_tv_value,
    direction = ifelse(raw_tv_value >= 0, "Positive", "Negative"),
    dominant_predictor = raw_var,
    driver_label = unname(labels[raw_var]),
    dominant_driver = ifelse(max_tv < 1.96, "No significant driver", driver_label),
    local_r2 = sdf$Local_R2
  ) |>
    bind_cols(metadata[rep(1, nrow(sdf)), , drop = FALSE])
}

signed_driver_bars <- function(dominant_drivers, spec, plot_dir) {
  signed_driver_plot <- dominant_drivers |>
    filter(dominant_driver != "No significant driver") |>
    count(trait_group, period_or_comparison, dominant_driver, direction, name = "n_cells") |>
    group_by(trait_group, period_or_comparison) |>
    mutate(
      prop_cells = n_cells / sum(n_cells),
      signed_prop = ifelse(direction == "Positive", prop_cells, -prop_cells)
    ) |>
    ungroup()

  if (nrow(signed_driver_plot) == 0) {
    return(NULL)
  }

  n_groups <- n_distinct(signed_driver_plot$trait_group)
  ncol_facets <- ifelse(n_groups <= 6, min(3, n_groups), 4)
  plot_height <- max(5.5, ceiling(n_groups / ncol_facets) * 2.7 + 1.4)

  p <- ggplot(signed_driver_plot, aes(x = period_or_comparison, y = signed_prop, fill = dominant_driver)) +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.25) +
    geom_col(width = 0.72) +
    facet_wrap(~ trait_group, ncol = ncol_facets) +
    scale_fill_manual(values = driver_colours, name = "Dominant driver") +
    scale_y_continuous(
      labels = function(x) paste0(abs(round(x * 100)), "%"),
      breaks = seq(-1, 1, by = 0.25)
    ) +
    labs(
      x = NULL,
      y = "Share of significant cells",
      title = paste("Dominant local drivers of", str_to_lower(spec$response_label)),
      subtitle = "Positive effects are above zero; negative effects are below zero"
    ) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")

  ggsave(
    file.path(plot_dir, paste0(pattern, "-local-", spec$output_stub, "-period-dominant-driver-signed-proportions.png")),
    p,
    width = 13,
    height = plot_height,
    dpi = 300
  )

  p
}

map_theme <- theme_minimal(base_size = 9) +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

make_probability_map <- function(df, period_name, limits, legend_title, show_legend = FALSE) {
  if (!is.null(shape_sf)) {
    plot_data <- shape_sf |>
      left_join(df, by = c("kvadratkod" = "site"))

    map_layers <- function(data) {
      list(
        geom_sf(data = data, aes(fill = group_probability), color = "grey30", linewidth = 0.08),
        scale_fill_gradientn(
          colors = probability_palette(10),
          limits = limits,
          na.value = "transparent",
          name = legend_title
        )
      )
    }

    p_main <- ggplot() +
      map_layers(plot_data) +
      coord_sf(
        xlim = c(mainland_bbox["xmin"], mainland_bbox["xmax"]),
        ylim = c(mainland_bbox["ymin"], mainland_bbox["ymax"]),
        expand = FALSE
      ) +
      labs(title = period_name) +
      map_theme +
      theme(legend.position = if (show_legend) "right" else "none")

    p_inset <- ggplot() +
      map_layers(plot_data) +
      coord_sf(
        xlim = c(bornholm_bbox["xmin"], bornholm_bbox["xmax"]),
        ylim = c(bornholm_bbox["ymin"], bornholm_bbox["ymax"]),
        expand = FALSE
      ) +
      theme_void() +
      theme(
        legend.position = "none",
        panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.45)
      )

    ggdraw(p_main) +
      draw_plot(
        p_inset,
        x = 1 - inset_w - 0.2,
        y = 1 - inset_h - 0.2,
        width = inset_w,
        height = inset_h
      )
  } else {
    ggplot(df, aes(X, Y, colour = group_probability)) +
      geom_point(size = 0.7) +
      coord_equal() +
      scale_colour_gradientn(colors = probability_palette(10), limits = limits, name = legend_title) +
      labs(title = period_name) +
      map_theme +
      theme(legend.position = if (show_legend) "right" else "none")
  }
}

make_group_probability_panel <- function(period_frames, group_name, spec) {
  group_df <- period_frames |>
    filter(trait_group == group_name)
  limits <- range(group_df$group_probability, na.rm = TRUE)

  maps <- map(year_lookup, function(period_name) {
    df <- group_df |> filter(period == period_name)
    make_probability_map(df, period_name, limits = limits, legend_title = spec$response_label)
  })

  legend_plot <- make_probability_map(
    group_df |> filter(period == year_lookup[[1]]),
    year_lookup[[1]],
    limits = limits,
    legend_title = spec$response_label,
    show_legend = TRUE
  )
  legend <- get_legend(legend_plot)

  title <- ggdraw() +
    draw_label(group_name, fontface = "bold", size = 12, hjust = 0.5)

  plot_grid(
    title,
    plot_grid(plotlist = maps, nrow = 1, align = "none"),
    legend,
    ncol = 1,
    rel_heights = c(0.08, 1, 0.14)
  )
}

write_probability_maps <- function(period_frames, group_species, spec, plot_dir) {
  group_names <- group_species$trait_group

  walk(group_names, function(group_name) {
    p <- make_group_probability_panel(period_frames, group_name, spec)
    safe_name <- str_to_lower(str_replace_all(group_name, "[^A-Za-z0-9]+", "-"))
    ggsave(
      file.path(plot_dir, paste0(pattern, "-", spec$output_stub, "-probability-maps-", safe_name, ".png")),
      p,
      width = 8.5,
      height = 3.8,
      dpi = 300
    )
  })

  if (length(group_names) <= 8) {
    probability_panels <- map(group_names, ~ make_group_probability_panel(period_frames, .x, spec))
    probability_map_panel <- plot_grid(plotlist = probability_panels, ncol = 1)
    ggsave(
      file.path(plot_dir, paste0(pattern, "-", spec$output_stub, "-probability-maps.png")),
      probability_map_panel,
      width = 8.5,
      height = max(3.8, 3.45 * length(group_names)),
      dpi = 300,
      limitsize = FALSE
    )
  }
}

run_trait_probability_analysis <- function(spec) {
  model_dir <- file.path(spec$out_dir, "models")
  plot_dir <- file.path(spec$out_dir, "plots")
  dir.create(spec$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  group_species <- make_group_species(traits, spec$trait_col)
  period_frames <- make_period_frames(predsY, designs, mods, traits, group_species, spec)

  period_splits <- period_frames |>
    filter(expected_total_richness >= min_expected_total_richness) |>
    group_by(trait_group, atlas, period) |>
    group_split()
  names(period_splits) <- map_chr(period_splits, ~ paste(unique(.x$trait_group), unique(.x$atlas), sep = "__"))

  period_models <- period_splits |>
    map(function(df) {
      vars <- predictor_names(df)
      fit <- safe_fit_gwr_model(df, response = "group_probability", vars = vars, bandwidth_m = bandwidth_m)
      list(data = df, fit = fit, predictor_set = "full_period")
    })

  saveRDS(
    list(
      period_models = period_models,
      group_species = group_species,
      min_species_per_group = min_species_per_group,
      min_expected_total_richness = min_expected_total_richness,
      bandwidth_m = bandwidth_m
    ),
    file.path(model_dir, paste0(pattern, "-local-", spec$output_stub, "-period-gwr-models-100km.rds"))
  )

  gwr_summary <- imap_dfr(period_models, function(x, key) {
    parts <- strsplit(key, "__", fixed = TRUE)[[1]]
    summarise_gwr_model(
      x$fit,
      response = "group_probability",
      trait_label = spec$trait_label,
      group_name = parts[1],
      period_or_comparison = unique(x$data$period),
      n_rows = nrow(x$data),
      label = "period",
      predictor_set = x$predictor_set
    )
  })

  dominant_drivers <- imap_dfr(period_models, function(x, key) {
    parts <- strsplit(key, "__", fixed = TRUE)[[1]]
    metadata <- tibble(
      response = "group_probability",
      trait = spec$trait_label,
      trait_group = parts[1],
      atlas = parts[2],
      period_or_comparison = unique(x$data$period),
      model_set = "period",
      predictor_set = x$predictor_set
    )
    dominant_driver_table(x$fit, metadata, predictor_labels)
  })

  write_csv(group_species |> select(trait_group, n_species), file.path(spec$out_dir, paste0(pattern, "-", spec$output_stub, "-groups-modelled.csv")))
  write_csv(period_frames, file.path(spec$out_dir, paste0(pattern, "-", spec$output_stub, "-period-frames.csv")))
  write_csv(gwr_summary, file.path(spec$out_dir, paste0(pattern, "-local-", spec$output_stub, "-period-gwr-summary-100km.csv")))
  write_csv(dominant_drivers, file.path(spec$out_dir, paste0(pattern, "-local-", spec$output_stub, "-period-dominant-drivers-100km.csv")))

  r2_plot <- gwr_summary |>
    ggplot(aes(x = period_or_comparison, y = local_r2_mean)) +
    geom_col(width = 0.65, na.rm = TRUE, fill = "#4575b4") +
    facet_wrap(~ trait_group, scales = "free_x", ncol = 4) +
    labs(x = NULL, y = "Mean local R2", title = paste("Local", str_to_lower(spec$response_label), "GWR fit at 100 km")) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(
    file.path(plot_dir, paste0(pattern, "-local-", spec$output_stub, "-period-gwr-local-r2-summary.png")),
    r2_plot,
    width = 13,
    height = max(4, ceiling(nrow(group_species) / 4) * 2.4 + 1.2),
    dpi = 300
  )

  signed_driver_bars(dominant_drivers, spec, plot_dir)
  write_probability_maps(period_frames, group_species, spec, plot_dir)

  message("Finished period-only local ", str_to_lower(spec$trait_label), " analysis.")
  message("Outputs written to: ", spec$out_dir)
}

walk(analysis_specs, run_trait_probability_analysis)
