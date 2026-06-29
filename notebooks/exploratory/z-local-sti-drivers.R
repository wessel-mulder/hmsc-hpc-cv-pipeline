rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  GWmodel, sf, sp, tidyverse, Hmsc,
  spdep, RColorBrewer, cowplot
)

source(file.path("support_scripts", "figure_data_helpers.R"))

base_dir <- "./HmscOutputs"
pattern <- "2026-03-13"
bandwidth_m <- 100000
min_guild_expected_richness <- 0.25

out_dir <- file.path("notebooks", "exploratory", "outputs", "local-sti-drivers")
model_dir <- file.path(out_dir, "models")
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

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
  select(species, foraging_guild_consensus, species_thermal_index)

taxonomy <- read.csv("Data/data/1_preprocessing/Taxonomy/taxonomy.csv", check.names = FALSE) |>
  as_tibble(.name_repair = "minimal")
names(taxonomy)[1:6] <- c("idx", "species", "latin", "family", "order", "genus")

guild_sets <- list(
  all_species = colnames(predsY[[1]]),
  woodpeckers = taxonomy |> filter(family == "Picidae") |> pull(species),
  dabbling_ducks = traits |> filter(foraging_guild_consensus == "Dabbling ducks") |> pull(species),
  aerial_insectivores = traits |> filter(foraging_guild_consensus == "Aerial insectivores") |> pull(species)
)

predictor_names <- function(df) {
  c(
    "tmean_breeding",
    "prec_breeding",
    "hh",
    grep("^perc_", colnames(df), value = TRUE)
  )
}

period_predictor_labels <- c(
  "tmean_breeding" = "Temperature",
  "prec_breeding" = "Precipitation",
  "hh" = "Land-use heterogeneity",
  "perc_urban" = "Urban (% coverage)",
  "perc_cropland" = "Cropland (% coverage)",
  "perc_pasture" = "Pasture (% coverage)",
  "perc_forest" = "Forest (% coverage)",
  "perc_grass_shrub" = "Grass/Shrubland (% coverage)"
)

delta_predictor_labels <- paste0("Delta ", period_predictor_labels)
names(delta_predictor_labels) <- paste0("delta_", names(period_predictor_labels))

driver_colours <- c(
  "Temperature" = "firebrick3",
  "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3",
  "Urban (% coverage)" = "#4d4d4d",
  "Cropland (% coverage)" = "goldenrod1",
  "Pasture (% coverage)" = "darkorange",
  "Forest (% coverage)" = "springgreen4",
  "Grass/Shrubland (% coverage)" = "springgreen2",
  "Delta Temperature" = "firebrick3",
  "Delta Precipitation" = "dodgerblue3",
  "Delta Land-use heterogeneity" = "orchid3",
  "Delta Urban (% coverage)" = "#4d4d4d",
  "Delta Cropland (% coverage)" = "goldenrod1",
  "Delta Pasture (% coverage)" = "darkorange",
  "Delta Forest (% coverage)" = "springgreen4",
  "Delta Grass/Shrubland (% coverage)" = "springgreen2",
  "No significant driver" = "#d9d9d9"
)

cwm_sti_for_species <- function(pred_y, sti, species_subset = NULL) {
  species <- colnames(pred_y)
  if (!is.null(species_subset)) {
    species <- intersect(species, species_subset)
  }
  species <- species[species %in% names(sti)]

  pred_sub <- pred_y[, species, drop = FALSE]
  richness <- rowSums(pred_sub, na.rm = TRUE)
  numerator <- as.numeric(pred_sub %*% sti[species])

  tibble(
    survey = rownames(pred_y),
    cwm_sti = ifelse(richness > 0, numerator / richness, NA_real_),
    expected_richness = richness,
    n_species = length(species)
  )
}

add_environment <- function(response_df, model) {
  response_df |>
    left_join(model$XData |> rownames_to_column(var = "survey"), by = "survey")
}

make_period_frames <- function(preds, designs, models, guild_sets) {
  sti <- models[[1]]$Tr[, "species_thermal_index"]

  imap_dfr(preds, function(pred_y, atlas) {
    imap_dfr(guild_sets, function(species_subset, guild) {
      cwm_sti_for_species(pred_y, sti, species_subset) |>
        left_join(designs[[atlas]], by = "survey") |>
        add_environment(models[[atlas]]) |>
        mutate(
          response = "cwm_sti",
          guild = guild,
          atlas = .env$atlas,
          period = year_lookup[[as.character(.env$atlas)]],
          .before = 1
        )
    })
  })
}

period_frames <- make_period_frames(predsY, designs, mods, guild_sets)

get_base_site <- function(x) sub("_[123]$", "", as.character(x))

make_delta_frames <- function(period_frames) {
  vars <- predictor_names(period_frames)

  period_frames |>
    mutate(base_site = get_base_site(site)) |>
    select(guild, atlas, period, base_site, survey, site, X, Y, cwm_sti,
           expected_richness, all_of(vars)) |>
    group_split(guild) |>
    map_dfr(function(guild_df) {
      guild_name <- unique(guild_df$guild)

      map_dfr(
        list(`1_2` = c("1", "2"), `2_3` = c("2", "3"), `1_3` = c("1", "3")),
        function(pair) {
          start <- guild_df |> filter(atlas == pair[1])
          end <- guild_df |> filter(atlas == pair[2])

          joined <- inner_join(
            start,
            end,
            by = "base_site",
            suffix = c("_start", "_end")
          )

          delta_predictors <- map_dfc(vars, function(v) {
            tibble(!!paste0("delta_", v) := joined[[paste0(v, "_end")]] - joined[[paste0(v, "_start")]])
          })

          bind_cols(
            tibble(
              response = "delta_cwm_sti",
              guild = guild_name,
              comparison = paste0(year_lookup[[pair[2]]], "_minus_", year_lookup[[pair[1]]]),
              atlas_start = pair[1],
              atlas_end = pair[2],
              period_start = year_lookup[[pair[1]]],
              period_end = year_lookup[[pair[2]]],
              base_site = joined$base_site,
              survey = joined$survey_end,
              site = joined$site_end,
              X = joined$X_end,
              Y = joined$Y_end,
              cwm_sti_start = joined$cwm_sti_start,
              cwm_sti_end = joined$cwm_sti_end,
              delta_cwm_sti = joined$cwm_sti_end - joined$cwm_sti_start,
              expected_richness_start = joined$expected_richness_start,
              expected_richness_end = joined$expected_richness_end,
              expected_richness_min = pmin(joined$expected_richness_start, joined$expected_richness_end)
            ),
            delta_predictors
          )
        }
      )
    })
}

delta_frames <- make_delta_frames(period_frames)

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

summarise_gwr_model <- function(fit, response, group_name, period_or_comparison, n_rows, label, predictor_set = NA_character_) {
  if (is.null(fit)) {
    return(tibble(
      response = response,
      group = group_name,
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
    group = group_name,
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

period_splits <- period_frames |>
  group_by(guild, atlas, period) |>
  group_split()
names(period_splits) <- map_chr(period_splits, ~ paste(unique(.x$guild), unique(.x$atlas), sep = "__"))

period_models <- period_splits |>
  map(function(df) {
    guild <- unique(df$guild)
    filtered <- if (guild == "all_species") {
      df
    } else {
      df |> filter(expected_richness >= min_guild_expected_richness)
    }
    vars <- predictor_names(filtered)
    fit <- safe_fit_gwr_model(filtered, response = "cwm_sti", vars = vars, bandwidth_m = bandwidth_m)
    list(data = filtered, fit = fit, predictor_set = "full_period")
  })

delta_splits <- delta_frames |>
  group_by(guild, comparison) |>
  group_split()
names(delta_splits) <- map_chr(delta_splits, ~ paste(unique(.x$guild), unique(.x$comparison), sep = "__"))

delta_models <- delta_splits |>
  map(function(df) {
    guild <- unique(df$guild)
    filtered <- if (guild == "all_species") {
      df
    } else {
      df |> filter(expected_richness_min >= min_guild_expected_richness)
    }
    vars <- setdiff(grep("^delta_", names(filtered), value = TRUE), "delta_cwm_sti")
    fit <- safe_fit_gwr_model(filtered, response = "delta_cwm_sti", vars = vars, bandwidth_m = bandwidth_m)
    predictor_set <- "full_delta"

    if (is.null(fit)) {
      vars <- c("delta_tmean_breeding", "delta_prec_breeding", "delta_hh")
      fit <- safe_fit_gwr_model(filtered, response = "delta_cwm_sti", vars = vars, bandwidth_m = bandwidth_m)
      predictor_set <- "delta_climate_hh_fallback"
    }

    list(data = filtered, fit = fit, predictor_set = predictor_set)
  })

saveRDS(
  list(
    period_models = period_models,
    delta_models = delta_models,
    guild_sets = guild_sets,
    min_guild_expected_richness = min_guild_expected_richness,
    bandwidth_m = bandwidth_m
  ),
  file.path(model_dir, paste0(pattern, "-local-sti-gwr-models-100km.rds"))
)

period_summary <- imap_dfr(period_models, function(x, key) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  period <- unique(x$data$period)
  summarise_gwr_model(
    x$fit,
    response = "cwm_sti",
    group_name = parts[1],
    period_or_comparison = period,
    n_rows = nrow(x$data),
    label = "period",
    predictor_set = x$predictor_set
  )
})

delta_summary <- imap_dfr(delta_models, function(x, key) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  summarise_gwr_model(
    x$fit,
    response = "delta_cwm_sti",
    group_name = parts[1],
    period_or_comparison = parts[2],
    n_rows = nrow(x$data),
    label = "change",
    predictor_set = x$predictor_set
  )
})

gwr_summary <- bind_rows(period_summary, delta_summary)

period_dominant <- imap_dfr(period_models, function(x, key) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  metadata <- tibble(
    response = "cwm_sti",
    group = parts[1],
    atlas = parts[2],
    period_or_comparison = unique(x$data$period),
    model_set = "period",
    predictor_set = x$predictor_set
  )
  dominant_driver_table(x$fit, metadata, period_predictor_labels)
})

delta_dominant <- imap_dfr(delta_models, function(x, key) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  metadata <- tibble(
    response = "delta_cwm_sti",
    group = parts[1],
    atlas = NA_character_,
    period_or_comparison = parts[2],
    model_set = "change",
    predictor_set = x$predictor_set
  )
  dominant_driver_table(x$fit, metadata, delta_predictor_labels)
})

dominant_drivers <- bind_rows(period_dominant, delta_dominant)

thermal_groups <- function(sti) {
  breaks <- quantile(sti, c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE)
  cut(sti, breaks = breaks, include.lowest = TRUE,
      labels = c("very_cold", "cold", "medium", "warm", "very_warm"))
}

thermal_contribution_change <- function(preds, models) {
  sti <- models[[1]]$Tr[, "species_thermal_index"]
  groups <- thermal_groups(sti)
  common_sites <- Reduce(intersect, map(preds, ~ get_base_site(rownames(.x))))

  rel <- map(preds, function(pred_y) {
    base <- get_base_site(rownames(pred_y))
    out <- pred_y[base %in% common_sites, names(sti), drop = FALSE]
    rownames(out) <- base[base %in% common_sites]
    out <- out[order(rownames(out)), , drop = FALSE]
    out / rowSums(out)
  })

  map_dfr(
    list(`1990s_minus_1970s` = c("1", "2"), `2010s_minus_1990s` = c("2", "3"), `2010s_minus_1970s` = c("1", "3")),
    function(pair) {
      delta_rel <- rel[[pair[2]]] - rel[[pair[1]]]
      map_dfr(levels(groups), function(g) {
        sp <- names(groups)[groups == g]
        tibble(
          base_site = rownames(delta_rel),
          thermal_group = g,
          contribution = rowSums(sweep(delta_rel[, sp, drop = FALSE], 2, sti[sp], "*"))
        )
      })
    },
    .id = "comparison"
  )
}

thermal_change_contributions <- thermal_contribution_change(predsY, mods)

write_csv(period_frames, file.path(out_dir, paste0(pattern, "-cwm-sti-period-frames.csv")))
write_csv(delta_frames, file.path(out_dir, paste0(pattern, "-cwm-sti-delta-frames.csv")))
write_csv(gwr_summary, file.path(out_dir, paste0(pattern, "-local-sti-gwr-summary-100km.csv")))
write_csv(dominant_drivers, file.path(out_dir, paste0(pattern, "-local-sti-dominant-drivers-100km.csv")))
write_csv(thermal_change_contributions, file.path(out_dir, paste0(pattern, "-thermal-group-change-contributions.csv")))

summary_plot <- gwr_summary |>
  mutate(
    group = str_replace_all(group, "_", " "),
    response_label = recode(response, cwm_sti = "CWM STI", delta_cwm_sti = "Delta CWM STI")
  ) |>
  ggplot(aes(x = period_or_comparison, y = local_r2_mean, fill = response_label)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, na.rm = TRUE) +
  facet_wrap(~ group, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c("CWM STI" = "#4575b4", "Delta CWM STI" = "#d73027"), name = NULL) +
  labs(x = NULL, y = "Mean local R2", title = "Local STI-driver GWR fit at 100 km") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")

ggsave(
  file.path(plot_dir, paste0(pattern, "-local-sti-gwr-local-r2-summary.png")),
  summary_plot,
  width = 9,
  height = 7,
  dpi = 300
)

driver_plot <- dominant_drivers |>
  count(response, group, model_set, period_or_comparison, dominant_driver, direction, name = "n_cells") |>
  group_by(response, group, model_set, period_or_comparison) |>
  mutate(prop_cells = n_cells / sum(n_cells)) |>
  ungroup() |>
  filter(dominant_driver != "No significant driver") |>
  mutate(group = str_replace_all(group, "_", " "))

if (nrow(driver_plot) > 0) {
  p <- ggplot(driver_plot, aes(x = period_or_comparison, y = prop_cells, fill = dominant_driver)) +
    geom_col(width = 0.75) +
    facet_grid(response + model_set ~ group, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = driver_colours, name = "Dominant driver") +
    labs(x = NULL, y = "Proportion of significant cells", title = "Dominant local STI drivers") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")

  ggsave(
    file.path(plot_dir, paste0(pattern, "-local-sti-dominant-driver-proportions.png")),
    p,
    width = 12,
    height = 8,
    dpi = 300
  )
}

message("Finished local STI-driver analysis.")
message("Outputs written to: ", out_dir)
