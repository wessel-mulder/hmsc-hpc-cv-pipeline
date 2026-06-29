rm(list = ls())

# GETTING STARTED  --------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, Hmsc, cowplot, ggbeeswarm, patchwork, scales)

source(file.path("support_scripts", "figure_data_helpers.R"))

base_dir <- "./HmscOutputs"
pattern <- "2026-03-13"
out_dir <- file.path("misc-figures", "outputs", "main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
trait_metadata_path <- file.path("Data", "preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData")
min_species_per_trait_level <- 5

matching_folders <- figure_model_folders(pattern = pattern, base_dir = base_dir)
mods <- load_hmsc_posteriors(matching_folders, base_dir = base_dir)
designs <- load_hmsc_study_designs(mods)
preds_y <- load_or_compute_site_predictions(mods, matching_folders, base_dir = base_dir)

period_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")
period_levels <- unname(period_lookup)
thermal_levels <- c("very cold", "cold", "medium", "warm", "very warm")
migration_levels <- c(
  "long-distance",
  "short-and long-distance",
  "short-distance",
  "sedentary and short-distance",
  "sedentary"
)

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

base_site <- function(x) sub("_[123]$", "", x)

sti_sp <- mods[[1]]$Tr[, "species_thermal_index"]
species <- Reduce(intersect, c(list(names(sti_sp)), map(preds_y, colnames)))
sti_sp <- sti_sp[species]
preds_y <- map(preds_y, ~ .x[, species, drop = FALSE])

trait_env <- new.env(parent = emptyenv())
load(trait_metadata_path, envir = trait_env)
trait_metadata <- trait_env$Tr

thermal_breaks <- quantile(sti_sp, seq(0, 1, 0.2), na.rm = TRUE)
sp_group <- cut(
  sti_sp,
  breaks = thermal_breaks,
  include.lowest = TRUE,
  labels = thermal_levels
)
names(sp_group) <- names(sti_sp)

# The regression helper below only needs named species-group vectors. Keep the
# thermal grouping from STI, then add the two categorical trait groupings from Tr.
migration_group <- trait_metadata[, "Migration_a3_DOF"]
names(migration_group) <- rownames(trait_metadata)

guild_group <- trait_metadata[, "foraging_guild_consensus"]
names(guild_group) <- rownames(trait_metadata)

# Add one reference group that ignores trait membership, so the trend panel can
# show the whole-community VP trend beside the trait-specific summaries.
all_species_group <- setNames(rep("All species", length(species)), species)

#### LOAD VP DATA ####
pattern <- "2026-03-13"
matching_folders <- figure_model_folders(pattern = pattern)
model <- matching_folders[1]
models_nums <- atlas_numbers(matching_folders)
scaled = T # set to NULL to not scale 
VPs <- load_vp_estimates(matching_folders, scaled = scaled)


# PREPARE FOR PLOT  -------------------------------------------------------

# Keep the environmental predictors in a stable manuscript-facing order. Any
# extra VP rows are appended after these known variables, but the random site
# effect is excluded from the temporal slope panels.
preferred_variable_order <- c(
  "tmean_breeding",
  "prec_breeding",
  "hh",
  "perc_urban",
  "perc_cropland",
  "perc_pasture",
  "perc_forest",
  "perc_grass_shrub"
)
variables <- setdiff(rownames(VPs[[1]]), "Random: site")
variables <- c(
  intersect(preferred_variable_order, variables),
  setdiff(variables, preferred_variable_order)
)

variable_label_lookup <- c(
  "tmean_breeding" = "Temperature",
  "prec_breeding" = "Precipitation",
  "hh" = "Habitat\nheterogeneity",
  "perc_urban" = "Urban",
  "perc_cropland" = "Cropland",
  "perc_pasture" = "Pasture",
  "perc_forest" = "Forest",
  "perc_grass_shrub" = "Grass/\nshrubland"
)

clean_variable_label <- function(variable) {
  dplyr::recode(
    variable,
    !!!variable_label_lookup,
    .default = stringr::str_to_sentence(stringr::str_replace_all(variable, "_", " "))
  )
}

variable_labels <- clean_variable_label(variables)

# Fit one temporal regression per environmental variable and species group.
# The response is the Tjur-R2-scaled VP value, because `scaled = TRUE` above.
fit_group_vp_regressions <- function(group_vector, group_family, group_levels) {
  group_vector <- group_vector[!is.na(group_vector) & group_vector != ""]
  group_levels <- intersect(group_levels, unique(as.character(group_vector)))
  
  purrr::map_dfr(variables, function(variable) {
    purrr::map_dfr(group_levels, function(group_name) {
      group_species <- names(group_vector)[as.character(group_vector) == group_name]
      
      VPs_g <- purrr::imap_dfr(VPs, function(vp_obj, atlas_name) {
        vp <- if (is.list(vp_obj) && "vals" %in% names(vp_obj)) {
          vp_obj$vals
        } else {
          vp_obj
        }
        
        species_keep <- intersect(group_species, colnames(vp))
        
        if (!variable %in% rownames(vp) || length(species_keep) == 0) {
          return(tibble::tibble())
        }
        
        vals <- as.numeric(as.matrix(vp[variable, species_keep, drop = FALSE]))
        
        tibble::tibble(
          atlas = atlas_name,
          atlas_numeric = as.numeric(stringr::str_extract(as.character(atlas_name), "\\d+")),
          variable = variable,
          trait_family = group_family,
          trait_value = group_name,
          species = species_keep,
          vp = vals
        )
      }) |>
        dplyr::filter(is.finite(.data$vp), !is.na(.data$atlas_numeric))
      
      if (nrow(VPs_g) < 3 || length(unique(VPs_g$atlas_numeric)) < 2) {
        return(tibble::tibble(
          variable = variable,
          trait_family = group_family,
          trait_value = group_name,
          n = nrow(VPs_g),
          n_species = dplyr::n_distinct(VPs_g$species),
          r_squared = NA_real_,
          slope = NA_real_,
          direction = NA_character_,
          p_value = NA_real_,
          significance = NA_character_
        ))
      }
      
      m <- lm(vp ~ atlas_numeric, data = VPs_g)
      sm <- summary(m)
      coef_table <- coef(sm)
      
      slope <- coef_table["atlas_numeric", "Estimate"]
      p_value <- coef_table["atlas_numeric", "Pr(>|t|)"]
      
      tibble::tibble(
        variable = variable,
        trait_family = group_family,
        trait_value = group_name,
        n = nrow(VPs_g),
        n_species = dplyr::n_distinct(VPs_g$species),
        r_squared = sm$r.squared,
        slope = slope,
        direction = dplyr::case_when(
          slope > 0 ~ "increasing",
          slope < 0 ~ "decreasing",
          TRUE ~ "flat"
        ),
        p_value = p_value,
        significance = dplyr::case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          p_value < 0.1 ~ ".",
          TRUE ~ "ns"
        )
      )
    })
  })
}

guild_counts <- sort(table(guild_group), decreasing = TRUE)
guild_levels <- names(guild_counts[guild_counts >= min_species_per_trait_level])

regression_results <- dplyr::bind_rows(
  fit_group_vp_regressions(all_species_group, "All species", "All species"),
  fit_group_vp_regressions(sp_group, "Thermal affinity", thermal_levels),
  fit_group_vp_regressions(migration_group, "Migratory strategy", migration_levels),
  fit_group_vp_regressions(guild_group, "Foraging guild", guild_levels)
) |>
  dplyr::mutate(
    variable_clean = factor(clean_variable_label(.data$variable), levels = variable_labels),
    sig_label = dplyr::case_when(
      is.na(.data$significance) | .data$significance == "ns" ~ "",
      TRUE ~ .data$significance
    )
  )

regression_csv <- file.path(out_dir, paste0(pattern, "-fig9-temporal-vp-regression-results.csv"))
readr::write_csv(regression_results, regression_csv)

slope_limit <- max(abs(regression_results$slope), na.rm = TRUE)
slope_limit <- ifelse(is.finite(slope_limit) && slope_limit > 0, slope_limit * 1.05, 1)

vp_panel_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 8.4, angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8.6),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),
    legend.title = element_text(face = "bold")
  )

make_vp_slope_panel <- function(plot_data, title, y_levels, show_x_labels = FALSE) {
  plot_data |>
    dplyr::mutate(
      trait_value = factor(as.character(.data$trait_value), levels = y_levels)
    ) |>
    ggplot(aes(x = .data$variable_clean, y = .data$trait_value, fill = .data$slope)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(
      aes(label = .data$sig_label),
      colour = "grey10",
      fontface = "bold",
      size = 3.1
    ) +
    scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      limits = c(-slope_limit, slope_limit),
      oob = scales::squish,
      breaks = scales::pretty_breaks(n = 5),
      labels = scales::label_number(accuracy = 0.01),
      name = "Temporal slope",
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(3.2, "in"),
        barheight = grid::unit(0.18, "in")
      )
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(title = title) +
    vp_panel_theme +
    theme(
      axis.text.x = if (show_x_labels) {
        element_text(size = 8.4, angle = 35, hjust = 1, vjust = 1)
      } else {
        element_blank()
      }
    )
}

all_species_panel <- regression_results |>
  dplyr::filter(.data$trait_family == "All species") |>
  make_vp_slope_panel(
    title = "All species",
    y_levels = "All species",
    show_x_labels = FALSE
  )

thermal_panel <- regression_results |>
  dplyr::filter(.data$trait_family == "Thermal affinity") |>
  make_vp_slope_panel(
    title = "Thermal affinity groups",
    y_levels = thermal_levels,
    show_x_labels = FALSE
  )

migration_panel <- regression_results |>
  dplyr::filter(.data$trait_family == "Migratory strategy") |>
  make_vp_slope_panel(
    title = "Migratory strategies",
    y_levels = rev(migration_levels),
    show_x_labels = FALSE
  )

guild_panel <- regression_results |>
  dplyr::filter(.data$trait_family == "Foraging guild") |>
  make_vp_slope_panel(
    title = "Foraging guilds",
    y_levels = rev(guild_levels),
    show_x_labels = TRUE
  )

vp_slope_panels <- all_species_panel / thermal_panel / migration_panel / guild_panel +
  patchwork::plot_layout(heights = c(0.42, 0.85, 0.9, 1.45), guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 14),
    plot.background = element_rect(fill = "white", colour = NA)
  )

vp_slope_png <- file.path(out_dir, paste0(pattern, "-fig9-temporal-vp-slope-panels.png"))

ggsave(vp_slope_png, vp_slope_panels, width = 11, height = 12.2, units = "in", dpi = 300, bg = "white")

message("Saved: ", regression_csv)
message("Saved: ", vp_slope_png)


# SPECIES-LEVEL DELTAS  ----------------------------------------------------

# This second section compares periods directly instead of fitting a temporal
# regression. Positive deltas mean the predictor explains more of a species'
# occurrence pattern in the later atlas period; negative deltas mean less.
vp_long <- purrr::imap_dfr(VPs, function(vp, atlas_name) {
  vp[variables, species, drop = FALSE] |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "variable") |>
    tidyr::pivot_longer(
      cols = -variable,
      names_to = "species",
      values_to = "vp"
    ) |>
    dplyr::mutate(
      atlas = as.character(atlas_name),
      period = unname(period_lookup[.data$atlas])
    )
})

delta_pairs <- tibble::tribble(
  ~from_atlas, ~to_atlas, ~comparison,
  "1", "3", "1970s to 2010s",
  "2", "3", "1990s to 2010s",
  "1", "2", "1970s to 1990s"
)

species_metadata <- tibble::tibble(
  species = species,
  species_label = stringr::str_replace_all(species, "_", " "),
  species_thermal_index = as.numeric(sti_sp[species]),
  thermal_group = as.character(sp_group[species]),
  migratory_strategy = as.character(migration_group[species]),
  foraging_guild = as.character(guild_group[species])
) |>
  dplyr::mutate(
    thermal_group = factor(.data$thermal_group, levels = thermal_levels)
  )

vp_species_deltas <- purrr::pmap_dfr(delta_pairs, function(from_atlas, to_atlas, comparison) {
  from_values <- vp_long |>
    dplyr::filter(.data$atlas == from_atlas) |>
    dplyr::select(variable, species, vp_from = vp)
  
  to_values <- vp_long |>
    dplyr::filter(.data$atlas == to_atlas) |>
    dplyr::select(variable, species, vp_to = vp)
  
  dplyr::left_join(
    from_values,
    to_values,
    by = c("variable", "species")
  ) |>
    dplyr::mutate(
      from_atlas = from_atlas,
      to_atlas = to_atlas,
      comparison = comparison,
      vp_delta = .data$vp_to - .data$vp_from,
      delta_direction = dplyr::case_when(
        .data$vp_delta > 0 ~ "increase",
        .data$vp_delta < 0 ~ "decrease",
        TRUE ~ "no change"
      )
    )
}) |>
  dplyr::left_join(species_metadata, by = "species") |>
  dplyr::mutate(
    comparison = factor(.data$comparison, levels = delta_pairs$comparison),
    variable_clean = factor(clean_variable_label(.data$variable), levels = variable_labels)
  )

species_order <- species_metadata |>
  dplyr::arrange(.data$thermal_group, .data$species_thermal_index, .data$species_label) |>
  dplyr::pull(.data$species_label)

vp_species_deltas <- vp_species_deltas |>
  dplyr::mutate(
    species_label = factor(.data$species_label, levels = species_order)
  )

delta_csv <- file.path(out_dir, paste0(pattern, "-fig9-species-vp-deltas.csv"))
readr::write_csv(vp_species_deltas, delta_csv)

delta_limit <- max(abs(vp_species_deltas$vp_delta), na.rm = TRUE)
delta_limit <- ifelse(is.finite(delta_limit) && delta_limit > 0, delta_limit * 1.05, 1)

vp_species_delta_plot <- ggplot(
  vp_species_deltas,
  aes(x = .data$variable_clean, y = .data$species_label, fill = .data$vp_delta)
) +
  geom_tile(colour = "white", linewidth = 0.08) +
  facet_grid(. ~ comparison) +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-delta_limit, delta_limit),
    oob = scales::squish,
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::label_number(accuracy = 0.01),
    name = "Delta variance\nexplained",
    guide = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(3.6, "in"),
      barheight = grid::unit(0.18, "in")
    )
  ) +
  labs(
    title = "Species-level change in variance explained",
    subtitle = "Tiles show later atlas minus earlier atlas for Tjur-R2-scaled variance partitioning.",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 7.8, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 4.2, lineheight = 0.82),
    strip.text = element_text(face = "bold", size = 9.5),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 8.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

#### THERMS 
thermal_cols <- c("#313695", "#74ADD1", "#FEE090", "#F46D43", "#9E0142")
                       
subset <- vp_species_deltas |> 
  filter(comparison == '1970s to 2010s',
         variable == 'tmean_breeding')

ggplot(subset,aes(x = thermal_group,
                  y = vp_delta,
                  col = thermal_group)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  scale_color_manual(values = thermal_cols)
#### MIGS 
migration_palette <- setNames(
  colorRampPalette(c("#542788", "#f2e5ff"))(length(migration_levels)),
  migration_levels
)

subset <- vp_species_deltas |> 
  filter(comparison == '1970s to 2010s',
         variable == 'tmean_breeding')

ggplot(subset,aes(x = migratory_strategy,
                  y = vp_delta)) +
  geom_violin()+
  geom_boxplot(width = 0.1)

  scale_color_manual(migration_palette)


vp_delta_png <- file.path(out_dir, paste0(pattern, "-fig9-species-vp-delta-panels.png"))
ggsave(vp_delta_png, vp_species_delta_plot, width = 13, height = 18, units = "in", dpi = 300, bg = "white")
message("Saved: ", delta_csv)
message("Saved: ", vp_delta_png)
