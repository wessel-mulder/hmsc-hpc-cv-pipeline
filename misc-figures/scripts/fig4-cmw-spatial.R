rm(list=ls())
#### GETTING STARTED ####
if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,Hmsc,RColorBrewer,ggplot2,
               rnaturalearth,rnaturalearthdata,
               gridExtra,patchwork,sf,cowplot,
               terra)
source(file.path("support_scripts", "figure_data_helpers.R"))

#### LOAD MODELS #### 
dir <- './HmscOutputs'
pattern <- '2026-03-13'

matching_folders <- figure_model_folders(pattern = pattern, base_dir = dir)
model <- matching_folders[1]
models_nums <- atlas_numbers(matching_folders)
mods <- load_hmsc_posteriors(matching_folders, base_dir = dir)
designs <- load_hmsc_study_designs(mods)


#### GET PREDS #### 
predsY <- load_or_compute_site_predictions(mods, matching_folders, base_dir = dir)

str(predsY)

CWMs <- community_weighted_means(predsY, mods)

# prepare for plotting 
dfs <- map2(CWMs, designs, ~ {
  # Calculate row sums and keep it as a data frame
  .x %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'survey') %>% 
    # Now join with .y (the design list element)
    left_join(.y, by = 'survey')
})

var <- 'species_thermal_index'

# get richness limits & aesthetics 
lims <- dfs %>%
  imap_dfr(~ data.frame(
    name = .y,
    min = min(.x %>% pull(var), na.rm = TRUE),
    max = max(.x %>% pull(var), na.rm = TRUE)
  ))
lim <- c(min(lims$min),max(lims$max))
my_breaks <- seq(lim[1], lim[2], length.out = 5)

# colors 
pal <- colorRampPalette(c("#313695", "#4575b4", "#74add1", "#abd9e9", 
                          "#ffffbf", 
                          "#fdae61", "#f46d43", "#d73027", "#a50026"))
#pal <- colorRampPalette(colorRamps::blue2red(n=10))

#background 
denmark <- ne_countries(scale = "large", country = "denmark", returnclass = "sf")
# positions for annotation 
inset <- 12
df <- dfs[[3]]
x_pos <- max(df$X) - (diff(range(df$X)) / inset)
y_pos <- max(df$Y) - (diff(range(df$Y)) / inset)
# years 
year_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")

#### READ ATLAS POLYGON ####
shape <- vect('~/box/PhD/logistics/data/distributions/DK5km_ED50grid_approx_kvadrkod_DOF/DK5km_ED50grid_approx_kvadrkod_DOF.shp')
shape_sf <- st_as_sf(shape)
str(dfs[[1]])

# Define a bounding box for Bornholm (approximate UTM 32N)
# You might need to tweak these numbers based on your specific grid extent
# Increase the range between min and max to "zoom out"
mainland_bbox <- st_bbox(c(xmin = 400000,
                           xmax = 750000,
                           ymin = 6000000,
                           ymax = 6450000),
                         crs = st_crs(25832))
bornholm_bbox <- st_bbox(c(xmin = 855000, xmax = 905000, 
                           ymin = 6100000, ymax = 6160000), 
                         crs = st_crs(25832))

# Compute inset size as fraction of main plot dimensions
mainland_width  <- as.numeric(mainland_bbox["xmax"] - mainland_bbox["xmin"])
mainland_height <- as.numeric(mainland_bbox["ymax"] - mainland_bbox["ymin"])
bornholm_width  <- as.numeric(bornholm_bbox["xmax"] - bornholm_bbox["xmin"])
bornholm_height <- as.numeric(bornholm_bbox["ymax"] - bornholm_bbox["ymin"])

inset_w <- bornholm_width  / mainland_width
inset_h <- bornholm_height / mainland_height

#### MAPS ####
plots <- imap(dfs, function(df, name) {
  
  year_label <- year_lookup[[name]]
  
  plot_data <- shape_sf %>%
    left_join(df, by = c("kvadratkod" = "site"))
  
  # --- shared geom builder to keep fill scale consistent ---
  make_geom <- function(data) {
    list(
      geom_sf(data = data, aes(fill = .data[[var]]), color = 'grey30', size = 0.1),
      scale_fill_gradientn(
        colors = pal(10),
        limits = lim,
        breaks = my_breaks,
        labels = round(my_breaks, 1),
        na.value = "transparent",
        name = expression(CWM[STI])
      )
    )
  }
  
  # --- Main plot (mainland) ---
  p_main <- ggplot() +
    make_geom(plot_data) +
    coord_sf(
      xlim = c(mainland_bbox["xmin"], mainland_bbox["xmax"]),
      ylim = c(mainland_bbox["ymin"], mainland_bbox["ymax"]),
      expand = FALSE
    ) +
    labs(x = NULL, y = NULL, title = year_label, fill = var) +
    theme_minimal() +
    theme(
      legend.position  = "none",
      plot.background  = element_rect(fill = "white", color = NA),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text        = element_blank(),
      panel.grid       = element_blank(),
      plot.margin      = margin(0, 0, 0, 0)  # remove all map margins
    )
  
  # --- Inset plot (Bornholm) ---
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
  
  # --- Combine with cowplot ---
  p_final <- ggdraw(p_main) +
    draw_plot(
      p_inset,
      x = 1 - inset_w - 0.2,
      y = 1 - inset_h - 0.2,
      width  = inset_w,
      height = inset_h
    )
  
  p_final <- as_grob(p_final)  # convert the complete object, not mid-pipe  
  
  list(
    grob = as_grob(p_final),
    ggplot = p_main  # keep original for legend extraction
  )
  
})

# WARMING DRIVEN BY?  -----------------------------------------------------

# ── 0. Extract species thermal indices ────────────────────────────────────────
# Assuming STI per species is a named vector from mods[[1]]$Tr
# ── 0. Species STI vector ─────────────────────────────────────────────────────
sti_sp <- mods[[1]]$Tr[, 'species_thermal_index']  # named numeric vector, length 157

# ── 1. Classify species into 5 thermal groups ─────────────────────────────────
breaks <- quantile(sti_sp, c(0, 0.2, 0.4, 0.6, 0.8, 1))
sp_group <- cut(sti_sp, breaks = breaks, include.lowest = TRUE,
                labels = c("very cold", "cold", "medium", "warm", "very warm"))
names(sp_group) <- names(sti_sp)

# ── 2. Shared sites only ──────────────────────────────────────────────────────
get_base_site <- function(x) sub("_[123]$", "", x)

preds_shared <- map(predsY, ~ {
  base <- get_base_site(rownames(.x))
  shared <- Reduce(intersect, map(predsY, ~ get_base_site(rownames(.x))))
  out <- .x[base %in% shared, ]
  rownames(out) <- get_base_site(rownames(out))
  out[order(rownames(out)), ]
})

# ── 3. Relative abundance (rows sum to 1) ─────────────────────────────────────
rel_abund <- map(preds_shared, ~ .x / rowSums(.x))

# ── 4. Contribution of each species to CWM = rel_abund * STI ─────────────────
# Result: cells x species matrix of STI contributions, per atlas
contrib <- map(rel_abund, ~ sweep(.x, 2, sti_sp, "*"))

# ── 5. Sum contributions by thermal group, per cell, per atlas ───────────────
groups <- c("very cold", "cold", "medium", "warm", "very warm")

contrib_by_group <- imap(contrib, ~ {
  map_dfc(groups, function(g) {
    sp <- names(sp_group)[sp_group == g]
    tibble(!!g := rowSums(.x[, sp]))
  }) %>%
    mutate(
      site  = rownames(.x),
      atlas = .y
    )
}) %>%
  bind_rows()

# ── 1. Find dominant thermal group per cell per atlas ─────────────────────────
dominant_group <- contrib_by_group %>%
  pivot_longer(all_of(groups), names_to = "group", values_to = "contribution") %>%
  group_by(site, atlas) %>%
  slice_max(contribution, n = 1) %>%
  ungroup() %>%
  mutate(
    group = factor(group, levels = groups),
    atlas = factor(atlas, labels = c("1970s", "1990s", "2010s"))
  )

# ── 2. Join to CWM STI values ─────────────────────────────────────────────────
cwm_sti <- imap(preds_shared, ~ {
  # recompute CWM STI per cell
  rel <- .x / rowSums(.x)
  cwm <- rel %*% sti_sp
  tibble(site = rownames(.x), sti = as.numeric(cwm), atlas = .y)
}) %>%
  bind_rows() %>%
  mutate(atlas = factor(atlas, labels = c("1970s", "1990s", "2010s")))

# ── 3. Join dominant group to CWM STI ─────────────────────────────────────────
plot_df <- cwm_sti %>%
  left_join(dominant_group %>% select(site, atlas, group), 
            by = c("site", "atlas"))


library(patchwork)
library(ggbeeswarm)



# ── BEES WITH SQUARES ────────────────────────────────────────────────────────
atlas_name <- '1970s'

pal <- c(
  "very cold" = "#313695",
  "cold"      = "#abd9e9",
  "medium"    = "#d9d9d9",
  "warm"      = "#fdae61",
  "very warm" = "#a50026"
)

make_bee_with_squares <- function(atlas_name, show_x = FALSE) {
  d <- plot_df %>% filter(atlas == atlas_name)
  
  # ensure all groups present, fill 0 if missing
  p <- tibble(group = factor(groups, levels = groups)) %>%
    left_join(
      pct_df %>% filter(atlas == atlas_name),
      by = "group"
    ) %>%
    mutate(pct = replace_na(pct, 0))
  
  sti_range <- range(plot_df$sti)
  sti_span  <- diff(sti_range)
  n_groups  <- nrow(p)
  
  # max box half-width — for 100%
  box_hw_max <- sti_span * 0.09
  
  # scale hw by sqrt(pct/100) so area proportional to pct
  # all boxes coloured regardless, even 0%
  p <- p %>%
    mutate(
      xcentre = seq(sti_range[1] + sti_span * 0.1,
                    sti_range[2] - sti_span * 0.1,
                    length.out = n_groups),
      hw      = box_hw_max * sqrt(pct / 100),
      # minimum visible size even for 0%
      hw      = pmax(hw, box_hw_max * 0.75),
      ycentre = 1.62 + box_hw_max
    )
  
  ggplot() +
    # ── scaled coloured squares, centred on fixed positions ───────────────────
    geom_rect(data = p,
              aes(xmin = xcentre - hw, xmax = xcentre + hw,
                  ymin = ycentre - hw, ymax = ycentre + hw,
                  fill = group),
              colour = "white", linewidth = 0.4) +
    # ── percentage label — always white ───────────────────────────────────────
    geom_text(data = p,
              aes(x = xcentre, y = ycentre,
                  label = paste0(round(pct), "%")),
              size = 2.8, fontface = "bold", colour = "white") +
    # ── group name below ──────────────────────────────────────────────────────
    # geom_text(data = p,
    #           aes(x = xcentre, y = 1.57, label = group),
    #           size = 2.2, colour = "grey40", hjust = 0.5) +
    # ── beeswarm coloured by group ────────────────────────────────────────────
    geom_beeswarm(data = d,
                  aes(x = sti, y = 1, colour = group),
                  size = 0.7, alpha = 0.6, cex = 0.4, orientation = "y") +
    geom_boxplot(data = d,
                 aes(x = sti, y = 1),
                 width = 0.1, outlier.shape = NA,
                 fill = NA, linewidth = 0.3) +
    #coord_fixed(ratio = diff(c(0.7, 2.2)) / diff(range(plot_df$sti))) +
    scale_fill_manual(values = pal) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(limits = c(sti_range[1] - 0.02, sti_range[2] + 0.02)) +
    scale_y_continuous(limits = c(0.7, 2.2)) +
    labs(x = if (show_x) expression(CWM[STI]) else NULL, y = NULL) +
    theme_minimal() +
    theme(
      legend.position    = "none",
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.x        = if (show_x) element_text() else element_blank(),
      axis.ticks.x       = if (show_x) element_line() else element_blank(),
      plot.margin        = margin(2, 5, 2, 5)
    )
    
}

# ── BUILD ROWS ────────────────────────────────────────────────────────────────
# ── BUILD ROWS ────────────────────────────────────────────────────────────────
# ── LEFT COLUMN: maps stacked ─────────────────────────────────────────────────
left_col <- plot_grid(
  plots[[1]]$grob,
  plots[[2]]$grob,
  plots[[3]]$grob,
  ncol        = 1,
  labels      = c("A", "B", "C"),
  label_size  = 12,
  align       = "v",
  axis        = "lr"
)

# ── RIGHT COLUMN: beeswarm plots stacked ─────────────────────────────────────
# bottom row shows x axis, others don't
bee_70 <- make_bee_with_squares("1970s", show_x = FALSE)
bee_90 <- make_bee_with_squares("1990s", show_x = FALSE)
bee_10 <- make_bee_with_squares("2010s", show_x = TRUE)

right_col <- plot_grid(
  bee_70, bee_90, bee_10,
  ncol        = 1,
  align       = "v",
  axis        = "lr",
)

# ── COMBINE COLUMNS ───────────────────────────────────────────────────────────
main <- plot_grid(
  left_col, right_col,
  ncol       = 2,
  rel_widths = c(0.5, 1),
  greedy = T
)

# ── LEGENDS ───────────────────────────────────────────────────────────────────
map_legend <- get_legend(
  plots[[1]]$ggplot + 
    guides(fill = guide_colourbar(barwidth = 8, barheight = 0.8)) +
    theme(legend.position = "bottom")
)

bee_legend <- get_legend(
  ggplot(pct_df, aes(x = "", y = pct, fill = group)) +
    geom_col() +
    scale_fill_manual(
      values = pal,
      breaks = groups,
      name   = "Dominant thermal affinity"
    ) +
    guides(fill = guide_legend(nrow = 2)) +
    theme_minimal() +
    theme(legend.position = "bottom")
)

both_legends <- plot_grid(map_legend, bee_legend, ncol = 2)

# ── 3. Reduce rel_widths gap and add negative margins ─────────────────────────
main <- plot_grid(
  left_col, right_col,
  ncol       = 2,
  rel_widths = c(0.5, 1)  # shrink map column
)

final <- plot_grid(
  main, both_legends,
  ncol        = 1,
  rel_heights = c(1, 0.1)
)

final_draw <- ggdraw(final) + 
  theme(
    plot.background = element_rect(fill = "white", colour = NA)
  )


ggsave(sprintf('misc-figures/%s-fig4-cwm-maps-dominance.png',pattern),width = 10,height= 7.5)
