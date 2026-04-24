rm(list=ls())
#### GETTING STARTED ####
if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,Hmsc,RColorBrewer,ggplot2,
               rnaturalearth,rnaturalearthdata,
               gridExtra,patchwork,sf,cowplot)
source(file.path("support_scripts", "figure_data_helpers.R"))

#### LOAD MODELS #### 
# GET PATTERNS + MODEL NAMES 
dir <- './HmscOutputs'
pattern <- '2026-03-13'
matching_folders <- figure_model_folders(pattern = pattern, base_dir = dir)
model <- matching_folders[1]
models_nums <- atlas_numbers(matching_folders)
mods <- load_hmsc_posteriors(matching_folders, base_dir = dir)
designs <- load_hmsc_study_designs(mods)

#### GET PREDS & EMPIRICAL #### 
# PREDS
predsY <- load_or_compute_site_predictions(mods, matching_folders, base_dir = dir)

# EMPIRICAL 
mod <- mods[[1]]
emps <- load_empirical_responses(mods)

#### PROBSORENSON INDEX ####
probsorensen <- function(preds.vector,emp.vector){
  # get all species 
  sp.present <- names(emp.vector)[emp.vector == 1]
  # probabilities of those species in predictions 
  sp.present.preds <- preds.vector[names(preds.vector) %in% sp.present]
  # minimum of those probabilities
  sp.present.min.preds <- min(sp.present.preds)
  # all probabilities higher than this 
  preds.greater.than.min.preds <- preds.vector[preds.vector >= sp.present.min.preds]
  # prob.sorensen
  prob.sorensen <- (2*sum(sp.present.preds)) / ((2 * sum(sp.present.preds)) + sum(preds.greater.than.min.preds))
  prob.sorensen
}

#### GET SORENSEN INDEX FOR SITES ####
sorensons <- map2(predsY,emps,function(pred.df,emp.df){
  
  # get row.names 
  row.names <- rownames(pred.df)
  # row.name <- row.names[[1]]
  # pred.df <- predsY[[1]]
  # emp.df <- emps[[1]]
  
  # get values 
  vals <- lapply(row.names,function(row.name){
    pred.vec <- pred.df[row.name,]
    emp.vec <- emp.df[row.name,]
    probsorensen(pred.vec,emp.vec)
  })
  
  data.frame(survey = row.names,
             sorensen = unlist(vals))
})


#### PREPARE MAPPING ####
dfs <- map2(sorensons,designs,function(sor,design){
  left_join(sor,design,by='survey')
})

#background 
denmark <- ne_countries(scale = "large", country = "denmark", returnclass = "sf")
# positions for annotation 
inset <- 12
df <- dfs[[1]]
x_pos <- max(df$X) - (diff(range(df$X)) / inset)
y_pos <- max(df$Y) - (diff(range(df$Y)) / inset)
# years 
year_lookup <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")
# color
pal <- colorRampPalette(brewer.pal(9,'RdYlGn'))

# get richness limits & aesthetics 
limits_df <- dfs %>%
  imap_dfr(~ data.frame(
    name = .y,
    min = min(.x$sorensen, na.rm = TRUE),
    max = max(.x$sorensen, na.rm = TRUE)
  ))
rich_limits <- c(min(limits_df$min),max(limits_df$max))
names(dfs) <- c("1", "2", "3")

### GET BOUNDING BOX AND SHAPEFILE FOR PLOTTING 
shape <- vect('~/box/PhD/logistics/data/distributions/DK5km_ED50grid_approx_kvadrkod_DOF/DK5km_ED50grid_approx_kvadrkod_DOF.shp')
shape_sf <- st_as_sf(shape)
str(dfs[[1]])

mainland_bbox <- st_bbox(c(xmin = 400000,
                           xmax = 750000,
                           ymin = 6000000,
                           ymax = 6450000),
                         crs = st_crs(25832))
bornholm_bbox <- st_bbox(c(xmin = 855000, xmax = 905000, 
                           ymin = 6100000, ymax = 6160000), 
                         crs = st_crs(25832))

# Compute inset fractions (reuse your existing bboxes)
mainland_width  <- as.numeric(mainland_bbox["xmax"] - mainland_bbox["xmin"])
mainland_height <- as.numeric(mainland_bbox["ymax"] - mainland_bbox["ymin"])
bornholm_width  <- as.numeric(bornholm_bbox["xmax"] - bornholm_bbox["xmin"])
bornholm_height <- as.numeric(bornholm_bbox["ymax"] - bornholm_bbox["ymin"])

inset_w <- bornholm_width  / mainland_width
inset_h <- bornholm_height / mainland_height

#### MAPS ####
plots <- imap(dfs, function(df, name) {
  
  year_label <- year_lookup[[name]]
  mean_sor   <- mean(df$sorensen, na.rm = TRUE)
  sd_sor     <- sd(df$sorensen,   na.rm = TRUE)
  
  # Join to shapefile
  plot_data <- shape_sf %>%
    left_join(df, by = c("kvadratkod" = "site"))
  
  make_geom <- function(data) {
    list(
      geom_sf(data = data, aes(fill = sorensen), color = 'grey30', size = 0.1),
      scale_fill_gradientn(
        colors   = pal(10),
        limits   = rich_limits,
        name     = 'Sørensen index',
        na.value = "transparent"
      )
    )
  }
  
  # --- Main plot (mainland) ---
  p_main <- ggplot() +
    make_geom(plot_data) +
    annotate("text", 
             x     = mainland_bbox["xmax"] - 0.05 * mainland_width,
             y     = mainland_bbox["ymin"] + 0.05 * mainland_height,
             label = paste('μ =', round(mean_sor, 2), '\nσ =', round(sd_sor, 2)),
             size  = 3, hjust = 1, alpha = 0.8) +
    coord_sf(
      xlim   = c(mainland_bbox["xmin"], mainland_bbox["xmax"]),
      ylim   = c(mainland_bbox["ymin"], mainland_bbox["ymax"]),
      expand = FALSE
    ) +
    labs(x = NULL, y = NULL, title = year_label, fill = 'Sørensen index') +
    theme_minimal() +
    theme(
      legend.position = 'none',
      plot.background = element_rect(fill = "white", color = NA),
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text       = element_blank(),
      panel.grid      = element_blank()
    )
  
  # --- Inset (Bornholm) ---
  p_inset <- ggplot() +
    make_geom(plot_data) +
    coord_sf(
      xlim   = c(bornholm_bbox["xmin"], bornholm_bbox["xmax"]),
      ylim   = c(bornholm_bbox["ymin"], bornholm_bbox["ymax"]),
      expand = FALSE
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      panel.border    = element_rect(color = "grey30", fill = NA, linewidth = 0.5)
    )
  
  p_final <- ggdraw(p_main) +
    draw_plot(
      p_inset,
      x      = 1 - inset_w - 0.02,
      y      = 1 - inset_h - 0.02,
      width  = inset_w,
      height = inset_h
    )
  
  list(
    grob   = as_grob(p_final),
    ggplot = p_main
  )
})


#### LAYOUT ####
leg_richness <- as_grob(get_legend(
  plots[[1]]$ggplot + theme(legend.position = "bottom")
))

layout <- "
AB
CD
"

final_layout <- wrap_plots(
  A = plots[[1]]$grob,
  B = plots[[2]]$grob,
  C = plots[[3]]$grob,
  D = leg_richness,
  design = layout
) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = list(c('A', 'B', 'C', ''))) &
  theme(plot.tag = element_text(face = "bold", size = 14))

final_layout
ggsave(sprintf('misc-figures/%s-sfig2-sorensen-similarity.png',pattern),final_layout,width=7.5,height=7.5)




