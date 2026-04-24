rm(list=ls())
#### GETTING STARTED ####
if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,Hmsc,RColorBrewer,ggplot2,
               rnaturalearth,rnaturalearthdata,
               gridExtra,patchwork,sf,cowplot)
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

#### GET RICHNESS ####
# prepare for plotting 
dfs <- predicted_richness_frames(predsY, designs)

# Atlas 1 to Atlas 2
diff_1_2 <- dfs[[1]] %>%
  select(site, X, Y, richness_1 = richness) %>%
  inner_join(select(dfs[[2]], site, richness_2 = richness), by = "site") %>%
  mutate(diff = richness_2 - richness_1)

# Atlas 2 to Atlas 3
diff_2_3 <- dfs[[2]] %>%
  select(site, X, Y, richness_2 = richness) %>%
  inner_join(select(dfs[[3]], site, richness_3 = richness), by = "site") %>%
  mutate(diff = richness_3 - richness_2)

# 2. Combine into a named list
diff_list <- list(
  atlas_1_to_2 = diff_1_2,
  atlas_2_to_3 = diff_2_3
)

# 3. Quick check of the results
str(diff_list, max.level = 1)


# get richness limits & aesthetics 
limits_df <- diff_list %>%
  imap_dfr(~ data.frame(
    name = .y,
    min_diff = min(.x$diff, na.rm = TRUE),
    max_diff = max(.x$diff, na.rm = TRUE)
  ))
rich_limits <- c(min(limits_df$min_diff),max(limits_df$max_diff))
lims <- c(-max(abs(rich_limits)),max(abs(rich_limits)))
my_breaks <- seq(rich_limits[1], rich_limits[2], length.out = 5)

# colors 
pal <- rev(brewer.pal(11,'RdYlBu'))

#background 
denmark <- ne_countries(scale = "large", country = "denmark", returnclass = "sf")
# positions for annotation 
inset <- 12
df <- dfs[[3]]
x_pos <- max(df$X) - (diff(range(df$X)) / inset)
y_pos <- max(df$Y) - (diff(range(df$Y)) / inset)
# years 
year_lookup <- c("atlas_1_to_2" = "1970s to 1990s", "atlas_2_to_3" = "1990s to 2010s")



#### RICHNESS MAPS ####
plots_richness <- imap(diff_list, function(df, name) {
  
  year_label <- year_lookup[[name]]
  
  ggplot() +
    # Background: The actual country shape
    geom_sf(data = denmark, fill='white', color = "grey80", size = 0.5) +
    
    # The Data Points: Using your ivory border (pch 21)
    geom_point(data = df, aes(x = X, y = Y, fill = diff),
               size = 1, pch = 21,color='white', alpha = 1) +
    
    # The Color Scale
    # Inside your plot loop:
    scale_fill_gradient2(
      low = pal[1],
      mid = 'white',
      high = pal[11],
      midpoint = 0,
      limits = lims,
      #breaks = my_breaks,             # Force these breaks
      #labels = round(my_breaks, 0)    # Force matching labels
    )+
    
    # coords 
    coord_sf(crs = st_crs(23032)) + # Forces the whole plot into UTM 32N +
    
    labs(x = NULL,
         y = NULL,
         title = year_label,
         fill = 'Δ Richness') +
    
    # annotate("text", x = x_pos, y = y_pos, 
    #         label = "Atlas 1 - 1970s", 
    #         size = 8,
    #         #hjust = 1.1, vjust = 1.5,
    #         fill = "white",         # Background color
    #         alpha = 0.8) +           # Slight transparency
    
    # theme 
    theme_minimal() + 
    theme(
      legend.position = 'none',
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_blank()
    ) 
})

#### PREPARE SORENSEN ####
#### GET SORENSEN INDEX FOR SITES ####
filter_common_sites <- function(preds_list, comparison_indices) {
  
  # 1. Identify common sites by stripping suffixes
  idx1 <- comparison_indices[1]
  idx2 <- comparison_indices[2]
  
  site_names_1 <- gsub("_[0-9]$", "", rownames(preds_list[[idx1]]))
  site_names_2 <- gsub("_[0-9]$", "", rownames(preds_list[[idx2]]))
  
  common_sites <- intersect(site_names_1, site_names_2)
  common_order <- sort(common_sites)
  
  # 2. Process first Atlas
  # Filter, Rename, and Sort
  mat1 <- preds_list[[idx1]][site_names_1 %in% common_sites, , drop = FALSE]
  rownames(mat1) <- gsub("_[0-9]$", "", rownames(mat1))
  mat1 <- mat1[common_order, , drop = FALSE]
  
  # 3. Process second Atlas
  # Filter, Rename, and Sort
  mat2 <- preds_list[[idx2]][site_names_2 %in% common_sites, , drop = FALSE]
  rownames(mat2) <- gsub("_[0-9]$", "", rownames(mat2))
  mat2 <- mat2[common_order, , drop = FALSE]
  
  # 4. Return as a list with two elements
  return(list(mat1, mat2))
}

prepped <- list(filter_common_sites(predsY,c(1,2)),
                filter_common_sites(predsY,c(2,3)))
dfs <- lapply(prepped,function(prep_choice){
  row.names <- rownames(prep_choice[[1]])
  row.name <- row.names[1]
  vals <- lapply(row.names,function(row.name){
    pred.vec <- prep_choice[[1]][row.name,]
    emp.vec <- prep_choice[[2]][row.name,]
    intersect <- sum(pmin(pred.vec,emp.vec))
    mass <- sum(pred.vec) + sum(emp.vec)
    prob_sim <- (2 * intersect) / mass
  }) 
  df <-   data.frame(survey = row.names,
                     sorensen = unlist(vals))
  return(df)
})

# prepare for plotting 
dfs_mapped <- map(dfs, ~ {
  # .x is the current element of sorensons_list
  # .y is the current element of designs
  .x %>% 
    as_tibble() %>% 
    rename(site = survey) %>% 
    left_join(designs[[3]], by = "site") %>% 
    select(site, X, Y, sorensen)
})


# get richness limits & aesthetics 
limits_df <- dfs_mapped %>%
  imap_dfr(~ data.frame(
    name = .y,
    min_diff = min(.x$sorensen, na.rm = TRUE),
    max_diff = max(.x$sorensen, na.rm = TRUE)
  ))
rich_limits <- c(min(limits_df$min_diff),max(limits_df$max_diff))
lims <- c(-max(abs(rich_limits)),max(abs(rich_limits)))
my_breaks <- seq(rich_limits[1], rich_limits[2], length.out = 5)

# colors 
pal <- colorRampPalette(rev(brewer.pal(9,'OrRd')))

#background 
denmark <- ne_countries(scale = "large", country = "denmark", returnclass = "sf")
# positions for annotation 
inset <- 12
df <- dfs_mapped[[2]]
x_pos <- max(df$X) - (diff(range(df$X)) / inset)
y_pos <- max(df$Y) - (diff(range(df$Y)) / inset)
# years 
year_lookup <- c("atlas_1_to_2" = "1970s to 1990s", "atlas_2_to_3" = "1990s to 2010s")

names(dfs_mapped) <- c('atlas_1_to_2','atlas_2_to_3')
name <- names(dfs_mapped)[[1]]

df <- dfs_mapped[[1]]
df
#### MAPS ####
plots_sorensen <- imap(dfs_mapped, function(df, name) {
  
  year_label <- year_lookup[[name]]
  
  ggplot() +
    # Background: The actual country shape
    geom_sf(data = denmark, fill='white', color = "grey80", size = 0.5) +
    
    # The Data Points: Using your ivory border (pch 21)
    geom_point(data = df, aes(x = X, y = Y, fill = sorensen),
               size = 1, pch = 21,color='white', alpha = 1) +
    
    # The Color Scale
    # Inside your plot loop:
    scale_fill_gradientn(
      colors = pal(10), 
      limits = rich_limits
    )+
    
    # coords 
    coord_sf(crs = st_crs(23032)) + # Forces the whole plot into UTM 32N +
    
    labs(x = NULL,
         y = NULL,
         #title = year_label,
         fill = 'Sorenson similarity') +
    
    # annotate("text", x = x_pos, y = y_pos, 
    #         label = "Atlas 1 - 1970s", 
    #         size = 8,
    #         #hjust = 1.1, vjust = 1.5,
    #         fill = "white",         # Background color
    #         alpha = 0.8) +           # Slight transparency
    
    # theme 
    theme_minimal() + 
    theme(
      legend.position = 'none',
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_blank()
    ) 
})

# ASSEMBLE ----------------------------------------------------------------
#### GET LEGENDS ####
# Extract the Richness legend from one of your map plots
# (Make sure plots[[1]] has legend.position = "right" temporarily to catch it)
leg_richness <- get_legend(
  plots_richness[[1]] + theme(legend.position = "right",
                              legend.direction = 'vertical')
)

leg_sor <- get_legend(
  plots_sorensen[[1]] + theme(legend.position = "right",
                              legend.direction = 'vertical')
)

layout <- "
ABE
CDF
"

# Combine your list of plots with the table
final_layout <- wrap_plots(
  A = plots_richness[[1]], 
  B = plots_richness[[2]], 
  C = plots_sorensen[[1]], 
  D = plots_sorensen[[2]], 
  E = leg_richness,
  F = leg_sor,
  design = layout
) + 
  plot_layout(widths = c(1, 1,0.5)) + # Make the legend row much shorter
  plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D', '', ''))) & # Don't label legends
  theme(plot.tag = element_text(face = "bold", size = 14))

final_layout

ggsave(sprintf('misc-figures/%s-sfig-richness-difference-temporal-sorenson.png',pattern),width = 7.5,height=7.5)
