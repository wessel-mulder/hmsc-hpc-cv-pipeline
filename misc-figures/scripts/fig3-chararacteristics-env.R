rm(list=ls())
#### GETTING STARTED ####
if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,Hmsc,RColorBrewer,ggplot2,
               rnaturalearth,rnaturalearthdata,
               gridExtra,patchwork,
               readxl,cowplot,abind,
               ggrepel)
source(file.path("support_scripts", "figure_data_helpers.R"))

devtools::install_github("davidsjoberg/ggbump")
devtools::install_github("davidsjoberg/ggsankey")
library(ggbump)
library(ggsankey)

scaled = F # set to NULL to not scale 

#### LOAD VP DATA ####
pattern <- "2026-03-13"
matching_folders <- figure_model_folders(pattern = pattern)
model <- matching_folders[1]
models_nums <- atlas_numbers(matching_folders)
preds <- load_gradient_predictions(matching_folders)
mods <- load_hmsc_posteriors(matching_folders)

#### HELPER FUNCTIONS ####
# Helper to extract category from trName
get_category <- function(trName) {
  dplyr::case_when(
    trName == "(Intercept)"                        ~ "Intercept",
    grepl("^Migration_a3_DOF", trName)             ~ "Migration_a3_DOF",
    grepl("^foraging_guild_consensus", trName)     ~ "foraging_guild_consensus",
    grepl("^species_thermal_index", trName)        ~ "species_thermal_index",
    TRUE                                           ~ trName  # fallback
  )
}

# Helper to extract the level label from trName
get_label <- function(trName) {
  dplyr::case_when(
    trName == "(Intercept)"                        ~ "Intercept",
    grepl("^Migration_a3_DOF", trName)             ~ sub("^Migration_a3_DOF", "", trName),
    grepl("^foraging_guild_consensus", trName)     ~ sub("^foraging_guild_consensus", "", trName),
    TRUE                                           ~ trName
  )
}

# modify source code to just get Pr
getPr <- function(hM, Gradient, predY, measure, index = 1) {
  
  if (!measure %in% c("S", "Y", "T")) {
    stop(paste("measure must be one of 'S', 'Y', 'T', got:", measure))
  }
  
  switch(class(hM$X)[1L],
         matrix = { xx = Gradient$XDataNew[, 1] },
         list   = {
           if (measure == "Y") {
             xx = Gradient$XDataNew[[index]][, 1]
           } else {
             xx = Gradient$XDataNew[[1]][, 1]
           }
         }
  )
  
  ngrid <- length(xx)
  
  Pr <- if (measure == "S") {
    predS <- abind(lapply(predY, rowSums), along = 2)
    mean(predS[ngrid, ] > predS[1, ])
    
  } else if (measure == "Y") {
    tmp <- abind(predY, along = 3)
    mean(tmp[ngrid, index, ] > tmp[1, index, ])
    
  } else if (measure == "T") {
    predT <- if (all(hM$distr[, 1] == 1)) {
      lapply(predY, function(a) (exp(a) %*% hM$Tr) / matrix(rep(rowSums(exp(a)), hM$nt), ncol = hM$nt))
    } else {
      lapply(predY, function(a) (a %*% hM$Tr) / matrix(rep(rowSums(a), hM$nt), ncol = hM$nt))
    }
    predT <- abind(predT, along = 3)
    mean(predT[ngrid, index, ] > predT[1, index, ])
  }
  Pr

}

#### EXTRACT PREDICTION INFO #####
pred <- preds[[1]]
mod <- mods[[1]]
k <- 1
df <- lapply(seq_along(preds),function(k){

    # set up my objects
  atlas <- names(preds)[[k]]
  pred <- preds[[k]]
  m <- mods[[k]]
  env_vars <- names(pred)
  
  # loop over traits and response 
  responses <- c('total','marginal')
  
  #### GET THE RICHNESS
  response <- responses[[1]]
  pls_richness <- lapply(responses,function(response){
    prediction <- ifelse(response == "total", "predY", "predY2")
    pls <- lapply(1:length(pred),function(k){
      pl = getPr(m, pred[[k]]$Gradient, pred=pred[[k]][[prediction]], measure="S")
      # if(inherits(pl, "ggplot")){
      #   print(pl + labs(title=paste0(modelnames,": community weighted mean trait (total effect)")))
      # }
      pl
    })
    category <- "community"
    label <- "richness"
    df <- data.frame(atlas = atlas,
                     category = category,
                     label = label,
                     variable = names(pred),
                     pr = unlist(pls),
                     response = response)
    
  }) %>% do.call(rbind,.)
  
  #### GET THE TRAITS
  response <- responses[[1]]
  pls_traits <- lapply(responses,function(response){ ### GET THE DIFFERENT RESPONSE
    prediction <- ifelse(response == "total", "predY", "predY2")
    index <- 2
    dfs_traits <- lapply(1:m$nt,function(index){ ### GET THE DIFFERENT GUILDS 
      
      pls <- lapply(1:length(pred),function(k){
        pl = getPr(m, pred[[k]]$Gradient, pred=pred[[k]][[prediction]], measure="T",index = index)
      })
      # GET CATEGORY AND LABEL INFO 
      name <- m$trNames[[index]]
      category <- get_category(name)
      label <- get_label(name)
      data.frame(atlas = atlas,
                       category = category,
                       label = label,
                       variable = names(pred),
                       pr = unlist(pls),
                       response = response)
    }) %>% do.call(rbind,.)
    
  }) %>% do.call(rbind,.)
  
  df <- rbind(pls_richness,pls_traits)
  
}) %>% do.call(rbind,.) %>% 
  mutate(support = case_when(
    pr > 0.95 ~ 'higher',
    pr < 0.05 ~ 'lower',
    TRUE ~ NA_character_
  ))

##### RICHNESS FOREST #####
trait <- 'species_thermal_index'
responses <- 'total'

subset <- df %>% 
  filter(category == trait,
         response == responses) %>%
  # Optional: Calculate a mean pr per label to ensure the reordering is consistent
  # if multiple atlases have different values for the same label.
  mutate(variable = reorder(variable, pr)) %>% 
  mutate(
    variable = case_match(variable,
                          "tmean_breeding"     ~ "Temperature",
                          "prec_breeding"      ~ "Precipitation",
                          "hh"                 ~ "Land-use heterogeneity",
                          "perc_urban"         ~ "Urban",
                          "perc_cropland"      ~ "Cropland",
                          "perc_pasture"       ~ "Pasture",
                          "perc_forest"        ~ "Forest",
                          "perc_grass_shrub"   ~ "Grass/Shrubland",
                          .default = variable
    ))
  

ggplot(subset, aes(x = factor(atlas, labels = c("1970s","1990s","2010s")), 
                   y = variable, 
                   fill = pr - 0.5)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low      = "firebrick",
    mid      = "white", 
    high     = "forestgreen",
    midpoint = 0,
    limits   = c(-0.5, 0.5),
    name     = "Support strength\n(centered probability)"
  ) +
  labs(
    x = "Atlas period",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.y     = element_text(size = 10),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    panel.grid      = element_blank()
  )

#### PCA OF TRAITS ##### --------------------------------------------------------------------
trait <- 'foraging_guild_consensus'
axes <- 'pca12'
n_species <- '5'
responses <- 'total'

m <- mods[[1]]
guild_counts <- m$TrData[trait] %>%
  # Pivot to long format to count species per trait label
  pivot_longer(everything(), names_to = "trait_column", values_to = "label") %>%
  # Filter out NAs or empty traits if necessary
  filter(!is.na(label)) %>%
  group_by(label) %>%
  summarise(n_species = n(), .groups = "drop")

target_guilds <- guild_counts$label[guild_counts$n_species>=n_species]

# 1. Prepare and Filter Data
subset_wide <- df %>% 
  filter(response == responses, 
         category == trait,
         label %in% target_guilds) %>% # <--- FILTERING HERE
  select(atlas, label, variable, pr) %>% 
  mutate(
    variable = case_match(variable,
                          "tmean_breeding"     ~ "Temperature",
                          "prec_breeding"      ~ "Precipitation",
                          "hh"                 ~ "Land-use heterogeneity",
                          "perc_urban"         ~ "Urban",
                          "perc_cropland"      ~ "Cropland",
                          "perc_pasture"       ~ "Pasture",
                          "perc_forest"        ~ "Forest",
                          "perc_grass_shrub"   ~ "Grass/Shrubland",
                          .default = variable
    )) %>% 
  pivot_wider(names_from = variable, values_from = pr)

# 2. Run PCA
pca_columns <- subset_wide %>% select(-atlas, -label)
pca_res <- prcomp(pca_columns, scale. = TRUE)

# 3. Extract Scores (Points) and Loadings (Arrows)
pca_scores <- as.data.frame(pca_res$x) %>%
  bind_cols(subset_wide %>% select(atlas, label))

# 1. PREPARE ARROWS
arrow_scale <- 1
pca_arrows <- as.data.frame(pca_res$rotation) %>%
  rownames_to_column("variable")

pca_arrows <- pca_arrows %>%
  mutate(x_coord = PC1, y_coord = PC2)

# Create the end points for the segments
pca_arrows <- pca_arrows %>%
  mutate(x_end = x_coord * arrow_scale, 
         y_end = y_coord * arrow_scale)

# 2. CLEAN ARROWS
pca_arrows_clean <- pca_arrows %>%
  filter(abs(x_coord) > 0.2 | abs(y_coord) > 0.2)

# 3. CALCULATE CENTROIDS
pca_centroids <- pca_scores %>%
  mutate(x_val = PC1, y_val = PC2) %>%
  group_by(label) %>%
  summarise(
    x_center = mean(x_val),
    y_center = mean(y_val),
    .groups = "drop"
  )

# --- CRITICAL ADDITION: JOIN CENTROIDS TO SCORES ---
# This allows geom_segment to see both the point (PC1) and its center (x_center)
pca_scores_joined <- pca_scores %>%
  left_join(pca_centroids, by = "label")

# 4. THE PLOT
ggplot(pca_scores_joined, aes(x = PC1, y = PC2)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey90") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey90") +
  
  # 1. Spider/Star Lines
  # Draws a line from every point to its group mean
  geom_segment(aes(xend = x_center, yend = y_center, color = label), 
               alpha = 0.2, linewidth = 1) +
  
  #2. Arrows (Variables/Loadings)
  geom_segment(data = pca_arrows_clean,
               aes(x = 0, y = 0, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.05, "cm")),
               color = "firebrick", alpha = 0.6, linewidth = 0.3) +

  geom_text_repel(data = pca_arrows_clean,
                  aes(x = x_end, y = y_end, label = variable),
                  color = "firebrick", fontface = "italic", size = 3,
                  box.padding = 0.3) +
  
  # 3. Points (Optional but recommended to see the 'end' of the lines)
  geom_point(aes(color = label, shape = atlas), size = 1.5, alpha = 0.6) +
  
  # 4. Centroid Labels (Placed exactly over the centerpoint)
  # Using shadowtext or bg.color from repel is better for visibility
  geom_text_repel(data = pca_centroids, 
                  aes(x = x_center, y = y_center, label = label, color = label),
                  fontface = "bold", size = 4,
                  bg.color = "white", bg.r = 0.1,
                  box.padding = 0) +
  
  theme_minimal() +
  labs(title = "Temporal stationarity of guild-responses to the environment",
       subtitle = "Lines connect individual survey results to the guild centroid") +
  theme(legend.position = 'right') +
  scale_color_viridis_d()

# 5. SAVE
ggsave(sprintf('misc-figures/%s-fig3-guilds-prob-of-presence-species-trheshold-%s.png', pattern, n_species), 
       width = 7.5, height = 7.5)



#### TRY A TILE PLOT ####
responses <- 'total'
trait <- 'foraging_guild_consensus'

subset_average <- df %>%
  # 1. Standardize Variable Names
  mutate(variable = case_match(variable,
                               "tmean_breeding"   ~ "Temperature",
                               "prec_breeding"    ~ "Precipitation",
                               "hh"               ~ "Land-use heterogeneity",
                               "perc_urban"       ~ "Urban",
                               "perc_cropland"    ~ "Cropland",
                               "perc_pasture"     ~ "Pasture",
                               "perc_forest"      ~ "Forest",
                               "perc_grass_shrub" ~ "Grass/Shrubland",
                               .default = variable)) %>% 
  
  # 2. Filter for specific JSDM outputs
  filter(response == responses, category == trait) %>%
  
  # 3. Transform PR (centering at 0)
  mutate(pr = pr - 0.5) %>%
  
  # 4. Handle NAs before averaging
  mutate(pr = replace_na(pr, 0)) %>%
  
  # 5. Group and Collapse
  group_by(label, variable) %>% 
  summarise(
    mean_pr = mean(pr, na.rm = TRUE),
    sd_pr   = sd(pr, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  filter(label %in% target_guilds)

mean_mat <- subset_average %>% 
  select(-sd_pr) %>% 
  pivot_wider(names_from = variable, values_from = mean_pr) %>% 
  column_to_rownames('label')

sd_mat <- subset_average %>% 
  select(-mean_pr) %>% 
  pivot_wider(names_from = variable, values_from = sd_pr) %>% 
  column_to_rownames('label') 

### MAKE MEAN MAT 
p1 <- pheatmap::pheatmap(mean_mat)  # silent = TRUE suppresses the plot
# extract row and column orders from the first clustering
row_order <- rownames(mean_mat)[p1$tree_row$order]
col_order <- colnames(mean_mat)[p1$tree_col$order]
### CHANGE ROW ORDER SD MAT 
sd_mat <- sd_mat[row_order, col_order]
p2 <- pheatmap::pheatmap(sd_mat,
                         cluster_rows=F,
                         cluster_cols=F)  # silent = TRUE suppresses the plot
### CLUSTER
# --- cluster on mean matrix ---
row_clust <- hclust(dist(mean_mat),      method = "ward.D2")
col_clust <- hclust(dist(t(mean_mat)),   method = "ward.D2")

row_order <- row_clust$labels[row_clust$order]
col_order <- col_clust$labels[col_clust$order]

# --- reshape with fixed order ---
mean_long <- mean_mat[row_order, col_order] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("guild") %>%
  pivot_longer(-guild, names_to = "variable", values_to = "mean") %>%
  mutate(
    guild    = factor(guild,    levels = rev(row_order)),
    variable = factor(variable, levels = col_order)
  )

sd_long <- sd_mat[row_order, col_order] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("guild") %>%
  pivot_longer(-guild, names_to = "variable", values_to = "sd")

library(scales)
plot_df <- left_join(mean_long, sd_long, by = c("guild", "variable")) %>%
  mutate(
    alpha_val = 1 - rescale(sd, to = c(0.2, 0.8)),
    guild    = factor(guild,    levels = rev(row_order)),
    variable = factor(variable, levels = col_order)
  )


ggplot(plot_df, aes(x = variable, y = guild)) +
  geom_tile(aes(fill = mean, alpha = alpha_val), color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "#D85A30", mid = "goldenrod", high = "#1D9E75",
    midpoint = 0, name = "mean"
  ) +
  scale_alpha_identity(name = "certainty\n(1 - scaled SD)",
                       guide = guide_legend()) +
  scale_x_discrete(position='top') +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 0, vjust = 0),
    panel.grid   = element_blank(),
    legend.position = "right"
  )

ggsave(sprintf('misc-figures/%s-fig3-guilds-prob-of-presence-species-trheshold-%s-tile.png', pattern, n_species), 
       width = 7.5, height = 7.5)


# VP TILE PLOT ------------------------------------------------------------



# test --------------------------------------------------------------------

var_rename <- function(df) {
  df %>% mutate(
    variable = case_match(variable,
                          "tmean_breeding"   ~ "Temperature",
                          "prec_breeding"    ~ "Precipitation",
                          "hh"               ~ "Land-use heterogeneity",
                          "perc_urban"       ~ "Urban",
                          "perc_cropland"    ~ "Cropland",
                          "perc_pasture"     ~ "Pasture",
                          "perc_forest"      ~ "Forest",
                          "perc_grass_shrub" ~ "Grass/Shrubland",
                          .default = variable
    ))
}

atlas_labels <- c("1" = "1970s", "2" = "1990s", "3" = "2010s")

# =========================================================
# HELPER: tile plot
# =========================================================
make_tile <- function(trait, responses = 'total', 
                      low = "blue4", high = "red2",
                      legend_title = "Support strength\n(centered probability)") {
  df %>%
    filter(category == trait, response == responses) %>%
    mutate(variable = reorder(variable, pr)) %>%
    var_rename() %>%
    ggplot(aes(x = factor(atlas, labels = c("1970s","1990s","2010s")),
               y = variable,
               fill = pr - 0.5)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = low, mid = "white", high = high,
      midpoint = 0, limits = c(-0.5, 0.5),
      name = legend_title
    ) +
    labs(x = "Atlas period", y = NULL) +
    theme_minimal() +
    theme(
      axis.text.y      = element_text(size = 10),
      legend.position  = "bottom",
      legend.key.width = unit(2, "cm"),
      panel.grid       = element_blank()
    )
}

# =========================================================
# HELPER: PCA plot with polygons
# =========================================================
make_pca <- function(trait, n_species_threshold = 5, responses = 'total',
                     arrow_scale = 1, arrow_threshold = 0.2) {
  
  m <- mods[[1]]
  guild_counts <- m$TrData[trait] %>%
    pivot_longer(everything(), names_to = "trait_column", values_to = "label") %>%
    filter(!is.na(label)) %>%
    group_by(label) %>%
    summarise(n_species = n(), .groups = "drop")
  
  target_guilds <- guild_counts$label[guild_counts$n_species >= n_species_threshold]
  
  subset_wide <- df %>%
    filter(response == responses,
           category == trait,
           label %in% target_guilds) %>%
    select(atlas, label, variable, pr) %>%
    var_rename() %>%
    pivot_wider(names_from = variable, values_from = pr)
  
  pca_res    <- prcomp(subset_wide %>% select(-atlas, -label), scale. = TRUE)
  pca_scores <- as.data.frame(pca_res$x) %>%
    bind_cols(subset_wide %>% select(atlas, label)) %>%
    mutate(atlas = factor(atlas, levels = c("1","2","3"), labels = c("1970s","1990s","2010s")))
  
  pca_arrows <- as.data.frame(pca_res$rotation) %>%
    rownames_to_column("variable") %>%
    mutate(x_end = PC1 * arrow_scale, y_end = PC2 * arrow_scale) %>%
    filter(abs(PC1) > arrow_threshold | abs(PC2) > arrow_threshold)
  
  pca_centroids <- pca_scores %>%
    group_by(label) %>%
    summarise(x_center = mean(PC1), y_center = mean(PC2), .groups = "drop")
  
  # Variance explained for axis labels
  var_exp <- summary(pca_res)$importance[2,] * 100
  
  ggplot(pca_scores, aes(x = PC1, y = PC2)) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey80") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey80") +
    
    # Polygons connecting the 3 atlas period points per guild
    geom_polygon(aes(fill = label, group = label), alpha = 0.15) +
    #geom_polygon(aes(color = label, group = label), fill = NA, linewidth = 0.4) +
    
    # Arrows
    geom_segment(data = pca_arrows,
                 aes(x = 0, y = 0, xend = x_end, yend = y_end),
                 arrow = arrow(length = unit(0.15, "cm")),
                 color = "firebrick", alpha = 0.7, linewidth = 0.4) +
    geom_text_repel(data = pca_arrows,
                    aes(x = x_end, y = y_end, label = variable),
                    color = "firebrick", fontface = "italic", size = 3,
                    box.padding = 0.3) +
    
    # Points
    #geom_point(aes(color = label, shape = atlas), size = 1.5, alpha = 0.8) +
    
    # Centroid labels
    geom_text_repel(data = pca_centroids,
                    aes(x = x_center, y = y_center, label = label, color = label),
                    fontface = "bold", size = 4,
                    bg.color = "white", bg.r = 0.1,
                    box.padding = 0) +
    
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    labs(
      x     = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y     = paste0("PC2 (", round(var_exp[2], 1), "%)"),
      shape = "Atlas period",
      color = NULL,
      fill  = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}

# =========================================================
# MAKE PLOTS
# =========================================================

# =========================================================
# ASSEMBLE
# =========================================================
layout <- "
AB
"

final <- wrap_plots(
  A = tile_richness,
  B = tile_sti
) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = list(c('A','B'))) &
  theme(plot.tag = element_text(face = "bold", size = 14))

final

ggsave(sprintf('misc-figures/%s-fig3-traits-summary.png', pattern),
       width = 14, height = 12)

# Tile plots
tile_richness <- make_tile(
  trait = 'community',       # adjust to match your df category name
  low   = "forestgreen",
  high  = "firebrick"
)

tile_sti <- make_tile(
  trait = 'species_thermal_index',
  low   = "firebrick",
  high  = "forestgreen"
)

# PCA plots
pca_foraging   <- make_pca(trait = 'foraging_guild_consensus',arrow_scale = 3)
pca_migratory  <- make_pca(trait = 'Migration_a3_DOF',arrow_scale = 3)        # adjust to match your trait name

