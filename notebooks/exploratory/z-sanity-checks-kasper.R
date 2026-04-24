### GETTING STARTED ###
rm(list=ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, tidyverse, cowplot)

# HELPER FUNCTIONS  -------------------------------------------------------
merge_by_rownames <- function(x,y) {
  merge <- merge(x,y,by='row.names')
  rownames(merge) <- merge$Row.names
  merge[,!colnames(merge) %in% c('Row.names')]
}
# MAKING DFs FUNCTIONS  ----------------------------------------------------
make_master_df <- function(dir,pattern){
  pattern2match <- pattern
  
  matching_folders <- list.dirs(dir, recursive = FALSE, full.names = F)
  matching_folders <- matching_folders[grepl(pattern2match,
                                             basename(matching_folders))]
  model <- matching_folders[1]
  models_nums <- as.numeric(str_extract(matching_folders, "(?<=Atlas)\\d+"))
  
  outputs <- lapply(matching_folders,function(model){
    # read objects
    MF <- readRDS(file.path(dir,
                            model,
                            'Results',
                            sprintf('%sMF.rds',
                                    model)))
    VP <- readRDS(file.path(dir,
                            model,
                            'Results',
                            sprintf('%sVP.rds',
                                    model)))
    load(file.path(dir,
                   model,
                   'Models/Fitted',
                   sprintf('HPC_samples_0250_thin_100_chains_4.Rdata')))
    m <- fitted_model$posteriors
    return(m)
  })
  names(outputs) <- models_nums
  # get enviro data
  Xs <- lapply(outputs, function(df) df$XData)
  Ys <- lapply(outputs, function(df) df$Y)
  Designs <- lapply(outputs, function(df) df$studyDesign)
  ranLevels <- lapply(outputs, function(df) df$ranLevels)
  s <- lapply(ranLevels, function(df) df$s)
  i<-1
  merged_list <- lapply(seq_along(outputs), function(i) {
    
    # Get the pieces for this specific model
    X_df <- as.data.frame(Xs[[i]])
    Y_df <- as.data.frame(Ys[[i]])
    Des  <- as.data.frame(Designs[[i]])
    ses <- as.data.frame(s[[i]]$s)
    
    head(X_df[,1:5])
    head(Y_df[,1:5])
    head(Des)
    head(ses)
    # Merge them together
    # We use '0' to merge by rownames in base R
    
    # 1. Merge the Sample-level data (X, Y, and Design)
    # merge(..., by = 0) merges by row names
    combined_samples <- merge(X_df, Y_df, by = 0) %>%
      column_to_rownames("Row.names") %>%
      merge(Des, by = 0) %>%
      column_to_rownames("Row.names")
    
    # 2. Add the Spatial data (ses)
    # Since 'ses' uses site names as rownames, we merge it with the 'site' column 
    # we just got from 'Des'
    final_df <- combined_samples %>%
      rownames_to_column("sample_id") %>%
      # Merge with ses: by.x is our column 'site', by.y is rownames of 'ses'
      merge(ses, by.x = "site", by.y = 0, all.x = TRUE) %>%
      # Put the sample ID back as the rowname
      column_to_rownames("sample_id")
    
    # 3. Clean up (Optional)
    # If you have an "(Intercept)" column from X_df that you don't need:
    final_df <- final_df %>% select(-any_of("(Intercept)"))%>%
      mutate(atlas = sapply(strsplit(rownames(.), "_"), "[", 2)) %>%
      
      # 3. REORDER: Move metadata to the front
      relocate(site, atlas, X, Y) 
    
    head(final_df)
    
    return(final_df)
  })
  
  # 4. Collapse the list into one large Master Dataframe
  master_df <- bind_rows(merged_list) %>%
    # Optional: Convert atlas to a factor for better plotting later
    mutate(atlas = factor(atlas, levels = c("1", "2", "3")))
  
}
make_merged <- function(dir,pattern){
  pattern2match <- pattern
  matching_folders <- list.dirs(dir, recursive = FALSE, full.names = F)
  matching_folders <- matching_folders[grepl(pattern2match, basename(matching_folders))]
  model <- matching_folders[1]
  
  outputs <- lapply(matching_folders,function(model){
    # read objects
    MF <- readRDS(file.path(dir,
                            model,
                            'Results',
                            sprintf('%sMF.rds',
                                    model)))
    VP <- readRDS(file.path(dir,
                            model,
                            'Results',
                            sprintf('%sVP.rds',
                                    model)))
    load(file.path(dir,
                   model,
                   'Models/Fitted',
                   sprintf('HPC_samples_0250_thin_100_chains_4.Rdata')))
    m <- fitted_model$posteriors
    # get the atlas used
    X <- m$XData
    atlas <- unique(sapply(strsplit(rownames(X),'_'),'[',2))
    # get number of occurrences in that atlas
    if(!length(atlas)==1){
      stop()
    }
    occs <- data.frame(occs = colSums(m$Y,na.rm=T),
                       atlas = atlas)
    row.names(occs) <- colnames(m$Y)
    # get the trait date
    Tr <- m$TrData
    #transpose VP
    df <-VP$vals # get vals
    df <- as.data.frame(t(df)) # transpose into dataframe
    # merge VP with species info
    Tjur <- data.frame(tjur = MF$TjurR2)
    df_rel <- Tjur$tjur*df
    spNames <- m$spNames
    row.names(Tjur) <- spNames
    TrOccs <- merge_by_rownames(Tr,occs)
    TrOccsTjur <- merge_by_rownames(TrOccs,Tjur)
    TrOccsTjurVP <- merge_by_rownames(TrOccsTjur,df_rel)
    TrOccsTjurVP
  })
  merged <- do.call(rbind,outputs)
  
  # count species per guild
  guild_counts <- merged %>%
    group_by(foraging_guild_consensus) %>%
    summarise(n_species = n()/3) %>%
    mutate(Guild_counts = paste0(foraging_guild_consensus,' (n=',n_species,')')) %>%
    ungroup()
  
  #str(merged)
  merged$Guild_counts <- merged$foraging_guild_consensus
  merged$Guild_counts <- guild_counts$Guild_counts[match(merged$foraging_guild_consensus,
                                                         guild_counts$foraging_guild_consensus)]
  return(merged)
}
make_VP_effect_sizes <- function(dir,pattern){
  pattern2match <- pattern
  
  matching_folders <- list.dirs(dir, recursive = FALSE, full.names = F)
  matching_folders <- matching_folders[grepl(pattern2match, basename(matching_folders))]
  model <- matching_folders[1]
  
  outputs <- lapply(matching_folders,function(model){
    # read objects
    MF <- readRDS(file.path(dir,
                            model,
                            'Results',
                            sprintf('%sMF.rds',
                                    model)))
    VP <- readRDS(file.path(dir,
                            model,
                            'Results',
                            sprintf('%sVP.rds',
                                    model)))
    effect_sizes <- read_excel(file.path(dir,
                                         model,
                                         'Results',
                                         sprintf('%sparameter_estimates_Beta_.xlsx',model)),
                               sheet = "Posterior mean")
    effect_sigs <- read_excel(file.path(dir,
                                        model,
                                        'Results',
                                        sprintf('%sparameter_estimates_Beta_.xlsx',model)),
                              sheet = "Pr(x>0)")
    
    load(file.path(dir,
                   model,
                   'Models/Fitted',
                   sprintf('HPC_samples_0250_thin_100_chains_4.Rdata')))
    m <- fitted_model$posteriors
    # get the atlas used
    X <- m$XData
    atlas <- unique(sapply(strsplit(rownames(X),'_'),'[',2))
    # get number of occurrences in that atlas
    if(!length(atlas)==1){
      stop()
    }
    occs <- data.frame(occs = colSums(m$Y,na.rm=T),
                       atlas = atlas)
    row.names(occs) <- colnames(m$Y)
    # get the trait date
    Tr <- m$TrData
    #transpose VP
    df <-VP$vals # get vals
    df <- as.data.frame(t(df)) # transpose into dataframe
    
    # merge VP with species info
    # 1. Prepare metrics and relative df first
    metrics <- data.frame(tjur = MF$TjurR2, auc = MF$AUC, row.names = m$spNames)
    df_rel <- metrics$tjur * df
    
    # 2. The main pipe, adding the scaled VPs 
    TrOccsTjurVP <- Tr %>%
      rownames_to_column("species") %>%
      inner_join(rownames_to_column(occs, "species"), by = "species") %>%
      inner_join(rownames_to_column(metrics, "species"), by = "species") %>%
      inner_join(rownames_to_column(as.data.frame(df_rel), "species"), by = "species") %>%
      pivot_longer(
        cols = c(tmean_year, prec_year, hh, `Random: site`,
                 starts_with("perc_")),   
        names_to = "environmental_variable",
        values_to = "VP_scaled_by_Tjur"
      )
    
    # 1. Pivot Effect Sizes to long format
    long_effects <- effect_sizes %>%
      select(-`(Intercept)`) %>%
      pivot_longer(cols = -Species, names_to = "environmental_variable", values_to = "effect_size") %>%
      rename(species = Species)
    
    # 2. Pivot Significance to long format
    long_sigs <- effect_sigs %>%
      select(-`(Intercept)`) %>%
      pivot_longer(cols = -Species, names_to = "environmental_variable", values_to = "sig_val") %>%
      rename(species = Species)
    
    # 3. Join, Mask, and Merge with your main data
    final_df <- long_effects %>%
      left_join(long_sigs, by = c("species", "environmental_variable")) %>%
      # Apply the mask: Keep only if sig <= 0.05 OR sig >= 0.95
      mutate(effect_size = ifelse(sig_val > 0.05 & sig_val < 0.95, NA, effect_size)) %>%
      select(-sig_val)
    # Join to your previous master table (assuming row names are Species)
    
    combined <- TrOccsTjurVP %>% 
      left_join(final_df,by=c("species","environmental_variable"))
    
    combined
    # adding effect sizes 
    
    
  })
  merged <- do.call(rbind,outputs)
  return(merged)
  
}
# PLOTTING FUNCTIONS  -----------------------------------------------------
# Guild vs environment VP across species 
guild_environment_vp_plot <- function(dir,env_var,plot_title){
  merged <- make_merged(dir,pattern = "2026-02-10")
  
  v2inspect <- sym(env_var)
  
  guild_means <- merged %>%
    group_by(Guild_counts, atlas) %>%
    summarise(!!v2inspect := mean(!!v2inspect, na.rm = TRUE), .groups = "drop")
  colors <- c('darkblue','darkgreen','darkred')
  
  ggplot(merged,
         aes(x = reorder(Guild_counts, !!v2inspect, FUN = mean, na.rm = TRUE),
             y=!!v2inspect,
             col = atlas, fill = atlas))+
    # violins for guilds with >1 species
    #geom_violin(drop = FALSE, alpha = 0.6) +
    geom_point(alpha = 0.5,size=0.5)+
    
    # Points for each species, colored by Atlas
    # Points for mean Tjur² per guild per Atlas
    geom_point(data = guild_means,
               aes(y = !!v2inspect, fill = atlas,color=atlas),
               size = 1, shape = 21, stroke = 1) +
    
    # set colors
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    
    #scale axes
    scale_y_continuous(limits=(c(0,0.175)))+
    
    # adjustements
    labs(title=plot_title)+
    xlab(label='Guilds (number of species)')+
    ylab(label='Raw variance explained')+
    
    coord_flip()+
    theme_minimal()+
    theme(
      axis.line = element_blank(),             # remove axis lines
      #panel.grid.major = element_blank(),      # remove major gridlines
      #panel.grid.minor = element_blank()       # remove minor gridlines
    )
}
# Single species VP vs effect sizes 
species_vp_effect_plots <- function(dir,species, english_name, remove_site = T){
  df <- make_VP_effect_sizes(dir,pattern = "2026-02-10")
  subset <- df[df$species == species,]
  
  if(remove_site == T) {
    subset <- subset[subset$environmental_variable != "Random: site",]
  }
  
  base_colors <- c(
    "Temperature"     = "red",
    "Precipitation"   = "blue",
    "Land-use heterogeneity"   = "pink",
    "Urban (% coverage)"           = "grey",
    "Cropland (% coverage)"        = "yellow3", # Slightly darker so it doesn't vanish
    "Pasture (% coverage)"         = "orange",
    "Forest (% coverage)"          = "darkgreen",
    "Grass/Shrubland (% coverage)" = "lightgreen"
  )
  
  # --- STEP 1: CALCULATE SCALING ---
  # We need to map effect_size into the R2 range (0 to max VP).
  # Adjust 'coeff' if your effect sizes are very large or small.
  max_vp <- max(subset$VP_scaled_by_Tjur, na.rm = TRUE)
  max_effect <- max(abs(subset$effect_size), na.rm = TRUE)
  coeff <- max_vp / max_effect 
  
  subset_clean <- subset %>%
    mutate(
      environmental_variable = case_match(environmental_variable,
                                          "tmean_year"       ~ "Temperature",
                                          "prec_year"        ~ "Precipitation",
                                          "hh"               ~ "Land-use heterogeneity",
                                          "perc_urban"       ~ "Urban (% coverage)",
                                          "perc_cropland"    ~ "Cropland (% coverage)",
                                          "perc_pasture"     ~ "Pasture (% coverage)",
                                          "perc_forest"      ~ "Forest (% coverage)",
                                          "perc_grass_shrub" ~ "Grass/Shrubland (% coverage)",
                                          .default = environmental_variable),
      environmental_variable = factor(environmental_variable, levels = c(
        "Temperature", "Precipitation", "Land-use heterogeneity", 
        "Urban (% coverage)", "Cropland (% coverage)", 
        "Pasture (% coverage)", "Forest (% coverage)", "Grass/Shrubland (% coverage)"
      )),
      atlas = factor(atlas),
      # Scale the effect size to fit the R2 axis
      effect_scaled = effect_size * coeff
    )
  
  ggplot(data = subset_clean, aes(x = environmental_variable)) +
    # Bars for Variance (Primary Axis)
    geom_bar(aes(y = VP_scaled_by_Tjur, alpha = atlas, fill = environmental_variable), 
             stat = "identity", position = "dodge") +
    
    # Points/Lines for Effect Sizes (Secondary Axis)
    geom_point(aes(y = effect_scaled, group = atlas), 
               position = position_dodge(width = 0.9), color = "black", shape = 18, size = 3) +
    
    # --- STEP 2: DUAL AXIS ---
    scale_y_continuous(
      name = "Raw variance explained (R²)",
      labels = scales::percent,
      sec.axis = sec_axis(~ . / coeff, name = "Effect Size (Pr>0.95)")
    ) +
    
    scale_fill_manual(values = base_colors, guide = "none") +
    scale_alpha_manual(values = c("1" = 0.4, "2" = 0.7, "3" = 1.0), 
                       name = "Atlas Period",
                       labels = c("1970s", "1990s", "2010s")) +
    labs(title = english_name, subtitle = "Bars = Variance Explained | Diamonds = Effect Size", x = NULL) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          legend.position = "top")
}
# Species presences in environmental space 
species_environment_plots <- function(dir,species, english_name, env_vars) {
  
  master_df <- make_master_df(dir,pattern="2026-02-10")
  ### NOW WE GET TO PLOTTING 
  subs_names <- c('atlas', env_vars, species)
  subs <- master_df[, colnames(master_df) %in% subs_names]
  
  v1 <- sym(env_vars[1])
  v2 <- sym(env_vars[2])
  
  xlims <- master_df %>% pull(!!v1) %>% range(na.rm = TRUE)
  ylims <- master_df %>% pull(!!v2) %>% range(na.rm = TRUE)
  
  # A. Generate a "dummy" plot just to grab the shared legend
  # We make a simple ggplot with the same alpha mapping
  legend_plot <- ggplot(subs, aes(x = !!v1, y = !!v2, alpha = factor(!!sym(species)))) +
    geom_point() +
    scale_alpha_manual(values = c("0" = 0.1, "1" = 1.0), 
                       labels = c("Absence", "Presence"),
                       name = "Status") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  shared_legend <- get_legend(legend_plot)
  
  # B. Generate the 3 actual plots (with legend = FALSE)
  colors <- c('darkblue','darkgreen','darkred')
  plots <- lapply(c(1, 2, 3), function(x) {
    single_plot(x, subs, color = colors[x], env_vars, xlims, ylims, show_legend = FALSE)
  })
  
  # C. Combine into a grid
  plot_row <- plot_grid(plotlist = plots, ncol = 3, labels = "AUTO", align = 'h')
  
  # D. Add the legend at the bottom
  final_combined <- plot_grid(plot_row, shared_legend, ncol = 1, rel_heights = c(1, 0.1))
  
  # Add a super title (optional)
  title_theme <- ggdraw() + draw_label(paste(english_name), fontface='bold')
  
  return(plot_grid(title_theme, final_combined, ncol = 1, rel_heights = c(0.1, 1)))
}
single_plot <- function(atlas_nr, subs, color, env_vars, xlims, ylims, show_legend = FALSE) {
  
  subset_data <- subs[subs$atlas == atlas_nr, ]
  colnames(subset_data)[ncol(subset_data)] <- "presenceabsence"
  subset_data$presenceabsence <- factor(subset_data$presenceabsence, levels = c("0", "1"))
  
  v1 <- sym(env_vars[1])
  v2 <- sym(env_vars[2])
  
  # MAIN PLOT
  pmain <- ggplot(data = subset_data, aes(x = !!v1, y = !!v2, alpha = presenceabsence)) +
    geom_point(col = color) +
    scale_alpha_manual(values = c("0" = 0.1, "1" = 1.0), 
                       labels = c("Absence", "Presence"),
                       name = "Status") +
    scale_x_continuous(limits = xlims) + 
    scale_y_continuous(limits = ylims) +
    theme_minimal() +
    labs(title = paste0("Atlas Period: ", atlas_nr)) +
    # Toggle legend based on function call
    theme(legend.position = if(show_legend) "bottom" else "none")
  
  # MARGINAL DENSITIES
  xdens <- axis_canvas(pmain, axis = "x") +
    geom_density(data = subset_data, aes(x = !!v1), color = color, linetype = 3) +
    geom_density(data = subset_data[subset_data$presenceabsence == 1, ], aes(x = !!v1), linetype = 1, color = color)
  
  ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
    geom_density(data = subset_data, aes(x = !!v2), linetype = 3, color = color) +
    geom_density(data = subset_data[subset_data$presenceabsence == 1, ], aes(x = !!v2), linetype = 1, color = color) +
    coord_flip()
  
  # COMBINE
  p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
  
  return(p2)
}

# Species presences & environmental variables in spatial space
species_spatial_comparison <- function(dir,species, english_name, env_var,
                                       col_low,col_high,legend_name) {
  
  master_df <- make_master_df(dir,pattern="2026-02-10")
  # 1. Prep data
  # Ensure these match your master_df exactly (X/Y vs coord_x/coord_y)
  subs_cols <- c('atlas', 'X', 'Y', env_var, species) 
  subs <- master_df[, colnames(master_df) %in% subs_cols]
  
  # 2. Global limits for consistent axes and colors
  xlims <- range(subs$X, na.rm = TRUE)
  ylims <- range(subs$Y, na.rm = TRUE)
  v_env <- sym(env_var)
  env_lims <- subs %>% pull(!!v_env) %>% range(na.rm = TRUE)
  
  # 3. Create a list for each Atlas containing TOP and BOTTOM plots
  all_plots <- lapply(c(1, 2, 3), function(x) {
    
    # Filter data for this atlas
    subset_data <- subs[subs$atlas == x, ]
    colnames(subset_data)[ncol(subset_data)] <- "presenceabsence"
    
    # TOP PLOT: Species Focus (with transparency)
    p_top <- single_spatial_plot(subset_data, v_env, xlims, ylims, env_lims, 
                                 title = paste0("Atlas ", x, ": Presence"), 
                                 use_alpha = TRUE,
                                 col_low,col_high,legend_name)
    
    # BOTTOM PLOT: Climate Focus (All points opaque)
    p_bottom <- single_spatial_plot(subset_data, v_env, xlims, ylims, env_lims, 
                                    title = paste0("Atlas ", x, ": Full Gradient"), 
                                    use_alpha = FALSE,
                                    col_low,col_high,legend_name)
    
    list(top = p_top, bottom = p_bottom)
  })
  
  # 4. Extract plots from the nested list for the grid
  top_row <- lapply(all_plots, `[[`, "top")
  bottom_row <- lapply(all_plots, `[[`, "bottom")
  
  # Combine into a 2-row grid
  final_grid <- plot_grid(
    plot_grid(plotlist = top_row, ncol = 3, align = 'h'),
    plot_grid(plotlist = bottom_row, ncol = 3, align = 'h'),
    ncol = 1,
    rel_heights = c(1, 1)
  )
  
  dummy_plot <- ggplot(master_df, aes(x=X, y=Y, color=!!v_env)) + 
    geom_point() +
    scale_color_gradient(low = col_low, high = col_high, name = legend_name)
  
  leg <- get_legend(dummy_plot + theme(legend.position = "bottom"))
  # Final assembly
  # 1. Create the title object
  title <- ggdraw() + 
    draw_label(
      english_name,
      fontface = 'bold',
      x = 0.5, # Center horizontally
      hjust = 0.5,
      size = 18
    )
  
  # 2. Stack everything together
  layout <- plot_grid(
    title, 
    final_grid, 
    leg, 
    ncol = 1, 
    # Adjust rel_heights so the title is small and the plot is large
    rel_heights = c(0.1, 1, 0.1) 
  )  
  return(layout)
}
single_spatial_plot <- function(subset_data, v_env, xlims, ylims, env_lims, title, use_alpha,
                                col_low, col_high, legend_name) {
  
  # Sort so presences are on top if using alpha
  if(use_alpha) {
    subset_data <- subset_data %>% arrange(presenceabsence)
    alpha_scale <- scale_alpha_manual(values = c("0" = 0.05, "1" = 1.0), guide = "none")
  } else {
    alpha_scale <- scale_alpha_manual(values = c("0" = 1.0, "1" = 1.0), guide = "none")
  }
  
  plot <- ggplot(subset_data, aes(x = X, y = Y)) +
    geom_point(aes(color = !!v_env, alpha = as.factor(presenceabsence)), size = 0.8) +
    alpha_scale +
    scale_color_gradient(limits = env_lims, 
                         low = col_low,
                         high = col_high) +
    coord_fixed(xlim = xlims, ylim = ylims) +
    theme_minimal() +
    theme(
      axis.text = element_blank(), axis.title = element_blank(),
      panel.grid = element_blank(), legend.position = "none"
    ) +
    labs(subtitle = title) 
  
}
# PLOTTING ----------------------------------------------------------------
### Set directory to end up in 'HmscOutputs' folder. 
path2HmscOutputs <- './HmscOutputs'

## environmental variables 
# env_var is one of tmean_year,prec_year,hh, perc_urban, perc_cropland,perc_forest,perc_pasture, perc_grass_shrub
# tmean_year is temperature
# prec_year is precipitation 
# hh is land-use heterogeneity
# perc_ are the % coverages of different land-use types 

## species 
# the species are all formatted as "Picus_viridis"
# to get an overview of all species and environmental variables run 
test <- make_master_df(dir=path2HmscOutputs,pattern="2026-02-10")
colnames(test)


### LOOK INTO GUILDS / ENVIRONMENTS BASED ON INDIVIDUAL SPECIES VARIANCE EXPLAINED BY ENVIRONMENT VARIABLES 
# this only takes an environmental variables 
guild_environment_vp_plot(dir = path2HmscOutputs,
                          env_var = 'perc_forest',
                          plot_title = 'Forest coverage %')
### LOOK INTO VARIANCE EXPLAINED BY DIFFERENT VARIABLES AND ESTIMATES OF EFFECT SIZES ACROSS PERIODS FOR A SINGLE SPECIES 
species_vp_effect_plots(dir = path2HmscOutputs,
                        species = 'Picus_viridis',
                        english_name = 'Green woodpecker',
                        remove_site=T)
# you'll get some errors with this one, that's correct, there are some insignificant effect sizes plotted as NAs and that 
# throws and error

### LOOK AT ENVIRONMENTAL SPACE OCCUPIED BY A SPECIES 
# here you can choose two variables to look at that environmental space 
species_environment_plots(dir = path2HmscOutputs,
                          species='Picus_viridis',
                          english_name='Green woodpecker',
                          env_vars=c('perc_forest','perc_urban'))
### LOOK AT THIS OCCUPIED ENVIRONMENTAL SPACE ON A MAP
# you you choose two colors, i like purple yellow for the stark contrast 
species_spatial_comparison(dir = path2HmscOutputs,
                           species='Picus_viridis',
                           english_name='Green woodpecker', 
                           env_var='perc_forest',
                           col_low = 'purple',col_high = 'yellow',
                           legend_name = 'Forest coverage %')






 
