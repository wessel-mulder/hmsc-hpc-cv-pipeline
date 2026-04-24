rm(list=ls())
#### GETTING STARTED ####
if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,Hmsc,RColorBrewer,ggplot2,
               rnaturalearth,rnaturalearthdata,
               gridExtra,patchwork,sf,cowplot,parallel)

#### LOAD MODELS #### 
# GET PATTERNS + MODEL NAMES 
dir <- './HmscOutputs'
pattern <- '2026-03-13'
matching_folders <- list.dirs(dir, recursive = FALSE, full.names = F)
matching_folders <- matching_folders[grepl(pattern,
                                           basename(matching_folders))]
model <- matching_folders[1]
models_nums <- as.numeric(str_extract(matching_folders, "(?<=Atlas)\\d+"))

# LOAD MODELS
mods <- lapply(matching_folders,function(model){
  load(file.path(dir,model,'Models','Fitted','HPC_samples_0250_thin_100_chains_4.Rdata'))
  m <- fitted_model$posteriors
})
names(mods) <- models_nums

# LOAD DESIGNS
designs <- lapply(mods,function(mod){
  mod$studyDesign %>% 
    rownames_to_column(.,var='survey') %>% 
    left_join(.,
              mod$ranLevels$site$s %>% rownames_to_column(.,var='site'),
              by = 'site')
})

# EMPIRICAL 
mod <- mods[[1]]
emps <- lapply(mods,function(mod){
  mod$Y
})

#### SIMULATIONS ####
emp_matrix <- emps[[1]]
simulate_prob_null <- function(emp_matrix, max_iter = 1000, tol = 1e-5) {
  
  # 1. Extract the empirical constraints
  target_richness <- rowSums(emp_matrix)
  target_prevalence <- colSums(emp_matrix)
  
  # shuffle starting
  random_prior <- matrix(runif(length(emp_matrix)), 
                         nrow = nrow(emp_matrix), 
                         ncol = ncol(emp_matrix))
  
  sim_mat <- random_prior
  # 3. Iterative Proportional Fitting (IPF) loop
  for (i in 1:max_iter) {
    mat_old <- sim_mat
    
    # Step A: Scale Rows to match Site Richness
    r_sums <- rowSums(sim_mat)
    # Avoid division by zero for empty sites
    r_multiplier <- target_richness / ifelse(r_sums == 0, 1, r_sums)
    sim_mat <- sim_mat * r_multiplier
    
    # Cap probabilities at 1.0 (ecologically necessary)
    sim_mat[sim_mat > 1] <- 1
    
    # Step B: Scale Columns to match Species Prevalence
    c_sums <- colSums(sim_mat)
    c_multiplier <- target_prevalence / ifelse(c_sums == 0, 1, c_sums)
    sim_mat <- sweep(sim_mat, 2, c_multiplier, FUN = "*")
    
    # Cap probabilities at 1.0 again
    sim_mat[sim_mat > 1] <- 1
    
    # 4. Check for convergence (Has the matrix stopped changing?)
    if (max(abs(sim_mat - mat_old)) < tol) {
      message(sprintf("Algorithm converged perfectly in %d iterations.", i))
      break
    }
  }
  
  # Clean up row and column names to match the input
  rownames(sim_mat) <- rownames(emp_matrix)
  colnames(sim_mat) <- colnames(emp_matrix)
  
  return(sim_mat)
}

#### SORENSON FUNCTION ####
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

#### GET SORENSEN INDEX FOR SIMULATIONS ####
n_sims <- 100
simu_sorensons <- lapply(emps,function(emp.df){
  
  lapply(1:n_sims,function(id){
  # get row.names 
  row.names <- rownames(emp.df)
  # get simulated community 
  simu.df <- simulate_prob_null(emp.df, max_iter = 1000, tol = 1e-5)
  
  row.name <- row.names[[1]]
  # pred.df <- predsY[[1]]
  # emp.df <- emps[[1]]
  
  # get values 
  vals <- lapply(row.names,function(row.name){
    pred.vec <- simu.df[row.name,]
    emp.vec <- emp.df[row.name,]
    probsorensen(pred.vec,emp.vec)
  })
  
  data.frame(survey = row.names,
             sorensen = unlist(vals))
  })
})

# average the simulations 
simu_sorensons_means <- lapply(simu_sorensons, function(sub_list) {
  # Combine dataframes within the sub-list and calculate mean
  sub_list %>%
    bind_rows() %>%
    group_by(survey) %>%
    summarize(sorensen = mean(sorensen, na.rm = TRUE))
})

#### PREPARE MAPPING ####
dfs <- map2(simu_sorensons_means,designs,function(sor,design){
  left_join(sor,design,by='survey')
})

#### MAP ####
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
# limits_df <- dfs %>%
#   imap_dfr(~ data.frame(
#     name = .y,
#     min = min(.x$sorensen, na.rm = TRUE),
#     max = max(.x$sorensen, na.rm = TRUE)
#   ))
# rich_limits <- c(min(limits_df$min),max(limits_df$max))

# copy limits from empirical
rich_limits <- c(0.3912802,0.6579209)

#### MAKE FAKE ATLAS 3 
dfs[[3]] <- dfs[[2]]
dfs[[3]]$sorensen <- rich_limits[1]
names(dfs) <- c("1", "2", "3")

#### MAPS ####
plots <- imap(dfs, function(df, name) {
  
  year_label <- year_lookup[[name]]
  
  mean_sor <- mean(df$sorensen)
  sd_sor <- sd(df$sorensen)
  
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
      limits = rich_limits,
      oob = scales::squish
    )+
    
    # coords 
    coord_sf(crs = st_crs(23032)) + # Forces the whole plot into UTM 32N +
    
    labs(x = NULL,
         y = NULL,
         title = year_label,
         fill = 'Sørenson index') +
    
    annotate("text", x = x_pos, y = y_pos,
            label = paste('μ =',round(mean_sor,2),'\n',
                          'σ =',round(sd_sor,2)),
            size = 4,
            #hjust = 1.1, vjust = 1.5,
            alpha = 0.8) +           # Slight transparency
    
    # theme 
    theme_minimal() + 
    theme(
      legend.position = 'none',
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_blank()
    ) 
})
plots[[1]]

#### LAYOUT ####
leg_richness <- get_legend(
  plots[[1]] + theme(legend.position = "bottom")
)



layout <- "
A
B
C
D
"

# Combine your list of plots with the table
final_layout <- wrap_plots(
  A = plots[[1]], 
  B = plots[[2]], 
  C = plots[[3]], 
  D = leg_richness,
  design = layout
) + 
  plot_layout(heights = c(1, 1, 1, 0.05)) + # Make the legend row much shorter
  plot_annotation(tag_levels = list(c('A', 'B', 'C',''))) & # Don't label legends
  theme(plot.tag = element_text(face = "bold", size = 14))

final_layout

ggsave(sprintf('misc-figures/%s-sfig3-sorensen-similarity-simu-nullmodels.png',pattern),final_layout,width=5,height=10)


