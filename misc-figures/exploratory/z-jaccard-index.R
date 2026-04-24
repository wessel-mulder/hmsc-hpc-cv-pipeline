# ==============================================================================
# z-spatial-predictions.R
# Spatial predictions, threshold optimization, and Null Model validation 
# for HMSC community composition.
# ==============================================================================

# 1. SETUP AND PACKAGE LOADING -------------------------------------------------
if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, Hmsc, PresenceAbsence, vegan)

# Global configuration variables
dir <- './HmscOutputs'
pattern <- '2026-03-13'
methods <- c("MaxSens+Spec", "Sens=Spec", "ObsPrev")
n_simulations <- 999 # Number of null model simulations

# 2. HELPER FUNCTIONS ----------------------------------------------------------

clean_name <- function(x) {
  x <- gsub("\\+", "_plus_", x)
  x <- gsub("=", "_equal_", x)
  return(tolower(x))
}

# Calculates the mean Jaccard similarity across all sites
calc_jaccard_mean <- function(m1, m2) {
  intersect <- rowSums((m1 == 1) & (m2 == 1))
  union <- rowSums((m1 == 1) | (m2 == 1))
  return(mean(ifelse(union == 0, 1, intersect / union)))
}

# 3. PREDICTIONS AND THRESHOLDS ------------------------------------------------

get_preds <- function(dir, pattern) {
  matching_folders <- list.dirs(dir, recursive = FALSE, full.names = FALSE)
  matching_folders <- matching_folders[grepl(pattern, basename(matching_folders))]
  
  lapply(matching_folders, function(model) {
    cat("Generating predictions for:", model, "\n")
    load(file.path(dir, model, 'Models', 'Fitted', 'HPC_samples_0250_thin_100_chains_4.Rdata'))
    m <- fitted_model$posteriors
    
    # Note: corrected to standard Hmsc computePredictedValues
    preds <- computePredictedValues(m) 
    
    EpredY <- rowMeans(preds, dims = 2)
    empirical <- m$Y
    comparison_list <- list(predictions = EpredY, empirical = empirical)
    
    test_dir <- file.path(dir, model, 'Tests')
    if(!dir.exists(test_dir)) dir.create(test_dir)
    
    saveRDS(comparison_list, file = file.path(test_dir, 'comparisonlist.RData'))
  })
}

get_species_thresholds <- function(preds, obs, method = "MaxSens+Spec") {
  thresh_vector <- setNames(numeric(ncol(preds)), colnames(preds))
  for (sp in colnames(preds)) {
    data_prep <- data.frame(
      ID = 1:nrow(preds), 
      obs = obs[, sp], 
      pred = preds[, sp]
    )
    opt <- optimal.thresholds(data_prep, opt.methods = method)
    thresh_vector[sp] <- opt[1, 2]
  }
  return(thresh_vector)
}

get_thresholds <- function(dir, pattern, methods) {
  matching_folders <- list.dirs(dir, recursive = FALSE, full.names = FALSE)
  matching_folders <- matching_folders[grepl(pattern, basename(matching_folders))]
  
  lapply(matching_folders, function(model) {
    cat("Calculating thresholds for:", model, "\n")
    comparison_list <- readRDS(file.path(dir, model, 'Tests', 'comparisonlist.RData'))
    
    vects <- lapply(methods, function(method) {
      get_species_thresholds(
        preds = comparison_list$predictions,
        obs = comparison_list$empirical,
        method = method
      )
    })
    
    names(vects) <- methods
    saveRDS(vects, file = file.path(dir, model, 'Tests', 'optimized-species-thresholds.Rdata'))
  })
}

# 4. NULL MODEL VALIDATION -----------------------------------------------------

evaluate_jaccard_null_models_fast <- function(dir, pattern, methods, n_sim = 999) {
  matching_folders <- list.dirs(dir, recursive = FALSE, full.names = FALSE)
  matching_folders <- matching_folders[grepl(pattern, basename(matching_folders))]
  
  all_results <- list()
  
  for (model_folder in matching_folders) {
    cat("\n--------------------------------------------------\n")
    cat("Processing Model:", model_folder, "\n")
    
    comparison <- readRDS(file.path(dir, model_folder, 'Tests', 'comparisonlist.RData'))
    thresholds <- readRDS(file.path(dir, model_folder, 'Tests', 'optimized-species-thresholds.Rdata'))
    
    emp_mat <- comparison$empirical
    pred_mat <- comparison$predictions
    
    # --- STEP 1: PARALLEL NULL GENERATION ---
    # Generating quasiswap nulls is the slowest part. 
    # We split the n_sim into smaller chunks to run on different cores.
    cat("Generating", n_sim, "null communities in parallel...\n")
    
    nm <- nullmodel(emp_mat, method = "quasiswap")
    
    # We generate the nulls in parallel batches
    null_list <- mclapply(1:n_sim, function(i) {
      simulate(nm, nsim = 1)[,,1]
    }, mc.cores=detectCores())
    
    # --- STEP 2: EVALUATE THRESHOLDS ---
    for (t_name in names(thresholds)) {
      cat("  -> Calculating statistics for:", t_name, "\n")
      
      # Binarize Predictions
      thresh_vals <- thresholds[[t_name]]
      bin_mat <- sweep(pred_mat, 2, thresh_vals[colnames(pred_mat)], ">") * 1
      
      # Observed Jaccard
      obs_jaccard <- calc_jaccard_mean(emp_mat, bin_mat)
      
      # Null Jaccards (This part is now fast as nulls are already generated)
      null_jaccards <- sapply(null_list, function(n_mat) {
        calc_jaccard_mean(emp_mat, n_mat)
      })
      
      # Final Stats
      p_val <- (sum(null_jaccards >= obs_jaccard) + 1) / (n_sim + 1)
      ses <- (obs_jaccard - mean(null_jaccards)) / sd(null_jaccards)
      
      all_results[[length(all_results) + 1]] <- data.frame(
        Model = model_folder,
        Threshold = t_name,
        Obs_Jaccard = round(obs_jaccard, 4),
        SES = round(ses, 3),
        P_Value = round(p_val, 4)
      )
    }
    # Clear memory for next model
    rm(null_list)
    gc()
  }
  return(bind_rows(all_results))
}

# ==============================================================================
# EXECUTION PIPELINE
# ==============================================================================

# STEP 1: Generate basic predictions and thresholds (Uncomment if you need to rerun)
# get_preds(dir, pattern)
# get_thresholds(dir, pattern, methods)
library(parallel)
# STEP 2: Run Null Model Validations
set.seed(42) # Ensure reproducible null communities
pattern <- "2026-03-13_06-58-56_Atlas1_MinOccs5_CoverageGoodAverage"
null_model_results <- evaluate_jaccard_null_models_fast(
  dir = dir, 
  pattern = pattern, 
  methods = methods, 
  n_sim = 100
)

# View the final statistical output
print(null_model_results)

# Optional: Save results to CSV
# write.csv(null_model_results, file = "jaccard_null_model_results.csv", row.names = FALSE)