

# GAMMA SHIFTS ------------------------------------------------------------
library(Hmsc)
library(ggplot2)
library(tidyr)
library(dplyr)
library(abind)

dir <- './HmscOutputs'
pattern <- '2026-03-13'
matching_folders <- list.dirs(dir, recursive = FALSE, full.names = FALSE)
matching_folders <- matching_folders[grepl(pattern, basename(matching_folders))]

# Load Model 1
load(file.path(dir, matching_folders[1], 'Models/Fitted', 'HPC_samples_0250_thin_100_chains_4.Rdata'))
m1 <- fitted_model$posteriors
# Load Model 2
load(file.path(dir, matching_folders[2], 'Models/Fitted', 'HPC_samples_0250_thin_100_chains_4.Rdata'))
m2 <- fitted_model$posteriors

analyze_hmsc_shifts <- function(model1 = m1, model2 = m2, parameter = "Beta", title_label = "Niche") {
  library(Hmsc)
  library(dplyr)
  library(tidyr)
  library(abind)
  
  # 1. Prepare Scaling Factors (SD of Environmental Variables)
  # HMSC stores X in different places depending on how it was loaded
  x1 <- if(!is.null(model1$X)) model1$X else model1$XData
  x2 <- if(!is.null(model2$X)) model2$X else model2$XData
  
  # Extract SDs excluding intercept
  sd_x1 <- apply(as.matrix(x1)[, -1], 2, sd)
  sd_x2 <- apply(as.matrix(x2)[, -1], 2, sd)
  
  # 2. Pool Chains
  post1 <- poolMcmcChains(model1$postList)
  post2 <- poolMcmcChains(model2$postList)
  n_samples <- length(post1)
  
  # 3. Identify Names and Dimensions
  cov_names <- model1$covNames[-1] # Exclude Intercept
  # Column names depend on whether we look at Species (Beta) or Traits (Gamma)
  col_names <- if(parameter == "Beta") model1$spNames else model1$trNames
  
  n_covs <- length(cov_names)
  n_cols <- length(col_names)
  
  # 4. Standardized Posterior Subtraction & Individual Significance
  # Function to get standardized samples for a specific parameter
  get_std_samples <- function(post_list, sd_vec, param) {
    lapply(post_list, function(x) {
      # Extract matrix, remove intercept row
      mat <- x[[param]][-1, , drop = FALSE]
      # Multiply each row by its corresponding covariate SD
      return(mat * sd_vec)
    })
  }
  
  std1_list <- get_std_samples(post1, sd_x1, parameter)
  std2_list <- get_std_samples(post2, sd_x2, parameter)
  
  # Convert to 3D arrays [Cov, Col, Iteration]
  arr1 <- abind(std1_list, along = 3)
  arr2 <- abind(std2_list, along = 3)
  arr_delta <- arr2 - arr1
  
  # 5. Statistical Testing (Credible Intervals)
  get_sig_mask <- function(arr) {
    low <- apply(arr, c(1, 2), quantile, 0.025)
    upp <- apply(arr, c(1, 2), quantile, 0.975)
    return((low > 0) | (upp < 0))
  }
  
  sig1 <- get_sig_mask(arr1)
  sig2 <- get_sig_mask(arr2)
  sig_delta <- get_sig_mask(arr_delta)
  
  # 6. Categorization Logic
  status_matrix <- matrix("No Change", nrow = n_covs, ncol = n_cols)
  status_matrix[!sig1 & sig2] <- "Gained Importance"
  status_matrix[sig1 & !sig2] <- "Lost Importance"
  status_matrix[sig1 & sig2 & sig_delta] <- "Reinforced/Shifted"
  status_matrix[sig1 & sig2 & !sig_delta] <- "Stable Significance"
  
  rownames(status_matrix) <- cov_names
  colnames(status_matrix) <- col_names
  
  # 7. Format for Plotting
  plot_df <- as.data.frame(status_matrix) %>%
    mutate(Covariate = rownames(.)) %>%
    pivot_longer(cols = -Covariate, names_to = "Unit", values_to = "ShiftType")
  
  # 8. Generate Plot
  p <- ggplot(plot_df, aes(x = Unit, y = Covariate, fill = ShiftType)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c(
      "Gained Importance"   = "forestgreen",
      "Lost Importance"     = "dodgerblue3",
      "Reinforced/Shifted"  = "firebrick3",
      "Stable Significance" = "#f0f0f0",
      "No Change"           = "#f0f0f0"
    )) +
    theme_minimal() +
    labs(title = paste("Standardized", parameter, title_label, "Dynamics"),
         subtitle = "Comparing Model 1 vs Model 2 (Intercept Removed)",
         x = ifelse(parameter == "Beta", "Species", "Traits/Guilds"),
         y = "Environmental Variable",
         fill = "Type of Shift") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))
  
  return(list(plot = p, data = plot_df, matrix = status_matrix))
}
t <- analyze_hmsc_shifts(m1,m2,'Beta','test')
print(t$plot)
