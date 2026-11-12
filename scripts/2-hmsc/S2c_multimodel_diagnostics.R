#!/usr/bin/env Rscript
# Two comparison heatmaps for all selected models. Run from the repository root.
# Edit this block; no fitted model or existing fit-cache file is changed.
.libPaths(c("~/Rlibs", .libPaths()))
patterns2match <- c("2026-11-11")
model_root <- Sys.getenv("CH1_MODEL_OUTPUT_DIR", "HmscOutputs")
samples <- 250L
thin <- 100L
nChains <- 4L
targets <- list(rhat = 1.01, ess_per_chain = 100, mcse_sd_ratio = 0.05)
output_root <- file.path(model_root, "multimodel-diagnostics")
atlas_labels <- c(Atlas1 = "1970s", Atlas2 = "1990s", Atlas3 = "2010s")
spatial_interval <- 0.95

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else
  "scripts/2-hmsc/S2c_multimodel_diagnostics.R"
source(file.path(dirname(normalizePath(script_path)), "multimodel-diagnostics", "helpers.R"))
source(file.path(dirname(normalizePath(script_path)), "multimodel-diagnostics", "spatial-factor-correlations.R"))
comparison_output <- run_multimodel_diagnostics(patterns2match, model_root, samples, thin, nChains,
                                               targets, output_root, atlas_labels)
run_spatial_factor_correlations(comparison_output, spatial_interval, atlas_labels)
