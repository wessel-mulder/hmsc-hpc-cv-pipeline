pattern <- "2026-07-27"
model_dir <- file.path("intermediate-outputs", "hmsc-output")
# get folders
matching_folders <- list.dirs(model_dir, recursive = FALSE, full.names = F)
matching_folders <- matching_folders[grepl(pattern, basename(matching_folders))]

model<-matching_folders[1]
models_nums <- as.numeric(str_extract(matching_folders, "(?<=Atlas)\\d+"))

model1 <- model[[1]]
load(file.path(model_dir,
                           model,
                           'Models',
                           'Fitted',
                           'HPC_samples_0250_thin_100_chains_4.rdata'))
atlas1<-fitted_model$posteriors
atlas1 <-poolMcmcChains(fitted_model)
draw <- atlas1$postList[[1]][[251]]
draw$Beta
