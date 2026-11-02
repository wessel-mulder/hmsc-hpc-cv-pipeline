# Run from the repository root:
# Rscript scripts/0-pre/preparing-environment.R --plan
# Rscript scripts/0-pre/preparing-environment.R --workers=3
# Options: --landuse-only, --refresh-landuse, --workers=N, --plan.
library(terra)
library(sf)
if (!requireNamespace('exactextractr', quietly = TRUE)) stop('Install exactextractr for climate extraction.')
args <- commandArgs(trailingOnly = TRUE)
plan_only <- '--plan' %in% args
run_climate <- !('--landuse-only' %in% args)
refresh_landuse <- '--refresh-landuse' %in% args
worker_arg <- grep('^--workers=', args, value = TRUE)
workers <- if (length(worker_arg)) as.integer(sub('^--workers=', '', worker_arg[1])) else 3L
stopifnot(length(workers) == 1L, is.finite(workers), workers >= 1L)

# Spatial grid is the only required shapefile. Climate comes directly from source rasters.
grid_path <- 'data-new/0-pre/raw/DOF_Shapefiles_/DK5km_ED50grid_approx_kvadrkod_DOF.shp'
grid <- st_read(grid_path, quiet = TRUE)
atlas <- list(one = 1971:1974, two = 1993:1996, three = 2014:2017)
# CHELSAcruts ends in 2016. Keep the climate/atlas coverage difference explicit.
climate_atlas <- list(one = 1971:1974, two = 1993:1996, three = 2014:2016)
output_path <- file.path(Sys.getenv('CH1_PREPROCESSING_DERIVED', 'data-new/0-pre/derived'), 'environmental')
hilda_path <- file.path(output_path, 'per-atlas', 'hilda-v2')
climate_path <- file.path(output_path, 'per-atlas', 'chelsacruts')
climate_cache <- 'data-new/0-pre/raw/chelsacruts-v1-grid-crops'
source('scripts/0-pre/environmental/process-hilda.R')
source('scripts/0-pre/environmental/process-chelsacruts.R')

# PART 1: HILDA+ v2 annual classes and period summaries ----------------------
hilda_files <- file.path(hilda_path, paste0('hilda_atlas', seq_along(atlas), '.csv'))
if (!plan_only && (refresh_landuse || !all(file.exists(hilda_files)))) {
  hilda_results <- process_hilda_v2(
    grid, atlas, 'data-new/0-pre/raw/hildap_vGLOB-2.0_geotiff_wgs84/states', hilda_path)
} else if (!plan_only) {
  # Recalculate heterogeneity from separate original classes even when reusing annual extractions.
  hilda_results <- summarise_hilda_classes(
    read.csv(file.path(hilda_path, 'hilda_v2_annual.csv')), atlas, hilda_path)
}

# PART 2: CHELSAcruts monthly rasters, grid means, annual/period metrics ------
# Breeding only: March-June mean (monthly tmin + monthly tmax)/2 and total precipitation.
if (run_climate) {
  climate_results <- process_chelsacruts(grid, climate_atlas, climate_cache,
                                       climate_path, plan_only = plan_only, workers = workers, offline = '--offline' %in% args)
  if (plan_only) {
    message('Plan: ', nrow(climate_results), ' monthly files; no raster reads or combined CSV written.')
  } else {
    # PART 3: one row per site and atlas, all environmental variables ---------
    x_environment <- write_x_environment(grid, atlas, climate_results$periods,
                                         hilda_path, file.path(output_path, 'X_environment.csv'))
    qa <- data.frame(variable = names(x_environment),
                     missing = vapply(x_environment, function(x) sum(is.na(x)), integer(1)))
    write.csv(qa, file.path(output_path, 'X_environment_missingness.csv'), row.names = FALSE)
    message('Saved ', nrow(x_environment), ' rows and ', ncol(x_environment),
            ' columns to ', file.path(output_path, 'X_environment.csv'))
  }
}

if (run_climate && !plan_only) {
  file.copy('scripts/0-pre/environmental/environment-description.txt',
            file.path(output_path,'X_environment_description.txt'),overwrite=TRUE)
}
