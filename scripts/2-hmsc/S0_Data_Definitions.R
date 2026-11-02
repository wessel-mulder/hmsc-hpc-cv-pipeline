# Prepare Hmsc inputs from the current preprocessing outputs.
# Run from repository root: Rscript scripts/2-hmsc/S0_Data_Definitions.R
# This creates a new forest-split dataset; the source CSVs remain unchanged.
suppressPackageStartupMessages({library(dplyr); library(readxl)})
coveragelevels <- 'good_and_average'
min_occs <- 5L
raw_dir <- 'data-new/0-pre/raw'
derived_dir <- Sys.getenv('CH1_PREPROCESSING_DERIVED','data-new/0-pre/derived')
output_dir <- Sys.getenv('CH1_HMSC_DERIVED','data-new/2-hmsc/derived')
paths <- c(environment = file.path(derived_dir,'environmental/X_environment.csv'),
           traits = file.path(derived_dir,'Tr/Trait_table.csv'),
           occurrences = file.path(raw_dir,'occurrences/Y_occurrences.csv'),
           design = file.path(raw_dir,'design/studyDesign.csv'),
           effort = file.path(raw_dir,'effort/dækning atlas1+2.xlsx'))
if (any(!file.exists(paths))) stop('Missing inputs:\n',paste(paths[!file.exists(paths)],collapse='\n'))

# COVERAGE: original Good/Average rule; Atlas 2 uses the corrected field.
get_effort <- function(path) {
  a1 <- read_excel(path,sheet='kvadrater atlas1')
  a2 <- read_excel(path,sheet='kvadrater atlas2')
  data.frame(survey=c(paste0(a1$KVADRATKOD,'_1'),paste0(a2$KVADRATKOD,'_2')),
             coverage_code=c(a1$UNDS_KODE_,a2$UNDS_KOD_1))
}
effort <- get_effort(paths[['effort']])
stopifnot(!anyDuplicated(effort$survey))

# ENVIRONMENT: user-defined grouping, unstandardised percentage points (0-100).
e <- read.csv(paths[['environment']],check.names=FALSE)
required <- c('survey','kvadratkod','atlas','tmean_breeding','prec_breeding','hh',
              paste0('pct_class_',c(0,11,22,33,40:45,55,77)))
stopifnot(all(required %in% names(e)),!anyDuplicated(e$survey),
          all(e$survey==paste(e$kvadratkod,e$atlas,sep='_')))
water_pct <- e$pct_class_0 + e$pct_class_77
# Exactly 75% water is retained; missing water coverage is excluded explicitly.
water_ok <- is.finite(water_pct) & water_pct <= 75
X_all <- data.frame(tmean_breeding=e$tmean_breeding,prec_breeding=e$prec_breeding,
                    hh=e$hh,
                    perc_cropland=e$pct_class_22,  # class 22 ONLY; not 23 or 24
                    perc_urban=e$pct_class_11,
                    perc_pasture=e$pct_class_33,
                    perc_grass_shrub=e$pct_class_55,
                    perc_forest_deciduous=e$pct_class_43+e$pct_class_44,
                    perc_forest_evergreen=e$pct_class_42+e$pct_class_41,
                    perc_forest_mixed_other=e$pct_class_40+e$pct_class_45,
                    row.names=e$survey)
complete_x <- apply(X_all,1,function(z) all(is.finite(z)))
allowed_codes <- switch(coveragelevels,good_and_average=c(1,2),good=1,
                        stop('Unknown coverage setting'))
coverage_ok <- e$survey %in% effort$survey[effort$coverage_code %in% allowed_codes]
# Preserve the original script's explicit Atlas 3 Good-coverage assumption.
coverage_ok[e$atlas==3] <- TRUE
eligible <- water_ok & complete_x & coverage_ok

# OCCURRENCES: retain matching surveys and species with >=5 occurrences in EACH atlas.
Y_all <- read.csv(paths[['occurrences']],row.names=1,check.names=FALSE)
stopifnot(!anyDuplicated(rownames(Y_all)),!anyDuplicated(names(Y_all)),
          !anyNA(Y_all),all(as.matrix(Y_all) %in% c(0,1)))
ids <- sort(intersect(e$survey[eligible],rownames(Y_all)))
Y <- Y_all[ids,,drop=FALSE]
X <- X_all[ids,,drop=FALSE]
atlas_by_row <- e$atlas[match(ids,e$survey)]
stopifnot(all(1:3 %in% atlas_by_row))
counts <- sapply(1:3,function(a) colSums(Y[atlas_by_row==a,,drop=FALSE]))
colnames(counts) <- paste0('atlas',1:3)
species_ok <- apply(counts,1,function(z) all(z>=min_occs))
Y <- Y[,species_ok,drop=FALSE]
stopifnot(ncol(Y)>0)

# TRAITS: consume the prepared species-level table without new proxies or species filters.
traits <- read.csv(paths[['traits']],stringsAsFactors=FALSE)
stopifnot(all(c('species','foraging_guild','migratory_strategy','STI') %in% names(traits)),
          !anyDuplicated(traits$species))
if(any(!names(Y) %in% traits$species)) stop('Species lack trait rows: ',paste(setdiff(names(Y),traits$species),collapse=', '))
traits <- traits[match(names(Y),traits$species),]
Tr <- data.frame(Migration_a3_DOF=traits$migratory_strategy,
                 foraging_guild_consensus=traits$foraging_guild,
                 species_thermal_index=traits$STI,row.names=traits$species)
missing_traits <- data.frame(species=rownames(Tr),Tr,check.names=FALSE)
missing_traits <- missing_traits[!complete.cases(Tr),,drop=FALSE]
# Do not silently change the >=5-occurrence species set to obtain complete traits.
if(nrow(missing_traits)) warning('Missing traits retained for: ',paste(missing_traits$species,collapse=', '),
                               '. These require resolution before a model uses those traits.')

# DESIGN: same projected coordinate convention as the previous script.
design <- read.csv(paths[['design']],check.names=FALSE)
stopifnot(all(c('survey','site','atlas','lon','lat') %in% names(design)),!anyDuplicated(design$survey))
if(any(!ids %in% design$survey)) stop('Some selected surveys lack design rows')
Design <- design[match(ids,design$survey),setdiff(names(design),'survey'),drop=FALSE]
rownames(Design) <- ids
stopifnot(all(Design$atlas==atlas_by_row),all(paste(Design$site,Design$atlas,sep='_')==ids),
          all(is.finite(Design$lon)),all(is.finite(Design$lat)))
Design$site <- factor(Design$site); Design$atlas <- factor(Design$atlas)
Design$year <- c('1'=1971,'2'=1992,'3'=2014)[as.character(Design$atlas)]
v <- terra::vect(Design[,c('lon','lat')],geom=c('lon','lat'),crs='EPSG:4326')
xy <- terra::crds(terra::project(v,'EPSG:23032'))
Design$lon <- xy[,1]; Design$lat <- xy[,2]

# Validate alignment, units and grouping before saving.
stopifnot(identical(rownames(X),rownames(Y)),identical(rownames(X),rownames(Design)),
          identical(colnames(Y),rownames(Tr)),ncol(X)==10L,!anyNA(X),
          all(as.matrix(X[grepl('^perc_',names(X))])>=0),
          all(as.matrix(X[grepl('^perc_',names(X))])<=100+1e-8),
          all(water_pct[match(ids,e$survey)]<=75),
          isTRUE(all.equal(X$perc_cropland,e$pct_class_22[match(ids,e$survey)])))
# Bundle the full species-level taxonomy tree; S1 prunes it to its selected species.
tree_file <- file.path(derived_dir,'taxonomy/tree_fromPD.tre')
if(!file.exists(tree_file)) stop('Taxonomy tree missing; run scripts/0-pre/preparing-taxonomy.R first.')
phy <- ape::read.tree(tree_file)
stopifnot(inherits(phy,'phylo'),!anyDuplicated(phy$tip.label),all(colnames(Y) %in% phy$tip.label))
# The original per-atlas coverage filtering is retained; no new common-cell restriction.
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
output_file <- file.path(output_dir,sprintf('preprocessed_data_min_occs_is_%s_coverage_is_%s_forest_split.RData',min_occs,coveragelevels))
save(X,Y,Tr,Design,phy,file=output_file)

# Read back the saved artifact independently.
check <- new.env(parent=emptyenv()); load(output_file,envir=check)
stopifnot(setequal(ls(check),c('X','Y','Tr','Design','phy')),identical(check$X,X),identical(check$Y,Y),
          identical(check$Tr,Tr),identical(check$Design,Design),identical(check$phy,phy))

# Human-readable audit files accompany the RData.
selection <- data.frame(survey=e$survey,atlas=e$atlas,ocean_freshwater_pct=water_pct,
                        water_ok,complete_predictors=complete_x,coverage_ok,
                        has_occurrences=e$survey %in% rownames(Y_all),retained=e$survey %in% ids)
write.csv(selection,file.path(output_dir,'forest_split_survey_filter_audit.csv'),row.names=FALSE)
write.csv(data.frame(species=rownames(counts),counts,retained=species_ok),
          file.path(output_dir,'forest_split_species_filter_audit.csv'),row.names=FALSE)
write.csv(missing_traits,file.path(output_dir,'forest_split_missing_traits.csv'),row.names=FALSE)
writeLines(c('S0 HMSC INPUTS: FOREST-SPLIT DATASET',
             paste('RData:',output_file),
             'Objects: X, Y, Tr, Design, phy (full tree; pruned in S1). X/Y/Design rows and Y columns/Tr rows explicitly aligned.',
             'Filter out ocean class 0 + freshwater class 77 >75%; exactly 75% retained. Applied per site-atlas row.',
             'Coverage: Good/Average for Atlases 1-2; corrected Atlas 2 UNDS_KOD_1. Atlas 3 retains original assumed-Good rule.',
             'Require complete selected predictors and occurrences. Species require >=5 presences in every atlas.',
             'No additional rare-species, trait-completeness or common-cell filters are applied.',
             'X: breeding temperature (degrees C), breeding precipitation (mm/season), hh, and seven cover percentages (0-100).',
             'hh: mean annual Simpson diversity across all 15 original HILDA classes, including ocean and freshwater.',
             'Cropland 22 ONLY; urban 11; pasture 33; grass/shrub 55;',
             'deciduous forest 43+44; evergreen forest 42+41; mixed/other forest 40+45.',
             'Tree crops (23), agroforestry (24), water, ocean and sparse vegetation are not separate cover predictors; all contribute to hh.',
             'Trait columns preserve S1 names. Missing trait values retained and reported; no new related-taxon substitution.',
             'Coordinates are projected EPSG:23032; lon/lat column names retained for compatibility.',
             'Design year labels preserve previous S0: 1971, 1992, 2014 (not literal climate-period definitions).',
             'Atlas 3 climate: 2014-2016; HILDA: 2014-2017.',
             paste('Surveys:',nrow(X)),paste('Species:',ncol(Y)),paste('Predictors:',ncol(X)),
             capture.output(table(Design$atlas)),
             'Shared sites across all atlases:',
             as.character(length(Reduce(intersect,lapply(1:3,function(a)as.character(Design$site[Design$atlas==a]))))),
             'S1 reads the portable RData from data-shared and prunes phy to the selected model species.'),
           file.path(output_dir,'forest_split_description.txt'))
message('Saved and read-back verified: ',output_file)
message(nrow(X),' surveys; ',ncol(Y),' species; ',ncol(X),' predictors.')
print(table(Design$atlas))

# Portable HPC bundle. Local raw and derived data remain outside Git.
shared_dir <- Sys.getenv('CH1_SHARED_DATA_DIR','data-shared')
dir.create(shared_dir,recursive=TRUE,showWarnings=FALSE)
destination <- file.path(shared_dir,basename(output_file))
if(!file.copy(output_file,destination,overwrite=TRUE))stop('Cannot export shared input: ',destination)
stopifnot(unname(tools::md5sum(output_file))==unname(tools::md5sum(destination)))
message('Shared RData (including full phy tree) exported to ',shared_dir)
