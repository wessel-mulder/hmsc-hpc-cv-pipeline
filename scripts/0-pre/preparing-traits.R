# Merge raw Danish bird traits and Devictor STI. Run from repository root:
# Rscript scripts/0-pre/preparing-traits.R
# No occurrence/status filters. IOC subspecies are collapsed to their parent species.
prepare_trait_table <- function(
    traits_file = 'data-new/0-pre/raw/traits/birdTraits-merged-all.xlsx',
    sti_file = 'data-new/0-pre/raw/sti/STI_Devictor.csv',
    output_dir = 'data-new/0-pre/derived/Tr',
    ioc_file = 'data-new/0-pre/raw/taxonomy/IOC_Names_File_Plus-15.1_wide.csv') {
  if (!requireNamespace('readxl', quietly = TRUE)) stop('Install readxl to read the trait workbook.')
  stopifnot(file.exists(traits_file), file.exists(sti_file))
  traits <- readxl::read_excel(traits_file, sheet = 'birdTraits-merged', .name_repair = 'unique_quiet')
  sti <- read.csv(sti_file, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c('latin_DOF', 'foraging_guild_consensus', 'Migration_a3_DOF', 'Migration_a3_DOF_description')
  stopifnot(all(required %in% names(traits)), all(c('SPECIES','STI') %in% names(sti)))
  normalise <- function(x) gsub('[[:space:]]+', '_', trimws(x))
  species <- normalise(traits$latin_DOF)
  source_taxa <- species
  ioc <- read.csv(ioc_file, stringsAsFactors = FALSE)
  ioc_subspecies <- normalise(ioc$ssp)
  ioc_species <- normalise(ioc$species)
  for (i in seq_along(species)) {
    if (length(strsplit(species[i], '_', fixed = TRUE)[[1]]) > 2L) {
      parent <- unique(ioc_species[which(!is.na(ioc_subspecies) & ioc_subspecies == species[i])])
      if (length(parent) != 1L || is.na(parent)) stop('Missing/ambiguous IOC parent for ',species[i])
      species[i] <- parent
    }
  }
  sti_species <- normalise(sti$SPECIES)
  if (anyNA(species) || any(!nzchar(species)) || anyDuplicated(source_taxa)) stop('Missing or duplicate trait species identifiers.')
  if (anyNA(sti_species) || anyDuplicated(sti_species)) stop('Missing or duplicate STI species identifiers.')
  # Carry forward only the six explicit taxonomic name substitutions already in S0.
  # The Corvus corone -> Corvus cornix proxy used by S0 is deliberately not a synonym.
  aliases <- c(Astur_gentilis = 'Accipiter_gentilis',
               Anarhynchus_alexandrinus = 'Charadrius_alexandrinus',
               Chroicocephalus_ridibundus = 'Larus_ridibundus',
               Coloeus_monedula = 'Corvus_monedula',
               Curruca_communis = 'Sylvia_communis',
               Curruca_curruca = 'Sylvia_curruca')
  match_name <- species
  use_alias <- !species %in% sti_species & species %in% names(aliases)
  match_name[use_alias] <- aliases[species[use_alias]]
  # User-authorised STI proxy; retain Corvus cornix as a separate species.
  use_proxy <- species == 'Corvus_cornix'
  if (any(use_proxy) && !('Corvus_corone' %in% sti_species)) stop('Corvus corone STI proxy source is missing.')
  match_name[use_proxy] <- 'Corvus_corone'
  index <- match(match_name, sti_species)
  if (any(use_proxy) && anyNA(sti$STI[index[use_proxy]])) stop('Corvus corone STI proxy value is missing.')
  result <- data.frame(species = species,
                       foraging_guild = traits$foraging_guild_consensus,
                       migratory_strategy = traits$Migration_a3_DOF,
                       STI = sti$STI[index], stringsAsFactors = FALSE)
  # Merge duplicate species only when nonmissing trait values agree; never average categories.
  result <- do.call(rbind, lapply(unique(species), function(name) {
    z <- result[result$species == name, , drop=FALSE]
    row <- z[1, , drop=FALSE]
    for (column in c('foraging_guild','migratory_strategy','STI')) {
      values <- unique(z[[column]][!is.na(z[[column]])])
      if (length(values)>1L) stop('Conflicting ',column,' values when collapsing ',name)
      row[[column]] <- if(length(values)) values[1] else z[[column]][1]
    }
    row
  }))
  rownames(result) <- NULL
  stopifnot(nrow(result) == length(unique(species)), ncol(result) == 4L,
            !anyDuplicated(result$species), all(lengths(strsplit(result$species,'_',fixed=TRUE))==2L))
  audit <- data.frame(species = species, source_taxon = source_taxa,
                      collapsed_to_species = source_taxa != species, source_latin_DOF = traits$latin_DOF,
                      source_foraging_guild_consensus = traits$foraging_guild_consensus,
                      source_Migration_a3_DOF = traits$Migration_a3_DOF,
                      migration_description = traits$Migration_a3_DOF_description,
                      sti_lookup_name = match_name,
                      sti_source_species = sti$SPECIES[index],
                      sti_match = ifelse(use_proxy, 'user_authorised_STI_proxy', ifelse(is.na(index), 'unmatched', ifelse(use_alias, 'S0_name_crosswalk', 'exact'))),
                      missing_guild = is.na(result$foraging_guild[match(species,result$species)]),
                      missing_migration = is.na(result$migratory_strategy[match(species,result$species)]),
                      missing_STI = is.na(sti$STI[index]), stringsAsFactors = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output <- file.path(output_dir, 'Trait_table.csv')
  write.csv(result, output, row.names = FALSE, na = 'NA')
  write.csv(audit, file.path(output_dir, 'Trait_table_match_audit.csv'), row.names = FALSE, na = 'NA')
  issues <- audit[audit$missing_guild | audit$missing_migration | audit$missing_STI, ]
  write.csv(issues, file.path(output_dir, 'Trait_table_missing_values.csv'), row.names = FALSE, na = 'NA')
  writeLines(c(
    'TRAIT_TABLE.CSV: DATA AND PROCESSING DESCRIPTION',
    'Reproduce: Rscript scripts/0-pre/preparing-traits.R (from repository root).',
    paste('Trait source:', traits_file, '; sheet birdTraits-merged.'),
    paste('STI source:', sti_file, '; columns SPECIES and STI.'),
    'Output columns:',
    'species: latin_DOF with spaces replaced by underscores; IOC subspecies collapsed to their species names.',
    'foraging_guild: unchanged foraging_guild_consensus from the workbook.',
    'migratory_strategy: unchanged numeric Migration_a3_DOF coding:',
    '  1 = resident',
    '  1.5 = resident and short-distance migrant',
    '  2 = short-distance migrant',
    '  2.5 = short- and long-distance migrant',
    '  3 = long-distance migrant',
    'STI: unchanged numeric value from the supplied Devictor file; not standardised.',
    'All workbook taxa retained at species level; duplicate species collapsed only if nonmissing trait values agree.',
    'Match species names after whitespace/underscore normalisation; apply only the six S0 name substitutions when no exact match exists.',
    'Subspecies are explicitly collapsed to IOC parent species BEFORE matching species-level STI. User-authorised exception: Corvus cornix is assigned the Corvus corone STI, explicitly labelled as a proxy in the audit.',
    'NA means unavailable, not zero. No missing values are filled except the explicit Corvus cornix STI proxy.',
    'Trait_table_match_audit.csv preserves original names, lookup names, match method and migration descriptions.',
    'Trait_table_missing_values.csv lists any rows missing guild, migration or STI.',
    paste('Rows:', nrow(result)),
    paste('STI matched:', sum(!is.na(result$STI))),
    paste('STI unmatched/missing:', sum(is.na(result$STI))),
    paste('Missing foraging guild:', sum(is.na(result$foraging_guild))),
    paste('Missing migratory strategy:', sum(is.na(result$migratory_strategy)))
  ), file.path(output_dir, 'Trait_table_description.txt'))
  # Read back the actual delivered CSV and validate its keys and values.
  saved <- read.csv(output, stringsAsFactors = FALSE)
  stopifnot(identical(saved$species, result$species), identical(names(saved), names(result)),
            isTRUE(all.equal(saved$STI, result$STI)),
            isTRUE(all.equal(saved$migratory_strategy, result$migratory_strategy)))
  message('Saved ', nrow(result), ' taxa to ', output, '; ', sum(!is.na(result$STI)),
          ' STI values matched, ', sum(is.na(result$STI)), ' missing.')
  invisible(list(traits = result, audit = audit, missing = issues))
}

if (sys.nframe() == 0L) prepare_trait_table()
