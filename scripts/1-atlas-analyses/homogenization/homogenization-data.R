# Run from the repository root or this analysis subfolder.
root <- normalizePath(if (dir.exists('data-shared')) '.' else '../../..')
out <- file.path(root, 'scripts/1-atlas-analyses/homogenization')
design <- read.csv(file.path(root, 'data-new/0-pre/raw/design/studyDesign.csv'), check.names=FALSE)
environment <- read.csv(file.path(root, 'data-new/0-pre/derived/environmental/X_environment.csv'), check.names=FALSE)
stopifnot(all(c('site','atlas','survey') %in% names(design)),
          all(c('kvadratkod','atlas','pct_class_0','pct_class_77') %in% names(environment)),
          !anyDuplicated(design$survey), !anyDuplicated(environment[c('kvadratkod','atlas')]),
          setequal(design$survey, environment$survey),
          all(environment$survey == paste0(environment$kvadratkod, '_', environment$atlas)))
environment$water_percent <- environment$pct_class_0 + environment$pct_class_77
stopifnot(!anyNA(environment$water_percent), all(environment$water_percent >= 0 & environment$water_percent <= 100 + 1e-6))
eligible_sites <- lapply(1:3, function(a) intersect(
  as.character(environment$kvadratkod[environment$atlas == a & environment$water_percent <= 75]),
  as.character(design$site[design$atlas == a])))
sites <- sort(Reduce(intersect, eligible_sites))
effort_path <- file.path(root, 'data-new/0-pre/raw/effort/dækning atlas1+2.xlsx')
effort1 <- readxl::read_excel(effort_path, sheet='kvadrater atlas1')
effort2 <- readxl::read_excel(effort_path, sheet='kvadrater atlas2')
good1 <- as.character(effort1$KVADRATKOD[effort1[['UNDS_KODE_']] %in% c(1, 2)])
good2 <- as.character(effort2$KVADRATKOD[effort2[['UNDS_KOD_1']] %in% c(1, 2)])
sites <- sort(Reduce(intersect, list(sites, good1, good2)))
y <- read.csv(file.path(root,'data-new/0-pre/raw/occurrences/Y_occurrences.csv'), row.names=1, check.names=FALSE)
matrices <- lapply(1:3,function(a) {
 ids <- paste(sites,a,sep='_'); stopifnot(all(ids %in% rownames(y)))
 z <- as.matrix(y[ids,,drop=FALSE]); rownames(z)<-sites; z
})
N <- length(sites); S <- ncol(y); total <- N*S
stopifnot(!anyDuplicated(rownames(y)),!anyDuplicated(names(y)),all(vapply(matrices,function(z) !anyNA(z)&&all(z%in%c(0,1)),logical(1))))
stopifnot(N >= 2L, S == 201L,
          !anyDuplicated(design[c('site','atlas')]),
          all(vapply(matrices, function(z) identical(dimnames(z), dimnames(matrices[[1]])), logical(1))))
# Validate original values above, before converting storage mode.
matrices <- lapply(matrices, function(z) { storage.mode(z) <- 'double'; z })
empty <- vapply(matrices, function(z) sum(rowSums(z) == 0), integer(1))
if (any(empty > 0)) stop('Empty communities found (counts by atlas: ',
  paste(empty, collapse=', '), '). Resolve the empty-community policy without silently changing the common grid.')

