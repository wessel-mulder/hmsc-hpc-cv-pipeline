# Taxonomic tree from IOC 15.1 ranks; not an evolutionary/time-calibrated phylogeny.
# Run from repository root: Rscript scripts/0-pre/preparing-taxonomy.R
prepare_taxonomy <- function(
    ioc_file = 'data-new/0-pre/raw/taxonomy/IOC_Names_File_Plus-15.1_wide.csv',
    taxa_file = 'data-new/0-pre/derived/Tr/Trait_table.csv',
    output_dir = 'data-new/0-pre/derived/taxonomy') {
  if (!requireNamespace('ape', quietly = TRUE)) stop('Install ape to build/write the tree.')
  stopifnot(file.exists(ioc_file), file.exists(taxa_file))
  ioc <- read.csv(ioc_file, stringsAsFactors = FALSE, check.names = FALSE)
  taxa <- read.csv(taxa_file, stringsAsFactors = FALSE)
  stopifnot(all(c('order','family','genus','species','ssp') %in% names(ioc)), 'species' %in% names(taxa))
  normalise <- function(z) gsub('[[:space:]]+', '_', trimws(z))
  tips <- normalise(taxa$species)
  if (length(tips) < 2L || anyNA(tips) || any(!nzchar(tips)) || anyDuplicated(tips)) stop('Need at least two distinct, nonmissing target taxa.')
  if (any(lengths(strsplit(tips, '_', fixed=TRUE)) != 2L)) stop('Species-level tips required; rerun preparing-traits.R to collapse subspecies.')
  ioc$species_key <- normalise(ioc$species)
  ioc$ssp_key <- normalise(ioc$ssp)
  # Multiple subspecies rows may repeat a species hierarchy; require a unique mapping.
  mapped <- lapply(tips, function(tip) {
    species_match <- which(!is.na(ioc$species_key) & ioc$species_key == tip)
    subspecies_match <- which(!is.na(ioc$ssp_key) & ioc$ssp_key == tip)
    idx <- species_match
    if (!length(idx)) stop('No IOC match for ', tip, '; no taxon is silently removed or replaced.')
    ranks <- unique(ioc[idx, c('order','family','genus','species'), drop = FALSE])
    if (nrow(ranks) != 1L || anyNA(ranks) || any(!nzchar(as.matrix(ranks)))) stop('Missing/ambiguous IOC ranks for ', tip)
    data.frame(tip = tip, class = 'Aves', ranks,
               subspecies = if (length(species_match)) NA_character_ else gsub('_',' ',tip),
               match_rank = if (length(species_match)) 'species' else 'subspecies',
               stringsAsFactors = FALSE)
  })
  taxonomy <- do.call(rbind, mapped)
  n <- nrow(taxonomy)
  distance <- matrix(0, n, n, dimnames = list(tips,tips))
  for (i in seq_len(n-1L)) for (j in (i+1L):n) {
    # Species-level rank distances; subspecies are collapsed upstream.
    value <- if (taxonomy$order[i] != taxonomy$order[j]) 8 else
      if (taxonomy$family[i] != taxonomy$family[j]) 6 else
      if (taxonomy$genus[i] != taxonomy$genus[j]) 4 else
      if (taxonomy$species[i] != taxonomy$species[j]) 2 else 1
    distance[i,j] <- distance[j,i] <- value
  }
  tree <- ape::as.phylo(hclust(as.dist(distance), method = 'average'))
  # UPGMA creates arbitrary binary resolutions of tied rank distances. Collapse
  # zero-length internal edges so these remain honest taxonomic polytomies.
  tree <- ape::di2multi(tree, tol = 1e-8)
  observed <- ape::cophenetic.phylo(tree)[tips,tips]
  stopifnot(setequal(tree$tip.label, tips), !anyDuplicated(tree$tip.label),
            isTRUE(all.equal(unname(observed), unname(distance), tolerance = 1e-8)),
            ape::is.ultrametric(tree), all(tree$edge.length >= 0))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output <- file.path(output_dir,'tree_fromPD.tre')
  ape::write.tree(tree, file = output)
  write.csv(taxonomy, file.path(output_dir,'IOC15.1_tip_taxonomy.csv'), row.names = FALSE, na = 'NA')
  write.csv(distance, file.path(output_dir,'taxonomic_distances.csv'), row.names = TRUE)
  writeLines(c(
    'IOC 15.1 TAXONOMIC TREE',
    'Reproduce: Rscript scripts/0-pre/preparing-taxonomy.R (repository root).',
    paste('IOC source:',ioc_file), paste('Target taxa:',taxa_file),
    'Target set: all taxa in Trait_table.csv, one tip per species, with subspecies collapsed upstream.',
    'Tip labels use underscores and match Trait_table.csv. No rarity/occurrence or status filters.',
    'Taxonomy ranks: Aves, order, family, genus, species; subspecies are collapsed upstream in preparing-traits.R.',
    'Pairwise distances: different orders 8; same order/different families 6;',
    'same family/different genera 4; same genus/different species 2;',
    'identical species 0; no subspecies tips.',
    'Build: UPGMA (average-linkage clustering), then collapse zero-length internal edges.',
    'These are arbitrary rank-distance units, not divergence times or estimated evolutionary rates.',
    'This is a taxonomy-based surrogate tree, not an independently inferred phylogeny.',
    'For Hmsc, prune to the final retained species set; phylogenetic interpretations require this limitation.',
    'Missing/ambiguous taxonomy or duplicate target taxa cause a clear error; no parent-species substitution.',
    paste('Tips:',length(tree$tip.label)),
    paste('Species-level tips:',sum(taxonomy$match_rank=='species')),
    paste('Subspecies tips:',sum(taxonomy$match_rank=='subspecies')),
    'Validation: exact target-tip agreement, patristic distances match the rank matrix, and nonnegative ultrametric edges.',
    paste('ape version:',as.character(utils::packageVersion('ape')))
  ),file.path(output_dir,'taxonomy_description.txt'))
  occurrence_file <- 'data-new/0-pre/raw/occurrences/Y_occurrences.csv'
  if (file.exists(occurrence_file)) {
    occurrence_species <- setdiff(names(read.csv(occurrence_file, nrows = 0L)), 'survey')
    absent <- setdiff(occurrence_species, tips)
    write.csv(data.frame(species = absent), file.path(output_dir,'occurrence_taxa_without_tree_tips.csv'), row.names = FALSE)
    if (length(absent)) warning('Occurrence taxa without exact tree tips: ', paste(absent,collapse=', '),
                               '. Reconcile taxonomy explicitly before model preparation; no substitution was applied.')
  }
  saved <- ape::read.tree(output)
  stopifnot(setequal(saved$tip.label,tips), length(saved$tip.label)==n,
            isTRUE(all.equal(unname(ape::cophenetic.phylo(saved)[tips,tips]),unname(distance),tolerance=1e-8)))
  plot_taxonomic_tree(saved, taxonomy, file.path(output_dir, 'tree_fromPD.png'))
  message('Saved ',n,' tips (',sum(taxonomy$match_rank=='subspecies'),' subspecies) to ',output)
  invisible(list(tree = saved, taxonomy = taxonomy, distance = distance))
}

# PNG only, preserving the stored tree; stars identify subspecies in the display.
plot_taxonomic_tree <- function(tree, taxonomy, filename) {
  taxonomy <- taxonomy[match(tree$tip.label, taxonomy$tip), ]
  stopifnot(!anyNA(taxonomy$tip), identical(taxonomy$tip, tree$tip.label))
  orders <- sort(unique(taxonomy$order))
  palette <- setNames(grDevices::hcl.colors(length(orders), palette = 'Dark 3'), orders)
  display_tree <- tree
  display_tree$tip.label <- tree$tip.label
  grDevices::png(filename, width = 4320, height = 4560, res = 240, pointsize = 12, bg = 'white')
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(c(1,2), ncol=1), heights=c(10,1.4))
  graphics::par(mar=c(0,1,4,1))
  ape::plot.phylo(display_tree, type='fan', cex=0.68, label.offset=0.08,
                 tip.color=unname(palette[taxonomy$order]), edge.color='#696969',
                 edge.width=0.9, no.margin=FALSE)
  graphics::title(main='IOC 15.1 taxonomic tree', cex.main=1.7, line=2.2,
                 sub=NULL)
  graphics::mtext(paste(length(tree$tip.label),'taxa | rank-based distances | colours indicate order'),
                  side=3, line=0.6, cex=0.95)
  graphics::par(mar=c(3,1,0,1))
  graphics::plot.new()
  graphics::legend('center', legend=orders, col=palette, pch=16,
                   ncol=4, bty='n', cex=0.9, x.intersp=0.7, y.intersp=1.25)
  graphics::mtext('Taxonomy-based surrogate; branch lengths are not evolutionary time. Species-level tips only.',
                  side=1, line=0.5, cex=0.9)
  invisible(filename)
}

if (sys.nframe() == 0L) prepare_taxonomy()
