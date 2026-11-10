# Run from the repository root or this analysis subfolder.
root <- normalizePath(if (dir.exists('data-new')) '.' else '../../..')
source(file.path(root,'scripts/1-atlas-analyses/homogenization/homogenization-data.R'),local=TRUE)

# Every unordered pair of distinct cells receives equal weight. No pair sampling.
spatial_beta <- function(x) {
  shared <- tcrossprod(x)
  richness <- rowSums(x)
  union <- outer(richness, richness, '+') - shared
  distance <- 1 - shared / union
  values <- distance[lower.tri(distance)]
  stopifnot(length(values) == choose(nrow(x), 2), all(is.finite(values)),
            all(values >= 0 & values <= 1), all(diag(distance) == 0),
            isTRUE(all.equal(distance, t(distance), check.attributes=FALSE)))
  # Independently check selected pairs using logical set membership.
  selected <- unique(round(seq(1, nrow(x), length.out=min(12L, nrow(x)))))
  checks <- combn(selected, 2)
  direct <- apply(checks, 2, function(ij) {
    a <- x[ij[1], ] == 1; b <- x[ij[2], ] == 1
    sum(xor(a, b)) / sum(a | b)
  })
  stopifnot(max(abs(direct - distance[t(checks)])) < 1e-12)
  mean(values)
}
beta <- vapply(matrices, spatial_beta, numeric(1))
atlas_summary <- data.frame(atlas=1:3, cells=N, species=S,
  cell_pairs=choose(N, 2), empty_cells=empty, mean_jaccard=beta)
pairs <- rbind(c(1,2), c(2,3), c(1,3))
changes <- data.frame(pair=c('1-2','2-3','1-3'),
  delta=beta[pairs[,2]] - beta[pairs[,1]])
changes$direction <- ifelse(changes$delta < 0, 'Homogenisation',
  ifelse(changes$delta > 0, 'Heterogenisation', 'No change'))
stopifnot(abs(sum(changes$delta[1:2]) - changes$delta[3]) < 1e-12)

library(ggplot2)
atlas_plot <- ggplot(atlas_summary, aes(atlas, mean_jaccard)) +
  geom_line(colour='#287D8E', linewidth=1) +
  geom_point(colour='#287D8E', size=3) +
  geom_text(aes(label=sprintf('%.5f', mean_jaccard)), vjust=-1, size=4) +
  scale_x_continuous(breaks=1:3, labels=paste('Atlas',1:3)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2)) +
  labs(x=NULL, y='Mean spatial Jaccard dissimilarity',
    title='Spatial differentiation of Danish bird assemblages',
    subtitle=sprintf('%s common cells | %s species | all unordered cell pairs',
      format(N,big.mark=','), S),
    caption='Lower values indicate more similar species lists among cells.\nEach cell pair has equal weight; self-comparisons are excluded.') +
  theme_bw(base_size=12) + theme(plot.caption=element_text(hjust=0))
changes$pair <- factor(changes$pair, levels=c('1-2','2-3','1-3'))
delta_limit <- max(abs(changes$delta), .001) * 1.45
change_plot <- ggplot(changes, aes(pair, delta, colour=direction)) +
  geom_hline(yintercept=0, colour='grey45', linetype='dashed') +
  geom_segment(aes(xend=pair, y=0, yend=delta), linewidth=1) +
  geom_point(size=3) +
  geom_text(aes(label=sprintf('%+.5f',delta)), vjust=-1, size=4, show.legend=FALSE) +
  scale_colour_manual(values=c(Homogenisation='#287D8E',
    Heterogenisation='#C86642', 'No change'='grey40')) +
  scale_x_discrete(labels=c('1-2'='Atlas 1 to 2','2-3'='Atlas 2 to 3','1-3'='Atlas 1 to 3')) +
  scale_y_continuous(limits=c(-delta_limit,delta_limit)) +
  labs(x=NULL, y='Change in mean spatial Jaccard (later minus earlier)',
    colour=NULL, title='Observed change in spatial differentiation',
    subtitle='Negative: homogenisation | Positive: heterogenisation',
    caption='Same cells and species in all comparisons. Descriptive changes; no uncertainty intervals or tests.') +
  theme_bw(base_size=12) + theme(legend.position='none', plot.caption=element_text(hjust=0))
ggsave(file.path(out,'homogenization-atlas-means.png'),atlas_plot,
  width=9,height=6,dpi=200,bg='white')
ggsave(file.path(out,'homogenization-changes.png'),change_plot,
  width=10,height=6,dpi=200,bg='white')

# Regenerate the numerical report with every run; retain only R, Markdown and PNG outputs.
report <- c('# Observed homogenisation of Danish bird assemblages', '',
  '## Question and method', '',
  'Following question five and its cookbook in Carsten’s `material/supp-analyses/Possibilities for Danish BirdAtlas Analysis of changes over time..pdf` (pp. 13–15; appendix code p. 26).', '',
  'Within each atlas, compare every unordered pair of different cells. Jaccard dissimilarity is the number of species present in only one of the two cells divided by the number present in either cell. Shared absences do not enter this denominator. Average these distances with equal weight per pair, then subtract the earlier atlas mean from the later atlas mean.', '',
  'This compares different cells within each period; it is not temporal turnover of a cell’s own species list.', '',
  '## Shared data and inclusion rules', '',
  '- Occurrences: `data-new/0-pre/raw/occurrences/Y_occurrences.csv`.',
  '- Cell/atlas identities: `data-new/0-pre/raw/design/studyDesign.csv`.',
  '- Land cover: `data-new/0-pre/derived/environmental/X_environment.csv`. Retain cells with `pct_class_0 + pct_class_77 <= 75` in every atlas: at least 25% land, excluding ocean and inland water.',
  '- Coverage: `data-new/0-pre/raw/effort/dækning atlas1+2.xlsx`. Retain codes 1 and 2 (Good/Average) in both Atlas 1 (`UNDS_KODE_`) and Atlas 2 (`UNDS_KOD_1`). Atlas 3 coverage is assumed adequate, not measured by this workbook.',
  '- Use the intersection across all three atlases, matching the current occupancy-change and unusual-temporal-turnover scripts. Keep all upstream species, including genuine zero-occurrence species within an atlas; apply no additional rarity filter.',
  '- Missing or nonbinary occurrence values stop execution. Empty communities also stop execution for an explicit policy decision; no cells are silently removed and no undefined distances are recoded.', '',
  sprintf('Validated sample: **%s cells × %s species**, with **%s unordered pairs per atlas**. Empty communities: %s in Atlases 1, 2 and 3.',
    format(N,big.mark=','), S, format(choose(N,2),big.mark=',',scientific=FALSE), paste(empty,collapse=', ')), '',
  '## Results', '', '| Atlas | Mean spatial Jaccard |', '| --- | ---: |',
  sprintf('| %d | %.8f |',1:3,beta), '',
  '| Atlas comparison | Change (later minus earlier) | Observed direction |',
  '| --- | ---: | --- |',
  sprintf('| %s | %+.8f | %s |',changes$pair,changes$delta,changes$direction), '',
  sprintf('From Atlas 1 to Atlas 3, mean spatial Jaccard changed from %.5f to %.5f (%+.5f), indicating observed %s on this fixed grid.',beta[1],beta[3],changes$delta[3],tolower(changes$direction[3])), '',
  '![Atlas means](homogenization-atlas-means.png)', '',
  '![Signed changes](homogenization-changes.png)', '',
  '## Interpretation and limits', '',
  'Negative changes mean that species lists became more similar on average (homogenisation); positive changes mean they became less similar (heterogenisation). These are complete descriptive summaries of the defined recorded-occurrence dataset. Pairwise distances are dependent and are not treated as independent replicates.', '',
  'The national all-pairs mean mixes nearby and distant comparisons and can conceal regional differences. Changes in richness and species occupancy may accompany the pattern. These results do not establish causes, correct for detection differences, or show departures from a random-allocation expectation. The descriptive script calculates no p-values. The separate paired spatial-block bootstrap runner supplies exploratory intervals in the R Markdown report. The separate independent-change null runner compares observed spatial change with baseline-fixed allocation of species gains and losses. Coverage validation, replacement/nestedness partitions and local analyses remain deferred.', '',
  '## Reproduce and validate', '',
  'From the repository root:', '', '```sh',
  'Rscript scripts/1-atlas-analyses/homogenization/homogenization.R', '```', '',
  'Or run `Rscript homogenization.R` from this folder. Requires R packages `readxl` and `ggplot2`. The script regenerates this report and the two PNGs; calculations remain available as `atlas_summary`, `changes`, `sites` and `matrices` when sourced into R.', '',
  'Checks cover binary values before conversion, aligned cell/species identities, duplicate identifiers, empty communities, pair counts, distance bounds and symmetry, zero self-distance, direct logical intersection/union checks on selected cell pairs, and the identity that changes 1–2 plus 2–3 equal change 1–3.', '',
  paste('Report generated:',format(Sys.Date())),
  paste('R version:',getRversion()))
writeLines(report, file.path(out,'README.md'))
print(atlas_summary, row.names=FALSE)
print(changes, row.names=FALSE)
message('Validated ',N,' common cells × ',S,' species; ',choose(N,2),' pairs per atlas.')

# Stage one: descriptive spatial diagnostics, not bootstrap inference.
jaccard_distances <- lapply(matrices, function(x) {
  shared <- tcrossprod(x)
  d <- 1 - shared / (outer(rowSums(x), rowSums(x), '+') - shared)
  diag(d) <- 0
  d
})
cell_change <- sapply(seq_len(nrow(pairs)), function(k) {
  rowSums(jaccard_distances[[pairs[k,2]]] - jaccard_distances[[pairs[k,1]]]) / (N-1)
})
colnames(cell_change) <- as.character(changes$pair)
rownames(cell_change) <- sites
stopifnot(max(abs(colMeans(cell_change) - changes$delta)) < 1e-12,
          max(abs(cell_change[,1] + cell_change[,2] - cell_change[,3])) < 1e-12)
coords <- design[design$atlas == 1, c('site','lon','lat')]
coords <- coords[match(sites, coords$site), ]
stopifnot(identical(as.character(coords$site), sites),
          !anyNA(coords), !anyDuplicated(coords$site),
          all(is.finite(coords$lon)), all(is.finite(coords$lat)))
# Coordinates describe the same cells in all three periods.
for (a in 2:3) {
  other <- design[match(paste0(sites,'_',a),design$survey),c('lon','lat')]
  stopifnot(max(abs(as.matrix(other)-as.matrix(coords[c('lon','lat')]))) < 1e-6)
}
projected <- sf::st_transform(sf::st_as_sf(coords,coords=c('lon','lat'),crs=4326),25832)
xy <- sf::st_coordinates(projected) / 1000
stopifnot(!anyDuplicated(as.data.frame(xy)),all(is.finite(xy)))
zxy <- scale(xy)
spatial_data <- data.frame(east=zxy[,1],north=zxy[,2])
trend_models <- lapply(seq_len(3),function(k) {
  lm(value ~ east + north + I(east^2) + I(north^2) + east:north,
     data=cbind(spatial_data,value=cell_change[,k]))
})
fitted_change <- sapply(trend_models,fitted)
residual_change <- sapply(trend_models,residuals)
stopifnot(all(is.finite(residual_change)),
          max(abs(fitted_change+residual_change-cell_change)) < 1e-12)
trend_summary <- data.frame(pair=as.character(changes$pair),
  spatial_R_squared=vapply(trend_models,function(m) summary(m)$r.squared,numeric(1)),
  raw_sd=apply(cell_change,2,sd),residual_sd=apply(residual_change,2,sd))
cell_diagnostics <- do.call(rbind,lapply(1:3,function(k) data.frame(
  site=sites,pair=as.character(changes$pair[k]),east_km=xy[,1],north_km=xy[,2],
  raw=cell_change[,k],trend=fitted_change[,k],residual=residual_change[,k])))
map_data <- do.call(rbind,lapply(c('raw','trend','residual'),function(type) {
  data.frame(cell_diagnostics[c('site','pair','east_km','north_km')],
    field=type,value=cell_diagnostics[[type]])
}))
map_data$field <- factor(map_data$field,levels=c('raw','trend','residual'),
  labels=c('Observed contribution','Quadratic spatial trend','Residual'))
map_data$pair <- factor(map_data$pair,levels=as.character(changes$pair),
  labels=c('Atlas 1 to 2','Atlas 2 to 3','Atlas 1 to 3'))
map_limit <- max(abs(map_data$value))
# Match the temporal-turnover maps: circular points, no axes/grid, grey
# headers, blue/cream/orange colours, and a framed Bornholm inset.
map_data$lon <- coords$lon[match(map_data$site,coords$site)]
map_data$lat <- coords$lat[match(map_data$site,coords$site)]
mainland_bbox <- c(xmin=7.95,xmax=13.0,ymin=54.4,ymax=57.85)
bornholm_bbox <- c(xmin=14.65,xmax=15.2,ymin=54.95,ymax=55.32)
map_fill <- function() scale_fill_gradient2(low='#2C7FB8',mid='#FFF7E6',
  high='#D95F0E',midpoint=0,limits=c(-map_limit,map_limit),name='Jaccard\nchange')
map_strip <- function(label) cowplot::ggdraw() + cowplot::draw_label(
  label,fontface='bold',size=11) + theme(
  plot.background=element_rect(fill='grey92',colour='white'),
  plot.margin=ggplot2::margin(0,0,1,0))
point_map <- function(z,bbox,border=FALSE,point_size=1.55) {
  ggplot(z,aes(lon,lat,fill=value)) +
    geom_point(shape=21,colour='white',stroke=.18,size=point_size) + map_fill() +
    coord_quickmap(xlim=bbox[c('xmin','xmax')],ylim=bbox[c('ymin','ymax')],expand=FALSE) +
    theme_minimal(base_size=10) + theme(panel.grid=element_blank(),
      panel.background=element_blank(),axis.text=element_blank(),
      axis.title=element_blank(),axis.ticks=element_blank(),legend.position='none',
      plot.background=element_rect(fill='white',colour=NA),
      panel.border=if(border) element_rect(colour='grey35',fill=NA,linewidth=.55) else element_blank(),
      plot.margin=ggplot2::margin(0,0,0,0))
}
map_rows <- lapply(levels(map_data$field),function(field_name) {
  panels <- lapply(levels(map_data$pair),function(pair_name) {
    z <- map_data[map_data$field==field_name & map_data$pair==pair_name,]
    main <- point_map(z,mainland_bbox)
    inset <- point_map(z,bornholm_bbox,border=TRUE,point_size=1.7)
    map <- cowplot::ggdraw(main) + cowplot::draw_plot(inset,x=.79,y=.68,width=.18,height=.22)
    cowplot::plot_grid(map_strip(pair_name),map,ncol=1,rel_heights=c(.07,1))
  })
  cowplot::plot_grid(map_strip(field_name),
    cowplot::plot_grid(plotlist=panels,nrow=1),ncol=1,rel_heights=c(.07,1))
})
legend_plot <- ggplot(map_data,aes(lon,lat,fill=value)) + geom_point(shape=21) +
  map_fill() + theme_void() + theme(legend.position='right',
    legend.title=element_text(face='bold'))
map_body <- cowplot::plot_grid(cowplot::plot_grid(plotlist=map_rows,ncol=1),
  cowplot::get_legend(legend_plot),nrow=1,rel_widths=c(1,.11))
map_title <- cowplot::ggdraw() + cowplot::draw_label(
  'Spatial pattern of contributions to national Jaccard change',
  x=0,y=.75,hjust=0,fontface='bold',size=14) + cowplot::draw_label(
  'Shared zero-centred scale: blue = negative, orange = positive; Bornholm shown in framed insets',
  x=0,y=.25,hjust=0,size=10)
map_caption <- cowplot::ggdraw() + cowplot::draw_label(
  'Negative observed values: a cell became more similar to other cells on average.\nResiduals remove an exploratory quadratic trend; points mark cell centres.',
  x=0,hjust=0,size=10)
map_plot <- cowplot::plot_grid(map_title,map_body,map_caption,ncol=1,
  rel_heights=c(.06,1,.05))
ggsave(file.path(out,'homogenization-spatial-diagnostics.png'),map_plot,
  width=13,height=15,dpi=250,bg='white')

# Unique unordered spatial pairs; distance bins are (0,10], (10,20], ..., (190,200].
ij <- which(lower.tri(jaccard_distances[[1]]),arr.ind=TRUE)
stopifnot(nrow(ij)==choose(N,2),all(ij[,1]>ij[,2]))
dx <- xy[ij[,1],1]-xy[ij[,2],1]
dy <- xy[ij[,1],2]-xy[ij[,2],2]
distance_km <- sqrt(dx^2+dy^2)
angle <- (atan2(dy,dx)*180/pi) %% 180
pair_bin <- ceiling(distance_km/10)
in_range <- distance_km>0 & distance_km<=200
sectors <- list(Omnidirectional=rep(TRUE,length(angle)),
  'East-west'=angle<=22.5 | angle>=157.5,
  'North-south'=angle>=67.5 & angle<=112.5)
stopifnot(!any(sectors[['East-west']] & sectors[['North-south']]),
          all(pair_bin[in_range] %in% 1:20))
variogram_list <- list(); v <- 0L
for (k in 1:3) for (field in c('Raw','Detrended')) {
  values <- if(field=='Raw') cell_change[,k] else residual_change[,k]
  field_variance <- var(values)
  semivariance <- .5*(values[ij[,1]]-values[ij[,2]])^2
  # Check an independent direct extraction of selected pair differences.
  check <- unique(round(seq(1,nrow(ij),length.out=20)))
  stopifnot(max(abs(semivariance[check]-vapply(check,function(q)
    (values[ij[q,1]]-values[ij[q,2]])^2/2,numeric(1)))) < 1e-12)
  for (direction in names(sectors)) {
    use <- in_range & sectors[[direction]]
    counts <- tabulate(pair_bin[use],nbins=20)
    stopifnot(sum(counts)==sum(use))
    means <- vapply(1:20,function(b) {
      ix <- use & pair_bin==b
      if(any(ix)) mean(semivariance[ix]) else NA_real_
    },numeric(1))
    v <- v+1L
    variogram_list[[v]] <- data.frame(pair=as.character(changes$pair[k]),field=field,
      direction=direction,bin=1:20,midpoint_km=seq(5,195,10),pairs=counts,
      semivariance=means,field_variance=field_variance,
      relative_semivariance=if(field_variance>0) means/field_variance else NA_real_,
      sparse=counts<100)
  }
}
variogram_summary <- do.call(rbind,variogram_list)
variogram_summary$direction <- factor(variogram_summary$direction,levels=names(sectors))
variogram_summary$field <- factor(variogram_summary$field,levels=c('Raw','Detrended'))
variogram_summary$pair <- factor(variogram_summary$pair,levels=as.character(changes$pair))
stopifnot(all(variogram_summary$semivariance[variogram_summary$pairs>0]>=0),
          all(is.finite(variogram_summary$relative_semivariance[variogram_summary$pairs>0])))
for (field_name in c('Raw','Detrended')) {
  z <- variogram_summary[variogram_summary$field==field_name,]
  upper <- data.frame(z,panel='Semivariance / field variance',value=z$relative_semivariance)
  lower <- data.frame(z,panel='Cell pairs per distance bin',value=z$pairs)
  plot_data <- rbind(upper,lower)
  plot_data$panel <- factor(plot_data$panel,
    levels=c('Semivariance / field variance','Cell pairs per distance bin'))
  vg <- ggplot(plot_data,aes(midpoint_km,value,colour=direction,group=direction)) +
    geom_hline(data=data.frame(panel=factor('Semivariance / field variance',levels=levels(plot_data$panel))),
      aes(yintercept=1),inherit.aes=FALSE,linetype='dashed',colour='grey55') +
    geom_line(linewidth=.6,na.rm=TRUE) +
    geom_point(aes(shape=sparse),size=1.6,na.rm=TRUE) +
    facet_grid(panel~pair,scales='free_y') +
    scale_shape_manual(values=c('FALSE'=16,'TRUE'=1),
      labels=c('FALSE'='At least 100 pairs','TRUE'='Fewer than 100 pairs'),name=NULL) +
    scale_colour_manual(values=c(Omnidirectional='#444444','East-west'='#287D8E','North-south'='#C86642'),name=NULL) +
    scale_x_continuous(breaks=seq(0,200,50),limits=c(0,200)) +
    labs(x='Separation between cell centres (km)',y=NULL,
      title=paste(field_name,'spatial dependence of cell contributions'),
      subtitle='Atlas comparisons 1-2, 2-3 and 1-3; 10-km bins; directional sectors +/-22.5 degrees',
      caption='Dashed line: semivariance equals the variance of the field. It is not an independence threshold.\nOpen symbols flag sparse bins; pair counts are not independent sample sizes. No range is fitted.') +
    theme_bw(base_size=11) + theme(legend.position='bottom',
      plot.caption=element_text(hjust=0),panel.grid.minor=element_blank())
  ggsave(file.path(out,paste0('homogenization-variogram-',tolower(field_name),'.png')),
    vg,width=13,height=8,dpi=200,bg='white')
}
# Concise binned summary for the R Markdown; no automated block-size selection.
regional_diagnostic_summary <- do.call(rbind,lapply(split(variogram_summary,
  interaction(variogram_summary$pair,variogram_summary$field,variogram_summary$direction,drop=TRUE)),
  function(z) data.frame(pair=as.character(z$pair[1]),field=as.character(z$field[1]),
    direction=as.character(z$direction[1]),populated_bins=sum(z$pairs>0),
    sparse_populated_bins=sum(z$pairs>0 & z$sparse),
    min_pairs=min(z$pairs),max_pairs=max(z$pairs))))
writeLines(c(report,'','## Stage-one spatial diagnostics','',
  'The script also produces cell-contribution, quadratic-trend and residual maps, and raw/detrended directional semivariograms with pair counts. See `homogenization.Rmd` and its rendered HTML for definitions, checks and interpretation. Requires `sf`; report rendering additionally requires `rmarkdown`, `knitr` and Pandoc.', '',
  'From the repository root:', '', '```sh',
  'Rscript -e \'rmarkdown::render("scripts/1-atlas-analyses/homogenization/homogenization.Rmd")\'',
  '```', '',
  'Additional in-memory objects: `cell_diagnostics`, `trend_summary`, `variogram_summary` and `regional_diagnostic_summary`. All figures are PNGs. Run `Rscript scripts/1-atlas-analyses/homogenization/run-spatial-bootstrap.R --replicates=2000 --seed=20261105` explicitly to save bootstrap results. Rendering reads completed compatible results without rerunning the bootstrap. No block width is selected automatically.'),
  file.path(out,'README.md'))
message('Stage-one spatial diagnostics validated for all three atlas comparisons.')

write(c('', '## Conditional allocation null', '',
  'Run `Rscript scripts/1-atlas-analyses/homogenization/run-homogenization-null.R --simulations=9999 --seed=20261105 --chunk-size=250 --workers=6` explicitly. Completed draws and provenance are saved under `null-results/independent-change/`; the R Markdown reads the latest compatible completed run without launching simulations.', '',
  'The null fixes the earlier occurrence matrix and each species’ observed gains and losses, sampling changes independently among eligible baseline cells. It preserves endpoint species occupancies while allowing their geographical arrangement to vary. The report shows observed and null changes separately, inclusive Monte Carlo tails, two-sided comparisons and Holm adjustment across three national comparisons.', '',
  'Validate the compiled statistic and sampler with `Rscript scripts/1-atlas-analyses/homogenization/validate-homogenization-null.R`. Use `--resume=/absolute/path/to/incomplete-run` to resume matching checkpoints; worker count may change. The runner retains no simulated matrices or distance arrays.'),file=file.path(out,'README.md'),append=TRUE)
