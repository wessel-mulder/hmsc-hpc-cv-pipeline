# Run from repository root: Rscript scripts/5-vis/env-viz/dominant-land-use.R
# Highest percentage among the seven land-cover predictors retained in S0.
suppressPackageStartupMessages({library(sf); library(ggplot2); library(cowplot)})
input <- 'data-shared/preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData'
out_dir <- 'scripts/5-vis/env-viz'
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
e <- new.env(parent=emptyenv()); load(input,envir=e)
X <- as.data.frame(e$X); d <- e$Design
stopifnot(identical(rownames(X),rownames(d)),!anyDuplicated(d[c('site','atlas')]))
shared <- sort(Reduce(intersect,lapply(1:3,function(a)
  as.character(d$site[as.character(d$atlas)==as.character(a)]))))
keep <- as.character(d$site) %in% shared
X <- X[keep,,drop=FALSE]; d <- droplevels(d[keep,,drop=FALSE])
stopifnot(length(shared)>0,all(table(d$site)==3),all(table(d$atlas)==length(shared)))
labels <- c(perc_cropland='Cropland',perc_urban='Urban',perc_pasture='Pasture',
  perc_grass_shrub='Grass / shrub',perc_forest_deciduous='Deciduous forest',
  perc_forest_evergreen='Evergreen forest',perc_forest_mixed_other='Mixed / other forest')
stopifnot(setequal(names(X)[grepl('^perc_',names(X))],names(labels)))
cover <- as.matrix(X[,names(labels),drop=FALSE])
stopifnot(all(is.finite(cover)),all(cover>=0 & cover<=100+1e-8))
highest <- apply(cover,1,max)
# Detect ties explicitly, rather than allowing max.col to choose one at random.
ties <- rowSums(abs(cover-highest)<=1e-8)
winner <- unname(labels[max.col(cover,ties.method='first')])
winner[ties>1] <- 'Tied maximum'
winner[highest==0] <- 'No modelled cover'
# Match variable_colours in scripts/5-vis/vp-distributions/vp-geom-sankey-bump.R.
palette <- c('Cropland'='goldenrod1','Urban'='snow3','Pasture'='darkorange',
  'Grass / shrub'='#B39DDB','Deciduous forest'='#A6D96A',
  'Evergreen forest'='#006400','Mixed / other forest'='#87CEEB',
  'Tied maximum'='#454545','No modelled cover'='#D4D4D4')
present <- names(palette)[names(palette) %in% winner]
d$dominant <- factor(winner,levels=present)
xy <- st_coordinates(st_transform(st_as_sf(d,coords=c('lon','lat'),crs=23032),4326))
d$lon <- xy[,1]; d$lat <- xy[,2]
stopifnot(all(d$lon>7 & d$lon<16),all(d$lat>54 & d$lat<59))
mainland <- c(xmin=7.95,xmax=13,ymin=54.4,ymax=57.85)
bornholm <- c(xmin=14.65,xmax=15.2,ymin=54.95,ymax=55.32)
point_map <- function(a,bbox=mainland,border=FALSE) {
  z <- d[as.character(d$atlas)==as.character(a),,drop=FALSE]
  ggplot(z,aes(lon,lat,fill=dominant)) +
    geom_point(shape=21,colour='black',stroke=.2,size=1.8,show.legend=TRUE) +
    scale_fill_manual(values=palette,limits=present,drop=FALSE,name='Largest modelled\nland-cover category') +
    coord_quickmap(xlim=bbox[c('xmin','xmax')],ylim=bbox[c('ymin','ymax')],expand=FALSE) +
    theme_void(base_size=11) +
    theme(legend.position='none',plot.background=element_rect(fill='white',colour=NA),
      panel.border=if(border) element_rect(colour='grey35',fill=NA,linewidth=.55)
        else element_blank(),plot.margin=margin(0,0,0,0))
}
panels <- lapply(1:3,function(a) {
  map <- ggdraw(point_map(a)) +
    draw_plot(point_map(a,bornholm,TRUE),x=.79,y=.68,width=.18,height=.22)
  strip <- ggdraw()+draw_label(sprintf('Atlas %s — %s',a,c('1970s','1990s','2010s')[a]),
                               fontface='bold',size=11)+
    theme(plot.background=element_rect(fill='grey92',colour='white'))
  plot_grid(strip,map,ncol=1,rel_heights=c(.07,1))
})
legend <- get_legend(point_map(1)+theme(legend.position='right',
  legend.title=element_text(face='bold'),legend.text=element_text(size=10)))
body <- plot_grid(plot_grid(plotlist=panels,nrow=1),legend,nrow=1,rel_widths=c(1,.21))
heading <- ggdraw()+
  draw_label('Dominant modelled land use across atlas periods',x=.02,y=.77,hjust=0,fontface='bold',size=15)+
  draw_label(sprintf('%s shared cells; category with highest cover among the seven modelled land-use predictors',
                     length(shared)),x=.02,y=.25,hjust=0,size=9)
caption <- ggdraw()+draw_label(
  'Largest modelled category need not cover >50% of a cell. Water and other unmodelled classes are excluded.',
  x=.02,hjust=0,size=9)
p <- plot_grid(heading,body,caption,ncol=1,rel_heights=c(.13,1,.05))
ggsave(file.path(out_dir,'dominant-land-use-atlas1-2-3.png'),p,
       width=14,height=5.9,dpi=300,bg='white')
if(interactive()) print(p)
counts <- table(d$dominant,d$atlas)
stopifnot(all(colSums(counts)==length(shared)))
rows <- vapply(seq_len(nrow(counts)),function(i)
  paste0('| ',rownames(counts)[i],' | ',paste(counts[i,],collapse=' | '),' |'),character(1))
writeLines(c('# Dominant modelled land-use counts','',paste('Source:',input),'',
  sprintf('Same %s grid cells in each atlas. Counts below sum to this denominator per column.',length(shared)),
  'Each cell takes the category with highest original percentage among the seven modelled cover predictors.',
  'Ties within 1e-8 percentage points are flagged; all-zero rows would be labelled No modelled cover.',
  'Water and unmodelled classes are excluded. Dominance does not imply a majority (>50%) of the cell.','',
  '| Category | Atlas 1 | Atlas 2 | Atlas 3 |','|---|---:|---:|---:|',rows),
  file.path(out_dir,'dominant-land-use-summary.md'))
print(counts)
message('Saved dominant-land-use-atlas1-2-3.png and dominant-land-use-summary.md')
