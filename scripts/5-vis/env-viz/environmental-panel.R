# Run from repository root: Rscript scripts/5-vis/env-viz/environmental-panel.R
# S0 environmental inputs in original units, on the same cells in every period.
suppressPackageStartupMessages({library(sf); library(ggplot2); library(cowplot)})
input <- 'data-shared/preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData'
out_dir <- 'scripts/5-vis/env-viz'
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
e <- new.env(parent=emptyenv()); load(input,envir=e)
X <- as.data.frame(e$X); d <- e$Design
stopifnot(identical(rownames(X),rownames(d)),
          !anyDuplicated(d[c('site','atlas')]),all(is.finite(as.matrix(X))))
shared <- sort(Reduce(intersect,lapply(1:3,function(a)
  as.character(d$site[as.character(d$atlas)==as.character(a)]))))
keep <- as.character(d$site) %in% shared
X <- X[keep,,drop=FALSE]; d <- droplevels(d[keep,,drop=FALSE])
stopifnot(length(shared)>0,all(table(d$site)==3),all(table(d$atlas)==length(shared)))
# S0 stores projected ED50 / UTM 32N metres despite the lon/lat column names.
xy <- st_coordinates(st_transform(st_as_sf(d,coords=c('lon','lat'),crs=23032),4326))
d$lon <- xy[,1]; d$lat <- xy[,2]
stopifnot(all(d$lon>7 & d$lon<16),all(d$lat>54 & d$lat<59))
labels <- c(tmean_breeding='Breeding-season temperature',
  prec_breeding='Breeding-season precipitation',hh='Habitat heterogeneity',
  perc_cropland='Cropland cover',perc_urban='Urban cover',
  perc_pasture='Pasture cover',perc_grass_shrub='Grass / shrub cover',
  perc_forest_deciduous='Deciduous forest cover',
  perc_forest_evergreen='Evergreen forest cover',
  perc_forest_mixed_other='Mixed / other forest cover')
units <- c(tmean_breeding='°C',prec_breeding='mm / season',hh='Simpson index',
           setNames(rep('%',7),names(labels)[4:10]))
stopifnot(setequal(names(X),names(labels)))
mainland <- c(xmin=7.95,xmax=13,ymin=54.4,ymax=57.85)
bornholm <- c(xmin=14.65,xmax=15.2,ymin=54.95,ymax=55.32)
periods <- c('1970s','1990s','2010s')
strip <- function(label) ggdraw() + draw_label(label,fontface='bold',size=11) +
  theme(plot.background=element_rect(fill='grey92',colour='white'))
figures <- list()
for(v in names(labels)) {
  d$value <- X[[v]]
  # One scale for all three periods, in original units; no within-atlas scaling.
  limits <- range(d$value)
  point_map <- function(a,bbox=mainland,border=FALSE) {
    z <- d[as.character(d$atlas)==as.character(a),,drop=FALSE]
    ggplot(z,aes(lon,lat,fill=value)) +
      geom_point(shape=21,colour='white',stroke=.12,size=1.8) +
      scale_fill_gradientn(colours=c('#2166AC','#FFFFBF','#B2182B'),
                           limits=limits,name=units[[v]]) +
      coord_quickmap(xlim=bbox[c('xmin','xmax')],ylim=bbox[c('ymin','ymax')],expand=FALSE) +
      theme_void(base_size=11) +
      theme(legend.position='none',plot.background=element_rect(fill='white',colour=NA),
        panel.border=if(border) element_rect(colour='grey35',fill=NA,linewidth=.55)
          else element_blank(),plot.margin=margin(0,0,0,0))
  }
  panels <- lapply(1:3,function(a) {
    map <- ggdraw(point_map(a)) +
      draw_plot(point_map(a,bornholm,TRUE),x=.79,y=.68,width=.18,height=.22)
    plot_grid(strip(sprintf('Atlas %s — %s',a,periods[a])),map,
              ncol=1,rel_heights=c(.07,1))
  })
  legend <- get_legend(point_map(1)+theme(legend.position='right',
    legend.title=element_text(face='bold'),legend.text=element_text(size=10)))
  body <- plot_grid(plot_grid(plotlist=panels,nrow=1),legend,nrow=1,rel_widths=c(1,.13))
  heading <- ggdraw() + draw_label(labels[[v]],x=.02,y=.77,hjust=0,fontface='bold',size=14) +
    draw_label(sprintf('%s shared training cells; original units; common scale across periods',
                       length(shared)),x=.02,y=.25,hjust=0,size=9)
  p <- plot_grid(heading,body,ncol=1,rel_heights=c(.13,1))
  figures[[v]] <- p
  message(v,': ',paste(signif(limits,5),collapse=' to '),' ',units[[v]])
}
# A tall overview preserves readable labels.
overview <- plot_grid(plotlist=figures,ncol=1)
ggsave(file.path(out_dir,'environmental-panel.png'),overview,
       width=13,height=57,dpi=160,bg='white',limitsize=FALSE)
if(interactive()) print(overview)
message('Saved environmental-panel.png to ',out_dir)
