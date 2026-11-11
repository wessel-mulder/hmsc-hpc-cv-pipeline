# Run from repository root:
# Rscript scripts/5-vis/env-viz/atlas-richness-environment-panel.R
# Five variables (columns) across three atlases (rows); PNG only.
suppressPackageStartupMessages({library(sf); library(ggplot2); library(cowplot)})
out_dir <- 'scripts/5-vis/env-viz'
input <- 'data-shared/preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData'
e <- new.env(parent=emptyenv()); load(input,envir=e)
d <- e$Design; X <- as.data.frame(e$X); Y <- as.matrix(e$Y)
stopifnot(identical(rownames(d),rownames(X)),identical(rownames(d),rownames(Y)),
          !anyDuplicated(d[c('site','atlas')]),!anyNA(Y),all(Y %in% 0:1))
shared <- sort(Reduce(intersect,lapply(1:3,function(a)
  as.character(d$site[as.character(d$atlas)==as.character(a)]))))
keep <- as.character(d$site) %in% shared
d <- droplevels(d[keep,,drop=FALSE]); X <- X[keep,,drop=FALSE]; Y <- Y[keep,,drop=FALSE]
stopifnot(length(shared)>0,all(table(d$site)==3),all(table(d$atlas)==length(shared)))
xy <- st_coordinates(st_transform(st_as_sf(d,coords=c('lon','lat'),crs=23032),4326))
d$lon <- xy[,1]; d$lat <- xy[,2]
d$richness <- rowSums(Y)
d$temperature <- X$tmean_breeding; d$precipitation <- X$prec_breeding
d$heterogeneity <- X$hh
# Dominance is defined only among the seven land-cover predictors in X.
# Preserve the model's three separate forest groups and their exact values.
cover_labels <- c(perc_cropland='Cropland',perc_urban='Urban',
  perc_pasture='Pasture',perc_grass_shrub='Grass / shrub',
  perc_forest_deciduous='Deciduous forest',
  perc_forest_evergreen='Evergreen forest',
  perc_forest_mixed_other='Mixed / other forest')
stopifnot(setequal(grep('^perc_',names(X),value=TRUE),names(cover_labels)))
cover <- X[,names(cover_labels),drop=FALSE]
names(cover) <- unname(cover_labels)
stopifnot(all(is.finite(as.matrix(cover))),all(rowSums(cover)>0))
# Ties retain the first predictor in the documented order above.
d$landuse <- factor(names(cover)[max.col(as.matrix(cover),ties.method='first')],
                    levels=names(cover))
# Exact category colours from vp-distributions/vp-geom-sankey-bump.R.
palette <- c(Cropland='goldenrod1',Urban='snow3',Pasture='darkorange',
             'Grass / shrub'='#B39DDB','Deciduous forest'='#A6D96A',
             'Evergreen forest'='#006400','Mixed / other forest'='#87CEEB')
continuous_palettes <- list(
  richness=c('#6A00A8','#2166AC','#00BFC4','#33A02C','#FFFF33','#FF8C00','#D7191C'),
  temperature=c('#081D58','#1597A5','#25C9BC',
                '#E9B52F','#EE8A33','#E77BA5','#D73045','#67001F'),
  precipitation=c('#8C510A','#5AB4AC','#2166AC'),
  heterogeneity=c('#FFFFCC','#78C679','#006837','#8ECAE6','#084081'))
variables <- c('richness','temperature','precipitation','heterogeneity','landuse')
titles <- c('Species richness','Temperature','Precipitation','Landscape heterogeneity','Dominant land use')
units <- c(richness='Species count',temperature='°C · March–June',
           precipitation='mm · March–June',heterogeneity='Simpson diversity',landuse='')
mainland <- c(xmin=7.95,xmax=13,ymin=54.4,ymax=57.85)
bornholm <- c(xmin=14.65,xmax=15.2,ymin=54.95,ymax=55.32)
limits <- lapply(d[variables[1:4]],range)
model_groups <- levels(d$landuse)
point_map <- function(a,v,bbox=mainland,border=FALSE) {
  z <- d[as.character(d$atlas)==as.character(a),,drop=FALSE]
  p <- ggplot(z,aes(lon,lat,fill=.data[[v]])) +
    geom_point(shape=21,colour='white',stroke=.08,size=1.12) +
    coord_quickmap(xlim=bbox[c('xmin','xmax')],ylim=bbox[c('ymin','ymax')],expand=FALSE) +
    theme_void(base_size=10) + theme(legend.position='none',
      plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(0,0,0,0),
      panel.border=if(border) element_rect(fill=NA,colour='grey45',linewidth=.4) else element_blank())
  if(v=='landuse') p + scale_fill_manual(values=palette,breaks=model_groups,drop=FALSE,name=NULL)
  else p + scale_fill_gradientn(colours=continuous_palettes[[v]],
    limits=limits[[v]],name=units[[v]])
}
strip <- function(label) ggdraw() + draw_label(label,fontface='bold',size=12) +
  theme(plot.background=element_rect(fill='grey92',colour='white'))
rows <- lapply(1:3,function(a) {
  maps <- lapply(variables,function(v) ggdraw(point_map(a,v)) +
    draw_plot(point_map(a,v,bornholm,TRUE),x=.79,y=.68,width=.18,height=.22))
  label <- ggdraw() + draw_label(paste0('Atlas ',a,'\n',c('1970s','1990s','2010s')[a]),
                                fontface='bold',size=13,lineheight=1.2)
  plot_grid(plotlist=c(list(label),maps),nrow=1,rel_widths=c(.35,rep(1,5)))
})
legends <- lapply(variables,function(v) {
  p <- point_map(1,v) + theme(legend.position='bottom',
    legend.title=element_text(size=10),legend.text=element_text(size=9),
    legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0))
  if(v=='landuse') p <- p + guides(fill=guide_legend(ncol=2,byrow=TRUE,
    override.aes=list(size=3),keywidth=unit(.3,'cm'),keyheight=unit(.4,'cm'))) +
      theme(legend.text=element_text(size=8),legend.spacing.x=unit(.08,'cm'))
  else p <- p + guides(fill=guide_colourbar(title.position='top',title.hjust=.5,
    barwidth=unit(4,'cm'),barheight=unit(.25,'cm')))
  get_legend(p)
})
header <- plot_grid(plotlist=c(list(ggdraw()),lapply(titles,strip)),nrow=1,
                    rel_widths=c(.35,rep(1,5)))
legend_row <- plot_grid(plotlist=c(list(ggdraw()),legends),nrow=1,
                        rel_widths=c(.35,rep(1,5)))
foot <- ggdraw() + draw_label(sprintf(
  '%s shared cells · %s model species · common scales across atlases · Bornholm insets',
  format(length(shared),big.mark=','),ncol(Y)),x=.08,hjust=0,size=10)
figure <- plot_grid(plotlist=c(list(header),rows,list(legend_row,foot)),ncol=1,
                    rel_heights=c(.11,1,1,1,.37,.09))
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(out_dir,'atlas-richness-environment-panel.png'),figure,
       width=18,height=12.6,dpi=300,bg='white')
