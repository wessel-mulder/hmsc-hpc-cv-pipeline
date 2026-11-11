# Run from repository root:
# Rscript scripts/5-vis/env-viz/atlas-chiffchaff-richness-panels.R
# Recorded Chiffchaff and richness (A–B); all species; PNG only.
suppressPackageStartupMessages({library(sf); library(ggplot2); library(cowplot)})
out_dir <- 'scripts/5-vis/env-viz'
# Reuse the exact land/effort/common-cell selection from 1-atlas-analyses.
# This loader performs no minimum-occurrence species filtering and writes nothing.
input <- 'scripts/1-atlas-analyses/homogenization/homogenization-data.R'
atlas_data <- new.env(parent=globalenv())
sys.source(input,envir=atlas_data)
shared <- atlas_data$sites
ids <- unlist(lapply(1:3,function(a) paste0(shared,'_',a)))
d <- atlas_data$design[match(ids,atlas_data$design$survey),,drop=FALSE]
rownames(d) <- d$survey
Y <- as.matrix(atlas_data$y[ids,,drop=FALSE])
env <- atlas_data$environment[match(ids,atlas_data$environment$survey),,drop=FALSE]
# Preserve model cover definitions, but not the model input's complete-predictor
# or species-frequency exclusions. Missing climate remains explicitly missing.
X <- data.frame(tmean_breeding=env$tmean_breeding,prec_breeding=env$prec_breeding,
  hh=env$hh,perc_cropland=env$pct_class_22,perc_urban=env$pct_class_11,
  perc_pasture=env$pct_class_33,perc_grass_shrub=env$pct_class_55,
  perc_forest_deciduous=env$pct_class_43+env$pct_class_44,
  perc_forest_evergreen=env$pct_class_42+env$pct_class_41,
  perc_forest_mixed_other=env$pct_class_40+env$pct_class_45,row.names=ids)
stopifnot(identical(rownames(d),rownames(X)),identical(rownames(d),rownames(Y)),
  !anyDuplicated(d[c('site','atlas')]),!anyNA(Y),all(Y %in% 0:1),
  all(table(d$site)==3),all(table(d$atlas)==length(shared)),
  all(d$lon>7 & d$lon<16),all(d$lat>54 & d$lat<59))
# Raw studyDesign coordinates already use geographic degrees.
stopifnot('Phylloscopus_collybita' %in% colnames(Y))
d$chiffchaff <- factor(Y[,'Phylloscopus_collybita'],levels=0:1,
                       labels=c('Absent','Present'))
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
variables <- c('chiffchaff','richness','temperature','precipitation','heterogeneity','landuse')
titles <- c('Common Chiffchaff','Species richness','Temperature','Precipitation','Landscape heterogeneity','Dominant land use')
units <- c(richness='Species count',temperature='°C · March–June',
           precipitation='mm · March–June',heterogeneity='Simpson diversity',landuse='')
mainland <- c(xmin=7.95,xmax=13,ymin=54.4,ymax=57.85)
bornholm <- c(xmin=14.65,xmax=15.2,ymin=54.95,ymax=55.32)
limits <- lapply(d[c('richness','temperature','precipitation','heterogeneity')],range,na.rm=TRUE)
model_groups <- levels(d$landuse)
point_map <- function(a,v,bbox=mainland,border=FALSE) {
  z <- d[as.character(d$atlas)==as.character(a),,drop=FALSE]
  p <- ggplot(z,aes(lon,lat,fill=.data[[v]])) +
    geom_point(shape=21,colour='white',stroke=.08,size=1.12) +
    coord_quickmap(xlim=bbox[c('xmin','xmax')],ylim=bbox[c('ymin','ymax')],expand=FALSE) +
    theme_void(base_size=10) + theme(legend.position='none',
      plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(0,0,0,0),
      panel.border=if(border) element_rect(fill=NA,colour='grey45',linewidth=.4) else element_blank())
  if(v=='chiffchaff') p + scale_fill_manual(
    values=c(Absent='#D4D4D4',Present='#238443'),drop=FALSE,name=NULL)
  else if(v=='landuse') p + scale_fill_manual(values=palette,breaks=model_groups,drop=FALSE,name=NULL)
  else p + scale_fill_gradientn(colours=continuous_palettes[[v]],
    limits=limits[[v]],name=units[[v]])
}

heading <- function(label,size=13) ggdraw() + draw_label(label,x=.015,hjust=0,
  fontface='bold',size=size) + theme(plot.background=element_rect(fill='grey92',colour='white'))
inset_map <- function(main,inset) ggdraw(main) + draw_plot(inset,x=.72,y=.70,width=.17,height=.21)
# A–B share one atlas header row. Rotated titles and right-hand legends
# allow the maps to occupy the full height without repeated panel strips.
status_row <- function(v,title) {
  panels <- lapply(1:3,function(a) inset_map(
    point_map(a,v) + geom_point(aes(fill=.data[[v]]),shape=21,
      colour='white',stroke=.1,size=1.55),
    point_map(a,v,bornholm,TRUE)))
  lp <- point_map(1,v) + theme(legend.position='right',
    legend.text=element_text(size=10),legend.title=element_text(size=10),
    legend.margin=margin(0,0,0,0))
  if(v=='chiffchaff') lp <- lp + guides(fill=guide_legend(ncol=1,
    override.aes=list(size=3)))
  else lp <- lp + guides(fill=guide_colourbar(direction='vertical',
    title.position='top',barheight=unit(4,'cm'),barwidth=unit(.3,'cm')))
  label <- ggdraw() + draw_label(title,angle=90,fontface='bold',size=13)
  plot_grid(plotlist=c(list(label),panels,list(get_legend(lp))),nrow=1,
            rel_widths=c(.14,1,1,1,.34))
}
atlas_header <- plot_grid(plotlist=c(list(ggdraw()),lapply(1:3,function(a)
  heading(paste0('Atlas ',a,' · ',c('1970s','1990s','2010s')[a]),12)),list(ggdraw())),
  nrow=1,rel_widths=c(.14,1,1,1,.34))
status_foot <- ggdraw() + draw_label(sprintf('%s shared cells · All %s supplied species · Bornholm insets',
  format(length(shared),big.mark=','),ncol(Y)),size=9)
status_figure <- plot_grid(atlas_header,
  status_row('chiffchaff','A  Common Chiffchaff'),
  status_row('richness','B  Species richness'),status_foot,
  ncol=1,rel_heights=c(.07,1,1,.06))
ggsave(file.path(out_dir,'atlas-chiffchaff-richness-panels.png'),status_figure,
       width=12.5,height=8.6,dpi=300,bg='white')
