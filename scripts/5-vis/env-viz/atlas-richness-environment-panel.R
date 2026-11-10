# Run from repository root:
# Rscript scripts/5-vis/env-viz/atlas-richness-environment-panel.R
# Five variables (columns) across three atlases (rows); PNG only.
suppressPackageStartupMessages({library(sf); library(ggplot2); library(cowplot)})
out_dir <- 'scripts/5-vis/env-viz'
input <- 'data-shared/preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData'
environment_input <- 'data-new/0-pre/derived/environmental/X_environment.csv'
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
# Include all terrestrial HILDA classes, rather than only the selected model
# predictors. Ocean and freshwater are excluded; all forest types are combined.
env <- read.csv(environment_input,check.names=FALSE)
stopifnot(!anyDuplicated(env$survey))
ix <- match(rownames(d),env$survey); stopifnot(!anyNA(ix))
env <- env[ix,,drop=FALSE]
stopifnot(isTRUE(all.equal(d$temperature,env$tmean_breeding)),
          isTRUE(all.equal(d$precipitation,env$prec_breeding)),
          isTRUE(all.equal(d$heterogeneity,env$hh)))
cover <- data.frame(Cropland=env$pct_class_22,Forest=rowSums(env[paste0('pct_class_',40:45)]),
  Pasture=env$pct_class_33,Urban=env$pct_class_11,
  'Grass / shrub'=env$pct_class_55,'Tree crops'=env$pct_class_23,
  Agroforestry=env$pct_class_24,'Sparse / bare'=env$pct_class_66,check.names=FALSE)
stopifnot(all(is.finite(as.matrix(cover))),all(rowSums(cover)>0))
# Ties retain the first group in the documented order above.
d$landuse <- factor(names(cover)[max.col(as.matrix(cover),ties.method='first')],
                    levels=names(cover))
palette <- c(Cropland='#D5AD48',Forest='#28734B',Pasture='#99BF69',Urban='#BA5762',
             'Grass / shrub'='#7CAEA3','Tree crops'='#BE7F32',
             Agroforestry='#8463A7','Sparse / bare'='#A7A39B')
variables <- c('richness','temperature','precipitation','heterogeneity','landuse')
titles <- c('Species richness','Temperature','Precipitation','Heterogeneity','Dominant land use')
units <- c(richness='Species count',temperature='°C · March–June',
           precipitation='mm · March–June',heterogeneity='Simpson diversity',landuse='')
mainland <- c(xmin=7.95,xmax=13,ymin=54.4,ymax=57.85)
bornholm <- c(xmin=14.65,xmax=15.2,ymin=54.95,ymax=55.32)
limits <- lapply(d[variables[1:4]],range)
observed_groups <- levels(droplevels(d$landuse))
point_map <- function(a,v,bbox=mainland,border=FALSE) {
  z <- d[as.character(d$atlas)==as.character(a),,drop=FALSE]
  p <- ggplot(z,aes(lon,lat,fill=.data[[v]])) +
    geom_point(shape=21,colour='white',stroke=.08,size=1.12) +
    coord_quickmap(xlim=bbox[c('xmin','xmax')],ylim=bbox[c('ymin','ymax')],expand=FALSE) +
    theme_void(base_size=10) + theme(legend.position='none',
      plot.background=element_rect(fill='white',colour=NA),plot.margin=margin(0,0,0,0),
      panel.border=if(border) element_rect(fill=NA,colour='grey45',linewidth=.4) else element_blank())
  if(v=='landuse') p + scale_fill_manual(values=palette,breaks=observed_groups,drop=FALSE,name=NULL)
  else p + scale_fill_viridis_c(option=if(v=='temperature') 'C' else 'D',
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
    override.aes=list(size=3),keywidth=unit(.35,'cm'),keyheight=unit(.4,'cm')))
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
                    rel_heights=c(.11,1,1,1,.29,.09))
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(out_dir,'atlas-richness-environment-panel.png'),figure,
       width=18,height=12.6,dpi=300,bg='white')
writeLines(c('# Atlas richness and environment panel',
  paste('Run: `Rscript scripts/5-vis/env-viz/atlas-richness-environment-panel.R`'),
  paste('Training input:',input),paste('Land-cover input:',environment_input),
  sprintf('Same %s cells in each atlas; recorded richness of the same %s model species.',length(shared),ncol(Y)),
  'Five columns correspond to the five requested variables; three rows represent the atlas periods.',
  'Temperature and precipitation summarize March–June in original units. Each continuous variable has one scale across all three atlases.',
  'Heterogeneity is mean annual Simpson diversity across all 15 original HILDA classes, including ocean and freshwater.',
  'Dominant land use is the largest period-mean terrestrial cover group; ocean and freshwater excluded. Forest combines classes 40–45.',
  'Groups: cropland (22), forest (40–45), pasture (33), urban (11), grass/shrub (55), tree crops (23), agroforestry (24), sparse/bare (66).',
  'Exact ties use the first group in that order. Dominance need not exceed 50% of a cell. Legend shows groups dominant in at least one cell.',
  'Climate years: 1971–1974, 1993–1996, 2014–2016; HILDA years: 1971–1974, 1993–1996, 2014–2017.',
  capture.output(table(d$atlas,d$landuse))),
  file.path(out_dir,'atlas-richness-environment-panel.md'))
message('Saved panel: ',length(shared),' cells per atlas; ',ncol(Y),' species')
print(table(d$atlas,d$landuse))
