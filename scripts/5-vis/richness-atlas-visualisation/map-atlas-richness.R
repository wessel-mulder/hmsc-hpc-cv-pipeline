# Run from the repository root:
# Rscript scripts/5-vis/richness-atlas-visualisation/map-atlas-richness.R
# Recorded richness of the All-species training set, not predicted richness.
# Read the same portable S0 bundle as S1; retain each atlas's training rows.
suppressPackageStartupMessages({library(sf); library(ggplot2); library(cowplot)})
input <- 'data-shared/preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData'
out_dir <- 'scripts/5-vis/richness-atlas-visualisation'
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
e <- new.env(parent=emptyenv())
load(input, envir=e)
Y <- as.matrix(e$Y)
d <- e$Design
species <- c(chiffchaff='Phylloscopus_collybita', willow_warbler='Phylloscopus_trochilus')
species_labels <- c(chiffchaff='Chiffchaff', willow_warbler='Willow warbler')
stopifnot(identical(rownames(Y), rownames(d)), !anyDuplicated(rownames(Y)),
          !anyDuplicated(colnames(Y)), !anyNA(Y), all(Y %in% 0:1),
          all(species %in% colnames(Y)), !anyDuplicated(d[c('site','atlas')]))
# S0's lon/lat columns are ED50 / UTM 32N metres, NOT geographic degrees.
points <- st_as_sf(d, coords=c('lon','lat'), crs=23032)
xy <- st_coordinates(st_transform(points, 4326))
d$lon <- xy[,1]; d$lat <- xy[,2]
d$richness <- rowSums(Y)
for (key in names(species)) {
  d[[key]] <- factor(Y[,species[[key]]], levels=0:1,
                    labels=c('Not recorded','Recorded'))
}
stopifnot(all(is.finite(d$lon)), all(is.finite(d$lat)),
          all(d$lon > 7 & d$lon < 16), all(d$lat > 54 & d$lat < 59))
limits <- range(d$richness)
mainland <- c(xmin=7.95,xmax=13,ymin=54.4,ymax=57.85)
bornholm <- c(xmin=14.65,xmax=15.2,ymin=54.95,ymax=55.32)

point_map <- function(atlas, variable, bbox=mainland, border=FALSE) {
  z <- d[as.character(d$atlas)==as.character(atlas),,drop=FALSE]
  p <- ggplot(z,aes(lon,lat,fill=.data[[variable]])) +
    geom_point(shape=21,colour='white',stroke=0.12,size=1.8) +
    coord_quickmap(xlim=bbox[c('xmin','xmax')],
                   ylim=bbox[c('ymin','ymax')],expand=FALSE) +
    theme_void(base_size=11) +
    theme(legend.position='none',plot.background=element_rect(fill='white',colour=NA),
          panel.border=if(border) element_rect(colour='grey35',fill=NA,linewidth=.55)
            else element_blank(),plot.margin=margin(0,0,0,0))
  if(variable=='richness') {
    p + scale_fill_gradientn(colours=c('#2166AC','#FFFFBF','#B2182B'),
                            limits=limits,name='Recorded\nspecies count')
  } else {
    p + scale_fill_manual(values=c('Not recorded'='#D4D4D4','Recorded'='#285B47'),
                          drop=FALSE,name=species_labels[[variable]])
  }
}
panel <- function(atlas,variable,label) {
  main <- point_map(atlas,variable)
  inset <- point_map(atlas,variable,bornholm,TRUE)
  map <- ggdraw(main) + draw_plot(inset,x=.79,y=.68,width=.18,height=.22)
  strip <- ggdraw() + draw_label(label,fontface='bold',size=11) +
    theme(plot.background=element_rect(fill='grey92',colour='white'))
  plot_grid(strip,map,ncol=1,rel_heights=c(.07,1))
}
legend_for <- function(variable) get_legend(
  point_map(1,variable) + theme(legend.position='right',
    legend.title=element_text(face='bold'),legend.text=element_text(size=10)))
save_map <- function(panels,variable,title,subtitle,filename,width) {
  body <- plot_grid(plot_grid(plotlist=panels,nrow=1),legend_for(variable),
                    nrow=1,rel_widths=c(1,if(length(panels)==1) .23 else .13))
  heading <- ggdraw() +
    draw_label(title,x=.02,y=.77,hjust=0,fontface='bold',size=14) +
    draw_label(subtitle,x=.02,y=.25,hjust=0,size=9)
  p <- plot_grid(heading,body,ncol=1,rel_heights=c(.13,1))
  ggsave(file.path(out_dir,filename),p,width=width,height=5.7,dpi=300,bg='white')
  if(interactive()) print(p)
}
n <- table(d$atlas)
save_map(list(panel(1,'richness','Atlas 1 — 1970s')),'richness',
  'Recorded species richness — Atlas 1',
  sprintf('%s training cells; %s model species; Bornholm shown in inset',n[['1']],ncol(Y)),
  'atlas1-species-richness.png',7)
for (key in names(species)) {
  recorded <- sum(Y[as.character(d$atlas)=='1',species[[key]]])
  save_map(list(panel(1,key,'Atlas 1 — 1970s')),key,
    paste0(species_labels[[key]],' — Atlas 1'),
    sprintf('%s; recorded in %s of %s training cells',
            gsub('_',' ',species[[key]]),recorded,n[['1']]),
    paste0('atlas1-',gsub('_','-',key),'-distribution.png'),7)
  message(species_labels[[key]],': ',recorded,' / ',n[['1']],' Atlas 1 cells')
}
save_map(lapply(1:3,function(a) panel(a,'richness',
  sprintf('Atlas %s — %s (%s cells)',a,c('1970s','1990s','2010s')[a],n[[as.character(a)]]))),
  'richness','Recorded species richness — model-training atlas data',
  sprintf('%s model species; common count scale; each atlas retains its own training cells',ncol(Y)),
  'atlas1-2-3-species-richness.png',13)
print(aggregate(richness~atlas,d,function(x) c(cells=length(x),min=min(x),max=max(x))))
message('Input: ',input,'\nSaved PNGs in: ',out_dir)
