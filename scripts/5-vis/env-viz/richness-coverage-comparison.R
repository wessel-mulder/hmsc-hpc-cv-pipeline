# Run from repository root: Rscript scripts/5-vis/env-viz/richness-coverage-comparison.R
# All supplied species; three nested coverage selections; PNG only.
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
out_dir <- 'scripts/5-vis/env-viz'
e <- new.env(parent=globalenv())
sys.source('scripts/1-atlas-analyses/homogenization/homogenization-data.R',envir=e)
# Land threshold in EVERY atlas, before any coverage restriction.
land_sites <- sort(Reduce(intersect,e$eligible_sites))
good1 <- as.character(e$effort1$KVADRATKOD[e$effort1[['UNDS_KODE_']]==1 & !is.na(e$effort1[['UNDS_KODE_']])])
good2 <- as.character(e$effort2$KVADRATKOD[e$effort2[['UNDS_KOD_1']]==1 & !is.na(e$effort2[['UNDS_KOD_1']])])
selections <- list(land_only=land_sites,good_fair=e$sites,
                  good_only=sort(Reduce(intersect,list(land_sites,good1,good2))))
stopifnot(all(selections$good_only %in% selections$good_fair),
          all(selections$good_fair %in% selections$land_only))
ids <- unlist(lapply(1:3,function(a)paste0(land_sites,'_',a)))
stopifnot(all(ids %in% rownames(e$y)))
d <- e$design[match(ids,e$design$survey),,drop=FALSE]
y <- as.matrix(e$y[ids,,drop=FALSE])
stopifnot(!anyNA(y),all(y %in% 0:1),!anyDuplicated(d[c('site','atlas')]),
          all(table(d$site)==3),all(table(d$atlas)==length(land_sites)))
d$richness <- rowSums(y)
limits <- range(d$richness)
colours <- c('#6A00A8','#2166AC','#00BFC4','#33A02C','#FFFF33','#FF8C00','#D7191C')
mainland <- c(xmin=7.95,xmax=13,ymin=54.4,ymax=57.85)
bornholm <- c(xmin=14.65,xmax=15.2,ymin=54.95,ymax=55.32)
point_map <- function(a,sites,bbox=mainland,border=FALSE) {
  z <- d[d$atlas==a & d$site %in% sites,,drop=FALSE]
  ggplot(z,aes(lon,lat,fill=richness)) +
    geom_point(shape=21,colour='white',stroke=.1,size=1.45) +
    scale_fill_gradientn(colours=colours,limits=limits,name='Species\ncount') +
    coord_quickmap(xlim=bbox[c('xmin','xmax')],ylim=bbox[c('ymin','ymax')],expand=FALSE) +
    theme_void(base_size=11) + theme(legend.position='none',plot.margin=margin(0,0,0,0),
      plot.background=element_rect(fill='white',colour=NA),
      panel.border=if(border)element_rect(fill=NA,colour='grey45',linewidth=.4)else element_blank())
}
row_titles <- c(land_only='A  All coverage levels',good_fair='B  Good + Fair coverage',good_only='C  Good coverage only')
rows <- lapply(names(selections),function(key) {
  sites <- selections[[key]]
  maps <- lapply(1:3,function(a)ggdraw(point_map(a,sites)) +
    draw_plot(point_map(a,sites,bornholm,TRUE),x=.72,y=.70,width=.17,height=.21))
  title <- ggdraw() + draw_label(sprintf('%s · %s cells',row_titles[[key]],
    format(length(sites),big.mark=',')),angle=90,fontface='bold',size=11)
  plot_grid(plotlist=c(list(title),maps),nrow=1,rel_widths=c(.15,1,1,1))
})
header <- plot_grid(plotlist=c(list(ggdraw()),lapply(1:3,function(a)ggdraw() +
  draw_label(paste0('Atlas ',a,' · ',c('1970s','1990s','2010s')[a]),fontface='bold',size=12) +
  theme(plot.background=element_rect(fill='grey92',colour='white')))),
  nrow=1,rel_widths=c(.15,1,1,1))
body <- plot_grid(plotlist=c(list(header),rows),ncol=1,rel_heights=c(.07,1,1,1))
legend <- get_legend(point_map(1,land_sites) + theme(legend.position='right') +
  guides(fill=guide_colourbar(direction='vertical',barheight=unit(5,'cm'),barwidth=unit(.35,'cm'))))
foot <- ggdraw() + draw_label(sprintf(
  'All %s species · ≥25%% land in every atlas · Coverage restricted jointly in Atlases 1 and 2 · Same cells across atlases within each row',ncol(y)),size=9)
figure <- plot_grid(plot_grid(body,legend,nrow=1,rel_widths=c(1,.075)),foot,
                    ncol=1,rel_heights=c(1,.025))
ggsave(file.path(out_dir,'richness-coverage-comparison.png'),figure,
       width=13,height=12.8,dpi=300,bg='white')
summary <- do.call(rbind,lapply(names(selections),function(key)do.call(rbind,lapply(1:3,function(a) {
  z <- d[d$atlas==a & d$site %in% selections[[key]],]
  data.frame(selection=key,atlas=a,cells=nrow(z),species=ncol(y),mean_richness=mean(z$richness),
             median_richness=median(z$richness),min_richness=min(z$richness),max_richness=max(z$richness))
}))))
stopifnot(identical(selections$good_fair,e$sites),
          all(summary$cells==rep(lengths(selections),each=3)))
