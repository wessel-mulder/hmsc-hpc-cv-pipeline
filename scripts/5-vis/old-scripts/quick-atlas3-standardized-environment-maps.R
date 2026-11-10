# Map the exact standardized predictor matrix stored in the fitted Atlas 3 model.
suppressPackageStartupMessages({library(sf); library(ggplot2); library(cowplot)})
fit_path <- 'intermediate-outputs/hmsc-output/2026-07-27_13-28-41_Atlas3_MinOccs5_CoverageGoodAverage/Models/Fitted/HPC_samples_0250_thin_100_chains_4.Rdata'
e <- new.env(); load(fit_path, e); m <- e$fitted_model$posteriors
x <- m$XScaled[, colnames(m$XScaled) != '(Intercept)', drop = FALSE]
stopifnot(max(abs(colMeans(x))) < 1e-8, max(abs(apply(x, 2, sd)-1)) < 1e-8,
          identical(rownames(m$XData), rownames(m$studyDesign)))
sites <- as.character(m$studyDesign$site)
grid <- st_transform(st_read('data-new/derived/atlas-grid/dof-5km-grid/DK5km_ED50grid_approx_kvadrkod_DOF.shp', quiet=TRUE), 25832)
stopifnot(!anyDuplicated(sites), !anyDuplicated(grid$kvadratkod), all(sites %in% grid$kvadratkod))
grid <- grid[match(sites, grid$kvadratkod), ]
stopifnot(nrow(grid) == nrow(x), !anyNA(x))
labels <- c(tmean_breeding='Breeding-season temperature', prec_breeding='Breeding-season precipitation', hh='Land-use heterogeneity', perc_urban='Urban cover', perc_cropland='Cropland cover', perc_pasture='Pasture cover', perc_forest='Forest cover', perc_grass_shrub='Grass / shrub cover')
lims <- c(floor(min(x)), ceiling(max(x)))
scale_env <- function() scale_fill_gradient2(low='#2166AC', mid='#F7F7F7', high='#B2182B', midpoint=0, limits=lims, breaks=c(-6,-3,0,3,6,8), name='Standardized value (SD from Atlas 3 mean)')
map <- function(v, inset=FALSE) {
  grid$value <- x[,v]
  ggplot(grid) + geom_sf(aes(fill=value), colour=NA) + scale_env() +
    coord_sf(xlim=if(inset)c(855000,905000) else c(400000,750000), ylim=if(inset)c(6100000,6160000) else c(6000000,6450000), expand=FALSE) +
    theme_void() + theme(legend.position='none', plot.background=element_rect(fill='white',colour=NA))
}
panels <- lapply(colnames(x), function(v) {
  main <- map(v) + labs(title=labels[[v]]) + theme(plot.title=element_text(size=12, face='bold', hjust=.5), plot.margin=margin(8,4,4,4))
  inset <- map(v,TRUE) + theme(panel.border=element_rect(fill=NA,colour='grey50',linewidth=.35))
  ggdraw(main) + draw_plot(inset,x=.70,y=.70,width=.145,height=.15)
})
legend <- get_legend(map(colnames(x)[1]) + theme(legend.position='bottom',legend.title=element_text(size=11)) + guides(fill=guide_colourbar(title.position='top',title.hjust=.5,barwidth=unit(16,'cm'),barheight=unit(.35,'cm'))))
body <- plot_grid(plotlist=panels,ncol=4)
header <- ggdraw() + draw_label('Atlas 3 | Standardized environmental predictors',fontface='bold',size=19)
footer <- ggdraw() + draw_label('1,904 model cells • Exact fitted-model XScaled values • Bornholm shown in each inset',size=11)
p <- plot_grid(header,body,legend,footer,ncol=1,rel_heights=c(.07,1,.10,.045))
out <- 'manuscript-outputs/figures/exploratory/2026-09-09-atlas3-standardized-environment.png'
dir.create(dirname(out), recursive=TRUE,showWarnings=FALSE)
ggsave(out,p,width=15,height=11,dpi=200,bg='white')
message('Saved: ',out)
print(data.frame(variable=colnames(x),mean=colMeans(x),sd=apply(x,2,sd),min=apply(x,2,min),max=apply(x,2,max)))
