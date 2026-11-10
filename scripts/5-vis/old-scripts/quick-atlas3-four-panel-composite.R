# Four-panel Atlas 3 composite; PNG only. Run from repository root.
suppressPackageStartupMessages({library(sf); library(ggplot2); library(cowplot)})
e <- new.env()
load('intermediate-outputs/hmsc-output/2026-07-27_13-28-41_Atlas3_MinOccs5_CoverageGoodAverage/Models/Fitted/HPC_samples_0250_thin_100_chains_4.Rdata', e)
m <- e$fitted_model$posteriors
sites <- as.character(m$studyDesign$site)
stopifnot(identical(rownames(m$Y),rownames(m$studyDesign)), identical(rownames(m$XData),rownames(m$studyDesign)))
grid <- st_transform(st_read('data-new/derived/atlas-grid/dof-5km-grid/DK5km_ED50grid_approx_kvadrkod_DOF.shp',quiet=TRUE),25832)
stopifnot(!anyDuplicated(sites), !anyDuplicated(grid$kvadratkod), all(sites %in% grid$kvadratkod))
grid <- grid[match(sites,grid$kvadratkod),]
grid$presence <- factor(m$Y[,'Dryocopus_martius'],levels=c(0,1),labels=c('Not recorded','Recorded'))
grid$richness <- rowSums(m$Y)
grid$forest <- m$XScaled[,'perc_forest']
grid$temperature <- m$XScaled[,'tmean_breeding']
stopifnot(!anyNA(st_drop_geometry(grid)[,c('presence','richness','forest','temperature')]), all(m$Y %in% c(0,1)), abs(mean(grid$forest))<1e-8, abs(sd(grid$forest)-1)<1e-8, abs(mean(grid$temperature))<1e-8, abs(sd(grid$temperature)-1)<1e-8)
map <- function(variable,scale,inset=FALSE) {
 ggplot(grid)+geom_sf(aes(fill=.data[[variable]]),colour=NA)+scale+
 coord_sf(xlim=if(inset)c(855000,905000) else c(400000,750000),ylim=if(inset)c(6100000,6160000) else c(6000000,6450000),expand=FALSE)+
 theme_void(base_size=12)+theme(legend.position='none',plot.background=element_rect(fill='white',colour=NA))
}
relative_scale <- function() scale_fill_gradient2(low='#2878AC',mid='#F3F1EA',high='#B94936',midpoint=0,limits=c(-3,3),oob=scales::squish,breaks=c(-3,0,3),labels=c('-3 or less','Mean','+3 or more'),name='SD from Atlas 3 mean')
panel <- function(variable,title,subtitle,scale) {
 main <- map(variable,scale)+labs(title=title,subtitle=subtitle)+theme(plot.title=element_text(size=17,face='bold'),plot.subtitle=element_text(size=11,colour='#505050'),plot.margin=margin(10,15,0,15))
 inset <- map(variable,scale,TRUE)+theme(panel.border=element_rect(fill=NA,colour='#999999',linewidth=.4))
 body <- ggdraw(main)+draw_plot(inset,x=.73,y=.64,width=.17,height=.19)
 legend <- get_legend(map(variable,scale)+theme(legend.position='bottom',legend.title=element_text(size=11),legend.text=element_text(size=10))+guides(fill=if(variable=='presence')guide_legend(title=NULL) else guide_colourbar(title.position='top',title.hjust=.5,barwidth=unit(6,'cm'),barheight=unit(.3,'cm'))))
 plot_grid(body,legend,ncol=1,rel_heights=c(1,.13))
}
panels <- list(
 panel('presence','A  Black woodpecker','Dryocopus martius | recorded distribution',scale_fill_manual(values=c('Not recorded'='#E6E6E6','Recorded'='#285B47'))),
 panel('richness','B  Species richness','Recorded species per grid cell',scale_fill_viridis_c(option='C',name='Number of species')),
 panel('forest','C  Relative forest cover','Percentage forest cover, standardized',relative_scale()),
 panel('temperature','D  Relative temperature','Breeding-season mean temperature, standardized',relative_scale()))
body <- plot_grid(plotlist=panels,ncol=2)
header <- ggdraw()+draw_label('Atlas 3 | Distributions, richness and environment',fontface='bold',size=22)
footer <- ggdraw()+draw_label(sprintf('%s model cells • Richness across %s modelled species • Bornholm shown in insets',format(nrow(grid),big.mark=','),ncol(m$Y)),size=11)
p <- plot_grid(header,body,footer,ncol=1,rel_heights=c(.055,1,.035))
out <- 'manuscript-outputs/figures/exploratory/2026-10-13-atlas3-four-panel-composite.png'
ggsave(out,p,width=12,height=13,dpi=220,bg='white')
message('Saved: ',out)
print(data.frame(cells=nrow(grid),species=ncol(m$Y),woodpecker_records=sum(grid$presence=='Recorded'),richness_min=min(grid$richness),richness_max=max(grid$richness)))
