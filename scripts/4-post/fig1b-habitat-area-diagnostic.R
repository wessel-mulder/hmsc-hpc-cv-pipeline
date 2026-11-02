# Descriptive diagnostic of Figure 1B, restricted to exact fitted-model sites.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(readxl);library(sf);library(cowplot)})
e <- new.env(); load('data-new/derived/hmsc-inputs/preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData', e)
grid <- st_drop_geometry(st_read('data-new/derived/atlas-grid/grids-ocean-thresholds/grids_ocean_thresholds.shp',quiet=TRUE))
vars <- c(perc_urban='Urban',perc_cropland='Cropland',perc_pasture='Pasture',perc_forest='Forest',perc_grass_shrub='Grass/shrubland')
cols <- c(Urban='grey55',Cropland='goldenrod',Pasture='darkorange',Forest='#008b45','Grass/shrubland'='#23c776')
periods <- c('1970s','1990s','2010s')
folders <- list.dirs('intermediate-outputs/hmsc-output',recursive=FALSE,full.names=TRUE)
folders <- folders[grepl('2026-07-27.*Atlas[123]_MinOccs5_CoverageGoodAverage$',folders)]
stopifnot(length(folders)==3)
results <- list(); coverage <- list(); cell_tables <- list()
for(folder in folders){
 a <- as.integer(sub('.*Atlas([123])_.*','\\1',folder)); period <- periods[a]
 f <- new.env();load(file.path(folder,'Models/Fitted/HPC_samples_0250_thin_100_chains_4.Rdata'),f)
 m <- f$fitted_model$posteriors
 surveys <- rownames(m$XData); sites <- as.character(m$studyDesign$site)
 stopifnot(length(surveys)==length(sites),!anyDuplicated(sites),all(surveys %in% rownames(e$X)))
 x <- e$X[surveys,names(vars),drop=FALSE]
 stopifnot(!anyNA(x), all(as.matrix(x)>=0 & as.matrix(x)<=1))
 g <- grid[match(sites,grid$kvdrtkd),];stopifnot(!anyNA(g$grid_ar),!anyNA(g$lnd_r_s))
 nsp <- ncol(m$Y)
 cells <- data.frame(survey=surveys,site=sites,period=period,grid_km2=g$grid_ar/1e6,land_km2=g$lnd_r_s/1e6,x)
 cell_tables[[period]] <- cells
 coverage[[period]] <- tibble(period=period,n_cells=nrow(cells),n_species=nsp,grid_km2=sum(cells$grid_km2),land_km2=sum(cells$land_km2))
 path <- file.path(folder,'Results',paste0(basename(folder),'parameter_estimates_Beta_.xlsx'))
 means <- read_excel(path,sheet='Posterior mean'); probs <- read_excel(path,sheet='Pr(x>0)')
 stopifnot(nrow(means)==nsp,setequal(means$Species,colnames(m$Y)),identical(means$Species,probs$Species))
 results[[period]] <- bind_rows(lapply(names(vars),function(v){
 supported <- probs[[v]]>=.95 | probs[[v]]<=.05
 tibble(period=period,variable=vars[[v]],direction=c('Positive','Negative'),
   supported_species=c(sum(supported & means[[v]]>0),sum(supported & means[[v]]<0)),
   percent_species=100*supported_species/nsp,
   mean_cover_percent=100*mean(x[[v]]),sd_cover_pp=100*sd(x[[v]]),
   estimated_habitat_km2=sum(x[[v]]*cells$grid_km2))
 }))
 rm(f,m);gc(FALSE)
}
dat <- bind_rows(results); area <- bind_rows(coverage)
dat$period <- factor(dat$period,levels=periods)
dat$direction <- factor(dat$direction,levels=c('Positive','Negative'))
long <- dat |> pivot_longer(c(mean_cover_percent,estimated_habitat_km2,sd_cover_pp),names_to='metric',values_to='availability') |>
 mutate(metric=factor(metric,levels=c('mean_cover_percent','estimated_habitat_km2','sd_cover_pp'),labels=c('Mean cover (%)','Estimated habitat area (km²)','SD of cover (percentage points)')))
p <- ggplot(long,aes(availability,percent_species,colour=variable,group=variable))+
 geom_path(linewidth=.65,alpha=.65)+geom_point(aes(alpha=period),shape=16,size=3.5)+
 facet_grid(metric~direction,scales='free_x')+
 scale_colour_manual(values=cols)+scale_alpha_manual(values=c('1970s'=.4,'1990s'=.7,'2010s'=1))+
 scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25))+
 labs(x=NULL,y='Species with supported associations (%)',colour='Land cover',alpha='Atlas period')+
 guides(alpha=guide_legend(override.aes=list(colour='black')),colour=guide_legend(override.aes=list(alpha=1)))+
 theme_bw(base_size=12)+theme(legend.position='bottom',strip.text=element_text(face='bold'),panel.grid.minor=element_blank())
# facet_grid shares x scales within columns, so use wrap to give each metric its own range.
p <- p + facet_wrap(vars(metric,direction),ncol=2,scales='free_x',dir='h')
area_text <- paste(sprintf('%s: %s cells | %.0f km² grid area | %.0f km² land',area$period,format(area$n_cells,big.mark=','),area$grid_km2,area$land_km2),collapse='\n')
header <- ggdraw()+draw_label('Does habitat availability track species–environment associations?',size=17,fontface='bold')
foot <- ggdraw()+draw_label(paste0(area_text,'\nExact final-model cells; all ',unique(area$n_species),' modelled species form each percentage denominator.\nHabitat area estimated as cover fraction × full cell area. Lines join periods; no trend model fitted.'),size=10,hjust=.5)
fig <- plot_grid(header,p,foot,ncol=1,rel_heights=c(.06,1,.15))
out <- 'manuscript-outputs/figures/exploratory';dir.create(out,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(out,'fig1b-habitat-area-diagnostic.png'),fig,width=13,height=12,dpi=300,bg='white')
write.csv(dat,file.path(out,'fig1b-habitat-area-diagnostic.csv'),row.names=FALSE)
write.csv(area,file.path(out,'fig1b-modelled-area-by-atlas.csv'),row.names=FALSE)
write.csv(bind_rows(cell_tables),file.path(out,'fig1b-habitat-area-model-cells.csv'),row.names=FALSE)
print(area); print(dat)
