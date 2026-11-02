#!/usr/bin/env Rscript
# Rerun cached observed-response GWR with forest and cropland separately.
# Uses shared sites and the full fixed species list; recomputes communities. PNG only.
suppressPackageStartupMessages({library(sf); library(sp); library(GWmodel); library(ggplot2); library(cowplot); library(vegan)})
source_path <- '.codex-work/exploratory-figures/gwr-raw-observed/gwr/models/2026-03-13-observed-gwr-analysis-inputs-100km.rds'
# Raw inputs below supersede the earlier HMSC-derived cache.
out <- 'intermediate-outputs/gwr-output/observed-single-predictor-shared-full-species-2026-10-22'
fig <- 'manuscript-outputs/figures/exploratory/observed-single-predictor-gwr-shared-full-species-2026-10-22'
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(fig,recursive=TRUE,showWarnings=FALSE)
grid <- st_transform(st_read('data-new/derived/atlas-grid/dof-5km-grid/DK5km_ED50grid_approx_kvadrkod_DOF.shp',quiet=TRUE),25832)
stopifnot(!anyDuplicated(grid$kvadratkod))
periods <- c('Atlas 1','Atlas 2','Atlas 3')
manifest <- read.csv('manuscript-outputs/atlas-occurrence-change/shared-cell-manifest.csv')
cells <- sort(as.character(manifest$site))
y <- read.csv('data-new/derived/occurrences/Y_occurrences.csv',check.names=FALSE)
env <- read.csv('data-new/derived/environmental/X_Environmental.csv')
design <- read.csv('data-new/derived/design/studyDesign.csv')
species <- setdiff(names(y),'survey')
stopifnot(!anyDuplicated(y$survey),!anyDuplicated(env$survey),!anyDuplicated(design$survey))
frames <- lapply(1:3,function(a){
 ids<-paste(cells,a,sep='_')
 stopifnot(all(ids%in%y$survey),all(ids%in%env$survey),all(ids%in%design$survey))
 community<-as.matrix(y[match(ids,y$survey),species,drop=FALSE])
 stopifnot(!anyNA(community),all(community%in%c(0,1)),all(rowSums(community)>0))
 dissimilarity<-as.matrix(vegan::vegdist(community,method='bray'))
 diag(dissimilarity)<-NA_real_
 d<-design[match(ids,design$survey),]
 xy<-st_coordinates(st_transform(st_as_sf(d,coords=c('lon','lat'),crs=4326),25832))
 e<-env[match(ids,env$survey),]
 r<-data.frame(survey=ids,site=cells,X=xy[,1],Y=xy[,2],atlas=a,
   richness=rowSums(community),turnover=rowMeans(dissimilarity,na.rm=TRUE),
   perc_forest=e$LULC_44,perc_cropland=e$LULC_22)
 stopifnot(!anyNA(r),all(r$site%in%grid$kvadratkod))
 r
})
stopifnot(all(vapply(frames,function(d) identical(d$site,cells),logical(1))))
write.csv(manifest,file.path(out,'shared-cell-manifest.csv'),row.names=FALSE)
write.csv(data.frame(species=species),file.path(out,'species-manifest.csv'),row.names=FALSE)
inputs<-do.call(rbind,frames)
stopifnot(all(inputs$perc_forest>=0 & inputs$perc_forest<=1),all(inputs$perc_cropland>=0 & inputs$perc_cropland<=1))
write.csv(inputs,file.path(out,'raw-inputs.csv'),row.names=FALSE)
# Predictors in the cache are fractions. Fit in percentage points: slope is
# response units per one percentage-point increase in habitat coverage.
fits<-list(); rows<-list(); summaries<-list()
for(a in 1:3) for(v in c('forest','cropland')) for(response in c('richness','turnover')) {
 key<-paste(a,v,response,sep='_'); message('Fitting ',key)
 d<-frames[[a]];d$coverage_percent<-100*d[[paste0('perc_',v)]]
 stopifnot(sd(d$coverage_percent)>0)
 spatial<-d;coordinates(spatial)<-c('X','Y')
 fit<-GWmodel::gwr.basic(reformulate('coverage_percent',response=response),data=spatial,
                       bw=100000,kernel='bisquare',adaptive=FALSE)
 fits[[key]]<-fit;saveRDS(fit,file.path(out,paste0('fit-',key,'.rds')))
 z<-as.data.frame(fit$SDF)
 stopifnot(nrow(z)==nrow(d),all(is.finite(z$coverage_percent)))
 rows[[key]]<-data.frame(atlas=a,site=d$site,survey=d$survey,predictor=v,response=response,
   slope=z$coverage_percent,se=z$coverage_percent_SE,t=z$coverage_percent_TV,
   local_r2=z$Local_R2)
 summaries[[key]]<-data.frame(atlas=a,predictor=v,response=response,n=nrow(d),
   bandwidth_m=100000,aicc=fit$GW.diagnostic$AICc,
   slope_min=min(z$coverage_percent),slope_median=median(z$coverage_percent),
   slope_max=max(z$coverage_percent),mean_local_r2=mean(z$Local_R2))
}
local<-do.call(rbind,rows);summary<-do.call(rbind,summaries)
# One predictor SD pooled over all shared cells and atlases, held fixed through time.
scaling <- data.frame(predictor=c('forest','cropland'),
  sd_cover_pp=vapply(c('forest','cropland'),function(v) sd(100*inputs[[paste0('perc_',v)]]),numeric(1)))
local$sd_cover_pp <- scaling$sd_cover_pp[match(local$predictor,scaling$predictor)]
# Turnover is displayed on a 0–100 dissimilarity scale; differences are percentage points.
response_factor <- ifelse(local$response=='turnover',100,1)
local$effect_per_cover_pp <- local$slope*response_factor
local$effect_per_cover_sd <- local$slope*local$sd_cover_pp*response_factor
local$se_per_cover_sd <- local$se*local$sd_cover_pp*response_factor
stopifnot(isTRUE(all.equal(local$effect_per_cover_sd/local$se_per_cover_sd,
                          local$slope/local$se)))
write.csv(scaling,file.path(out,'predictor-scaling.csv'),row.names=FALSE)
write.csv(local,file.path(out,'local-coefficients.csv'),row.names=FALSE)
write.csv(summary,file.path(out,'model-summary.csv'),row.names=FALSE)
map_panel<-function(d,value,a,limits,diverging=FALSE,label='') {
 g<-grid[match(d$site,grid$kvadratkod),];g$value<-d[[value]]
 base<-ggplot(g)+geom_sf(aes(fill=value),colour=NA)+
   theme_void(base_size=11)+theme(plot.background=element_rect(fill='white',colour=NA))
 if(diverging) base<-base+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,limits=limits,name=label)
 else base<-base+scale_fill_viridis_c(limits=limits,name=label)
 main<-base+coord_sf(xlim=c(400000,750000),ylim=c(6000000,6450000),expand=FALSE)+
   labs(title=periods[a])+theme(legend.position='none',plot.title=element_text(hjust=.5))
 inset<-base+coord_sf(xlim=c(855000,905000),ylim=c(6100000,6160000),expand=FALSE)+
   theme(legend.position='none',panel.border=element_rect(fill=NA,colour='grey60',linewidth=.3))
 list(panel=ggdraw(main)+draw_plot(inset,x=.73,y=.70,width=.19,height=.2),
      legend=get_legend(base+theme(legend.position='bottom')+guides(fill=guide_colourbar(title.position='top',barwidth=unit(7,'cm')))))
}
map_three<-function(ds,value,limits,title,label,path,diverging=FALSE,subtitle=''){
 panels<-lapply(1:3,function(a) map_panel(ds[[a]],value,a,limits,diverging,label))
 header<-ggdraw()+draw_label(paste(title,subtitle,sep='\n'),size=13)
 p<-plot_grid(header,plot_grid(plotlist=lapply(panels,`[[`,'panel'),nrow=1),panels[[1]]$legend,ncol=1,rel_heights=c(.15,1,.15))
 ggsave(path,p,width=12,height=7,dpi=200,bg='white')
}
# Fixed symmetric limits shared by BOTH predictors and ALL atlases, per response/unit.
legend_scales <- do.call(rbind,lapply(c('richness','turnover'),function(response)
 do.call(rbind,lapply(c('effect_per_cover_sd','effect_per_cover_pp'),function(unit){
  lim<-max(abs(local[local$response==response,unit]))
  data.frame(response=response,unit=unit,min=-lim,max=lim)
 }))))
write.csv(legend_scales,file.path(out,'legend-scales.csv'),row.names=FALSE)
for(v in c('forest','cropland')) for(response in c('richness','turnover')){
 d<-local[local$predictor==v & local$response==response,]
 for(unit in c('effect_per_cover_sd','effect_per_cover_pp')) {
 lim<-legend_scales$max[legend_scales$response==response & legend_scales$unit==unit]
 per_sd<-unit=='effect_per_cover_sd'
 label<-paste(if(response=='richness') 'Species' else 'Turnover percentage points',
              if(per_sd) 'per +1 SD of cover' else 'per +1 percentage point of cover')
 subtitle<-if(per_sd) sprintf('1 pooled SD = %.2f percentage points of cover; shared scale across forest and cropland',unique(d$sd_cover_pp)) else 'Shared scale across forest and cropland; fixed 100-km bisquare kernel'
 map_three(lapply(1:3,function(a)d[d$atlas==a,]),unit,c(-lim,lim),
   paste(tools::toTitleCase(v),'only GWR:',response),label,
   file.path(fig,paste0(v,'-',response,if(per_sd) '-slope.png' else '-slope-per-cover-pp.png')),TRUE,subtitle)
 }
 map_three(lapply(1:3,function(a)d[d$atlas==a,]),'local_r2',c(0,1),
   paste(tools::toTitleCase(v),'only GWR:',response),'Local R² (common 0–1 scale)',
   file.path(fig,paste0(v,'-',response,'-local-r2.png')))
}
for(v in c('forest','cropland')){
 ds<-lapply(frames,function(d){d$coverage<-100*d[[paste0('perc_',v)]];d})
 map_three(ds,'coverage',c(0,100),paste('Input:',v,'coverage'),'Grid-cell cover (%)',file.path(fig,paste0('input-',v,'.png')))
}
for(response in c('richness','turnover')){
 map_three(frames,response,range(inputs[[response]]),paste('Input: observed',response),
   if(response=='richness') 'Species count' else 'Mean Bray–Curtis dissimilarity',
   file.path(fig,paste0('input-',response,'.png')),subtitle=if(response=='turnover') 'Dissimilarity to all other cells within the same atlas' else paste(length(species),'fixed species;',length(cells),'shared cells'))
}
writeLines(capture.output(sessionInfo()),file.path(out,'session-info.txt'))
writeLines(c('Observed single-predictor GWR, 22 October 2026','Sources: data-new/derived/occurrences/Y_occurrences.csv; environmental/X_Environmental.csv; design/studyDesign.csv; agreed shared-cell manifest',
 paste('Constant sample:',length(cells),'cells and',length(species),'species in every atlas. Good/Average effort in Atlas 1 and 2, >25% land; Atlas 3 effort assumed Good.'),
 'Observed richness = row sum of full occurrence table restricted to the shared cells and fixed species list.',
 'Turnover = mean Bray-Curtis dissimilarity to all other sites within each atlas on the shared grid (compositional distinctiveness); recomputed after filtering.',
 'Models fitted per percentage point of cover. Main slope maps multiply slopes by predictor-specific pooled SD across all cells and atlases; see predictor-scaling.csv. Supplemental maps show per-cover-percentage-point slopes. Turnover effects multiplied by 100 to express dissimilarity percentage points; richness remains species counts.',
 'Separate intercept + one predictor models; fixed 100 km bisquare kernel; coordinates projected from study design to EPSG:25832 in metres.',
 'No covariate adjustment or causal interpretation. Slopes mapped without significance filtering.',
 'Raw local t values saved for exploration, not multiplicity-adjusted significance claims.',
 'All maps include Bornholm. Slope limits are identical across forest/cropland and all atlases within response and scaling; see legend-scales.csv. All local R2 scales are 0–1; both habitat inputs use 0–100%.',
 paste('Figures:',fig)),file.path(out,'README.txt'))
message('Complete: ',out,'; ',fig)
