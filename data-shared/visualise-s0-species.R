# PNG summary using ONLY the portable S0 RData (no raw grid shapefile).
# Rscript data-shared/visualise-s0-species.R
# Optional: --species=Phylloscopus_collybita --output-dir=/path/to/output
visualise_s0_species <- function(
    rdata_file='data-shared/preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData',
    species='Phylloscopus_collybita',
    comparisons=c('Willow warbler'='Phylloscopus_trochilus',
                  'Blackcap'='Sylvia_atricapilla',
                  'Green woodpecker'='Picus_viridis'),
    output_dir='data-shared') {
  for(package in c('ggplot2','cowplot','ape'))
    if(!requireNamespace(package,quietly=TRUE))stop('Install ',package)
  e<-new.env(parent=emptyenv());load(rdata_file,envir=e)
  stopifnot(all(c('X','Y','Tr','Design','phy') %in% ls(e)),
            identical(rownames(e$Y),rownames(e$Design)),
            identical(rownames(e$X),rownames(e$Y)),identical(colnames(e$Y),rownames(e$Tr)))
  if(!species %in% colnames(e$Y))stop('Species not retained in S0: ',species)
  if(any(!c(species,unname(comparisons)) %in% e$phy$tip.label))stop('A comparison taxon is missing from the bundled tree')
  d<-e$Design
  stopifnot(all(is.finite(d$lon)),all(is.finite(d$lat)),all(e$Y[[species]] %in% c(0,1)))
  # These columns contain projected EPSG:23032 coordinates in metres, despite their legacy names.
  # Show nominal 5-km square footprints, not a reconstructed coastline or survey outside the model input.
  sites<-unique(data.frame(site=as.character(d$site),x=d$lon,y=d$lat))
  if(anyDuplicated(sites$site))stop('Sites have inconsistent coordinates across atlases')
  atlas_labels<-c('Atlas 1 · 1971–1974','Atlas 2 · 1993–1996','Atlas 3 · 2014–2017')
  maps<-lapply(1:3,function(a){
    ii<-which(as.character(d$atlas)==as.character(a))
    z<-sites;j<-match(z$site,as.character(d$site[ii]));valid<-!is.na(j)
    z$status<-'Not included'
    z$status[valid]<-ifelse(e$Y[[species]][ii[j[valid]]]==1,'Presence','Absence')
    z$status<-factor(z$status,levels=c('Presence','Absence','Not included'))
    z$atlas<-atlas_labels[a];z
  })
  map_data<-do.call(rbind,maps);map_data$atlas<-factor(map_data$atlas,levels=atlas_labels)
  counts<-data.frame(atlas=1:3,sites=sapply(1:3,function(a)sum(as.character(d$atlas)==a)),
                     presences=sapply(1:3,function(a)sum(e$Y[[species]][as.character(d$atlas)==a])))
  counts$absences<-counts$sites-counts$presences
  map_data$atlas<-factor(as.integer(map_data$atlas),levels=1:3,
    labels=paste0(atlas_labels,'\n',counts$presences,' presences / ',counts$sites,' model cells'))
  map<-ggplot2::ggplot(map_data,ggplot2::aes(x,y,fill=status))+
    ggplot2::geom_tile(width=5000,height=5000,colour=NA)+
    ggplot2::facet_wrap(~atlas,nrow=1)+ggplot2::coord_equal()+
    ggplot2::scale_fill_manual(values=c(Presence='#167a53',Absence='#fff3d4','Not included'='#d3d7dc'),drop=FALSE,name=NULL)+
    ggplot2::labs(title=paste0(if(species=='Phylloscopus_collybita')'Common chiffchaff'else gsub('_',' ',species),' — S0 presence data'),
                  subtitle=paste(gsub('_',' ',species),'| one 5-km cell per square | the filtered S0 inputs'))+
    ggplot2::theme_void(base_size=12)+ggplot2::theme(legend.position='bottom',strip.text=ggplot2::element_text(face='bold',size=12),
      plot.title=ggplot2::element_text(face='bold',size=19),plot.subtitle=ggplot2::element_text(size=12),
      plot.margin=ggplot2::margin(14,15,4,15))
  trait<-e$Tr[species,,drop=FALSE]
  migration_labels<-c('1'='Resident','1.5'='Resident and short-distance migrant',
                      '2'='Short-distance migrant','2.5'='Short- and long-distance migrant','3'='Long-distance migrant')
  code<-as.character(trait$Migration_a3_DOF)
  migration<-unname(migration_labels[code]);if(is.na(migration))migration<-'Unknown'
  trait_text<-paste0('Foraging guild\n',trait$foraging_guild_consensus,'\n\nMigratory strategy\n',migration,
                     ' (code ',code,')\n\nSpecies thermal index (STI)\n',format(trait$species_thermal_index,digits=8))
  traits_panel<-cowplot::ggdraw()+cowplot::draw_label('Species traits',x=.04,y=.97,hjust=0,vjust=1,fontface='bold',size=15)+
    cowplot::draw_label(trait_text,x=.04,y=.84,hjust=0,vjust=1,size=12,lineheight=1.05)
  pd<-ape::cophenetic.phylo(e$phy)
  distances<-data.frame(comparison=names(comparisons),species=unname(comparisons),
                        taxonomic_distance=as.numeric(pd[species,unname(comparisons)]))
  distances$label<-paste0(distances$comparison,'\n',gsub('_',' ',distances$species))
  distances$label<-factor(distances$label,levels=rev(distances$label))
  distance_panel<-ggplot2::ggplot(distances,ggplot2::aes(taxonomic_distance,label))+
    ggplot2::geom_col(fill='#3c6384',width=.55)+ggplot2::geom_text(ggplot2::aes(label=taxonomic_distance),hjust=-.4,size=4.5)+
    ggplot2::scale_x_continuous(limits=c(0,max(distances$taxonomic_distance)*1.12),breaks=seq(0,8,2),expand=c(0,0))+
    ggplot2::labs(title='Taxonomic distance from the focal species',x='IOC rank-distance units',y=NULL)+
    ggplot2::theme_minimal(base_size=12)+ggplot2::theme(panel.grid.major.y=ggplot2::element_blank(),panel.grid.minor=ggplot2::element_blank(),
      plot.title=ggplot2::element_text(face='bold',size=14),axis.text.y=ggplot2::element_text(size=11),plot.margin=ggplot2::margin(12,22,10,12))
  lower<-cowplot::plot_grid(traits_panel,distance_panel,nrow=1,rel_widths=c(1,1.5))
  footer<-cowplot::ggdraw()+cowplot::draw_label(paste(
    'Maps use the actual S0 samples: cell sets differ between atlases. Grey = no row in that atlas model; it is not an absence.',
    'Nominal squares from EPSG:23032 cell centres. Distances come from the bundled taxonomic surrogate, not evolutionary time.',
    'Rank distances: same genus/different species = 2; different families within an order = 6; different orders = 8.',sep='\n'),
    x=.015,y=.9,hjust=0,vjust=1,size=10,lineheight=1.15)
  figure<-cowplot::plot_grid(map,lower,footer,ncol=1,rel_heights=c(2.7,1.35,.4))
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  filename<-file.path(output_dir,paste0(species,'-s0-summary.png'))
  ggplot2::ggsave(filename,figure,width=15,height=12,dpi=220,bg='white')
  write.csv(counts,file.path(output_dir,paste0(species,'-presence-counts.csv')),row.names=FALSE)
  write.csv(distances[c('comparison','species','taxonomic_distance')],file.path(output_dir,paste0(species,'-taxonomic-distances.csv')),row.names=FALSE)
  write.csv(data.frame(species,foraging_guild=trait$foraging_guild_consensus,migratory_strategy=trait$Migration_a3_DOF,
                       migration_label=migration,STI=trait$species_thermal_index),file.path(output_dir,paste0(species,'-traits.csv')),row.names=FALSE)
  if(interactive())print(figure)
  message('Saved ',normalizePath(filename))
  invisible(list(plot=figure,counts=counts,traits=trait,distances=distances,map_data=map_data,file=filename))
}
if(sys.nframe()==0L){
  args<-commandArgs(trailingOnly=TRUE)
  value<-function(key,default){v<-grep(paste0('^--',key,'='),args,value=TRUE);if(length(v))sub(paste0('^--',key,'='),'',v[1])else default}
  visualise_s0_species(species=value('species','Phylloscopus_collybita'),
                       rdata_file=value('rdata','data-shared/preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData'),
                       output_dir=value('output-dir','data-shared'))
}
